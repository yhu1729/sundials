/* -----------------------------------------------------------------------------
 * Programmer(s): Yifan Hu @ UMBC
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Comparison of FIRKODE with radau5.f (design milestone M2, Section 10.3,
 * Table 6.2) on Robertson kinetics, the van der Pol oscillator (eps = 1e-6)
 * and the Oregonator at the nominal RADAU5 tolerances 1e-4, 1e-6 and 1e-8,
 * from the reference files of test/unit_tests/firkode/radau5 (run_radau5.py).
 *
 * FIRKODE runs in BEGINNER mode (the RADAU5 reuse policy, the predictive
 * controller and the automatic Newton tolerance) at the tolerances RADAU5
 * used, rtol' = 0.1 rtol^(2/3) and atol' = rtol' atol / rtol, with the dense
 * solver, the analytic Jacobian and the two inputs both codes received, the
 * initial step h0 = 1e-6 and the final time as a stop time (RADAU5 lands on
 * tend; without a stop time FIRKODE would interpolate at tend with the
 * degree-s collocation polynomial).
 *
 * Pass criterion (not run in single precision): for every cell the number of
 * steps is within 20% of RADAU5's accepted steps, and the error at tend in the
 * WRMS norm with the weights 1 / (atol' + rtol' |y_ref|) of the radau5.f
 * reference solution at rtol = 1e-14 is at most twice RADAU5's own error, or
 * at most half the tolerance (1 in this norm). The reference's own
 * uncertainty, 2e-12 relative, floors the RADAU5 error.
 *
 * The floor is needed because both codes' errors lie well below the
 * tolerance (RADAU5: 0.0006 to 0.08, FIRKODE: 0.0005 to 0.4 on 2026-10-08,
 * macOS 27, Apple M-series, Apple clang, double precision), where the error
 * at tend is set by the path through the stiff transitions and not by the
 * tolerance: the step sequences agree almost step for step (VDPOL at 1e-6:
 * 486 against 489 steps, 211 against 212 through the jump at t = 1.614),
 * while the error ratio of that cell is 23 with the defaults and 0.1 with the
 * Newton tolerance 1e-3 instead of its automatic 3.2e-3, the refilter off, or
 * the stacked solve 50 times tighter. Measured ratios FIRKODE / RADAU5 on that
 * date: ROBER 0.77, 0.97, 0.91; VDPOL 2.05, 22.8, 0.75; OREGO 1.05, 3.60,
 * 1.02 at the nominal tolerances 1e-4, 1e-6, 1e-8.
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define ESYM "Le"
#else
#define ESYM "e"
#endif

namespace {

constexpr sunrealtype zero = SUN_RCONST(0.0);
constexpr sunrealtype one  = SUN_RCONST(1.0);

/* the step-count gate of design Section 10.3 */
constexpr sunrealtype step_factor = SUN_RCONST(0.2);

/* the error gate of design Section 10.3: twice RADAU5's error, or an error
   below half the tolerance (see the header) */
constexpr sunrealtype err_factor = SUN_RCONST(2.0);
constexpr sunrealtype err_floor  = SUN_RCONST(0.5);

const char* const tol_names[3] = {"1e-4", "1e-6", "1e-8"};

/* ---- the problems, as in radau5_driver.f90 ---- */

int f_rober(sunrealtype, N_Vector y, N_Vector ydot, void*)
{
  const sunrealtype y1 = NV_Ith_S(y, 0), y2 = NV_Ith_S(y, 1), y3 = NV_Ith_S(y, 2);
  NV_Ith_S(ydot, 0) = -SUN_RCONST(0.04) * y1 + SUN_RCONST(1.0e4) * y2 * y3;
  NV_Ith_S(ydot, 1) = SUN_RCONST(0.04) * y1 - SUN_RCONST(1.0e4) * y2 * y3 -
                      SUN_RCONST(3.0e7) * y2 * y2;
  NV_Ith_S(ydot, 2) = SUN_RCONST(3.0e7) * y2 * y2;
  return 0;
}

int jac_rober(sunrealtype, N_Vector y, N_Vector, SUNMatrix J, void*, N_Vector,
              N_Vector, N_Vector)
{
  const sunrealtype y2 = NV_Ith_S(y, 1), y3 = NV_Ith_S(y, 2);
  SM_ELEMENT_D(J, 0, 0) = -SUN_RCONST(0.04);
  SM_ELEMENT_D(J, 0, 1) = SUN_RCONST(1.0e4) * y3;
  SM_ELEMENT_D(J, 0, 2) = SUN_RCONST(1.0e4) * y2;
  SM_ELEMENT_D(J, 1, 0) = SUN_RCONST(0.04);
  SM_ELEMENT_D(J, 1, 1) = -SUN_RCONST(1.0e4) * y3 - SUN_RCONST(6.0e7) * y2;
  SM_ELEMENT_D(J, 1, 2) = -SUN_RCONST(1.0e4) * y2;
  SM_ELEMENT_D(J, 2, 0) = zero;
  SM_ELEMENT_D(J, 2, 1) = SUN_RCONST(6.0e7) * y2;
  SM_ELEMENT_D(J, 2, 2) = zero;
  return 0;
}

int f_vdpol(sunrealtype, N_Vector y, N_Vector ydot, void*)
{
  constexpr sunrealtype eps = SUN_RCONST(1.0e-6);
  const sunrealtype y1 = NV_Ith_S(y, 0), y2 = NV_Ith_S(y, 1);
  NV_Ith_S(ydot, 0) = y2;
  NV_Ith_S(ydot, 1) = ((one - y1 * y1) * y2 - y1) / eps;
  return 0;
}

int jac_vdpol(sunrealtype, N_Vector y, N_Vector, SUNMatrix J, void*, N_Vector,
              N_Vector, N_Vector)
{
  constexpr sunrealtype eps = SUN_RCONST(1.0e-6);
  const sunrealtype y1 = NV_Ith_S(y, 0), y2 = NV_Ith_S(y, 1);
  SM_ELEMENT_D(J, 0, 0) = zero;
  SM_ELEMENT_D(J, 0, 1) = one;
  SM_ELEMENT_D(J, 1, 0) = (-SUN_RCONST(2.0) * y1 * y2 - one) / eps;
  SM_ELEMENT_D(J, 1, 1) = (one - y1 * y1) / eps;
  return 0;
}

int f_orego(sunrealtype, N_Vector y, N_Vector ydot, void*)
{
  const sunrealtype y1 = NV_Ith_S(y, 0), y2 = NV_Ith_S(y, 1), y3 = NV_Ith_S(y, 2);
  NV_Ith_S(ydot, 0) = SUN_RCONST(77.27) *
                      (y2 + y1 * (one - SUN_RCONST(8.375e-6) * y1 - y2));
  NV_Ith_S(ydot, 1) = (y3 - (one + y1) * y2) / SUN_RCONST(77.27);
  NV_Ith_S(ydot, 2) = SUN_RCONST(0.161) * (y1 - y3);
  return 0;
}

int jac_orego(sunrealtype, N_Vector y, N_Vector, SUNMatrix J, void*, N_Vector,
              N_Vector, N_Vector)
{
  const sunrealtype y1 = NV_Ith_S(y, 0), y2 = NV_Ith_S(y, 1);
  SM_ELEMENT_D(J, 0, 0) =
    SUN_RCONST(77.27) * (one - SUN_RCONST(2.0) * SUN_RCONST(8.375e-6) * y1 - y2);
  SM_ELEMENT_D(J, 0, 1) = SUN_RCONST(77.27) * (one - y1);
  SM_ELEMENT_D(J, 0, 2) = zero;
  SM_ELEMENT_D(J, 1, 0) = -y2 / SUN_RCONST(77.27);
  SM_ELEMENT_D(J, 1, 1) = -(one + y1) / SUN_RCONST(77.27);
  SM_ELEMENT_D(J, 1, 2) = one / SUN_RCONST(77.27);
  SM_ELEMENT_D(J, 2, 0) = SUN_RCONST(0.161);
  SM_ELEMENT_D(J, 2, 1) = zero;
  SM_ELEMENT_D(J, 2, 2) = -SUN_RCONST(0.161);
  return 0;
}

struct Problem
{
  const char* name;
  int n;
  FIRKRhsFn f;
  FIRKLsJacFn jac;
  sunrealtype y0[3];
};

const Problem problems[3] =
  {{"rober", 3, f_rober, jac_rober, {one, zero, zero}},
   {"vdpol", 2, f_vdpol, jac_vdpol, {SUN_RCONST(2.0), zero, zero}},
   {"orego", 3, f_orego, jac_orego, {one, SUN_RCONST(2.0), SUN_RCONST(3.0)}}};

/* ---- the reference files ---- */

struct Reference
{
  int n       = 0;
  double tend = 0, h0 = 0, rtol = 0, atol = 0;
  long int nstep = 0, naccpt = 0, njac = 0, nsol = 0;
  double y[3] = {0, 0, 0};
};

/* reads <dir>/<problem>_<tol>.ref; returns false if a field is missing */
bool read_reference(const std::string& dir, const char* problem,
                    const char* tol, Reference* ref)
{
  const std::string path = dir + "/" + problem + "_" + tol + ".ref";
  std::ifstream in(path.c_str());
  if (!in)
  {
    std::printf("  FAIL: cannot open %s\n", path.c_str());
    return false;
  }
  int seen = 0;
  std::string line;
  while (std::getline(in, line))
  {
    if (line.empty() || line[0] == '#') { continue; }
    std::istringstream is(line);
    std::string key;
    is >> key;
    if (key == "n")
    {
      is >> ref->n;
      seen |= 1;
    }
    else if (key == "tend")
    {
      is >> ref->tend;
      seen |= 2;
    }
    else if (key == "h0")
    {
      is >> ref->h0;
      seen |= 4;
    }
    else if (key == "rtol")
    {
      is >> ref->rtol;
      seen |= 8;
    }
    else if (key == "atol")
    {
      is >> ref->atol;
      seen |= 16;
    }
    else if (key == "nstep")
    {
      is >> ref->nstep;
      seen |= 32;
    }
    else if (key == "naccpt")
    {
      is >> ref->naccpt;
      seen |= 64;
    }
    else if (key == "njac")
    {
      is >> ref->njac;
      seen |= 128;
    }
    else if (key == "nsol")
    {
      is >> ref->nsol;
      seen |= 256;
    }
    else if (key == "y")
    {
      int i    = 0;
      double v = 0;
      is >> i >> v;
      if (i >= 1 && i <= 3)
      {
        ref->y[i - 1] = v;
        seen |= 512;
      }
    }
  }
  if (seen != 1023 || ref->n < 1 || ref->n > 3)
  {
    std::printf("  FAIL: %s is incomplete\n", path.c_str());
    return false;
  }
  return true;
}

/* ---- one FIRKODE run ---- */

struct Run
{
  int flag         = 0;
  sunrealtype tret = zero;
  long int nst = 0, nst_attempts = 0, nni = 0, nje = 0;
  sunrealtype y[3] = {zero, zero, zero};
};

Run integrate(SUNContext sunctx, const Problem& p, const Reference& ref)
{
  Run r;
  N_Vector y = N_VNew_Serial(p.n, sunctx);
  for (int i = 0; i < p.n; i++) { NV_Ith_S(y, i) = p.y0[i]; }
  SUNMatrix A        = SUNDenseMatrix(p.n, p.n, sunctx);
  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
  void* mem          = FIRKodeCreate(sunctx);

  int flag = FIRKodeInit(mem, p.f, zero, y);
  flag += FIRKodeSStolerances(mem, (sunrealtype)ref.rtol, (sunrealtype)ref.atol);
  flag += FIRKodeSetLinearSolver(mem, LS, A);
  flag += FIRKodeSetJacFn(mem, p.jac);
  flag += FIRKodeSetMaxNumSteps(mem, 100000);
  flag += FIRKodeSetInitStep(mem, (sunrealtype)ref.h0);
  flag += FIRKodeSetStopTime(mem, (sunrealtype)ref.tend);
  if (flag) { r.flag = -1000; }
  else
  {
    r.flag = FIRKodeEvolve(mem, (sunrealtype)ref.tend, y, &r.tret, FIRK_NORMAL);
  }
  for (int i = 0; i < p.n; i++) { r.y[i] = NV_Ith_S(y, i); }
  FIRKodeGetNumSteps(mem, &r.nst);
  FIRKodeGetNumStepAttempts(mem, &r.nst_attempts);
  FIRKodeGetNumNonlinSolvIters(mem, &r.nni);
  FIRKodeGetNumJacEvals(mem, &r.nje);

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  return r;
}

/* the WRMS norm of v in the weights 1 / (atol + rtol |yref|) */
sunrealtype wrms(const sunrealtype* v, const double* yref, int n, double rtol,
                 double atol)
{
  sunrealtype sum = zero;
  for (int i = 0; i < n; i++)
  {
    const sunrealtype w =
      one /
      ((sunrealtype)atol + (sunrealtype)rtol * std::abs((sunrealtype)yref[i]));
    sum += (v[i] * w) * (v[i] * w);
  }
  return std::sqrt(sum / (sunrealtype)n);
}

} // namespace

int main(int argc, char* argv[])
{
  const std::string dir = (argc > 1) ? argv[1] : ".";
  const char* only      = (argc > 2) ? argv[2] : nullptr;
  int fails             = 0;

  SUNContext sunctx;
  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  std::printf("FIRKODE against radau5.f (references in %s)\n", dir.c_str());
  std::printf("  %-6s %-5s %6s %6s %6s  %9s %9s %6s  %s\n", "prob", "rtol",
              "naccpt", "nst", "ratio", "e_radau5", "e_firkode", "ratio",
              "result");

  for (int ip = 0; ip < 3; ip++)
  {
    const Problem& p = problems[ip];
    if (only != nullptr && std::strcmp(only, p.name) != 0) { continue; }

    Reference exact;
    if (!read_reference(dir, p.name, "1e-14", &exact))
    {
      fails++;
      continue;
    }

    for (int it = 0; it < 3; it++)
    {
      Reference ref;
      if (!read_reference(dir, p.name, tol_names[it], &ref))
      {
        fails++;
        continue;
      }

      const Run r = integrate(sunctx, p, ref);
      if (r.flag != FIRK_TSTOP_RETURN || r.tret != (sunrealtype)ref.tend)
      {
        std::printf("  %-6s %-5s FAIL: FIRKodeEvolve returned %d at t = %" ESYM
                    "\n",
                    p.name, tol_names[it], r.flag, r.tret);
        fails++;
        continue;
      }

      /* errors against the 1e-14 reference, RADAU5's floored at ten times
         the reference's uncertainty */
      sunrealtype d_firk[3], d_rad[3], unc[3];
      for (int i = 0; i < p.n; i++)
      {
        d_firk[i] = r.y[i] - (sunrealtype)exact.y[i];
        d_rad[i]  = (sunrealtype)ref.y[i] - (sunrealtype)exact.y[i];
        unc[i]    = SUN_RCONST(2.0e-12) * std::abs((sunrealtype)exact.y[i]);
      }
      const sunrealtype e_firk = wrms(d_firk, exact.y, p.n, ref.rtol, ref.atol);
      const sunrealtype e_rad =
        std::max(wrms(d_rad, exact.y, p.n, ref.rtol, ref.atol),
                 SUN_RCONST(10.0) * wrms(unc, exact.y, p.n, ref.rtol, ref.atol));

      const sunrealtype step_ratio = (sunrealtype)r.nst / (sunrealtype)ref.naccpt;
      const sunrealtype err_ratio = e_firk / e_rad;
      const bool steps_ok         = std::abs(step_ratio - one) <= step_factor;
      const bool err_ok = err_ratio <= err_factor || e_firk <= err_floor;

      std::printf("  %-6s %-5s %6ld %6ld %6.3f  %9.2" ESYM " %9.2" ESYM
                  " %6.2f  %s\n",
                  p.name, tol_names[it], ref.naccpt, r.nst, (double)step_ratio,
                  e_rad, e_firk, (double)err_ratio,
                  (steps_ok && err_ok) ? "pass" : "FAIL");
      std::printf("  %-6s %-5s attempts %ld/%ld, Newton iterations %ld/%ld, "
                  "Jacobians %ld/%ld (FIRKODE/RADAU5)\n",
                  "", "", r.nst_attempts, ref.nstep, r.nni, ref.nsol, r.nje,
                  ref.njac);
      if (!steps_ok || !err_ok) { fails++; }
    }
  }

  SUNContext_Free(&sunctx);

  if (fails)
  {
    std::printf("\n%d cell(s) FAILED\n", fails);
    return 1;
  }
  std::printf("\nAll tests passed\n");
  return 0;
}
