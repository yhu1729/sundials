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
 * Mass matrix scaling test for FIRKODE with iterative linear solvers: the P1
 * finite element heat problem of firk_test_mass_fem,
 * u_t = u_xx + (pi^2 - 1) e^{-t} sin(pi x) on [0, 1] with homogeneous
 * Dirichlet conditions, c M u' = c (-K u + M f(t)), is integrated with the
 * s-stage Radau IIA method. The interior rows of M are of order
 * dx = 1/(N - 1) and its boundary rows are those of the identity. The factor
 * c (1 or 1e6) does not change the solution; c = 1e6 gives entries of order
 * 1e6 dx, as for a heat equation in physical units with the volumetric heat
 * capacity in M. The semi-discrete solution is known exactly (see
 * firk_test_mass_fem), so the error of the time integration alone is measured.
 *
 * The user's linear solver is the matrix-free SPGMR solver with the lumped
 * mass preconditioner P = c (M_L + gamma K) applied on the left or on the
 * right, and the mass solves use the matrix-free PCG solver. SPGMR solves the
 * error estimate and, for s = 1, the Newton systems (c M - gamma c J) x = b;
 * for s > 1 the Newton systems are solved by the flexible GMRES iteration of
 * FIRKODE on the coupled stage system, with P applied to each stage. The
 * right-hand sides and the Krylov residuals have the units of c M u, while the
 * tolerances of these solves are in the units of u. If they were compared
 * directly, the Newton iteration would accept the predictor on fine meshes
 * (c dx small) and the step size would collapse, while for c dx large the
 * tolerances would be tighter than intended. The conversion must also hold
 * for left preconditioning, where the Krylov solver tests the preconditioned
 * residual, which already has the units of u. The reference is the
 * integration with c = 1 and banded direct solvers, which have no
 * tolerances: for both preconditioning sides, both values of c and both
 * meshes, the number of steps must match that of the reference, the error
 * must not exceed that of the reference by more than the tolerance, and the
 * number of linear iterations must not depend on c.
 *
 * With the argument "step", the initial data is the step function u = 1 on
 * [0.25, 0.75], u = 0 elsewhere, and the errors are measured against an
 * integration with banded direct solvers at a tight tolerance. At t = 0,
 * f = c (-K u + M f(0)) then has entries of order c/dx at the
 * discontinuities, so y' = M^{-1} f, which the initial step size estimate
 * computes with mass solves, has entries of order 1/dx^2. These mass solves
 * need a tolerance relative to ||f||: PCG could not reach a tolerance in the
 * units of u applied to the residual in the units of f, which stopped the
 * integration at t = 0.
 *
 * Pass criterion (not run in single precision): every integration succeeds;
 * the step count of each configuration differs from the reference's by at
 * most 20% and its error is at most 1.25 times the reference error plus
 * rtol; with c = 1e6 the iteration counts of the user's solver and of the
 * stacked solver differ from those with c = 1 by at most 20%.
 *
 * Usage: firk_test_mass_fem_krylov [stages] [step]
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "firkode/firkode.h"
#include "firkode/firkode_impl.h" /* firkSetNumStagesAny */
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "sunlinsol/sunlinsol_band.h"
#include "sunlinsol/sunlinsol_pcg.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunmatrix/sunmatrix_band.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define ESYM "Le"
#else
#define ESYM "e"
#endif

namespace {

constexpr sunrealtype zero = SUN_RCONST(0.0);
constexpr sunrealtype half = SUN_RCONST(0.5);
constexpr sunrealtype one  = SUN_RCONST(1.0);
constexpr sunrealtype two  = SUN_RCONST(2.0);
constexpr sunrealtype four = SUN_RCONST(4.0);
constexpr sunrealtype six  = SUN_RCONST(6.0);
constexpr sunrealtype pi   = SUN_RCONST(3.141592653589793238462643383279502884);
constexpr sunrealtype tf   = SUN_RCONST(1.0);

struct Problem
{
  sunindextype N; // number of nodes
  sunrealtype dx; // mesh spacing
  sunrealtype c;  // factor multiplying both sides of the equations

  // LU factorization of the interior rows of P = c (M_L + gamma K) by the
  // Thomas algorithm: off-diagonal entry, multipliers of the subdiagonal and
  // inverse pivots, for the gamma of the last preconditioner setup
  sunrealtype off;
  std::vector<sunrealtype> mult, pivinv;

  Problem(sunindextype n, sunrealtype factor)
    : N(n),
      dx(one / static_cast<sunrealtype>(n - 1)),
      c(factor),
      off(zero),
      mult(static_cast<size_t>(n), zero),
      pivinv(static_cast<size_t>(n), zero)
  {}
};

// f = c (-K u + M f_src(t)) in the interior rows, zero in the boundary rows
int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  const Problem* prob  = static_cast<const Problem*>(user_data);
  const sunindextype N = prob->N;
  const sunrealtype dx = prob->dx;
  const sunrealtype* Y = N_VGetArrayPointer(y);
  sunrealtype* Ydot    = N_VGetArrayPointer(ydot);
  const sunrealtype a  = (pi * pi - one) * std::exp(-t);

  Ydot[0]     = zero;
  Ydot[N - 1] = zero;
  for (sunindextype i = 1; i < N - 1; i++)
  {
    const sunrealtype fl = a *
                           std::sin(pi * dx * static_cast<sunrealtype>(i - 1));
    const sunrealtype fc = a * std::sin(pi * dx * static_cast<sunrealtype>(i));
    const sunrealtype fr = a *
                           std::sin(pi * dx * static_cast<sunrealtype>(i + 1));
    Ydot[i] = prob->c * (-(-Y[i - 1] + two * Y[i] - Y[i + 1]) / dx +
                         dx * (fl + four * fc + fr) / six);
  }
  return 0;
}

// J = -c K (band, mu = ml = 1)
int jac(sunrealtype, N_Vector, N_Vector, SUNMatrix J, void* user_data, N_Vector,
        N_Vector, N_Vector)
{
  const Problem* prob = static_cast<const Problem*>(user_data);
  SUNMatZero(J);
  for (sunindextype i = 1; i < prob->N - 1; i++)
  {
    SM_ELEMENT_B(J, i, i - 1) = prob->c / prob->dx;
    SM_ELEMENT_B(J, i, i)     = -two * prob->c / prob->dx;
    SM_ELEMENT_B(J, i, i + 1) = prob->c / prob->dx;
  }
  return 0;
}

// J v = -c K v
int jtimes(N_Vector v, N_Vector Jv, sunrealtype, N_Vector, N_Vector,
           void* user_data, N_Vector)
{
  const Problem* prob  = static_cast<const Problem*>(user_data);
  const sunindextype N = prob->N;
  const sunrealtype* V = N_VGetArrayPointer(v);
  sunrealtype* JV      = N_VGetArrayPointer(Jv);

  JV[0]     = zero;
  JV[N - 1] = zero;
  for (sunindextype i = 1; i < N - 1; i++)
  {
    JV[i] = prob->c * (V[i - 1] - two * V[i] + V[i + 1]) / prob->dx;
  }
  return 0;
}

// The P1 mass matrix c M with boundary rows c I. The boundary values of the
// solution and of all corrections are zero, so the coupling of the first and
// last interior rows to the boundary values is dropped; this keeps M
// symmetric for PCG without changing the solution.

// c M (band, mu = ml = 1)
int mass(sunrealtype, SUNMatrix M, void* user_data, N_Vector, N_Vector, N_Vector)
{
  const Problem* prob  = static_cast<const Problem*>(user_data);
  const sunindextype N = prob->N;
  const sunrealtype m  = prob->c * prob->dx / six;
  SUNMatZero(M);
  SM_ELEMENT_B(M, 0, 0)         = prob->c;
  SM_ELEMENT_B(M, N - 1, N - 1) = prob->c;
  for (sunindextype i = 1; i < N - 1; i++)
  {
    if (i > 1) { SM_ELEMENT_B(M, i, i - 1) = m; }
    SM_ELEMENT_B(M, i, i) = four * m;
    if (i < N - 2) { SM_ELEMENT_B(M, i, i + 1) = m; }
  }
  return 0;
}

// c M v
int mass_times(N_Vector v, N_Vector Mv, sunrealtype, void* mtimes_data)
{
  const Problem* prob  = static_cast<const Problem*>(mtimes_data);
  const sunindextype N = prob->N;
  const sunrealtype m  = prob->c * prob->dx / six;
  const sunrealtype* V = N_VGetArrayPointer(v);
  sunrealtype* MV      = N_VGetArrayPointer(Mv);

  MV[0]     = prob->c * V[0];
  MV[N - 1] = prob->c * V[N - 1];
  for (sunindextype i = 1; i < N - 1; i++)
  {
    const sunrealtype vl = (i > 1) ? V[i - 1] : zero;
    const sunrealtype vr = (i < N - 2) ? V[i + 1] : zero;
    MV[i]                = m * (vl + four * V[i] + vr);
  }
  return 0;
}

// preconditioner setup: factor the interior rows of P = c (M_L + gamma K),
// where M_L = dx I is the lumped mass matrix; the boundary rows of P are those
// of c I. P is spectrally equivalent to c (M + gamma K) uniformly in dx.
int psetup(sunrealtype, N_Vector, N_Vector, sunbooleantype jok,
           sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data)
{
  Problem* prob         = static_cast<Problem*>(user_data);
  const sunindextype N  = prob->N;
  const sunrealtype dia = prob->c * (prob->dx + two * gamma / prob->dx);

  prob->off       = -prob->c * gamma / prob->dx;
  prob->pivinv[1] = one / dia;
  for (sunindextype i = 2; i < N - 1; i++)
  {
    prob->mult[i]   = prob->off * prob->pivinv[i - 1];
    prob->pivinv[i] = one / (dia - prob->mult[i] * prob->off);
  }
  *jcurPtr = !jok;
  return 0;
}

// preconditioner solve P z = r with the factorization of the last setup, whose
// gamma may lag behind the current one
int psolve(sunrealtype, N_Vector, N_Vector, N_Vector r, N_Vector z, sunrealtype,
           sunrealtype, int, void* user_data)
{
  const Problem* prob  = static_cast<const Problem*>(user_data);
  const sunindextype N = prob->N;
  const sunrealtype* R = N_VGetArrayPointer(r);
  sunrealtype* Z       = N_VGetArrayPointer(z);

  Z[0]     = R[0] / prob->c;
  Z[N - 1] = R[N - 1] / prob->c;
  Z[1]     = R[1];
  for (sunindextype i = 2; i < N - 1; i++)
  {
    Z[i] = R[i] - prob->mult[i] * Z[i - 1];
  }
  Z[N - 2] *= prob->pivinv[N - 2];
  for (sunindextype i = N - 3; i >= 1; i--)
  {
    Z[i] = (Z[i] - prob->off * Z[i + 1]) * prob->pivinv[i];
  }
  return 0;
}

// amplitude a(t) of the exact semi-discrete solution: K v = mu M v for
// v_i = sin(pi x_i), mu = 12 sin^2(pi dx / 2) / (dx^2 (2 + cos(pi dx)))
sunrealtype amplitude(sunrealtype t, sunrealtype dx)
{
  const sunrealtype sn = std::sin(half * pi * dx);
  const sunrealtype mu = SUN_RCONST(12.0) * sn * sn /
                         (dx * dx * (two + std::cos(pi * dx)));
  const sunrealtype c = (pi * pi - one) / (mu - one);
  return c * std::exp(-t) + (one - c) * std::exp(-mu * t);
}

struct Config
{
  const char* name;
  bool band;     // banded direct solvers, otherwise SPGMR and PCG
  int pretype;   // SPGMR preconditioning side
  sunrealtype c; // factor multiplying the equations
};

struct Result
{
  long int nst, nni, nli, nsli, nmsolves;
  sunrealtype err;            // max-norm error at tf
  std::vector<sunrealtype> y; // solution at tf
};

// integrate from the sine (step = false) or step-function initial data; the
// error is measured against the exact semi-discrete solution or, for step
// data, against yref, unless yref is empty
int run(SUNContext sunctx, const Config& cfg, int s, sunindextype N,
        sunrealtype rtol, bool step, const std::vector<sunrealtype>& yref,
        Result* res)
{
  Problem prob(N, cfg.c);

  N_Vector y      = N_VNew_Serial(N, sunctx);
  sunrealtype* yd = N_VGetArrayPointer(y);
  for (sunindextype i = 0; i < N; i++)
  {
    if (step) { yd[i] = (4 * i >= N - 1 && 4 * i <= 3 * (N - 1)) ? one : zero; }
    else { yd[i] = std::sin(pi * prob.dx * static_cast<sunrealtype>(i)); }
  }

  void* mem = FIRKodeCreate(sunctx);
  int flag  = FIRKodeInit(mem, rhs, zero, y);
  flag += FIRKodeSetUserData(mem, &prob);
  flag += firkSetNumStagesAny(mem, s);
  flag += FIRKodeSStolerances(mem, rtol, SUN_RCONST(1.0e-4) * rtol);

  SUNMatrix A         = nullptr;
  SUNMatrix M         = nullptr;
  SUNLinearSolver LS  = nullptr;
  SUNLinearSolver MLS = nullptr;
  if (cfg.band)
  {
    A   = SUNBandMatrix(N, 1, 1, sunctx);
    LS  = SUNLinSol_Band(y, A, sunctx);
    M   = SUNBandMatrix(N, 1, 1, sunctx);
    MLS = SUNLinSol_Band(y, M, sunctx);
    flag += FIRKodeSetLinearSolver(mem, LS, A);
    flag += FIRKodeSetJacFn(mem, jac);
    flag += FIRKodeSetMassLinearSolver(mem, MLS, M, SUNFALSE);
    flag += FIRKodeSetMassFn(mem, mass);
  }
  else
  {
    LS  = SUNLinSol_SPGMR(y, cfg.pretype, 0, sunctx);
    MLS = SUNLinSol_PCG(y, SUN_PREC_NONE, 20, sunctx);
    flag += FIRKodeSetLinearSolver(mem, LS, nullptr);
    flag += FIRKodeSetJacTimes(mem, nullptr, jtimes);
    flag += FIRKodeSetPreconditioner(mem, psetup, psolve);
    flag += FIRKodeSetMassLinearSolver(mem, MLS, nullptr, SUNFALSE);
    flag += FIRKodeSetMassTimes(mem, nullptr, mass_times, &prob);
  }
  flag += FIRKodeSetStopTime(mem, tf);
  flag += FIRKodeSetMaxNumSteps(mem, 20000);

  sunrealtype t = zero;
  if (flag == 0) { flag = FIRKodeEvolve(mem, tf, y, &t, FIRK_NORMAL); }

  if (flag == FIRK_SUCCESS || flag == FIRK_TSTOP_RETURN)
  {
    flag                = 0;
    const sunrealtype a = amplitude(t, prob.dx);
    res->err            = zero;
    for (sunindextype i = 0; i < N; i++)
    {
      if (step && yref.empty()) { break; }
      const sunrealtype exact =
        step ? yref[static_cast<size_t>(i)]
             : a * std::sin(pi * prob.dx * static_cast<sunrealtype>(i));
      res->err = std::max(res->err, std::abs(yd[i] - exact));
    }
    res->y.assign(yd, yd + N);
    FIRKodeGetNumSteps(mem, &res->nst);
    FIRKodeGetNumNonlinSolvIters(mem, &res->nni);
    FIRKodeGetNumLinIters(mem, &res->nli);
    FIRKodeGetNumStageLinIters(mem, &res->nsli);
    FIRKodeGetNumMassSolves(mem, &res->nmsolves);
  }

  FIRKodeFree(&mem);
  SUNLinSolFree(MLS);
  SUNLinSolFree(LS);
  if (M) { SUNMatDestroy(M); }
  if (A) { SUNMatDestroy(A); }
  N_VDestroy(y);
  return flag;
}

// true if a and b differ by more than 20% of b
bool differ(long int a, long int b) { return 5 * std::labs(a - b) > b; }

} // namespace

int main(int argc, char* argv[])
{
  const int s     = (argc > 1) ? std::atoi(argv[1]) : 1;
  const bool step = (argc > 2) && (std::strcmp(argv[2], "step") == 0);
  int fails       = 0;

  SUNContext sunctx;
  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  std::printf("Mass matrix scaling test, P1 FEM heat equation with Radau IIA "
              "s = %d and iterative linear solvers%s\n",
              s, step ? ", step-function initial data" : "");

  const sunindextype meshes[2] = {401, 1601};
  const sunrealtype rtols[2]   = {SUN_RCONST(1.0e-4), SUN_RCONST(1.0e-6)};

  // the reference first, then each SPGMR configuration with c = 1 followed by
  // the same configuration with c = 1e6
  const sunrealtype big = SUN_RCONST(1.0e6);
  const Config cfgs[5]  = {{"band", true, SUN_PREC_NONE, one},
                           {"SPGMR left", false, SUN_PREC_LEFT, one},
                           {"SPGMR left, c = 1e6", false, SUN_PREC_LEFT, big},
                           {"SPGMR right", false, SUN_PREC_RIGHT, one},
                           {"SPGMR right, c = 1e6", false, SUN_PREC_RIGHT, big}};
  constexpr int ncfgs   = 5;

  // for step data, reference solutions with band solvers, s = 5 and a tight
  // tolerance, accurate to about 1e-13
  std::vector<sunrealtype> yrefs[2];
  for (int m = 0; step && m < 2; m++)
  {
    Result ref     = {};
    const int flag = run(sunctx, cfgs[0], 5, meshes[m], SUN_RCONST(1.0e-9),
                         step, {}, &ref);
    if (flag != 0)
    {
      std::printf("FAIL: the reference integration failed for N = %li\n",
                  static_cast<long int>(meshes[m]));
      SUNContext_Free(&sunctx);
      return 1;
    }
    yrefs[m] = ref.y;
  }

  for (const sunrealtype rtol : rtols)
  {
    std::printf("rtol = %.0" ESYM ", atol = %.0" ESYM ":\n", rtol,
                SUN_RCONST(1.0e-4) * rtol);
    for (int m = 0; m < 2; m++)
    {
      const sunindextype N = meshes[m];
      Result res[ncfgs]    = {};
      bool ran[ncfgs];
      for (int k = 0; k < ncfgs; k++)
      {
        ran[k] = (run(sunctx, cfgs[k], s, N, rtol, step, yrefs[m], &res[k]) == 0);
        if (!ran[k])
        {
          std::printf("  FAIL: integration with %s failed for N = %li\n",
                      cfgs[k].name, static_cast<long int>(N));
          fails++;
          continue;
        }
        std::printf("  N = %4li, %-20s: nst = %li, nni = %li, linear iters = "
                    "%li, ",
                    static_cast<long int>(N), cfgs[k].name, res[k].nst,
                    res[k].nni, res[k].nli);
        if (s > 1) { std::printf("stage linear iters = %li, ", res[k].nsli); }
        std::printf("mass solves = %li, max error = %.1" ESYM "\n",
                    res[k].nmsolves, res[k].err);
      }
      if (!ran[0]) { continue; }

      for (int k = 1; k < ncfgs; k++)
      {
        if (!ran[k]) { continue; }

        // with linear solves to the requested tolerances, the integration
        // follows the reference; the inexact solves may change the step
        // sequence slightly, which for s > 1, with errors far below the
        // tolerance, changes the error by a large factor relative to the
        // reference but not relative to the tolerance (|u| <= 1)
        if (differ(res[k].nst, res[0].nst))
        {
          std::printf("  FAIL: %s takes %li steps against %li with band "
                      "solvers\n",
                      cfgs[k].name, res[k].nst, res[0].nst);
          fails++;
        }
        if (res[k].err > SUN_RCONST(1.25) * res[0].err + rtol)
        {
          std::printf("  FAIL: %s error %.1" ESYM " against %.1" ESYM
                      " with band solvers\n",
                      cfgs[k].name, res[k].err, res[0].err);
          fails++;
        }

        // multiplying the equations by c changes no tolerance
        if (cfgs[k].c != one && ran[k - 1] &&
            (differ(res[k].nli, res[k - 1].nli) ||
             differ(res[k].nsli, res[k - 1].nsli)))
        {
          std::printf("  FAIL: %s takes %li + %li linear iterations against "
                      "%li + %li with c = 1\n",
                      cfgs[k].name, res[k].nli, res[k].nsli, res[k - 1].nli,
                      res[k - 1].nsli);
          fails++;
        }
      }
    }
  }

  SUNContext_Free(&sunctx);

  if (fails)
  {
    std::printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  std::printf("\nAll tests passed\n");
  return 0;
}
