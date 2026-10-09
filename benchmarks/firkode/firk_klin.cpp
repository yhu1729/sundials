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
 * The stacked Krylov iterations per Newton iteration of FIRKODE (design
 * Section 11.4) on Robertson kinetics, the van der Pol oscillator, the
 * Oregonator and the 2-D heat equation; see README.md.
 * ---------------------------------------------------------------------------*/

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "firkode/firkode.h"
#include "nvector/nvector_serial.h"
#include "sunlinsol/sunlinsol_band.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunmatrix/sunmatrix_band.h"
#include "sunmatrix/sunmatrix_dense.h"

namespace {

constexpr sunrealtype zero = SUN_RCONST(0.0);
constexpr sunrealtype one  = SUN_RCONST(1.0);

/* ---- the small stiff problems (dense block solver) ---- */

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

/* ---- the 2-D heat equation on the unit square, interior nodes only ---- */

struct Heat
{
  sunindextype nx;
  sunrealtype dx2inv; /* 1 / dx^2 */
  sunrealtype gamma;  /* of the last preconditioner setup */
};

/* u_t = u_xx + u_yy, 5-point stencil, homogeneous Dirichlet data */
int f_heat(sunrealtype, N_Vector y, N_Vector ydot, void* data)
{
  const Heat* h        = static_cast<Heat*>(data);
  const sunindextype n = h->nx;
  const sunrealtype* u = N_VGetArrayPointer(y);
  sunrealtype* f       = N_VGetArrayPointer(ydot);
  for (sunindextype j = 0; j < n; j++)
  {
    for (sunindextype i = 0; i < n; i++)
    {
      const sunindextype k = j * n + i;
      sunrealtype s        = -SUN_RCONST(4.0) * u[k];
      if (i > 0) { s += u[k - 1]; }
      if (i < n - 1) { s += u[k + 1]; }
      if (j > 0) { s += u[k - n]; }
      if (j < n - 1) { s += u[k + n]; }
      f[k] = h->dx2inv * s;
    }
  }
  return 0;
}

int jac_heat(sunrealtype, N_Vector, N_Vector, SUNMatrix J, void* data, N_Vector,
             N_Vector, N_Vector)
{
  const Heat* h        = static_cast<Heat*>(data);
  const sunindextype n = h->nx;
  for (sunindextype j = 0; j < n; j++)
  {
    for (sunindextype i = 0; i < n; i++)
    {
      const sunindextype k  = j * n + i;
      SM_ELEMENT_B(J, k, k) = -SUN_RCONST(4.0) * h->dx2inv;
      if (i > 0) { SM_ELEMENT_B(J, k, k - 1) = h->dx2inv; }
      if (i < n - 1) { SM_ELEMENT_B(J, k, k + 1) = h->dx2inv; }
      if (j > 0) { SM_ELEMENT_B(J, k, k - n) = h->dx2inv; }
      if (j < n - 1) { SM_ELEMENT_B(J, k, k + n) = h->dx2inv; }
    }
  }
  return 0;
}

int jtimes_heat(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector, N_Vector,
                void* data, N_Vector)
{
  return f_heat(t, v, Jv, data);
}

/* Jacobi preconditioner of I - gamma J: the diagonal 1 + 4 gamma / dx^2 */
int psetup_heat(sunrealtype, N_Vector, N_Vector, sunbooleantype,
                sunbooleantype* jcur, sunrealtype gamma, void* data)
{
  static_cast<Heat*>(data)->gamma = gamma;
  *jcur                           = SUNTRUE;
  return 0;
}

int psolve_heat(sunrealtype, N_Vector, N_Vector, N_Vector r, N_Vector z,
                sunrealtype gamma, sunrealtype, int, void* data)
{
  const Heat* h = static_cast<Heat*>(data);
  N_VScale(one / (one + SUN_RCONST(4.0) * gamma * h->dx2inv), r, z);
  return 0;
}

/* ---- one run ---- */

struct Result
{
  int flag     = 0;
  long int nst = 0, nni = 0, nsli = 0, nscf = 0, nje = 0;
  double seconds = 0;
};

enum class Solver
{
  dense,
  band,
  spgmr_jacobi
};

Result run(SUNContext sunctx, const char* name, FIRKRhsFn f, FIRKLsJacFn jac,
           N_Vector y0, sunrealtype tend, sunrealtype rtol, sunrealtype atol,
           Solver solver, Heat* heat)
{
  Result r;
  const sunindextype N = N_VGetLength(y0);
  N_Vector y           = N_VClone(y0);
  N_VScale(one, y0, y);
  SUNMatrix A        = nullptr;
  SUNLinearSolver LS = nullptr;
  void* mem          = FIRKodeCreate(sunctx);
  int flag           = FIRKodeInit(mem, f, zero, y);
  flag += FIRKodeSStolerances(mem, rtol, atol);
  flag += FIRKodeSetUserData(mem, heat);
  flag += FIRKodeSetMaxNumSteps(mem, 1000000);
  flag += FIRKodeSetStopTime(mem, tend);
  switch (solver)
  {
  case Solver::dense:
    A  = SUNDenseMatrix(N, N, sunctx);
    LS = SUNLinSol_Dense(y, A, sunctx);
    flag += FIRKodeSetLinearSolver(mem, LS, A);
    flag += FIRKodeSetJacFn(mem, jac);
    break;
  case Solver::band:
    A  = SUNBandMatrix(N, heat->nx, heat->nx, sunctx);
    LS = SUNLinSol_Band(y, A, sunctx);
    flag += FIRKodeSetLinearSolver(mem, LS, A);
    flag += FIRKodeSetJacFn(mem, jac);
    break;
  case Solver::spgmr_jacobi:
    LS = SUNLinSol_SPGMR(y, SUN_PREC_LEFT, 20, sunctx);
    flag += FIRKodeSetLinearSolver(mem, LS, nullptr);
    flag += FIRKodeSetJacTimes(mem, nullptr, jtimes_heat);
    flag += FIRKodeSetPreconditioner(mem, psetup_heat, psolve_heat);
    break;
  }
  if (flag) { r.flag = -1000; }
  else
  {
    sunrealtype t;
    const auto start = std::chrono::steady_clock::now();
    r.flag           = FIRKodeEvolve(mem, tend, y, &t, FIRK_NORMAL);
    r.seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  }
  FIRKodeGetNumSteps(mem, &r.nst);
  FIRKodeGetNumNonlinSolvIters(mem, &r.nni);
  FIRKodeGetNumStageLinIters(mem, &r.nsli);
  FIRKodeGetNumStageLinConvFails(mem, &r.nscf);
  FIRKodeGetNumJacEvals(mem, &r.nje);

  const char* solver_name = (solver == Solver::dense)  ? "dense"
                            : (solver == Solver::band) ? "band"
                                                       : "spgmr+jacobi";
  std::printf("%s,%ld,%s,%.0e,%d,%ld,%ld,%ld,%.2f,%ld,%ld,%.3f\n", name,
              (long)N, solver_name, (double)rtol, r.flag, r.nst, r.nni, r.nsli,
              (r.nni > 0) ? (double)r.nsli / (double)r.nni : 0.0, r.nscf, r.nje,
              r.seconds);

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  return r;
}

} // namespace

int main(int argc, char* argv[])
{
  const sunindextype nx = (argc > 1) ? std::atol(argv[1]) : 100;
  SUNContext sunctx;
  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  std::printf(
    "problem,N,solver,rtol,flag,nst,nni,nsli,klin,nscf,nje,seconds\n");
  const sunrealtype rtols[2] = {SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-8)};
  int fails                  = 0;

  /* the small problems */
  struct Small
  {
    const char* name;
    int n;
    FIRKRhsFn f;
    FIRKLsJacFn jac;
    sunrealtype y0[3];
    sunrealtype tend;
  };

  const Small small[3] =
    {{"rober", 3, f_rober, jac_rober, {one, zero, zero}, SUN_RCONST(1.0e11)},
     {"vdpol", 2, f_vdpol, jac_vdpol, {SUN_RCONST(2.0), zero, zero}, SUN_RCONST(2.0)},
     {"orego",
      3,
      f_orego,
      jac_orego,
      {one, SUN_RCONST(2.0), SUN_RCONST(3.0)},
      SUN_RCONST(360.0)}};
  for (const Small& p : small)
  {
    N_Vector y0 = N_VNew_Serial(p.n, sunctx);
    for (int i = 0; i < p.n; i++) { NV_Ith_S(y0, i) = p.y0[i]; }
    for (sunrealtype rtol : rtols)
    {
      const Result r = run(sunctx, p.name, p.f, p.jac, y0, p.tend, rtol,
                           SUN_RCONST(1.0e-4) * rtol, Solver::dense, nullptr);
      if (r.flag != FIRK_TSTOP_RETURN) { fails++; }
    }
    N_VDestroy(y0);
  }

  /* the heat equation: u(0) = sin(pi x) sin(pi y) plus a rough component */
  Heat heat;
  heat.nx     = nx;
  heat.dx2inv = (sunrealtype)((nx + 1) * (nx + 1));
  heat.gamma  = zero;
  N_Vector y0 = N_VNew_Serial(nx * nx, sunctx);
  {
    sunrealtype* u       = N_VGetArrayPointer(y0);
    const sunrealtype pi = SUN_RCONST(3.14159265358979323846);
    for (sunindextype j = 0; j < nx; j++)
    {
      for (sunindextype i = 0; i < nx; i++)
      {
        const sunrealtype x  = (sunrealtype)(i + 1) / (sunrealtype)(nx + 1);
        const sunrealtype yy = (sunrealtype)(j + 1) / (sunrealtype)(nx + 1);
        u[j * nx + i]        = std::sin(pi * x) * std::sin(pi * yy) +
                        SUN_RCONST(0.1) * std::sin(SUN_RCONST(7.0) * pi * x) *
                          std::sin(SUN_RCONST(5.0) * pi * yy);
      }
    }
  }
  for (Solver solver : {Solver::band, Solver::spgmr_jacobi})
  {
    for (sunrealtype rtol : rtols)
    {
      const Result r = run(sunctx, "heat2d", f_heat, jac_heat, y0,
                           SUN_RCONST(0.1), rtol, SUN_RCONST(1.0e-10), solver,
                           &heat);
      if (r.flag != FIRK_TSTOP_RETURN) { fails++; }
    }
  }
  N_VDestroy(y0);

  SUNContext_Free(&sunctx);
  return fails ? 1 : 0;
}
