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
 * Mass matrix test for FIRKODE: the Kvaerno-Prothero-Robinson right-hand side
 * f is integrated as the system M y' = f(t, y) with the constant mass matrix
 * M = [[2, 1], [1, 3]], using (a) dense linear solvers for the Newton
 * matrices and the mass matrix with user-supplied Jacobian and mass matrix
 * routines and (b) matrix-free SPGMR (Newton) and PCG (mass) solvers with the
 * internal difference-quotient Jacobian-vector product and a user-supplied
 * mass-times routine. Both are compared against the equivalent system
 * y' = M^{-1} f(t, y) integrated without a mass matrix.
 * With fixed steps the solutions must agree to roundoff and solver
 * tolerances; an adaptive run additionally exercises the mass solves in the
 * initial step estimate and the error estimate.
 *
 * Usage: firk_test_mass [stages]
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "problems/kpr.hpp"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunlinsol/sunlinsol_pcg.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunmatrix/sunmatrix_dense.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define ESYM "Le"
#else
#define ESYM "e"
#endif

namespace {

constexpr sunrealtype zero = SUN_RCONST(0.0);

// mass matrix M and its inverse
constexpr sunrealtype m00 = SUN_RCONST(2.0);
constexpr sunrealtype m01 = SUN_RCONST(1.0);
constexpr sunrealtype m11 = SUN_RCONST(3.0);
constexpr sunrealtype det = m00 * m11 - m01 * m01;

sunrealtype kpr_data[4] = {SUN_RCONST(-2.0), SUN_RCONST(0.5), SUN_RCONST(0.5),
                           SUN_RCONST(-1.0)};

// mass matrix routine for the dense SUNMatrix
int mass_fn(sunrealtype, SUNMatrix M, void*, N_Vector, N_Vector, N_Vector)
{
  sunrealtype* Md = SUNDenseMatrix_Data(M);
  Md[0]           = m00; // column-major: (0,0)
  Md[1]           = m01; // (1,0)
  Md[2]           = m01; // (0,1)
  Md[3]           = m11; // (1,1)
  return 0;
}

// matrix-free mass-vector product Mv = M v
int mass_times(N_Vector v, N_Vector Mv, sunrealtype, void*)
{
  const sunrealtype* vd = N_VGetArrayPointer(v);
  sunrealtype* Mvd      = N_VGetArrayPointer(Mv);
  Mvd[0]                = m00 * vd[0] + m01 * vd[1];
  Mvd[1]                = m01 * vd[0] + m11 * vd[1];
  return 0;
}

// reference system: y' = M^{-1} f(t, y)
int rhs_minv(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  const int flag = problems::kpr::ode_rhs(t, y, ydot, user_data);
  if (flag) { return flag; }
  sunrealtype* fd      = N_VGetArrayPointer(ydot);
  const sunrealtype f0 = fd[0];
  const sunrealtype f1 = fd[1];
  fd[0]                = (m11 * f0 - m01 * f1) / det;
  fd[1]                = (-m01 * f0 + m00 * f1) / det;
  return 0;
}

int jac_minv(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
             void* user_data, N_Vector t1, N_Vector t2, N_Vector t3)
{
  const int flag = problems::kpr::ode_rhs_jac(t, y, fy, J, user_data, t1, t2, t3);
  if (flag) { return flag; }
  sunrealtype* Jd = SUNDenseMatrix_Data(J);
  for (int col = 0; col < 2; col++)
  {
    const sunrealtype j0 = Jd[2 * col];
    const sunrealtype j1 = Jd[2 * col + 1];
    Jd[2 * col]          = (m11 * j0 - m01 * j1) / det;
    Jd[2 * col + 1]      = (-m01 * j0 + m00 * j1) / det;
  }
  return 0;
}

struct Stats
{
  long int nst, nni, nmsetups, nmsolves, nmvevals;
};

// integrate KPR to t = 0.5; mode 0 = reference without mass matrix,
// 1 = dense mass solver, 2 = matrix-free PCG mass solver
int run(SUNContext sunctx, int s, sunrealtype h, int mode, sunrealtype* yout,
        Stats* st)
{
  N_Vector y      = N_VNew_Serial(2, sunctx);
  sunrealtype* yd = N_VGetArrayPointer(y);
  problems::kpr::true_sol(zero, &yd[0], &yd[1]);

  void* mem = FIRKodeCreate(sunctx);
  int flag  = FIRKodeInit(mem, (mode == 0) ? rhs_minv : problems::kpr::ode_rhs,
                          zero, y);
  flag += FIRKodeSetUserData(mem, static_cast<void*>(kpr_data));
  flag += FIRKodeSetNumStages(mem, s);
  if (h > zero)
  {
    flag += FIRKodeSStolerances(mem, SUN_RCONST(1.0e-10), SUN_RCONST(1.0e-10));
    flag += FIRKodeSetFixedStep(mem, h);
  }
  else
  {
    flag += FIRKodeSStolerances(mem, SUN_RCONST(1.0e-8), SUN_RCONST(1.0e-8));
  }
  flag += FIRKodeSetNonlinConvCoef(mem, SUN_RCONST(1.0e-3));
  flag += FIRKodeSetMaxNonlinIters(mem, 20);
  flag += FIRKodeSetMaxNumSteps(mem, 100000);

  SUNMatrix A         = nullptr;
  SUNLinearSolver LS  = nullptr;
  SUNMatrix M         = nullptr;
  SUNLinearSolver MLS = nullptr;
  if (mode < 2)
  {
    A  = SUNDenseMatrix(2, 2, sunctx);
    LS = SUNLinSol_Dense(y, A, sunctx);
    flag += FIRKodeSetLinearSolver(mem, LS, A);
    flag += FIRKodeSetJacFn(mem, (mode == 0) ? jac_minv
                                             : problems::kpr::ode_rhs_jac);
  }
  else
  {
    LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 10, sunctx);
    flag += FIRKodeSetLinearSolver(mem, LS, nullptr);
    flag += FIRKodeSetEpsLin(mem, SUN_RCONST(1.0e-3));
  }
  flag += FIRKodeSetStageEpsLin(mem, SUN_RCONST(1.0e-3));

  if (mode == 1)
  {
    M   = SUNDenseMatrix(2, 2, sunctx);
    MLS = SUNLinSol_Dense(y, M, sunctx);
    flag += FIRKodeSetMassLinearSolver(mem, MLS, M, SUNFALSE);
    flag += FIRKodeSetMassFn(mem, mass_fn);
  }
  else if (mode == 2)
  {
    MLS = SUNLinSol_PCG(y, SUN_PREC_NONE, 10, sunctx);
    flag += FIRKodeSetMassLinearSolver(mem, MLS, nullptr, SUNFALSE);
    flag += FIRKodeSetMassTimes(mem, nullptr, mass_times, nullptr);
    flag += FIRKodeSetMassEpsLin(mem, SUN_RCONST(1.0e-4));
  }
  if (flag) { return flag; }

  sunrealtype t;
  flag = FIRKodeEvolve(mem, SUN_RCONST(0.5), y, &t, FIRK_NORMAL);
  if (flag != FIRK_SUCCESS) { return flag; }

  yout[0] = yd[0];
  yout[1] = yd[1];
  FIRKodeGetNumSteps(mem, &st->nst);
  FIRKodeGetNumNonlinSolvIters(mem, &st->nni);
  st->nmsetups = st->nmsolves = st->nmvevals = 0;
  if (mode > 0)
  {
    FIRKodeGetNumMassSetups(mem, &st->nmsetups);
    FIRKodeGetNumMassSolves(mem, &st->nmsolves);
    FIRKodeGetNumMassMult(mem, &st->nmvevals);
  }

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  if (A) { SUNMatDestroy(A); }
  if (MLS) { SUNLinSolFree(MLS); }
  if (M) { SUNMatDestroy(M); }
  N_VDestroy(y);
  return 0;
}

const char* mode_name(int mode)
{
  return (mode == 0) ? "no mass" : (mode == 1) ? "dense mass" : "matrix-free";
}

int compare(SUNContext sunctx, int s, sunrealtype h, const sunrealtype* tol,
            const char* label)
{
  sunrealtype y[3][2];
  Stats st[3];
  int fails = 0;

  std::printf("%s\n", label);
  for (int mode = 0; mode < 3; mode++)
  {
    if (run(sunctx, s, h, mode, y[mode], &st[mode]))
    {
      std::printf("  FAIL: %s run failed\n", mode_name(mode));
      return 1;
    }
    std::printf("  %-11s y = (%.10" ESYM ", %.10" ESYM "), nst = %li, nni = "
                "%li, mass setups/solves/mults = %li/%li/%li\n",
                mode_name(mode), y[mode][0], y[mode][1], st[mode].nst,
                st[mode].nni, st[mode].nmsetups, st[mode].nmsolves,
                st[mode].nmvevals);
  }
  for (int mode = 1; mode < 3; mode++)
  {
    const sunrealtype err = std::max(std::abs(y[mode][0] - y[0][0]),
                                     std::abs(y[mode][1] - y[0][1]));
    std::printf("  %s vs no mass: max difference = %.3" ESYM "\n",
                mode_name(mode), err);
    if (err > tol[mode])
    {
      std::printf("  FAIL: %s solution differs from the reference\n",
                  mode_name(mode));
      fails++;
    }
    if (st[mode].nmvevals == 0)
    {
      std::printf("  FAIL: no mass matrix-vector products were counted\n");
      fails++;
    }
  }
  if (st[1].nmsetups != 1)
  {
    std::printf("  FAIL: the constant mass matrix should be set up once\n");
    fails++;
  }
  if (h <= zero && (st[1].nmsolves == 0 || st[2].nmsolves == 0))
  {
    std::printf("  FAIL: the adaptive run should perform mass solves\n");
    fails++;
  }
  return fails;
}

} // namespace

int main(int argc, char* argv[])
{
  const int s = (argc > 1) ? std::atoi(argv[1]) : 3;
  int fails   = 0;

  SUNContext sunctx;
  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  /* tolerances per mode: the matrix-free runs use inexact Krylov solves */
  const sunrealtype tol_fixed[3] = {zero, SUN_RCONST(1.0e-9), SUN_RCONST(1.0e-8)};
  const sunrealtype tol_adapt[3] = {zero, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-6)};

  std::printf("Mass matrix test, KPR with Radau IIA s = %d\n", s);
  fails += compare(sunctx, s, SUN_RCONST(0.01), tol_fixed,
                   "Fixed step h = 0.01:");
  fails += compare(sunctx, s, zero, tol_adapt, "Adaptive, rtol = atol = 1e-8:");

  SUNContext_Free(&sunctx);

  if (fails)
  {
    std::printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  std::printf("\nAll tests passed\n");
  return 0;
}
