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
 * mass-times routine, and (c) the solvers of (b) with both sides of the
 * system multiplied by 1e-4, which leaves the solution unchanged but would
 * make the tolerances of the linear solves 1e4 times too loose if they were
 * applied to the residual, in the units of M y, instead of the units of y
 * (SPGMR has no preconditioner to convert them). All are compared against
 * the equivalent system y' = M^{-1} f(t, y) integrated without a mass matrix.
 * With fixed steps the solutions must agree to roundoff and solver
 * tolerances; an adaptive run additionally exercises the mass solves in the
 * initial step estimate and the error estimate. A last adaptive run with
 * s = 1 and dense Newton and mass matrices, where the mass linear solver is
 * only used for y' = M^{-1} f in the initial step estimate, limits PCG to
 * one iteration: its approximation of y', with the residual reduced but
 * above the tolerance, must be accepted since y' only feeds the heuristic,
 * and the integration must match that with a dense mass solver. Finally, a
 * time-dependent mass matrix, and a matrix-free mass solver with a dense
 * Newton solver, must be rejected (design hazard H28).
 *
 * Pass criterion (not run in single precision): each mass matrix run differs
 * from the run without a mass matrix by at most 1e-9 (dense) or 1e-8 (Krylov)
 * with h = 0.01, and 1e-6 adaptively at rtol = atol = 1e-8, in the max norm;
 * every mass matrix run counts mass products; the dense M is set up once; the
 * dense fixed-step run counts s (nni + nsli) products and s nni solves (nni
 * and no solves for s = 1); the adaptive dense and unscaled Krylov runs
 * perform mass solves; in the last run PCG misses its tolerance at least once
 * and the step count differs from that with a dense mass solver by at most
 * 20%, with an error at most 1.25 times its error plus rtol; time_dep makes
 * FIRKodeSetMassLinearSolver return FIRKLS_ILL_INPUT, and the mixed solvers
 * make FIRKodeEvolve return FIRK_LINIT_FAIL.
 *
 * Usage: firk_test_mass [stages]
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "firkode/firkode.h"
#include "firkode/firkode_impl.h" /* firkSetNumStagesAny */
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

// both sides of M y' = f multiplied by mfac in mode 3
constexpr sunrealtype mfac = SUN_RCONST(1.0e-4);

int rhs_scaled(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  const int flag = problems::kpr::ode_rhs(t, y, ydot, user_data);
  if (flag) { return flag; }
  N_VScale(mfac, ydot, ydot);
  return 0;
}

int mass_times_scaled(N_Vector v, N_Vector Mv, sunrealtype t, void* mtimes_data)
{
  const int flag = mass_times(v, Mv, t, mtimes_data);
  if (flag) { return flag; }
  N_VScale(mfac, Mv, Mv);
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
  long int nst, nni, nsli, nmsetups, nmsolves, nmvevals;
};

// integrate KPR to t = 0.5; mode 0 = reference without mass matrix,
// 1 = dense mass solver, 2 = matrix-free PCG mass solver, 3 = mode 2 with
// both sides multiplied by mfac
int run(SUNContext sunctx, int s, sunrealtype h, int mode, sunrealtype* yout,
        Stats* st)
{
  N_Vector y      = N_VNew_Serial(2, sunctx);
  sunrealtype* yd = N_VGetArrayPointer(y);
  problems::kpr::true_sol(zero, &yd[0], &yd[1]);

  void* mem = FIRKodeCreate(sunctx);
  int flag  = FIRKodeInit(mem,
                         (mode == 0)   ? rhs_minv
                          : (mode == 3) ? rhs_scaled
                                        : problems::kpr::ode_rhs,
                          zero, y);
  flag += FIRKodeSetUserData(mem, static_cast<void*>(kpr_data));
  flag += FIRKodeSetMode(mem, FIRK_EXPERT);
  flag += firkSetNumStagesAny(mem, s);
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
  else if (mode > 1)
  {
    MLS = SUNLinSol_PCG(y, SUN_PREC_NONE, 10, sunctx);
    flag += FIRKodeSetMassLinearSolver(mem, MLS, nullptr, SUNFALSE);
    flag += FIRKodeSetMassTimes(mem, nullptr,
                                (mode == 3) ? mass_times_scaled : mass_times,
                                nullptr);
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
  FIRKodeGetNumStageLinIters(mem, &st->nsli);
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
  const char* names[4] = {"no mass", "dense mass", "matrix-free", "M, f x 1e-4"};
  return names[mode];
}

int compare(SUNContext sunctx, int s, sunrealtype h, const sunrealtype* tol,
            const char* label)
{
  sunrealtype y[4][2];
  Stats st[4];
  int fails = 0;

  std::printf("%s\n", label);
  for (int mode = 0; mode < 4; mode++)
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
  for (int mode = 1; mode < 4; mode++)
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
  // fixed steps with a direct mass solver: s products per Newton iteration in
  // the residual and per FGMRES iteration in the stacked operator, and s mass
  // solves per stacked Newton solve for its units conversion; with s = 1 the
  // direct block solver needs no conversion
  if (h > zero)
  {
    const long int nmult_expected  = (s > 1) ? s * (st[1].nni + st[1].nsli)
                                             : st[1].nni;
    const long int nsolve_expected = (s > 1) ? s * st[1].nni : 0;
    if (st[1].nmvevals != nmult_expected || st[1].nmsolves != nsolve_expected)
    {
      std::printf("  FAIL: dense mass run has %li products and %li solves, "
                  "expected %li and %li\n",
                  st[1].nmvevals, st[1].nmsolves, nmult_expected,
                  nsolve_expected);
      fails++;
    }
  }
  if (h <= zero && (st[1].nmsolves == 0 || st[2].nmsolves == 0))
  {
    std::printf("  FAIL: the adaptive run should perform mass solves\n");
    fails++;
  }
  return fails;
}

// integrate M y' = f to t = 0.5 with s = 1 and dense Newton and mass
// matrices, solving with M by the dense solver or by PCG limited to one
// iteration; the mass linear solver is then only used for y' = M^{-1} f in
// the initial step size estimate
int run_one_stage(SUNContext sunctx, bool pcg, sunrealtype rtol,
                  sunrealtype* yout, Stats* st, long int* nmcfails)
{
  N_Vector y      = N_VNew_Serial(2, sunctx);
  sunrealtype* yd = N_VGetArrayPointer(y);
  problems::kpr::true_sol(zero, &yd[0], &yd[1]);

  void* mem = FIRKodeCreate(sunctx);
  int flag  = FIRKodeInit(mem, problems::kpr::ode_rhs, zero, y);
  flag += FIRKodeSetUserData(mem, static_cast<void*>(kpr_data));
  flag += FIRKodeSetMode(mem, FIRK_EXPERT);
  flag += firkSetNumStagesAny(mem, 1);
  flag += FIRKodeSStolerances(mem, rtol, rtol);
  flag += FIRKodeSetMaxNumSteps(mem, 100000);

  SUNMatrix A         = SUNDenseMatrix(2, 2, sunctx);
  SUNLinearSolver LS  = SUNLinSol_Dense(y, A, sunctx);
  SUNMatrix M         = SUNDenseMatrix(2, 2, sunctx);
  SUNLinearSolver MLS = pcg ? SUNLinSol_PCG(y, SUN_PREC_NONE, 1, sunctx)
                            : SUNLinSol_Dense(y, M, sunctx);
  flag += FIRKodeSetLinearSolver(mem, LS, A);
  flag += FIRKodeSetJacFn(mem, problems::kpr::ode_rhs_jac);
  flag += FIRKodeSetMassLinearSolver(mem, MLS, M, SUNFALSE);
  flag += FIRKodeSetMassFn(mem, mass_fn);

  sunrealtype t;
  if (flag == 0)
  {
    flag = FIRKodeEvolve(mem, SUN_RCONST(0.5), y, &t, FIRK_NORMAL);
  }
  if (flag == FIRK_SUCCESS)
  {
    yout[0] = yd[0];
    yout[1] = yd[1];
    FIRKodeGetNumSteps(mem, &st->nst);
    FIRKodeGetNumMassSolves(mem, &st->nmsolves);
    FIRKodeGetNumMassConvFails(mem, nmcfails);
  }

  FIRKodeFree(&mem);
  SUNLinSolFree(MLS);
  SUNMatDestroy(M);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  return flag;
}

// mass solves of the initial step size estimate that only reduce the residual
int inexact_initial_mass_solves(SUNContext sunctx)
{
  const sunrealtype rtol = SUN_RCONST(1.0e-5);
  const char* names[2]   = {"dense mass", "PCG, 1 iter"};
  sunrealtype yref[2], y[2][2], err[2];
  Stats stref          = {};
  Stats st[2]          = {};
  long int nmcfails[2] = {0, 0};
  int fails            = 0;

  std::printf("Adaptive, s = 1, rtol = atol = 1e-5, mass solves only in the "
              "initial step estimate:\n");

  // reference: the system without mass matrix, s = 3, rtol = atol = 1e-8
  if (run(sunctx, 3, zero, 0, yref, &stref))
  {
    std::printf("  FAIL: reference run failed\n");
    return 1;
  }
  for (int k = 0; k < 2; k++)
  {
    if (run_one_stage(sunctx, k == 1, rtol, y[k], &st[k], &nmcfails[k]))
    {
      std::printf("  FAIL: %s run failed\n", names[k]);
      return 1;
    }
    err[k] = std::max(std::abs(y[k][0] - yref[0]), std::abs(y[k][1] - yref[1]));
    std::printf("  %-11s y = (%.10" ESYM ", %.10" ESYM "), nst = %li, mass "
                "solves/conv fails = %li/%li, error = %.1" ESYM "\n",
                names[k], y[k][0], y[k][1], st[k].nst, st[k].nmsolves,
                nmcfails[k], err[k]);
  }
  if (nmcfails[1] == 0)
  {
    std::printf("  FAIL: PCG reached the tolerance of the mass solves\n");
    fails++;
  }
  if (5 * std::labs(st[1].nst - st[0].nst) > st[0].nst ||
      err[1] > SUN_RCONST(1.25) * err[0] + rtol)
  {
    std::printf("  FAIL: the PCG run differs from the dense mass run\n");
    fails++;
  }
  return fails;
}

// a time-dependent mass matrix is rejected, and the Newton and mass linear
// solvers must both be matrix-based or both matrix-free: here the Newton
// solver is dense and the mass solver matrix-free, so Evolve fails in the
// initialization of the linear solvers
int mass_input_checks(SUNContext sunctx)
{
  int fails       = 0;
  N_Vector y      = N_VNew_Serial(2, sunctx);
  sunrealtype* yd = N_VGetArrayPointer(y);
  problems::kpr::true_sol(zero, &yd[0], &yd[1]);

  SUNMatrix A          = SUNDenseMatrix(2, 2, sunctx);
  SUNMatrix M          = SUNDenseMatrix(2, 2, sunctx);
  SUNLinearSolver LS   = SUNLinSol_Dense(y, A, sunctx);
  SUNLinearSolver MLS  = SUNLinSol_Dense(y, M, sunctx);
  SUNLinearSolver MPCG = SUNLinSol_PCG(y, SUN_PREC_NONE, 10, sunctx);
  void* mem            = FIRKodeCreate(sunctx);
  int flag             = FIRKodeInit(mem, problems::kpr::ode_rhs, zero, y);
  flag += FIRKodeSetUserData(mem, static_cast<void*>(kpr_data));
  flag += FIRKodeSStolerances(mem, SUN_RCONST(1.0e-8), SUN_RCONST(1.0e-8));
  flag += FIRKodeSetLinearSolver(mem, LS, A);
  flag += FIRKodeSetJacFn(mem, problems::kpr::ode_rhs_jac);
  if (flag) { return 1; }

  std::printf("Mass matrix input checks:\n");
  flag = FIRKodeSetMassLinearSolver(mem, MLS, M, SUNTRUE);
  std::printf("  time-dependent M: %d\n", flag);
  if (flag != FIRKLS_ILL_INPUT)
  {
    std::printf("  FAIL: a time-dependent mass matrix is not rejected\n");
    fails++;
  }

  flag = FIRKodeSetMassLinearSolver(mem, MPCG, nullptr, SUNFALSE);
  flag += FIRKodeSetMassTimes(mem, nullptr, mass_times, nullptr);
  if (flag) { return 1; }
  sunrealtype t;
  flag = FIRKodeEvolve(mem, SUN_RCONST(0.5), y, &t, FIRK_ONE_STEP);
  std::printf("  dense Newton and matrix-free mass solvers: %d\n", flag);
  if (flag != FIRK_LINIT_FAIL)
  {
    std::printf("  FAIL: a matrix-free mass solver with a matrix-based "
                "Newton solver is not rejected\n");
    fails++;
  }

  FIRKodeFree(&mem);
  SUNLinSolFree(MPCG);
  SUNLinSolFree(MLS);
  SUNLinSolFree(LS);
  SUNMatDestroy(M);
  SUNMatDestroy(A);
  N_VDestroy(y);
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
  const sunrealtype tol_fixed[4] = {zero, SUN_RCONST(1.0e-9),
                                    SUN_RCONST(1.0e-8), SUN_RCONST(1.0e-8)};
  const sunrealtype tol_adapt[4] = {zero, SUN_RCONST(1.0e-6),
                                    SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-6)};

  std::printf("Mass matrix test, KPR with Radau IIA s = %d\n", s);
  fails += compare(sunctx, s, SUN_RCONST(0.01), tol_fixed,
                   "Fixed step h = 0.01:");
  fails += compare(sunctx, s, zero, tol_adapt, "Adaptive, rtol = atol = 1e-8:");
  fails += inexact_initial_mass_solves(sunctx);
  fails += mass_input_checks(sunctx);

  SUNContext_Free(&sunctx);

  if (fails)
  {
    std::printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  std::printf("\nAll tests passed\n");
  return 0;
}
