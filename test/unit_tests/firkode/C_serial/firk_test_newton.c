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
 * Unit test for the Newton driver (design Listing 8.3, decision D3).
 *
 * Part 1 calls firkNewton directly. For y' = lambda y, the one-stage method and
 * a Jacobian mu != lambda, the simplified Newton iteration on the linear stage
 * equation reduces the error by exactly rho = h (lambda - mu) / (1 - h mu) per
 * iteration, so that from the trivial predictor it ends with
 *   - FIRK_NLS_SUCCESS after 3 iterations and the rate rho, for rho = 0.2;
 *   - FIRK_NLS_DIVERGED at k = 2 and the factor 0.5, for rho = 1.5;
 *   - FIRK_NLS_PREDICTED_FAIL at k = 2 and the factor 0.8 q^(-1/(4+m-1-k)),
 *     q clamped to 20, for rho = 0.9;
 *   - FIRK_NLS_MAXITER and the factor 0.5 with m = 2, the rate left at
 *     theta_jac since no rate is estimated at k = m;
 *   - FIRK_NLS_RHS_RECVR and the factor etacf for a recoverable RHS failure,
 *     and FIRK_RHSFUNC_FAIL for an unrecoverable one, before any solve;
 *   - FIRK_NLS_DIVERGED and the factor 0.5 at k = 1 for a non-finite update.
 * Every linear solve is one iteration in nni.
 *
 * Part 2 (hazard H19) uses the public API with the three-stage method: a
 * recoverable RHS failure in a stage reduces h by 0.25 and the step succeeds;
 * maxncf of them return FIRK_REPTD_RHSFUNC_ERR; an unrecoverable failure returns
 * FIRK_RHSFUNC_FAIL.
 *
 * Pass criterion: every status, step factor, rate and count matches.
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>

#include "firkode/firkode.h"
#include "firkode/firkode_impl.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ONE SUN_RCONST(1.0)

/* the problem, its Jacobian, and injected RHS failures */
static sunrealtype lambda = SUN_RCONST(-1.0);
static sunrealtype mu     = SUN_RCONST(-1.0);
static int rhs_mode = 0; /* 0 normal, 1 recoverable, 2 unrecoverable, 3 NaN */
static sunrealtype fail_after_t = SUN_RCONST(0.0); /* fail only at t > this */
static int fail_count           = 0; /* number of failures left to inject */

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)user_data;
  NV_Ith_S(ydot, 0) = lambda * NV_Ith_S(y, 0);
  if (rhs_mode == 3) { NV_Ith_S(ydot, 0) = (sunrealtype)NAN; }
  if (rhs_mode != 0 && rhs_mode != 3 && t > fail_after_t && fail_count > 0)
  {
    fail_count--;
    return (rhs_mode == 1) ? 1 : -1;
  }
  return 0;
}

static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  (void)t;
  (void)y;
  (void)fy;
  (void)user_data;
  (void)tmp1;
  (void)tmp2;
  (void)tmp3;
  SM_ELEMENT_D(J, 0, 0) = mu;
  return 0;
}

static int check(int ok, const char* what)
{
  if (!ok)
  {
    printf("  FAIL: %s\n", what);
    return 1;
  }
  return 0;
}

static int close_to(sunrealtype a, sunrealtype b)
{
  return SUNRabs(a - b) <= SUN_RCONST(1000.0) * SUN_UNIT_ROUNDOFF * SUNRabs(b);
}

/* the Jacobian with rate rho for h = 0.1 and lambda = -1 */
static sunrealtype mu_for_rate(sunrealtype h, sunrealtype rho)
{
  /* rho (1 - h mu) = h (lambda - mu)  =>  mu = (h lambda - rho) / (h - h rho) */
  return (h * lambda - rho) / (h * (ONE - rho));
}

/* one direct call of firkNewton from the trivial predictor with a fresh
   setup at the Jacobian mu_case */
static int newton(void* mem, sunrealtype h, sunrealtype mu_case, int maxcor,
                  int mode, long int* nsolves)
{
  FIRKodeMem fm = (FIRKodeMem)mem;
  sunbooleantype jcur;
  long int nni0;
  int status;

  mu          = mu_case;
  fm->maxcor  = maxcor;
  fm->h       = h;
  fm->gamma   = h * fm->tbl->gamma0;
  fm->eta_nls = ONE;
  if (fm->lsetup(fm, SUNTRUE, &jcur) != 0) { return -999; }
  N_VConst(SUN_RCONST(0.0), fm->Z[0]);

  rhs_mode     = mode;
  fail_after_t = -SUN_BIG_REAL;
  fail_count   = 1;
  nni0         = fm->nni;
  status       = firkNewton(fm);
  rhs_mode     = 0;
  fail_count   = 0;
  *nsolves     = fm->nni - nni0;
  return status;
}

static int part1(SUNContext sunctx)
{
  const sunrealtype h = SUN_RCONST(0.1);
  int fails           = 0, status;
  long int nsolves;
  sunrealtype tret, theta, hfac;
  void* mem          = NULL;
  N_Vector y         = N_VNew_Serial(1, sunctx);
  SUNMatrix A        = SUNDenseMatrix(1, 1, sunctx);
  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
  FIRKodeMem fm;

  NV_Ith_S(y, 0) = ONE;
  mem            = FIRKodeCreate(sunctx);
  if (FIRKodeInit(mem, f, SUN_RCONST(0.0), y) ||
      FIRKodeSStolerances(mem, SUN_RCONST(1.0e-2), SUN_RCONST(1.0e-2)) ||
      FIRKodeSetLinearSolver(mem, LS, A) || FIRKodeSetJacFn(mem, jac) ||
      FIRKodeSetMode(mem, FIRK_EXPERT) || firkSetNumStagesAny(mem, 1) ||
      FIRKodeSetFixedStep(mem, h) ||
      FIRKodeSetNonlinConvCoef(mem, SUN_RCONST(0.1)))
  {
    return 1;
  }
  fm = (FIRKodeMem)mem;

  /* one step to set up the integrator */
  mu = mu_for_rate(h, SUN_RCONST(0.2));
  fails += check(FIRKodeEvolve(mem, ONE, y, &tret, FIRK_ONE_STEP) == FIRK_SUCCESS,
                 "the first step succeeds");

  /* convergence: the rate is rho */
  status = newton(mem, h, mu_for_rate(h, SUN_RCONST(0.2)), 7, 0, &nsolves);
  FIRKodeGetLastNonlinConvRate(mem, &theta);
  fails += check(status == FIRK_NLS_SUCCESS && nsolves == 3 &&
                   fm->nni_last == 3 && close_to(theta, SUN_RCONST(0.2)),
                 "rho = 0.2: success after 3 iterations with theta = 0.2");

  /* divergence */
  status = newton(mem, h, mu_for_rate(h, SUN_RCONST(1.5)), 7, 0, &nsolves);
  fails += check(status == FIRK_NLS_DIVERGED && nsolves == 2 &&
                   fm->nls_hfac == SUN_RCONST(0.5) &&
                   close_to(fm->theta, SUN_RCONST(1.5)),
                 "rho = 1.5: divergence at k = 2, factor 0.5");

  /* predicted failure: q = 20 at k = 2 of m = 7 */
  status = newton(mem, h, mu_for_rate(h, SUN_RCONST(0.9)), 7, 0, &nsolves);
  hfac = SUN_RCONST(0.8) * SUNRpowerR(SUN_RCONST(20.0), -ONE / SUN_RCONST(8.0));
  fails += check(status == FIRK_NLS_PREDICTED_FAIL && nsolves == 2 &&
                   close_to(fm->nls_hfac, hfac),
                 "rho = 0.9: predicted failure at k = 2, factor 0.8 20^(-1/8)");

  /* iteration limit */
  status = newton(mem, h, mu_for_rate(h, SUN_RCONST(0.2)), 2, 0, &nsolves);
  fails += check(status == FIRK_NLS_MAXITER && nsolves == 2 &&
                   fm->nls_hfac == SUN_RCONST(0.5) && fm->theta == fm->theta_jac,
                 "m = 2: iteration limit, factor 0.5, rate theta_jac");

  /* recoverable and unrecoverable RHS failures */
  status = newton(mem, h, mu_for_rate(h, SUN_RCONST(0.2)), 7, 1, &nsolves);
  fails += check(status == FIRK_NLS_RHS_RECVR && nsolves == 0 &&
                   fm->nls_hfac == fm->hadapt_mem->etacf,
                 "recoverable RHS failure: factor etacf");
  status = newton(mem, h, mu_for_rate(h, SUN_RCONST(0.2)), 7, 2, &nsolves);
  fails += check(status == FIRK_RHSFUNC_FAIL && nsolves == 0,
                 "unrecoverable RHS failure: FIRK_RHSFUNC_FAIL");

  /* a non-finite update is a divergence */
  status = newton(mem, h, mu_for_rate(h, SUN_RCONST(0.2)), 7, 3, &nsolves);
  fails += check(status == FIRK_NLS_DIVERGED && nsolves == 1 &&
                   fm->nls_hfac == SUN_RCONST(0.5),
                 "non-finite update: divergence at k = 1");

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  return fails;
}

/* H19: recoverable and unrecoverable RHS failures in the stages of a step */
static int part2_case(SUNContext sunctx, int mode, int count, int expect,
                      const char* what)
{
  int fails = 0, flag;
  long int ncfn0, ncfn1;
  sunrealtype tret, hnext, hlast;
  void* mem          = NULL;
  N_Vector y         = N_VNew_Serial(1, sunctx);
  SUNMatrix A        = SUNDenseMatrix(1, 1, sunctx);
  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);

  mu             = lambda;
  NV_Ith_S(y, 0) = ONE;
  mem            = FIRKodeCreate(sunctx);
  if (FIRKodeInit(mem, f, SUN_RCONST(0.0), y) ||
      FIRKodeSStolerances(mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-8)) ||
      FIRKodeSetLinearSolver(mem, LS, A) || FIRKodeSetJacFn(mem, jac))
  {
    return 1;
  }

  /* one step, then failures at the stages of the next one (t > tn) */
  fails += check(FIRKodeEvolve(mem, SUN_RCONST(10.0), y, &tret, FIRK_ONE_STEP) ==
                   FIRK_SUCCESS,
                 "the first step succeeds");
  FIRKodeGetCurrentStep(mem, &hnext);
  FIRKodeGetNumStepSolveFails(mem, &ncfn0);
  rhs_mode     = mode;
  fail_after_t = tret;
  fail_count   = count;
  flag         = FIRKodeEvolve(mem, SUN_RCONST(10.0), y, &tret, FIRK_ONE_STEP);
  rhs_mode     = 0;
  fail_count   = 0;
  FIRKodeGetNumStepSolveFails(mem, &ncfn1);
  FIRKodeGetLastStep(mem, &hlast);

  fails += check(flag == expect, what);
  if (expect == FIRK_SUCCESS)
  {
    fails += check(ncfn1 - ncfn0 == 1 && hlast == SUN_RCONST(0.25) * hnext,
                   "one failure, the step reduced by 0.25");
  }
  else if (expect == FIRK_REPTD_RHSFUNC_ERR)
  {
    fails += check(ncfn1 - ncfn0 == MXNCF, "maxncf failures");
  }

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  return fails;
}

int main(void)
{
  SUNContext sunctx = NULL;
  int fails         = 0;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  fails += part1(sunctx);
  fails += part2_case(sunctx, 1, 1, FIRK_SUCCESS,
                      "a recoverable stage failure is recovered");
  fails += part2_case(sunctx, 1, 100, FIRK_REPTD_RHSFUNC_ERR,
                      "repeated recoverable failures: FIRK_REPTD_RHSFUNC_ERR");
  fails += part2_case(sunctx, 2, 1, FIRK_RHSFUNC_FAIL,
                      "an unrecoverable stage failure: FIRK_RHSFUNC_FAIL");

  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("All tests passed\n");
  return 0;
}
