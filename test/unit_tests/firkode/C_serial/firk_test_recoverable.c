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
 * Unit test for failures of the right-hand side function in a stage (hazard
 * H19, design Listing 8.2) and in the error estimate (hazard H20, design
 * Listing 8.5).
 *
 * H19, on y' = -y with a first step of 0.1, failures injected by call count:
 *   - one recoverable failure in the first stage evaluation: the attempt
 *     fails, and the step is retried with h reduced by etacf = 0.25;
 *   - recoverable failures from then on: FIRK_REPTD_RHSFUNC_ERR after
 *     maxncf = 10 failed attempts;
 *   - one unrecoverable failure: FIRK_RHSFUNC_FAIL.
 *
 * H20:
 *
 * y' = -y from y(0) = 1 with a first step of 1 fails the error test, so the
 * estimate of the first step, and of every rejected attempt after it, is
 * filtered again with f evaluated at (t0, y0 + err). The right-hand side
 * function fails there, and only there (at t0 with y != y0; the Jacobian is
 * analytic):
 *   - recoverably: the first estimate is kept, the step is rejected by the
 *     error test alone (no solve failure), and the integration to t = 1
 *     succeeds with the accuracy of the tolerances;
 *   - unrecoverably: FIRKodeEvolve returns FIRK_RHSFUNC_FAIL.
 *
 * Pass criterion: the failures are injected, the return values match, no step
 * fails in its solve, and the solution error is below 10 times the tolerance.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

static int fail_value = 0; /* +1 recoverable, -1 unrecoverable */
static int ninjected  = 0;

/* H19: f fails with stage_fail_value on its calls fail_from to fail_to
   (0: forever) */
static long int ncalls = 0, fail_from = 0, fail_to = 0;
static int stage_fail_value = 0;

static int f_stage(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)t;
  (void)user_data;
  ncalls++;
  if (fail_from > 0 && ncalls >= fail_from && (fail_to == 0 || ncalls <= fail_to))
  {
    return stage_fail_value;
  }
  NV_Ith_S(ydot, 0) = -NV_Ith_S(y, 0);
  return 0;
}

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)user_data;
  if (t == ZERO && NV_Ith_S(y, 0) != ONE)
  {
    ninjected++;
    return fail_value;
  }
  NV_Ith_S(ydot, 0) = -NV_Ith_S(y, 0);
  return 0;
}

/* an analytic Jacobian, so that only the refilter evaluates f at t0 with
   y != y0 (a difference quotient would too) */
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
  SM_ELEMENT_D(J, 0, 0) = -ONE;
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

static int run(SUNContext sunctx, int value, int expect, const char* what)
{
  const sunrealtype rtol = SUN_RCONST(1.0e-6);
  int fails              = 0, flag;
  long int ncfn          = -1;
  sunrealtype tret, err;
  void* mem          = NULL;
  N_Vector y         = N_VNew_Serial(1, sunctx);
  SUNMatrix A        = SUNDenseMatrix(1, 1, sunctx);
  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);

  fail_value     = value;
  ninjected      = 0;
  NV_Ith_S(y, 0) = ONE;
  mem            = FIRKodeCreate(sunctx);
  if (FIRKodeInit(mem, f, ZERO, y) || FIRKodeSStolerances(mem, rtol, rtol) ||
      FIRKodeSetLinearSolver(mem, LS, A) || FIRKodeSetJacFn(mem, jac) ||
      FIRKodeSetMode(mem, FIRK_EXPERT) || FIRKodeSetInitStep(mem, ONE))
  {
    return 1;
  }

  flag = FIRKodeEvolve(mem, ONE, y, &tret, FIRK_NORMAL);
  FIRKodeGetNumStepSolveFails(mem, &ncfn);
  err = SUNRabs(NV_Ith_S(y, 0) - SUNRexp(-ONE));

  fails += check(ninjected > 0, "the refilter evaluates f at (t0, y0 + err)");
  fails += check(flag == expect, what);
  if (expect == FIRK_SUCCESS)
  {
    fails += check(ncfn == 0, "no step fails in its solve");
    fails += check(err <= SUN_RCONST(10.0) * rtol, "the solution is accurate");
  }

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  return fails;
}

/* H19: the first step of 0.1 with failures from the second call (the first
   is f(t0, y0)); returns the flag of FIRK_ONE_STEP, the step taken and the
   failed attempts */
static int run_stage(SUNContext sunctx, int value, long int last,
                     sunrealtype* hlast, long int* ncfn)
{
  int flag;
  sunrealtype tret;
  void* mem          = NULL;
  N_Vector y         = N_VNew_Serial(1, sunctx);
  SUNMatrix A        = SUNDenseMatrix(1, 1, sunctx);
  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);

  ncalls           = 0;
  fail_from        = 2;
  fail_to          = last;
  stage_fail_value = value;
  NV_Ith_S(y, 0)   = ONE;
  mem              = FIRKodeCreate(sunctx);
  if (FIRKodeInit(mem, f_stage, ZERO, y) ||
      FIRKodeSStolerances(mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-6)) ||
      FIRKodeSetLinearSolver(mem, LS, A) || FIRKodeSetJacFn(mem, jac) ||
      FIRKodeSetInitStep(mem, SUN_RCONST(0.1)))
  {
    return -999;
  }
  flag = FIRKodeEvolve(mem, ONE, y, &tret, FIRK_ONE_STEP);
  FIRKodeGetLastStep(mem, hlast);
  FIRKodeGetNumStepSolveFails(mem, ncfn);

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  fail_from = 0;
  return flag;
}

static int check_stage_failures(SUNContext sunctx)
{
  int fails         = 0, flag;
  sunrealtype hlast = ZERO;
  long int ncfn     = -1;
  char what[160];

  flag = run_stage(sunctx, 1, 2, &hlast, &ncfn);
  snprintf(what, sizeof(what),
           "H19: one recoverable stage failure: flag %d, step " SUN_FORMAT_G
           " (0.25 x 0.1), failed attempts %ld",
           flag, hlast, ncfn);
  fails += check(flag == FIRK_SUCCESS && ncfn == 1 &&
                   SUNRabs(hlast - SUN_RCONST(0.025)) <=
                     SUN_RCONST(10.0) * SUN_UNIT_ROUNDOFF,
                 what);

  flag = run_stage(sunctx, 1, 0, &hlast, &ncfn);
  snprintf(what, sizeof(what),
           "H19: repeated recoverable failures: flag %d, failed attempts %ld",
           flag, ncfn);
  fails += check(flag == FIRK_REPTD_RHSFUNC_ERR && ncfn == 10, what);

  flag = run_stage(sunctx, -1, 2, &hlast, &ncfn);
  snprintf(what, sizeof(what), "H19: an unrecoverable failure: flag %d", flag);
  fails += check(flag == FIRK_RHSFUNC_FAIL, what);
  return fails;
}

int main(void)
{
  SUNContext sunctx = NULL;
  int fails         = 0;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  fails += check_stage_failures(sunctx);
  fails += run(sunctx, 1, FIRK_SUCCESS,
               "a recoverable refilter failure keeps the first estimate");
  fails += run(sunctx, -1, FIRK_RHSFUNC_FAIL,
               "an unrecoverable refilter failure: FIRK_RHSFUNC_FAIL");

  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("All tests passed\n");
  return 0;
}
