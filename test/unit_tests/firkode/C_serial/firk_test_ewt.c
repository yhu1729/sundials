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
 * Unit test of the error weights (design hazard H35) on y' = (-1, -y_2):
 *
 * - a zero weight is FIRK_ILL_INPUT: FIRKodeSStolerances(rtol, 0) with
 *   y_1(0) = 0 at the first step; with y_1(0) = 1, fixed steps of 1/4 reach
 *   y_1 = 1 - t = 0 exactly at t = 1, where the next step is refused;
 * - the weights are evaluated at y_n once per step, whatever the number of
 *   attempts: a user weight function (FIRKodeWFtolerances) is called once
 *   per FIRKodeEvolve call in FIRK_ONE_STEP mode after the first, with the
 *   solution the previous call returned;
 * - a second FIRKodeInit discards the tolerances (its vectors, among them
 *   the vector absolute tolerances, are freed): FIRKodeEvolve then returns
 *   FIRK_ILL_INPUT until they are set again, and succeeds afterwards;
 * - the tolerances may be changed between FIRKodeEvolve calls: each of the
 *   three tolerance functions takes effect at the next call, and a user
 *   weight function set then receives the user data pointer.
 *
 * Pass criterion: the return values, the call counts and the arguments match.
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
#define RTOL SUN_RCONST(1.0e-6)

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)t;
  (void)user_data;
  NV_Ith_S(ydot, 0) = -ONE;
  NV_Ith_S(ydot, 1) = -NV_Ith_S(y, 1);
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
  SM_ELEMENT_D(J, 1, 1) = -ONE;
  return 0;
}

/* user weights, counting their calls and recording their argument */
static long int nefun = 0;
static sunrealtype efun_y[2];
static void* efun_data = NULL;

static int efun(N_Vector y, N_Vector ewt, void* user_data)
{
  int i;
  efun_data = user_data;
  nefun++;
  for (i = 0; i < 2; i++)
  {
    efun_y[i]        = NV_Ith_S(y, i);
    NV_Ith_S(ewt, i) = ONE / (RTOL * SUNRabs(NV_Ith_S(y, i)) + RTOL);
  }
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

struct run
{
  N_Vector y;
  SUNMatrix A;
  SUNLinearSolver LS;
  void* mem;
};

static int setup(SUNContext sunctx, struct run* R, sunrealtype y1)
{
  R->y              = N_VNew_Serial(2, sunctx);
  NV_Ith_S(R->y, 0) = y1;
  NV_Ith_S(R->y, 1) = ONE;
  R->A              = SUNDenseMatrix(2, 2, sunctx);
  R->LS             = SUNLinSol_Dense(R->y, R->A, sunctx);
  R->mem            = FIRKodeCreate(sunctx);
  return FIRKodeInit(R->mem, f, ZERO, R->y) ||
         FIRKodeSetLinearSolver(R->mem, R->LS, R->A) ||
         FIRKodeSetJacFn(R->mem, jac);
}

static void teardown(struct run* R)
{
  FIRKodeFree(&R->mem);
  SUNLinSolFree(R->LS);
  SUNMatDestroy(R->A);
  N_VDestroy(R->y);
}

int main(void)
{
  SUNContext sunctx = NULL;
  int fails         = 0, flag, k;
  long int n0;
  sunrealtype tret, yprev[2];
  char what[200];
  struct run R;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  /* a zero weight at the first step */
  if (setup(sunctx, &R, ZERO) || FIRKodeSStolerances(R.mem, RTOL, ZERO))
  {
    return 1;
  }
  flag = FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_NORMAL);
  snprintf(what, sizeof(what), "a zero weight at t0 is FIRK_ILL_INPUT (%d)",
           flag);
  fails += check(flag == FIRK_ILL_INPUT, what);
  teardown(&R);

  /* a weight that becomes zero: y_1 = 1 - t, exact at the fixed steps */
  if (setup(sunctx, &R, ONE) || FIRKodeSStolerances(R.mem, RTOL, ZERO) ||
      FIRKodeSetFixedStep(R.mem, SUN_RCONST(0.25)))
  {
    return 1;
  }
  flag = FIRKodeEvolve(R.mem, SUN_RCONST(2.0), R.y, &tret, FIRK_NORMAL);
  snprintf(what, sizeof(what),
           "a weight that becomes zero is FIRK_ILL_INPUT (%d) at t = 1 "
           "(" SUN_FORMAT_G "), y_1 = " SUN_FORMAT_G,
           flag, tret, NV_Ith_S(R.y, 0));
  fails += check(flag == FIRK_ILL_INPUT && tret == ONE && NV_Ith_S(R.y, 0) == ZERO,
                 what);
  teardown(&R);

  /* weights once per step, at y_n */
  if (setup(sunctx, &R, ONE) || FIRKodeWFtolerances(R.mem, efun)) { return 1; }
  nefun = 0;
  flag  = FIRKodeEvolve(R.mem, SUN_RCONST(0.5), R.y, &tret, FIRK_ONE_STEP);
  fails += check(flag == FIRK_SUCCESS, "the first step succeeds");
  for (k = 0; k < 10 && flag == FIRK_SUCCESS; k++)
  {
    yprev[0] = NV_Ith_S(R.y, 0);
    yprev[1] = NV_Ith_S(R.y, 1);
    n0       = nefun;
    flag     = FIRKodeEvolve(R.mem, SUN_RCONST(0.5), R.y, &tret, FIRK_ONE_STEP);
    snprintf(what, sizeof(what),
             "step %d: %ld weight evaluation(s), at the solution of the "
             "previous step",
             k + 2, nefun - n0);
    fails += check(flag >= 0 && nefun - n0 == 1 && efun_y[0] == yprev[0] &&
                     efun_y[1] == yprev[1],
                   what);
  }
  teardown(&R);

  /* a second Init discards the tolerances, vector ones included */
  if (setup(sunctx, &R, ONE)) { return 1; }
  {
    N_Vector abstol = N_VClone(R.y);
    N_VConst(SUN_RCONST(1.0e-8), abstol);
    flag = FIRKodeSVtolerances(R.mem, RTOL, abstol);
    N_VDestroy(abstol);
  }
  fails += check(flag == FIRK_SUCCESS, "SVtolerances before the second Init");
  fails += check(FIRKodeInit(R.mem, f, ZERO, R.y) == FIRK_SUCCESS,
                 "a second Init succeeds");
  flag = FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP);
  fails += check(flag == FIRK_ILL_INPUT,
                 "Evolve after a second Init without tolerances is "
                 "FIRK_ILL_INPUT");
  fails += check(FIRKodeSStolerances(R.mem, RTOL, RTOL) == FIRK_SUCCESS,
                 "tolerances set again after the second Init");
  flag = FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP);
  fails += check(flag == FIRK_SUCCESS, "Evolve succeeds with new tolerances");
  teardown(&R);

  /* tolerances changed between Evolve calls: scalar, vector, then a user
     weight function, which receives the user data pointer */
  if (setup(sunctx, &R, ONE) || FIRKodeSStolerances(R.mem, RTOL, RTOL) ||
      FIRKodeSetUserData(R.mem, &R))
  {
    return 1;
  }
  flag = FIRKodeEvolve(R.mem, SUN_RCONST(0.5), R.y, &tret, FIRK_ONE_STEP);
  fails += check(flag == FIRK_SUCCESS, "the first step succeeds");
  fails += check(FIRKodeSStolerances(R.mem, RTOL, RTOL * RTOL) == FIRK_SUCCESS,
                 "scalar tolerances between Evolve calls");
  flag = FIRKodeEvolve(R.mem, SUN_RCONST(0.5), R.y, &tret, FIRK_ONE_STEP);
  fails += check(flag == FIRK_SUCCESS, "a step after new scalar tolerances");
  {
    N_Vector abstol = N_VClone(R.y);
    N_VConst(RTOL * RTOL, abstol);
    flag = FIRKodeSVtolerances(R.mem, RTOL, abstol);
    N_VDestroy(abstol);
  }
  fails += check(flag == FIRK_SUCCESS,
                 "vector tolerances between Evolve calls");
  flag = FIRKodeEvolve(R.mem, SUN_RCONST(0.5), R.y, &tret, FIRK_ONE_STEP);
  fails += check(flag == FIRK_SUCCESS, "a step after new vector tolerances");
  fails += check(FIRKodeWFtolerances(R.mem, efun) == FIRK_SUCCESS,
                 "a weight function between Evolve calls");
  nefun     = 0;
  efun_data = NULL;
  flag      = FIRKodeEvolve(R.mem, SUN_RCONST(0.5), R.y, &tret, FIRK_ONE_STEP);
  snprintf(what, sizeof(what),
           "a step after a new weight function: %ld evaluation(s), user data "
           "%s",
           nefun, (efun_data == &R) ? "passed" : "missing");
  fails += check(flag == FIRK_SUCCESS && nefun == 1 && efun_data == &R, what);
  teardown(&R);

  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("All tests passed\n");
  return 0;
}
