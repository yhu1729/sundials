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
 * Unit test for the Jacobian and setup reuse policies (design decision D9),
 * under FIRK_POLICY_RADAU5 and FIRK_POLICY_PERIODIC, with a counting analytic
 * Jacobian:
 *
 * - H8: the first step fails twice in its stage evaluations (a recoverable RHS
 *   failure); its Jacobian, evaluated at the first attempt, is current, so the
 *   two retries set up again (with the reduced step) but never evaluate it
 *   again: one Jacobian, three setups;
 * - H1: with SetJacBadConvRate 1e-12, each step after a Newton solve whose
 *   rate exceeded it evaluates the Jacobian and sets up, also under
 *   PERIODIC, which otherwise reuses the Jacobian for up to 51 steps, and
 *   under RADAU5 after rates of at most theta_jac = 1e-3, which RADAU5
 *   otherwise reuses it for;
 * - H30: a step of one attempt whose size the dead band kept, and whose
 *   Jacobian was not re-evaluated, does no setup; under RADAU5 a step of one
 *   attempt after a Newton rate of at most theta_jac reuses the Jacobian.
 *
 * Pass criterion: every count matches.
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

/* y' = lambda y + eps y^3, a mildly nonlinear decay */
static const sunrealtype lambda = SUN_RCONST(-10.0);
static const sunrealtype eps    = SUN_RCONST(0.5);

static int fail_count = 0; /* recoverable failures left to inject at t > 0 */
static long int njac  = 0;

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype u = NV_Ith_S(y, 0);
  (void)user_data;
  if (t > ZERO && fail_count > 0)
  {
    fail_count--;
    return 1;
  }
  NV_Ith_S(ydot, 0) = lambda * u + eps * u * u * u;
  return 0;
}

static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype u = NV_Ith_S(y, 0);
  (void)t;
  (void)fy;
  (void)user_data;
  (void)tmp1;
  (void)tmp2;
  (void)tmp3;
  njac++;
  SM_ELEMENT_D(J, 0, 0) = lambda + SUN_RCONST(3.0) * eps * u * u;
  return 0;
}

static int check(int ok, const char* what, int policy)
{
  if (!ok)
  {
    printf("  FAIL (%s): %s\n",
           (policy == FIRK_POLICY_RADAU5) ? "RADAU5" : "PERIODIC", what);
    return 1;
  }
  return 0;
}

struct problem
{
  void* mem;
  N_Vector y;
  SUNMatrix A;
  SUNLinearSolver LS;
};

static int setup(SUNContext sunctx, struct problem* P, int policy)
{
  P->y              = N_VNew_Serial(1, sunctx);
  P->A              = SUNDenseMatrix(1, 1, sunctx);
  P->LS             = SUNLinSol_Dense(P->y, P->A, sunctx);
  NV_Ith_S(P->y, 0) = ONE;
  P->mem            = FIRKodeCreate(sunctx);
  njac              = 0;
  fail_count        = 0;
  return FIRKodeInit(P->mem, f, ZERO, P->y) ||
         FIRKodeSStolerances(P->mem, SUN_RCONST(1.0e-5), SUN_RCONST(1.0e-5)) ||
         FIRKodeSetLinearSolver(P->mem, P->LS, P->A) ||
         FIRKodeSetJacFn(P->mem, jac) || FIRKodeSetMode(P->mem, FIRK_EXPERT) ||
         FIRKodeSetReusePolicy(P->mem, policy);
}

static void teardown(struct problem* P)
{
  FIRKodeFree(&P->mem);
  SUNLinSolFree(P->LS);
  SUNMatDestroy(P->A);
  N_VDestroy(P->y);
}

/* H8: two failed attempts with a current Jacobian */
static int check_h8(SUNContext sunctx, int policy)
{
  int fails        = 0;
  long int nsetups = 0, nje = 0, ncfn = 0;
  sunrealtype tret;
  struct problem P;

  if (setup(sunctx, &P, policy)) { return 1; }
  if (FIRKodeSetInitStep(P.mem, SUN_RCONST(0.01))) { return 1; }
  fail_count = 2;
  fails +=
    check(FIRKodeEvolve(P.mem, ONE, P.y, &tret, FIRK_ONE_STEP) == FIRK_SUCCESS,
          "the first step succeeds", policy);
  FIRKodeGetNumLinSolvSetups(P.mem, &nsetups);
  FIRKodeGetNumJacEvals(P.mem, &nje);
  FIRKodeGetNumStepSolveFails(P.mem, &ncfn);
  if (ncfn != 2 || njac != 1 || nje != 1 || nsetups != 3)
  {
    printf("  ncfn = %ld, jac calls = %ld, nje = %ld, nsetups = %ld\n", ncfn,
           njac, nje, nsetups);
  }
  fails += check(ncfn == 2 && njac == 1 && nje == 1 && nsetups == 3,
                 "H8: one Jacobian and three setups across two failures", policy);
  teardown(&P);
  return fails;
}

/* H1: a Newton rate above the bad-convergence rate evaluates J at the next
   step; H30: steps kept in the dead band reuse the setup */
static int check_h1_h30(SUNContext sunctx, int policy, sunbooleantype h1)
{
  int k, fails = 0, ok_h1 = 1, ok_h30 = 1, ok_reuse = 1, nkept = 0;
  int nslow_fast = 0; /* H1 steps after a rate in (1e-12, theta_jac] */
  long int nsetups0, nsetups1, nje0, nje1, natt0, natt1;
  sunrealtype tret, h0, h1v, theta, theta_prev = ZERO;
  struct problem P;

  if (setup(sunctx, &P, policy)) { return 1; }
  if (h1 && FIRKodeSetJacBadConvRate(P.mem, SUN_RCONST(1.0e-12))) { return 1; }

  /* the first step */
  if (FIRKodeEvolve(P.mem, SUN_RCONST(5.0), P.y, &tret, FIRK_ONE_STEP))
  {
    return 1;
  }
  FIRKodeGetLastStep(P.mem, &h0);
  FIRKodeGetLastNonlinConvRate(P.mem, &theta_prev);

  for (k = 0; k < 40 && tret < SUN_RCONST(5.0); k++)
  {
    FIRKodeGetNumLinSolvSetups(P.mem, &nsetups0);
    FIRKodeGetNumJacEvals(P.mem, &nje0);
    FIRKodeGetNumStepAttempts(P.mem, &natt0);
    if (FIRKodeEvolve(P.mem, SUN_RCONST(5.0), P.y, &tret, FIRK_ONE_STEP))
    {
      return 1;
    }
    FIRKodeGetNumLinSolvSetups(P.mem, &nsetups1);
    FIRKodeGetNumJacEvals(P.mem, &nje1);
    FIRKodeGetNumStepAttempts(P.mem, &natt1);
    FIRKodeGetLastStep(P.mem, &h1v);
    FIRKodeGetLastNonlinConvRate(P.mem, &theta);

    /* H1: a Newton rate above 1e-12 evaluates J and sets up at the next
       step */
    if (h1 && theta_prev > SUN_RCONST(1.0e-12))
    {
      if (nje1 == nje0 || nsetups1 == nsetups0)
      {
        printf("  step %d: theta_prev = %g, nje %ld -> %ld, nsetups %ld -> "
               "%ld\n",
               k, (double)theta_prev, nje0, nje1, nsetups0, nsetups1);
        ok_h1 = 0;
      }
      if (theta_prev <= SUN_RCONST(1.0e-3)) { nslow_fast++; }
    }

    /* H30: one attempt of the same step size and no new Jacobian: no setup */
    if (natt1 - natt0 == 1 && h1v == h0 && nje1 == nje0)
    {
      nkept++;
      if (nsetups1 != nsetups0) { ok_h30 = 0; }
    }

    /* RADAU5: a converged rate of at most theta_jac reuses the Jacobian on
       the next step, unless an attempt fails */
    if (!h1 && policy == FIRK_POLICY_RADAU5 && natt1 - natt0 == 1 &&
        theta_prev <= SUN_RCONST(1.0e-3) && nje1 != nje0)
    {
      printf("  step %d: theta_prev = %g, nje %ld -> %ld\n", k,
             (double)theta_prev, nje0, nje1);
      ok_reuse = 0;
    }

    h0         = h1v;
    theta_prev = theta;
  }

  if (h1)
  {
    fails += check(ok_h1, "H1: every step evaluates J and sets up", policy);
    fails += check(nslow_fast > 0,
                   "H1: some steps follow a rate in (1e-12, 1e-3]", policy);
  }
  else
  {
    fails += check(nkept > 0, "some steps keep their size", policy);
    fails += check(ok_h30, "H30: a kept step does no setup", policy);
    fails += check(ok_reuse, "RADAU5: a fast rate reuses J", policy);
  }
  teardown(&P);
  return fails;
}

int main(void)
{
  SUNContext sunctx     = NULL;
  int fails             = 0, i;
  const int policies[2] = {FIRK_POLICY_RADAU5, FIRK_POLICY_PERIODIC};

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  for (i = 0; i < 2; i++)
  {
    fails += check_h8(sunctx, policies[i]);
    fails += check_h1_h30(sunctx, policies[i], SUNTRUE);
    fails += check_h1_h30(sunctx, policies[i], SUNFALSE);
  }

  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("All tests passed\n");
  return 0;
}
