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
 * Unit test for setting the stop time (adapted from the CVODE test), with
 * adaptive steps and with a fixed step size, forward and backward in time
 * (design hazard H21): stop returns, a stop time behind the current time
 * refused, the stop time cleared by its return, full fixed steps resumed
 * after the step shortened to reach it, and the shortened step itself,
 * (tstop - tn)(1 - 4u), after which the current time is the stop time.
 * y' = 1, so y = t.
 *
 * Pass criterion: every return is FIRK_SUCCESS or FIRK_TSTOP_RETURN as
 * expected, at the expected time to 100 u; a stop time behind the current
 * time is FIRK_ILL_INPUT, forward and backward; after the shortened step the
 * forward fixed-step run ends with 11 steps, the last of the fixed size, and
 * y = t to 1e-4; in one-step mode the steps before the stop end before it,
 * the step that reaches it is (tstop - tn)(1 - 4u) and leaves the current
 * time at tstop exactly.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_matrix.h"
#include "sundials/sundials_nvector.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#else
#define GSYM "g"
#endif

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

static int ode_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* ydot_data = N_VGetArrayPointer(ydot);
  ydot_data[0]           = ONE;
  return 0;
}

static int ode_jac(sunrealtype t, N_Vector y, N_Vector f, SUNMatrix J,
                   void* user_data, N_Vector tempv1, N_Vector tempv2,
                   N_Vector tempv3)
{
  sunrealtype* J_data = SUNDenseMatrix_Data(J);
  J_data[0]           = ZERO;
  return 0;
}

/* In fixed-step mode the step shortened to reach the stop time must not
   shrink the steps that follow it */
static int test_fixed_step(SUNContext sunctx)
{
  N_Vector y         = NULL;
  SUNMatrix A        = NULL;
  SUNLinearSolver LS = NULL;
  void* firkode_mem  = NULL;
  int flag, result = 1;
  long int nst      = 0;
  sunrealtype tret  = ZERO;
  sunrealtype hlast = ZERO;

  const sunrealtype hfixed = SUN_RCONST(0.1);
  const sunrealtype tstop  = SUN_RCONST(0.25);
  const sunrealtype tout   = ONE;
  const sunrealtype tol    = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  y = N_VNew_Serial(1, sunctx);
  if (!y) { goto cleanup; }
  N_VConst(ZERO, y);

  firkode_mem = FIRKodeCreate(sunctx);
  if (!firkode_mem) { goto cleanup; }
  if (FIRKodeInit(firkode_mem, ode_rhs, ZERO, y)) { goto cleanup; }
  if (FIRKodeSStolerances(firkode_mem, SUN_RCONST(1.0e-4), SUN_RCONST(1.0e-8)))
  {
    goto cleanup;
  }

  A = SUNDenseMatrix(1, 1, sunctx);
  if (!A) { goto cleanup; }
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (!LS) { goto cleanup; }
  if (FIRKodeSetLinearSolver(firkode_mem, LS, A)) { goto cleanup; }
  if (FIRKodeSetJacFn(firkode_mem, ode_jac)) { goto cleanup; }

  if (FIRKodeSetFixedStep(firkode_mem, hfixed)) { goto cleanup; }
  if (FIRKodeSetStopTime(firkode_mem, tstop)) { goto cleanup; }

  /* steps 0.1, 0.1 and 0.05 reach the stop time */
  flag = FIRKodeEvolve(firkode_mem, tout, y, &tret, FIRK_NORMAL);
  printf("fixed step: return = %i, tret = %" GSYM "\n", flag, tret);
  if (flag != FIRK_TSTOP_RETURN || SUNRabs(tret - tstop) > tol)
  {
    printf("ERROR: Expected stop return at %" GSYM "!\n", tstop);
    goto cleanup;
  }

  /* eight more steps of 0.1 pass the output time */
  flag = FIRKodeEvolve(firkode_mem, tout, y, &tret, FIRK_NORMAL);
  if (flag < 0) { goto cleanup; }
  if (FIRKodeGetLastStep(firkode_mem, &hlast)) { goto cleanup; }
  if (FIRKodeGetNumSteps(firkode_mem, &nst)) { goto cleanup; }
  printf("fixed step: return = %i, tret = %" GSYM ", last step = %" GSYM
         ", steps = %li, y = %" GSYM "\n",
         flag, tret, hlast, nst, N_VGetArrayPointer(y)[0]);
  /* y = t up to the integration tolerance: when the residual of the predicted
     stages is below the tolerance of the stacked Krylov solve, the solve
     returns the zero correction and the Newton iteration accepts them */
  if (flag != FIRK_SUCCESS || SUNRabs(tret - tout) > tol ||
      SUNRabs(hlast - hfixed) > tol || nst != 11 ||
      SUNRabs(N_VGetArrayPointer(y)[0] - tout) > SUN_RCONST(1.0e-4))
  {
    printf("ERROR: Expected steps of %" GSYM " after the stop time!\n", hfixed);
    goto cleanup;
  }

  result = 0;

cleanup:
  FIRKodeFree(&firkode_mem);
  N_VDestroy(y);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  return result;
}

/* In one-step mode with adaptive steps: the steps before the stop time end
   before it; the step that reaches it is (tstop - tn)(1 - 4u), after which
   the current time is tstop; once past it, a stop time behind the current
   time is refused */
static int test_clamp(SUNContext sunctx)
{
  N_Vector y         = NULL;
  SUNMatrix A        = NULL;
  SUNLinearSolver LS = NULL;
  void* firkode_mem  = NULL;
  int k, flag, result = 1;
  sunrealtype tret = ZERO, tprev = ZERO, tcur = ZERO, hlast = ZERO, hclamp;

  const sunrealtype tstop = SUN_RCONST(0.25);
  const sunrealtype tout  = ONE;

  y = N_VNew_Serial(1, sunctx);
  if (!y) { goto cleanup; }
  N_VConst(ZERO, y);

  firkode_mem = FIRKodeCreate(sunctx);
  if (!firkode_mem) { goto cleanup; }
  if (FIRKodeInit(firkode_mem, ode_rhs, ZERO, y)) { goto cleanup; }
  if (FIRKodeSStolerances(firkode_mem, SUN_RCONST(1.0e-4), SUN_RCONST(1.0e-8)))
  {
    goto cleanup;
  }
  A = SUNDenseMatrix(1, 1, sunctx);
  if (!A) { goto cleanup; }
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (!LS) { goto cleanup; }
  if (FIRKodeSetLinearSolver(firkode_mem, LS, A)) { goto cleanup; }
  if (FIRKodeSetJacFn(firkode_mem, ode_jac)) { goto cleanup; }
  if (FIRKodeSetStopTime(firkode_mem, tstop)) { goto cleanup; }

  for (k = 0; k < 100; k++)
  {
    if (FIRKodeGetCurrentTime(firkode_mem, &tprev)) { goto cleanup; }
    flag = FIRKodeEvolve(firkode_mem, tout, y, &tret, FIRK_ONE_STEP);
    if (FIRKodeGetCurrentTime(firkode_mem, &tcur)) { goto cleanup; }
    if (flag == FIRK_SUCCESS && tcur < tstop) { continue; }
    if (FIRKodeGetLastStep(firkode_mem, &hlast)) { goto cleanup; }
    hclamp = (tstop - tprev) * (ONE - SUN_RCONST(4.0) * SUN_UNIT_ROUNDOFF);
    printf("clamp: return = %i at step %i, tret = %" GSYM ", tcur = %" GSYM
           ", last step = %" GSYM "\n",
           flag, k, tret, tcur, hlast);
    if (flag != FIRK_TSTOP_RETURN || tret != tstop || tcur != tstop ||
        hlast != hclamp)
    {
      printf("ERROR: Expected a stop return at %" GSYM
             " after the step (tstop - tn)(1 - 4u) = %" GSYM "!\n",
             tstop, hclamp);
      goto cleanup;
    }
    break;
  }
  if (k == 100) { goto cleanup; }

  /* past the stop time, a stop time behind the current time is refused */
  flag = FIRKodeEvolve(firkode_mem, SUN_RCONST(0.5), y, &tret, FIRK_NORMAL);
  if (flag != FIRK_SUCCESS) { goto cleanup; }
  if (FIRKodeSetStopTime(firkode_mem, SUN_RCONST(0.3)) != FIRK_ILL_INPUT)
  {
    printf("ERROR: Expected a stop time behind the current time refused!\n");
    goto cleanup;
  }

  result = 0;

cleanup:
  FIRKodeFree(&firkode_mem);
  N_VDestroy(y);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  return result;
}

/* Backward in time from t0 = 0: adaptive steps reach tout = -0.1, stop at
   tstop = -0.3, refuse a stop time of -0.2 (passed), integrate past -0.3 once
   the stop time is cleared; fixed steps of -0.1 stop at -0.25 and then take
   full steps to -1 */
static int test_backward(SUNContext sunctx)
{
  N_Vector y         = NULL;
  SUNMatrix A        = NULL;
  SUNLinearSolver LS = NULL;
  void* mem          = NULL;
  int flag, k, result = 1;
  long int nst          = 0;
  sunrealtype tret      = ZERO;
  sunrealtype hlast     = ZERO;
  const sunrealtype tol = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  for (k = 0; k < 2; k++)
  {
    y = N_VNew_Serial(1, sunctx);
    N_VConst(ZERO, y);
    A   = SUNDenseMatrix(1, 1, sunctx);
    LS  = SUNLinSol_Dense(y, A, sunctx);
    mem = FIRKodeCreate(sunctx);
    if (FIRKodeInit(mem, ode_rhs, ZERO, y) ||
        FIRKodeSStolerances(mem, SUN_RCONST(1.0e-4), SUN_RCONST(1.0e-8)) ||
        FIRKodeSetLinearSolver(mem, LS, A) || FIRKodeSetJacFn(mem, ode_jac))
    {
      goto cleanup;
    }

    if (k == 0)
    {
      /* adaptive */
      if (FIRKodeSetMaxStep(mem, SUN_RCONST(0.05)) ||
          FIRKodeSetStopTime(mem, -SUN_RCONST(0.3)))
      {
        goto cleanup;
      }
      flag = FIRKodeEvolve(mem, -SUN_RCONST(0.1), y, &tret, FIRK_NORMAL);
      printf("backward: return = %i, tret = %" GSYM "\n", flag, tret);
      if (flag != FIRK_SUCCESS || SUNRabs(tret + SUN_RCONST(0.1)) > tol)
      {
        printf("ERROR: Expected an output return at -0.1!\n");
        goto cleanup;
      }
      flag = FIRKodeEvolve(mem, -SUN_RCONST(0.35), y, &tret, FIRK_NORMAL);
      printf("backward: return = %i, tret = %" GSYM ", y = %" GSYM "\n", flag,
             tret, N_VGetArrayPointer(y)[0]);
      if (flag != FIRK_TSTOP_RETURN || SUNRabs(tret + SUN_RCONST(0.3)) > tol ||
          SUNRabs(N_VGetArrayPointer(y)[0] - tret) > SUN_RCONST(1.0e-4))
      {
        printf("ERROR: Expected a stop return at -0.3!\n");
        goto cleanup;
      }
      if (FIRKodeSetStopTime(mem, -SUN_RCONST(0.2)) != FIRK_ILL_INPUT)
      {
        printf("ERROR: Expected a passed stop time to be refused!\n");
        goto cleanup;
      }
      flag = FIRKodeEvolve(mem, -SUN_RCONST(0.6), y, &tret, FIRK_NORMAL);
      printf("backward: return = %i, tret = %" GSYM "\n", flag, tret);
      if (flag != FIRK_SUCCESS || SUNRabs(tret + SUN_RCONST(0.6)) > tol)
      {
        printf("ERROR: Expected the stop time to be cleared!\n");
        goto cleanup;
      }
    }
    else
    {
      /* fixed steps */
      if (FIRKodeSetFixedStep(mem, -SUN_RCONST(0.1)) ||
          FIRKodeSetStopTime(mem, -SUN_RCONST(0.25)))
      {
        goto cleanup;
      }
      flag = FIRKodeEvolve(mem, -ONE, y, &tret, FIRK_NORMAL);
      printf("backward fixed step: return = %i, tret = %" GSYM "\n", flag, tret);
      if (flag != FIRK_TSTOP_RETURN || SUNRabs(tret + SUN_RCONST(0.25)) > tol)
      {
        printf("ERROR: Expected a stop return at -0.25!\n");
        goto cleanup;
      }
      flag = FIRKodeEvolve(mem, -ONE, y, &tret, FIRK_NORMAL);
      FIRKodeGetLastStep(mem, &hlast);
      FIRKodeGetNumSteps(mem, &nst);
      printf("backward fixed step: return = %i, tret = %" GSYM
             ", last step = %" GSYM ", steps = %li, y = %" GSYM "\n",
             flag, tret, hlast, nst, N_VGetArrayPointer(y)[0]);
      if (flag != FIRK_SUCCESS || SUNRabs(tret + ONE) > tol ||
          SUNRabs(hlast + SUN_RCONST(0.1)) > tol || nst != 11 ||
          SUNRabs(N_VGetArrayPointer(y)[0] + ONE) > SUN_RCONST(1.0e-4))
      {
        printf("ERROR: Expected steps of -0.1 after the stop time!\n");
        goto cleanup;
      }
    }

    FIRKodeFree(&mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    N_VDestroy(y);
    mem = NULL;
    LS  = NULL;
    A   = NULL;
    y   = NULL;
  }
  result = 0;

cleanup:
  FIRKodeFree(&mem);
  if (LS) { SUNLinSolFree(LS); }
  if (A) { SUNMatDestroy(A); }
  if (y) { N_VDestroy(y); }
  return result;
}

int main(int argc, char* argv[])
{
  SUNContext sunctx  = NULL;
  N_Vector y         = NULL;
  SUNMatrix A        = NULL;
  SUNLinearSolver LS = NULL;
  void* firkode_mem  = NULL;

  int flag             = 0;
  int firk_flag        = 0;
  int i                = 0;
  sunrealtype tout     = SUN_RCONST(0.10);
  sunrealtype dt_tout  = SUN_RCONST(0.25);
  sunrealtype tstop    = SUN_RCONST(0.30);
  sunrealtype dt_tstop = SUN_RCONST(0.30);
  sunrealtype tret     = ZERO;
  sunrealtype tcur     = ZERO;

  flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (flag)
  {
    fprintf(stderr, "SUNContext_Create returned %i\n", flag);
    return 1;
  }

  y = N_VNew_Serial(1, sunctx);
  if (!y) { return 1; }
  N_VConst(ZERO, y);

  firkode_mem = FIRKodeCreate(sunctx);
  if (!firkode_mem) { return 1; }

  flag = FIRKodeInit(firkode_mem, ode_rhs, ZERO, y);
  if (flag) { return 1; }

  flag = FIRKodeSStolerances(firkode_mem, SUN_RCONST(1.0e-4), SUN_RCONST(1.0e-8));
  if (flag) { return 1; }

  A = SUNDenseMatrix(1, 1, sunctx);
  if (!A) { return 1; }

  LS = SUNLinSol_Dense(y, A, sunctx);
  if (!LS) { return 1; }

  flag = FIRKodeSetLinearSolver(firkode_mem, LS, A);
  if (flag) { return 1; }

  flag = FIRKodeSetJacFn(firkode_mem, ode_jac);
  if (flag) { return 1; }

  /* use a modest maximum step so that tout and tstop are approached
     separately in the first returns */
  flag = FIRKodeSetMaxStep(firkode_mem, SUN_RCONST(0.05));
  if (flag) { return 1; }

  flag = FIRKodeSetStopTime(firkode_mem, tstop);
  if (flag) { return 1; }

  printf("0: tout = %" GSYM ", tstop = %" GSYM ", tret = %" GSYM
         ", tcur = %" GSYM "\n",
         tout, tstop, tret, tcur);

  for (i = 1; i <= 6; i++)
  {
    firk_flag = FIRKodeEvolve(firkode_mem, tout, y, &tret, FIRK_NORMAL);
    if (firk_flag < 0)
    {
      flag = 1;
      break;
    }

    flag = FIRKodeGetCurrentTime(firkode_mem, &tcur);
    if (flag) { break; }

    printf("%i: tout = %" GSYM ", tstop = %" GSYM ", tret = %" GSYM
           ", tcur = %" GSYM ", return = %i\n",
           i, tout, tstop, tret, tcur, firk_flag);

    /* First return: output time < stop time */
    if (i == 1 && firk_flag != FIRK_SUCCESS)
    {
      printf("ERROR: Expected output return!\n");
      flag = 1;
      break;
    }

    /* Second return: output time > stop time */
    if (i == 2)
    {
      if (firk_flag != FIRK_TSTOP_RETURN)
      {
        printf("ERROR: Expected stop return!\n");
        flag = 1;
        break;
      }
      tstop += dt_tstop;
      flag = FIRKodeSetStopTime(firkode_mem, tstop);
      if (flag) { break; }
    }

    /* Third return: output time = stop time */
    if (i == 3)
    {
      if (firk_flag != FIRK_TSTOP_RETURN)
      {
        printf("ERROR: Expected stop return!\n");
        flag = 1;
        break;
      }
      tstop += dt_tstop;
      flag = FIRKodeSetStopTime(firkode_mem, tstop);
      if (flag) { break; }
    }

    /* Fourth return: output time < stop time */
    if (i == 4)
    {
      if (firk_flag != FIRK_SUCCESS)
      {
        printf("ERROR: Expected output return!\n");
        flag = 1;
        break;
      }
    }

    /* Fifth return: output time > stop time */
    if (i == 5)
    {
      if (firk_flag != FIRK_TSTOP_RETURN)
      {
        printf("ERROR: Expected stop return!\n");
        flag = 1;
        break;
      }
    }

    /* Sixth return: output time < stop time (not updated) */
    if (i == 6 && firk_flag != FIRK_SUCCESS)
    {
      printf("ERROR: Expected output return!\n");
      flag = 1;
      break;
    }

    tout += dt_tout;
  }

  FIRKodeFree(&firkode_mem);
  N_VDestroy(y);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);

  if (!flag) { flag = test_fixed_step(sunctx); }
  if (!flag) { flag = test_clamp(sunctx); }
  if (!flag) { flag = test_backward(sunctx); }

  SUNContext_Free(&sunctx);

  if (!flag) { printf("SUCCESS\n"); }

  return flag;
}

/*---- end of file ----*/
