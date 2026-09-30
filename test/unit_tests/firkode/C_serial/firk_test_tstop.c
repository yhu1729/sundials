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
 * adaptive steps and with a fixed step size
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
  if (flag != FIRK_SUCCESS || SUNRabs(tret - tout) > tol ||
      SUNRabs(hlast - hfixed) > tol || nst != 11 ||
      SUNRabs(N_VGetArrayPointer(y)[0] - tout) > tol)
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

  SUNContext_Free(&sunctx);

  if (!flag) { printf("SUCCESS\n"); }

  return flag;
}

/*---- end of file ----*/
