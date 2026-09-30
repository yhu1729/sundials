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
 * Unit test for the error estimate with an iterative block linear solver.
 * After one step on y' = lambda y, the filtered local error estimate computed
 * with the matrix-free SPGMR block solver must agree with the one computed
 * with the dense block solver, including when the estimate is small compared
 * with the nonlinear solver tolerance.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* The estimate must lie well above roundoff but below the Newton tolerance of
   the block solve (0.005 in the weighted norm), where the old small right-hand
   side shortcut returned it as zero. */
#if defined(SUNDIALS_SINGLE_PRECISION)
#define RTOL SUN_RCONST(1.0e-2)
#define ATOL SUN_RCONST(1.0e-6)
#define H0   SUN_RCONST(0.1)
#else
#define RTOL SUN_RCONST(1.0e-6)
#define ATOL SUN_RCONST(1.0e-10)
#define H0   SUN_RCONST(0.01)
#endif

static const sunrealtype lambda = SUN_RCONST(-1.0);

static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VGetArrayPointer(ydot)[0] = lambda * N_VGetArrayPointer(y)[0];
  return 0;
}

static int jac(sunrealtype t, N_Vector y, N_Vector f, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  SUNDenseMatrix_Data(J)[0] = lambda;
  return 0;
}

/* take one step and return the local error estimate; iterative = 0 uses the
   dense block solver, 1 the matrix-free SPGMR block solver */
static int one_step(SUNContext sunctx, int iterative, sunrealtype* ele_out)
{
  N_Vector y         = NULL;
  N_Vector ele       = NULL;
  SUNMatrix A        = NULL;
  SUNLinearSolver LS = NULL;
  void* firkode_mem  = NULL;
  sunrealtype tret   = ZERO;
  int result         = 1;

  y   = N_VNew_Serial(1, sunctx);
  ele = N_VNew_Serial(1, sunctx);
  if (!y || !ele) { goto cleanup; }
  N_VConst(ONE, y);

  firkode_mem = FIRKodeCreate(sunctx);
  if (!firkode_mem) { goto cleanup; }
  if (FIRKodeInit(firkode_mem, rhs, ZERO, y)) { goto cleanup; }
  if (FIRKodeSStolerances(firkode_mem, RTOL, ATOL)) { goto cleanup; }
  if (FIRKodeSetInitStep(firkode_mem, H0)) { goto cleanup; }

  if (iterative)
  {
    LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 5, sunctx);
    if (!LS) { goto cleanup; }
    if (FIRKodeSetLinearSolver(firkode_mem, LS, NULL)) { goto cleanup; }
  }
  else
  {
    A = SUNDenseMatrix(1, 1, sunctx);
    if (!A) { goto cleanup; }
    LS = SUNLinSol_Dense(y, A, sunctx);
    if (!LS) { goto cleanup; }
    if (FIRKodeSetLinearSolver(firkode_mem, LS, A)) { goto cleanup; }
    if (FIRKodeSetJacFn(firkode_mem, jac)) { goto cleanup; }
  }

  if (FIRKodeEvolve(firkode_mem, ONE, y, &tret, FIRK_ONE_STEP) < 0)
  {
    goto cleanup;
  }
  if (FIRKodeGetEstLocalErrors(firkode_mem, ele)) { goto cleanup; }
  *ele_out = N_VGetArrayPointer(ele)[0];
  result   = 0;

cleanup:
  FIRKodeFree(&firkode_mem);
  N_VDestroy(y);
  N_VDestroy(ele);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  return result;
}

int main(void)
{
  SUNContext sunctx   = NULL;
  sunrealtype e_dense = ZERO, e_iter = ZERO;
  int fails = 0;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  if (one_step(sunctx, 0, &e_dense) || one_step(sunctx, 1, &e_iter))
  {
    printf("FAIL: integration failed\n");
    fails++;
  }
  else
  {
    printf("local error estimate: dense = " SUN_FORMAT_E
           ", SPGMR = " SUN_FORMAT_E "\n",
           e_dense, e_iter);
    if (e_dense == ZERO ||
        SUNRabs(e_iter - e_dense) > SUN_RCONST(0.5) * SUNRabs(e_dense))
    {
      printf("FAIL: the SPGMR estimate differs from the dense estimate\n");
      fails++;
    }
  }

  SUNContext_Free(&sunctx);

  if (!fails) { printf("SUCCESS\n"); }
  return fails ? 1 : 0;
}

/*---- end of file ----*/
