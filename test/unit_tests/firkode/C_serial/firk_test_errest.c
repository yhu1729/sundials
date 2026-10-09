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
 * Unit test for the error estimate on y' = -y (design milestone M2):
 *
 * - its order: the estimate of one accepted step of h = HORD / 2^k is
 *   O(h^(p+1)) with p = 3; the stages are converged far below the estimate
 *   (EXPERT nonlinear tolerance 1e-8), whose iteration error would otherwise
 *   swamp it at the loose tolerances that accept these steps;
 * - with an iterative block solver: after one step the filtered estimate
 *   computed with the matrix-free SPGMR block solver agrees with the one
 *   computed with the dense block solver, also when the estimate is small
 *   compared with the nonlinear solver tolerance.
 *
 * Pass criterion: each step is taken as requested, the slopes of the two
 * finest pairs are within 0.3 of 4, and the two estimates agree within 50%.
 * ---------------------------------------------------------------------------*/

#include <math.h>
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

/* the estimate order: steps HORD / 2^k, k < NH, from y(0) = 1, whose
   estimates stay far above roundoff */
#define NH 5
#if defined(SUNDIALS_SINGLE_PRECISION)
#define HORD SUN_RCONST(0.8)
#else
#define HORD SUN_RCONST(0.4)
#endif

static const sunrealtype lambda = SUN_RCONST(-1.0);
static int tight_nls = 0; /* converge the stages to 1e-8 of the tolerance */

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

/* take one step of h0 at the tolerances rtol, atol and return the local
   error estimate and the step taken; iterative = 0 uses the dense block
   solver, 1 the matrix-free SPGMR block solver */
static int one_step(SUNContext sunctx, int iterative, sunrealtype h0,
                    sunrealtype rtol, sunrealtype atol, sunrealtype* ele_out,
                    sunrealtype* h_out)
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
  if (FIRKodeSStolerances(firkode_mem, rtol, atol)) { goto cleanup; }
  if (FIRKodeSetInitStep(firkode_mem, h0)) { goto cleanup; }
  /* the order measurement converges the stages far below the estimate,
     which their iteration error would swamp */
  if (tight_nls && (FIRKodeSetMode(firkode_mem, FIRK_EXPERT) ||
                    FIRKodeSetNonlinConvCoef(firkode_mem, SUN_RCONST(1.0e-8))))
  {
    goto cleanup;
  }

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
  if (FIRKodeGetLastStep(firkode_mem, h_out)) { goto cleanup; }
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
  sunrealtype e_dense = ZERO, e_iter = ZERO, h = ZERO;
  sunrealtype e[NH], hk[NH], slope;
  int fails = 0, k;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  /* the estimate order: the estimate of one step is O(h^(p+1)), p = 3 */
  tight_nls = 1;
  for (k = 0; k < NH; k++)
  {
    hk[k] = HORD / SUNRpowerI(SUN_RCONST(2.0), k);
    if (one_step(sunctx, 0, hk[k], SUN_RCONST(1.0e-2), SUN_RCONST(1.0e-2),
                 &e[k], &h) ||
        h != hk[k])
    {
      printf("FAIL: a step of " SUN_FORMAT_G " was not taken as such\n", hk[k]);
      fails++;
      break;
    }
    e[k] = SUNRabs(e[k]);
    printf("h = " SUN_FORMAT_G ": estimate = " SUN_FORMAT_E, hk[k], e[k]);
    if (k > 0)
    {
      slope = (sunrealtype)(log((double)(e[k - 1] / e[k])) / log(2.0));
      printf(", slope = %.3f", (double)slope);
      /* the two finest pairs are in the asymptotic range */
      if (k >= NH - 2 && SUNRabs(slope - SUN_RCONST(4.0)) > SUN_RCONST(0.3))
      {
        printf("  FAIL: the slope is not p + 1 = 4");
        fails++;
      }
    }
    printf("\n");
  }

  tight_nls = 0;
  if (one_step(sunctx, 0, H0, RTOL, ATOL, &e_dense, &h) ||
      one_step(sunctx, 1, H0, RTOL, ATOL, &e_iter, &h))
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
