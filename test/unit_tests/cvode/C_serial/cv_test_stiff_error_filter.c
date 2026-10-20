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
 * Unit test for the Shampine-Baca stiff error filter (CVodeSetStiffErrorFilter)
 *
 * The test problem is the Prothero-Robinson equation
 *
 *   y' = lambda (y - cos(t)) - sin(t),  y(0) = 1,  y(t) = cos(t),
 *
 * with lambda = -1e6. For this problem the standard BDF local error estimate
 * C' (y_n - y_n(0)) ~ C h^{q+1} cos^{(q+1)}(t) does not see the stiff damping
 * of the implicit formula, whose actual local error in the stiff component is
 * smaller by the factor 1 / (1 - gamma lambda). The filtered estimate recovers
 * this factor, so the filtered integration should take far fewer steps while
 * still meeting the requested tolerances.
 *
 * The test checks:
 *   1. input validation of CVodeSetStiffErrorFilter and of the first CVode
 *      call (BDF only, after CVodeInit, Newton-type nonlinear solver with a
 *      linear solver),
 *   2. the solution error at the final time with and without the filter is
 *      within a small multiple of the requested tolerance,
 *   3. the filtered integration takes fewer steps than the unfiltered one.
 *
 * A stop time is used so that the final solution is a step-end value. The
 * filter controls the error of the implicit formula at step ends only;
 * interpolated (dense) output over the large steps the filter permits is not
 * stiffly damped and can be much less accurate.
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cvode/cvode.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_matrix.h"
#include "sundials/sundials_nvector.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM   "Lg"
#define COS(x) cosl(x)
#define SIN(x) sinl(x)
#elif defined(SUNDIALS_DOUBLE_PRECISION)
#define GSYM   "g"
#define COS(x) cos(x)
#define SIN(x) sin(x)
#else
#define GSYM   "g"
#define COS(x) cosf(x)
#define SIN(x) sinf(x)
#endif

#define ZERO   SUN_RCONST(0.0)
#define ONE    SUN_RCONST(1.0)
#define LAMBDA SUN_RCONST(-1.0e6)
#define TFINAL SUN_RCONST(10.0)
#define RTOL   SUN_RCONST(1.0e-6)
#define ATOL   SUN_RCONST(1.0e-10)

/* Allowed solution error relative to the tolerance rtol * |y| + atol. The
   local error test bounds each step's error by the tolerance, and for this
   stiff problem the global error is dominated by the last local errors. */
#define ERR_FACTOR SUN_RCONST(10.0)

static int ode_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* y_data    = N_VGetArrayPointer(y);
  sunrealtype* ydot_data = N_VGetArrayPointer(ydot);
  ydot_data[0]           = LAMBDA * (y_data[0] - COS(t)) - SIN(t);
  return 0;
}

static int ode_jac(sunrealtype t, N_Vector y, N_Vector f, SUNMatrix J,
                   void* user_data, N_Vector tempv1, N_Vector tempv2,
                   N_Vector tempv3)
{
  sunrealtype* J_data = SUNDenseMatrix_Data(J);
  J_data[0]           = LAMBDA;
  return 0;
}

/* Integrate the test problem with or without the filter. Returns 0 on
   success and sets the number of steps taken and the final solution error. */
static int integrate(SUNContext sunctx, sunbooleantype filter, long int* nst,
                     sunrealtype* err)
{
  N_Vector y         = NULL;
  SUNMatrix A        = NULL;
  SUNLinearSolver LS = NULL;
  void* cvode_mem    = NULL;
  sunrealtype tret   = ZERO;
  long int netf      = 0;
  int flag           = 0;
  int retval         = 1;

  y = N_VNew_Serial(1, sunctx);
  if (!y) { goto cleanup; }
  N_VConst(ONE, y);

  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (!cvode_mem) { goto cleanup; }

  flag = CVodeInit(cvode_mem, ode_rhs, ZERO, y);
  if (flag) { goto cleanup; }

  flag = CVodeSStolerances(cvode_mem, RTOL, ATOL);
  if (flag) { goto cleanup; }

  A = SUNDenseMatrix(1, 1, sunctx);
  if (!A) { goto cleanup; }

  LS = SUNLinSol_Dense(y, A, sunctx);
  if (!LS) { goto cleanup; }

  flag = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (flag) { goto cleanup; }

  flag = CVodeSetJacFn(cvode_mem, ode_jac);
  if (flag) { goto cleanup; }

  flag = CVodeSetStiffErrorFilter(cvode_mem, filter);
  if (flag)
  {
    fprintf(stderr, "CVodeSetStiffErrorFilter returned %i\n", flag);
    goto cleanup;
  }

  flag = CVodeSetMaxNumSteps(cvode_mem, 100000);
  if (flag) { goto cleanup; }

  flag = CVodeSetStopTime(cvode_mem, TFINAL);
  if (flag) { goto cleanup; }

  flag = CVode(cvode_mem, TFINAL, y, &tret, CV_NORMAL);
  if (flag < 0)
  {
    fprintf(stderr, "CVode returned %i\n", flag);
    goto cleanup;
  }

  flag = CVodeGetNumSteps(cvode_mem, nst);
  if (flag) { goto cleanup; }

  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  if (flag) { goto cleanup; }

  *err = SUNRabs(N_VGetArrayPointer(y)[0] - COS(tret));

  printf("filter = %i: steps = %li, error test fails = %li, y(%" GSYM
         ") = %" GSYM ", error = %" GSYM "\n",
         (int)filter, *nst, netf, tret, N_VGetArrayPointer(y)[0], *err);

  retval = 0;

cleanup:
  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);

  return retval;
}

/* Check that invalid configurations are rejected */
static int check_inputs(SUNContext sunctx)
{
  N_Vector y             = NULL;
  SUNNonlinearSolver NLS = NULL;
  void* cvode_mem        = NULL;
  sunrealtype tret       = ZERO;
  int flag               = 0;
  int retval             = 1;

  y = N_VNew_Serial(1, sunctx);
  if (!y) { goto cleanup; }
  N_VConst(ONE, y);

  /* Adams method is rejected */
  cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
  if (!cvode_mem) { goto cleanup; }

  flag = CVodeInit(cvode_mem, ode_rhs, ZERO, y);
  if (flag) { goto cleanup; }

  flag = CVodeSetStiffErrorFilter(cvode_mem, SUNTRUE);
  if (flag != CV_ILL_INPUT)
  {
    fprintf(stderr, "Expected CV_ILL_INPUT with CV_ADAMS, got %i\n", flag);
    goto cleanup;
  }
  CVodeFree(&cvode_mem);

  /* Enabling before CVodeInit is rejected, disabling is not */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (!cvode_mem) { goto cleanup; }

  flag = CVodeSetStiffErrorFilter(cvode_mem, SUNTRUE);
  if (flag != CV_NO_MALLOC)
  {
    fprintf(stderr, "Expected CV_NO_MALLOC before CVodeInit, got %i\n", flag);
    goto cleanup;
  }

  flag = CVodeSetStiffErrorFilter(cvode_mem, SUNFALSE);
  if (flag != CV_SUCCESS)
  {
    fprintf(stderr, "Expected CV_SUCCESS when disabling, got %i\n", flag);
    goto cleanup;
  }

  /* A fixed-point nonlinear solver (no linear solver) is rejected on the
     first call to CVode */
  flag = CVodeInit(cvode_mem, ode_rhs, ZERO, y);
  if (flag) { goto cleanup; }

  flag = CVodeSStolerances(cvode_mem, RTOL, ATOL);
  if (flag) { goto cleanup; }

  NLS = SUNNonlinSol_FixedPoint(y, 0, sunctx);
  if (!NLS) { goto cleanup; }

  flag = CVodeSetNonlinearSolver(cvode_mem, NLS);
  if (flag) { goto cleanup; }

  flag = CVodeSetStiffErrorFilter(cvode_mem, SUNTRUE);
  if (flag != CV_SUCCESS)
  {
    fprintf(stderr, "Expected CV_SUCCESS after CVodeInit, got %i\n", flag);
    goto cleanup;
  }

  flag = CVode(cvode_mem, ONE, y, &tret, CV_NORMAL);
  if (flag != CV_ILL_INPUT)
  {
    fprintf(stderr, "Expected CV_ILL_INPUT with fixed-point solver, got %i\n",
            flag);
    goto cleanup;
  }

  retval = 0;

cleanup:
  CVodeFree(&cvode_mem);
  SUNNonlinSolFree(NLS);
  N_VDestroy(y);

  return retval;
}

int main(int argc, char* argv[])
{
  SUNContext sunctx = NULL;
  long int nst_raw  = 0;
  long int nst_filt = 0;
  sunrealtype err   = ZERO;
  sunrealtype bound = ERR_FACTOR * (RTOL * SUNRabs(COS(TFINAL)) + ATOL);
  int flag          = 0;

  flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (flag)
  {
    fprintf(stderr, "SUNContext_Create returned %i\n", flag);
    return 1;
  }

  if (check_inputs(sunctx))
  {
    fprintf(stderr, "FAILED: input checks\n");
    SUNContext_Free(&sunctx);
    return 1;
  }

  if (integrate(sunctx, SUNFALSE, &nst_raw, &err))
  {
    fprintf(stderr, "FAILED: unfiltered integration\n");
    SUNContext_Free(&sunctx);
    return 1;
  }
  if (err > bound)
  {
    fprintf(stderr, "FAILED: unfiltered error %" GSYM " > %" GSYM "\n", err,
            bound);
    SUNContext_Free(&sunctx);
    return 1;
  }

  if (integrate(sunctx, SUNTRUE, &nst_filt, &err))
  {
    fprintf(stderr, "FAILED: filtered integration\n");
    SUNContext_Free(&sunctx);
    return 1;
  }
  if (err > bound)
  {
    fprintf(stderr, "FAILED: filtered error %" GSYM " > %" GSYM "\n", err, bound);
    SUNContext_Free(&sunctx);
    return 1;
  }

  if (nst_filt >= nst_raw)
  {
    fprintf(stderr, "FAILED: filtered steps %li >= unfiltered steps %li\n",
            nst_filt, nst_raw);
    SUNContext_Free(&sunctx);
    return 1;
  }

  SUNContext_Free(&sunctx);

  printf("SUCCESS\n");

  return 0;
}
