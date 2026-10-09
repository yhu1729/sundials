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
 * Mock objects shared by the FIRKODE unit tests.
 *
 * A recording SUNLinearSolver wraps another solver of any type and forwards
 * every operation the wrapped solver implements. For each solve it records
 * whether SUNLinSolSetZeroGuess(.., SUNTRUE) was called since the last solve
 * (the flag resets after a solve, as in the SUNDIALS solvers), the tolerance
 * delta, and the iterations of the wrapped solver, and whether a solve is
 * running (for callbacks of the wrapped solver); it can replace the return
 * value of a range of solves with an injected one.
 * ---------------------------------------------------------------------------*/

#ifndef _FIRK_TEST_MOCKS_H
#define _FIRK_TEST_MOCKS_H

#include <sundials/sundials_linearsolver.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  SUNLinearSolver inner;     /* the wrapped solver (owned)                 */
  long int nsolves;          /* number of solves                           */
  long int nsolves_no_zero;  /* solves without a zero guess set before it  */
  sunbooleantype zero_guess; /* the zero-guess flag since the last solve   */
  sunrealtype last_delta;    /* tolerance of the last solve                */
  long int niters;           /* iterations of the wrapped solver, summed   */
  sunbooleantype in_solve;   /* is a solve of the wrapped solver running?  */
  long int inject_first;     /* first solve (1-based) to override, 0 none  */
  long int inject_last;      /* last solve to override                     */
  int inject_value;          /* the value returned by overridden solves    */
} FirkMockLSContent;

/* Wraps inner (taking ownership); returns NULL on failure */
SUNLinearSolver FirkMockLS_New(SUNLinearSolver inner, SUNContext sunctx);

/* The recording content of a mock solver */
FirkMockLSContent* FirkMockLS_Content(SUNLinearSolver LS);

#ifdef __cplusplus
}
#endif

#endif
