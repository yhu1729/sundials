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
 * Mock objects shared by the FIRKODE unit tests (see firk_test_mocks.h).
 * ---------------------------------------------------------------------------*/

#include "firk_test_mocks.h"

#include <stdlib.h>

#define MOCK(LS)  ((FirkMockLSContent*)((LS)->content))
#define INNER(LS) (MOCK(LS)->inner)

static SUNLinearSolver_Type mock_gettype(SUNLinearSolver LS)
{
  return SUNLinSolGetType(INNER(LS));
}

static SUNLinearSolver_ID mock_getid(SUNLinearSolver LS)
{
  return SUNLinSolGetID(INNER(LS));
}

static SUNErrCode mock_setatimes(SUNLinearSolver LS, void* A_data,
                                 SUNATimesFn ATimes)
{
  return SUNLinSolSetATimes(INNER(LS), A_data, ATimes);
}

static SUNErrCode mock_setpreconditioner(SUNLinearSolver LS, void* P_data,
                                         SUNPSetupFn Pset, SUNPSolveFn Psol)
{
  return SUNLinSolSetPreconditioner(INNER(LS), P_data, Pset, Psol);
}

static SUNErrCode mock_setscalingvectors(SUNLinearSolver LS, N_Vector s1,
                                         N_Vector s2)
{
  return SUNLinSolSetScalingVectors(INNER(LS), s1, s2);
}

static SUNErrCode mock_setzeroguess(SUNLinearSolver LS, sunbooleantype onoff)
{
  MOCK(LS)->zero_guess = onoff;
  if (INNER(LS)->ops->setzeroguess)
  {
    return SUNLinSolSetZeroGuess(INNER(LS), onoff);
  }
  return SUN_SUCCESS;
}

static SUNErrCode mock_initialize(SUNLinearSolver LS)
{
  return SUNLinSolInitialize(INNER(LS));
}

static int mock_setup(SUNLinearSolver LS, SUNMatrix A)
{
  return SUNLinSolSetup(INNER(LS), A);
}

static int mock_solve(SUNLinearSolver LS, SUNMatrix A, N_Vector x, N_Vector b,
                      sunrealtype tol)
{
  FirkMockLSContent* c = MOCK(LS);
  int retval;

  c->nsolves++;
  if (!c->zero_guess) { c->nsolves_no_zero++; }
  c->zero_guess = SUNFALSE;
  c->last_delta = tol;

  c->in_solve = SUNTRUE;
  retval      = SUNLinSolSolve(c->inner, A, x, b, tol);
  c->in_solve = SUNFALSE;
  if (c->inner->ops->numiters) { c->niters += SUNLinSolNumIters(c->inner); }
  if (c->inject_first > 0 && c->nsolves >= c->inject_first &&
      c->nsolves <= c->inject_last)
  {
    retval = c->inject_value;
  }
  return retval;
}

static int mock_numiters(SUNLinearSolver LS)
{
  return SUNLinSolNumIters(INNER(LS));
}

static sunrealtype mock_resnorm(SUNLinearSolver LS)
{
  return SUNLinSolResNorm(INNER(LS));
}

static sunindextype mock_lastflag(SUNLinearSolver LS)
{
  return SUNLinSolLastFlag(INNER(LS));
}

static N_Vector mock_resid(SUNLinearSolver LS)
{
  return SUNLinSolResid(INNER(LS));
}

static SUNErrCode mock_free(SUNLinearSolver LS)
{
  if (LS == NULL) { return SUN_SUCCESS; }
  if (LS->content != NULL)
  {
    SUNLinSolFree(INNER(LS));
    free(LS->content);
    LS->content = NULL;
  }
  SUNLinSolFreeEmpty(LS);
  return SUN_SUCCESS;
}

SUNLinearSolver FirkMockLS_New(SUNLinearSolver inner, SUNContext sunctx)
{
  SUNLinearSolver LS;
  FirkMockLSContent* c;

  if (inner == NULL) { return NULL; }
  LS = SUNLinSolNewEmpty(sunctx);
  if (LS == NULL) { return NULL; }
  c = (FirkMockLSContent*)calloc(1, sizeof(FirkMockLSContent));
  if (c == NULL)
  {
    SUNLinSolFreeEmpty(LS);
    return NULL;
  }
  c->inner    = inner;
  LS->content = c;

  /* forward exactly the operations of the wrapped solver */
  LS->ops->gettype = mock_gettype;
  if (inner->ops->getid) { LS->ops->getid = mock_getid; }
  if (inner->ops->setatimes) { LS->ops->setatimes = mock_setatimes; }
  if (inner->ops->setpreconditioner)
  {
    LS->ops->setpreconditioner = mock_setpreconditioner;
  }
  if (inner->ops->setscalingvectors)
  {
    LS->ops->setscalingvectors = mock_setscalingvectors;
  }
  LS->ops->setzeroguess = mock_setzeroguess;
  if (inner->ops->initialize) { LS->ops->initialize = mock_initialize; }
  if (inner->ops->setup) { LS->ops->setup = mock_setup; }
  LS->ops->solve = mock_solve;
  if (inner->ops->numiters) { LS->ops->numiters = mock_numiters; }
  if (inner->ops->resnorm) { LS->ops->resnorm = mock_resnorm; }
  if (inner->ops->lastflag) { LS->ops->lastflag = mock_lastflag; }
  if (inner->ops->resid) { LS->ops->resid = mock_resid; }
  LS->ops->free = mock_free;
  return LS;
}

FirkMockLSContent* FirkMockLS_Content(SUNLinearSolver LS) { return MOCK(LS); }
