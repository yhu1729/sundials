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
 * Unit test for the selection of the number of stages (design decision D7):
 *
 * - the internal mapping firkOrderToStages gives the smallest s with
 *   2s - 1 >= ord, the default for ord <= 0, and -1 above 2 FIRK_MAX_STAGES - 1
 *   (including INT_MAX, which must not overflow);
 * - in EXPERT mode FIRKodeSetOrder accepts only the orders that select s = 3
 *   (4, 5 and ord <= 0) and FIRKodeSetNumStages only s = 3; every other value
 *   returns FIRK_ILL_INPUT and leaves the number of stages unchanged;
 * - in BEGINNER mode both return FIRK_ILL_INPUT;
 * - the internal firkSetNumStagesAny accepts 1 <= s <= FIRK_MAX_STAGES.
 *
 * Pass criterion: every return value and number of stages matches.
 * ---------------------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>

#include "firkode/firkode.h"
#include "firkode/firkode_impl.h"
#include "sundials/sundials_context.h"

static int check(int ok, const char* what, int value)
{
  if (!ok)
  {
    printf("  FAIL: %s (value %d)\n", what, value);
    return 1;
  }
  return 0;
}

int main(void)
{
  SUNContext sunctx = NULL;
  void* mem         = NULL;
  int i, s = 0, fails = 0;

  /* mapping of orders to stage counts */
  const int ords[] = {-1, 0, 1, 2, 3, 4, 5, 16, 17, 18, INT_MAX - 1, INT_MAX};
  const int stages[] =
    {FIRK_DEFAULT_STAGES, FIRK_DEFAULT_STAGES, 1,  2,  2, 3, 3,
     FIRK_MAX_STAGES,     FIRK_MAX_STAGES,     -1, -1, -1};

  /* public values that are accepted and rejected */
  const int good_ords[] = {-1, 0, 4, 5};
  const int bad_ords[]  = {1, 2, 3, 6, 17, 18, INT_MAX};
  const int bad_stages[] =
    {-1, 0, 1, 2, 4, FIRK_MAX_STAGES, FIRK_MAX_STAGES + 1};

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  for (i = 0; i < (int)(sizeof(ords) / sizeof(ords[0])); i++)
  {
    fails += check(firkOrderToStages(ords[i]) == stages[i], "firkOrderToStages",
                   ords[i]);
  }

  mem = FIRKodeCreate(sunctx);
  if (mem == NULL) { return 1; }

  /* BEGINNER mode rejects both setters */
  fails += check(FIRKodeSetOrder(mem, 5) == FIRK_ILL_INPUT,
                 "SetOrder in BEGINNER mode", 5);
  fails += check(FIRKodeSetNumStages(mem, 3) == FIRK_ILL_INPUT,
                 "SetNumStages in BEGINNER mode", 3);

  /* EXPERT mode accepts s = 3 only */
  if (FIRKodeSetMode(mem, FIRK_EXPERT) != FIRK_SUCCESS) { return 1; }
  for (i = 0; i < (int)(sizeof(good_ords) / sizeof(good_ords[0])); i++)
  {
    fails += check(FIRKodeSetOrder(mem, good_ords[i]) == FIRK_SUCCESS &&
                     FIRKodeGetNumStages(mem, &s) == FIRK_SUCCESS && s == 3,
                   "SetOrder accepts an order of s = 3", good_ords[i]);
  }
  fails += check(FIRKodeSetNumStages(mem, 3) == FIRK_SUCCESS,
                 "SetNumStages accepts s = 3", 3);
  for (i = 0; i < (int)(sizeof(bad_ords) / sizeof(bad_ords[0])); i++)
  {
    fails += check(FIRKodeSetOrder(mem, bad_ords[i]) == FIRK_ILL_INPUT &&
                     FIRKodeGetNumStages(mem, &s) == FIRK_SUCCESS && s == 3,
                   "SetOrder rejects an order of s != 3", bad_ords[i]);
  }
  for (i = 0; i < (int)(sizeof(bad_stages) / sizeof(bad_stages[0])); i++)
  {
    fails += check(FIRKodeSetNumStages(mem, bad_stages[i]) == FIRK_ILL_INPUT &&
                     FIRKodeGetNumStages(mem, &s) == FIRK_SUCCESS && s == 3,
                   "SetNumStages rejects s != 3", bad_stages[i]);
  }

  /* the internal setter keeps the general-s path available to tests */
  for (i = 1; i <= FIRK_MAX_STAGES; i++)
  {
    fails += check(firkSetNumStagesAny(mem, i) == FIRK_SUCCESS &&
                     FIRKodeGetNumStages(mem, &s) == FIRK_SUCCESS && s == i,
                   "firkSetNumStagesAny accepts 1 <= s <= FIRK_MAX_STAGES", i);
  }
  fails += check(firkSetNumStagesAny(mem, 0) == FIRK_ILL_INPUT,
                 "firkSetNumStagesAny rejects s = 0", 0);
  fails += check(firkSetNumStagesAny(mem, FIRK_MAX_STAGES + 1) == FIRK_ILL_INPUT,
                 "firkSetNumStagesAny rejects s > FIRK_MAX_STAGES",
                 FIRK_MAX_STAGES + 1);

  FIRKodeFree(&mem);
  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("All tests passed\n");
  return 0;
}
