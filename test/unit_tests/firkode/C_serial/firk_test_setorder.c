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
 * Unit test for FIRKodeSetOrder: the selected stage count is the smallest s
 * with 2s - 1 >= ord, a nonpositive order selects the default, and orders
 * above 2 FIRK_MAX_STAGES - 1 (including INT_MAX, which must not overflow)
 * are rejected without changing the stage count.
 * ---------------------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>

#include "firkode/firkode.h"
#include "firkode/firkode_impl.h"
#include "firkode/firkode_tables.h"
#include "sundials/sundials_context.h"

int main(void)
{
  SUNContext sunctx = NULL;
  void* firkode_mem = NULL;
  int flag, i, s = 0, fails = 0;

  const int ords[] =
    {-1, 0, 1, 2, 3, 4, 5, 2 * FIRK_MAX_STAGES - 2, 2 * FIRK_MAX_STAGES - 1};
  const int stages[] = {FIRK_DEFAULT_STAGES, FIRK_DEFAULT_STAGES, 1, 2, 2, 3, 3,
                        FIRK_MAX_STAGES,     FIRK_MAX_STAGES};
  const int bad[]    = {2 * FIRK_MAX_STAGES, INT_MAX - 1, INT_MAX};

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }
  firkode_mem = FIRKodeCreate(sunctx);
  if (!firkode_mem) { return 1; }

  for (i = 0; i < (int)(sizeof(ords) / sizeof(ords[0])); i++)
  {
    flag = FIRKodeSetOrder(firkode_mem, ords[i]);
    if (flag == FIRK_SUCCESS) { flag = FIRKodeGetNumStages(firkode_mem, &s); }
    if (flag != FIRK_SUCCESS || s != stages[i])
    {
      printf("FAIL: order %d gives flag %d, stages %d (expected %d)\n", ords[i],
             flag, s, stages[i]);
      fails++;
    }
  }

  for (i = 0; i < (int)(sizeof(bad) / sizeof(bad[0])); i++)
  {
    if (FIRKodeSetOrder(firkode_mem, 5) != FIRK_SUCCESS) { return 1; }
    flag = FIRKodeSetOrder(firkode_mem, bad[i]);
    if (FIRKodeGetNumStages(firkode_mem, &s) != FIRK_SUCCESS) { return 1; }
    if (flag != FIRK_ILL_INPUT || s != 3)
    {
      printf("FAIL: order %d gives flag %d, stages %d (expected rejection)\n",
             bad[i], flag, s);
      fails++;
    }
  }

  FIRKodeFree(&firkode_mem);
  SUNContext_Free(&sunctx);

  if (fails) { printf("%d failures\n", fails); }
  else { printf("SUCCESS\n"); }
  return fails ? 1 : 0;
}

/*---- end of file ----*/
