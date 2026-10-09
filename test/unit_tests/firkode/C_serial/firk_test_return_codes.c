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
 * Unit test for the FIRKODE return codes (design Table 8.2, decision D10) and
 * constants (design Section 9.1): every code has its value,
 * FIRKodeGetReturnFlagName and FIRKodeGetLinReturnFlagName map it to its name,
 * an unknown value maps to "NONE" (the FIRKLS function names integrator codes
 * such as FIRK_ILL_INPUT), and every constant has its value.
 *
 * Pass criterion: every value and name matches.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"

struct code
{
  long int value;
  long int expected;
  const char* name;
};

#define CODE(c, v) \
  {                \
    c, v, #c       \
  }

static const struct code firk_codes[] = {
  CODE(FIRK_SUCCESS, 0),
  CODE(FIRK_TSTOP_RETURN, 1),
  CODE(FIRK_ROOT_RETURN, 2),
  CODE(FIRK_WARNING, 99),
  CODE(FIRK_TOO_MUCH_WORK, -1),
  CODE(FIRK_TOO_MUCH_ACC, -2),
  CODE(FIRK_ERR_FAILURE, -3),
  CODE(FIRK_CONV_FAILURE, -4),
  CODE(FIRK_LINIT_FAIL, -5),
  CODE(FIRK_LSETUP_FAIL, -6),
  CODE(FIRK_LSOLVE_FAIL, -7),
  CODE(FIRK_RHSFUNC_FAIL, -8),
  CODE(FIRK_FIRST_RHSFUNC_ERR, -9),
  CODE(FIRK_REPTD_RHSFUNC_ERR, -10),
  CODE(FIRK_UNREC_RHSFUNC_ERR, -11),
  CODE(FIRK_RTFUNC_FAIL, -12),
  CODE(FIRK_MEM_FAIL, -20),
  CODE(FIRK_MEM_NULL, -21),
  CODE(FIRK_ILL_INPUT, -22),
  CODE(FIRK_NO_MALLOC, -23),
  CODE(FIRK_BAD_K, -24),
  CODE(FIRK_BAD_T, -25),
  CODE(FIRK_BAD_DKY, -26),
  CODE(FIRK_TOO_CLOSE, -27),
  CODE(FIRK_VECTOROP_ERR, -28),
  CODE(FIRK_CONTEXT_ERR, -32),
  CODE(FIRK_MASSINIT_FAIL, -40),
  CODE(FIRK_MASSSETUP_FAIL, -41),
  CODE(FIRK_MASSSOLVE_FAIL, -42),
  CODE(FIRK_MASSFREE_FAIL, -43),
  CODE(FIRK_MASSMULT_FAIL, -44),
  CODE(FIRK_INVALID_TABLE, -45),
  CODE(FIRK_TABLE_FAIL, -46),
  CODE(FIRK_CONTROLLER_ERR, -47),
  CODE(FIRK_STEPPER_ERR, -48),
  CODE(FIRK_ORDER_CHANGE_FAIL, -49),
  CODE(FIRK_UNRECOGNIZED_ERR, -99),
};

/* the FIRKLS codes equal the ARKLS values */
static const struct code firkls_codes[] = {
  CODE(FIRKLS_SUCCESS, 0),           CODE(FIRKLS_MEM_NULL, -1),
  CODE(FIRKLS_LMEM_NULL, -2),        CODE(FIRKLS_ILL_INPUT, -3),
  CODE(FIRKLS_MEM_FAIL, -4),         CODE(FIRKLS_MASSMEM_NULL, -6),
  CODE(FIRKLS_JACFUNC_UNRECVR, -7),  CODE(FIRKLS_JACFUNC_RECVR, -8),
  CODE(FIRKLS_MASSFUNC_UNRECVR, -9), CODE(FIRKLS_MASSFUNC_RECVR, -10),
  CODE(FIRKLS_SUNMAT_FAIL, -11),     CODE(FIRKLS_SUNLS_FAIL, -12),
};

static int check_codes(const struct code* codes, size_t n,
                       char* (*name_fn)(long int))
{
  int fails = 0;
  size_t i;
  char* name;

  for (i = 0; i < n; i++)
  {
    if (codes[i].value != codes[i].expected)
    {
      printf("  FAIL: %s = %ld, expected %ld\n", codes[i].name, codes[i].value,
             codes[i].expected);
      fails++;
    }
    name = name_fn(codes[i].value);
    if (name == NULL || strcmp(name, codes[i].name) != 0)
    {
      printf("  FAIL: name of %ld is %s, expected %s\n", codes[i].value,
             (name == NULL) ? "NULL" : name, codes[i].name);
      fails++;
    }
    free(name);
  }

  name = name_fn(12345);
  if (name == NULL || strcmp(name, "NONE") != 0)
  {
    printf("  FAIL: an unknown value is named %s\n",
           (name == NULL) ? "NULL" : name);
    fails++;
  }
  free(name);

  return fails;
}

/* the constants of design Section 9.1 */
static const struct code firk_constants[] = {
  CODE(FIRK_NORMAL, 1),
  CODE(FIRK_ONE_STEP, 2),
  CODE(FIRK_BEGINNER, 0),
  CODE(FIRK_EXPERT, 1),
  CODE(FIRK_PREDICT_EXTRAPOLATE, 0),
  CODE(FIRK_PREDICT_TRIVIAL, 1),
  CODE(FIRK_POLICY_RADAU5, 0),
  CODE(FIRK_POLICY_PERIODIC, 1),
  CODE(FIRK_MAX_STAGES, 9),
  CODE(FIRK_DEFAULT_STAGES, 3),
};

int main(void)
{
  int fails = 0;
  size_t i;

  for (i = 0; i < sizeof(firk_constants) / sizeof(firk_constants[0]); i++)
  {
    if (firk_constants[i].value != firk_constants[i].expected)
    {
      printf("  FAIL: %s = %ld, expected %ld\n", firk_constants[i].name,
             firk_constants[i].value, firk_constants[i].expected);
      fails++;
    }
  }

  fails += check_codes(firk_codes, sizeof(firk_codes) / sizeof(firk_codes[0]),
                       FIRKodeGetReturnFlagName);
  fails += check_codes(firkls_codes,
                       sizeof(firkls_codes) / sizeof(firkls_codes[0]),
                       FIRKodeGetLinReturnFlagName);

  /* FIRKLS functions can return FIRK_ILL_INPUT (EXPERT options in BEGINNER
     mode), which the FIRKLS name function names too */
  {
    char* name = FIRKodeGetLinReturnFlagName(FIRK_ILL_INPUT);
    if (name == NULL || strcmp(name, "FIRK_ILL_INPUT") != 0)
    {
      printf("  FAIL: FIRKodeGetLinReturnFlagName(FIRK_ILL_INPUT) is %s\n",
             (name == NULL) ? "NULL" : name);
      fails++;
    }
    free(name);
    /* a positive flag (an LU factorization failure) has no name */
    name = FIRKodeGetLinReturnFlagName(1);
    if (name == NULL || strcmp(name, "NONE") != 0)
    {
      printf("  FAIL: FIRKodeGetLinReturnFlagName(1) is %s\n",
             (name == NULL) ? "NULL" : name);
      fails++;
    }
    free(name);
  }

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("All tests passed\n");
  return 0;
}
