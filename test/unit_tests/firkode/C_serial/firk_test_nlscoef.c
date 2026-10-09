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
 * Unit test of the automatic Newton stopping tolerance (design Section 4.4,
 * Table 5.9): by default eps_nls is RADAU5's
 *   max(10 u / reltol, min(0.03, sqrt(reltol))),
 * evaluated from the relative tolerance whenever the tolerances change, and
 * 0.03 before the tolerances are set, with a zero relative tolerance or with a
 * user weight function; the EXPERT setter fixes a positive value, which the
 * tolerances no longer change, and a non-positive value or a return to BEGINNER
 * mode restores the rule; the option key "nonlin_conv_coef" does the same;
 * FIRKodeWriteParameters marks the automatic value.
 *
 * Pass criterion: every effective value, return value and marker matches.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "firkode/firkode.h"
#include "firkode/firkode_impl.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)t;
  (void)user_data;
  N_VScale(-ONE, y, ydot);
  return 0;
}

static int ewt(N_Vector y, N_Vector w, void* user_data)
{
  (void)y;
  (void)user_data;
  N_VConst(ONE, w);
  return 0;
}

static int check(int ok, const char* what)
{
  if (!ok)
  {
    printf("  FAIL: %s\n", what);
    return 1;
  }
  return 0;
}

/* RADAU5's rule, as the test computes it */
static sunrealtype rule(sunrealtype reltol)
{
  return SUNMAX(NLSCOEF_UROUND_FAC * SUN_UNIT_ROUNDOFF / reltol,
                SUNMIN(NLSCOEF_MAX, SUNRsqrt(reltol)));
}

static int same(sunrealtype a, sunrealtype b)
{
  return SUNRabs(a - b) <= SUN_RCONST(10.0) * SUN_UNIT_ROUNDOFF * SUNRabs(b);
}

/* the automatic value in use is v and the user value is unset */
static int automatic(FIRKodeMem fm, sunrealtype v)
{
  return fm->nlscoef_user <= ZERO && same(fm->nlscoef, v);
}

/* the fixed value v is in use */
static int fixed(FIRKodeMem fm, sunrealtype v)
{
  return fm->nlscoef_user == v && fm->nlscoef == v;
}

/* 1 if FIRKodeWriteParameters marks the coefficient as automatic, 0 if it
   prints the coefficient without the marker, -1 otherwise */
static int marker(void* mem)
{
  FILE* fp = tmpfile();
  char line[256];
  int found = -1;

  if (fp == NULL || FIRKodeWriteParameters(mem, fp) != FIRK_SUCCESS)
  {
    if (fp) { fclose(fp); }
    return -1;
  }
  rewind(fp);
  while (fgets(line, (int)sizeof(line), fp) != NULL)
  {
    if (strstr(line, "Nonlinear convergence coefficient = ") != NULL)
    {
      found = (strstr(line, "(automatic)") != NULL) ? 1 : 0;
    }
  }
  fclose(fp);
  return found;
}

static int set_option(void* mem, const char* key, const char* value)
{
  char* argv[3];
  int retval;

  argv[0] = (char*)malloc(strlen("prog") + 1);
  argv[1] = (char*)malloc(strlen(key) + 1);
  argv[2] = (char*)malloc(strlen(value) + 1);
  if (argv[0] == NULL || argv[1] == NULL || argv[2] == NULL) { return -999; }
  strcpy(argv[0], "prog");
  strcpy(argv[1], key);
  strcpy(argv[2], value);
  retval = FIRKodeSetOptions(mem, NULL, NULL, 3, argv);
  free(argv[0]);
  free(argv[1]);
  free(argv[2]);
  return retval;
}

int main(void)
{
  SUNContext sunctx  = NULL;
  void* mem          = NULL;
  N_Vector y         = NULL;
  N_Vector abstol    = NULL;
  SUNMatrix A        = NULL;
  SUNLinearSolver LS = NULL;
  FIRKodeMem fm;
  const sunrealtype atol       = SUN_RCONST(1.0e-10);
  const sunrealtype reltols[3] = {SUN_RCONST(1.0e-4), SUN_RCONST(1.0e-6),
                                  SUN_RCONST(1.0e-12)};
  char what[128];
  int i, fails = 0;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }
  y = N_VNew_Serial(1, sunctx);
  N_VConst(ONE, y);
  abstol = N_VClone(y);
  N_VConst(atol, abstol);
  A  = SUNDenseMatrix(1, 1, sunctx);
  LS = SUNLinSol_Dense(y, A, sunctx);

  printf("Automatic Newton tolerance: u = " SUN_FORMAT_G ", cap = " SUN_FORMAT_G
         "\n",
         (sunrealtype)SUN_UNIT_ROUNDOFF, (sunrealtype)NLSCOEF_MAX);

  mem = FIRKodeCreate(sunctx);
  if (mem == NULL) { return 1; }
  fm = (FIRKodeMem)mem;

  /* before any tolerance: the cap */
  fails += check(automatic(fm, NLSCOEF_MAX), "the cap 0.03 after Create");

  if (FIRKodeInit(mem, f, ZERO, y) || FIRKodeSetLinearSolver(mem, LS, A))
  {
    printf("  FAIL: setup failed\n");
    return 1;
  }
  fails += check(automatic(fm, NLSCOEF_MAX), "the cap 0.03 after Init");

  /* the rule follows every scalar relative tolerance */
  for (i = 0; i < 3; i++)
  {
    fails += check(FIRKodeSStolerances(mem, reltols[i], atol) == FIRK_SUCCESS,
                   "SStolerances");
    snprintf(what, sizeof(what), "the rule at reltol = %g", (double)reltols[i]);
    fails += check(automatic(fm, rule(reltols[i])), what);
  }
  /* the roundoff floor dominates the last one in double precision */
  fails += check(rule(reltols[2]) >=
                   NLSCOEF_UROUND_FAC * SUN_UNIT_ROUNDOFF / reltols[2],
                 "the floor 10 u / reltol is part of the rule");

  /* vector absolute tolerances: the same rule */
  fails += check(FIRKodeSVtolerances(mem, reltols[1], abstol) == FIRK_SUCCESS &&
                   automatic(fm, rule(reltols[1])),
                 "the rule after SVtolerances");

  /* no relative tolerance: the cap */
  fails += check(FIRKodeSStolerances(mem, ZERO, atol) == FIRK_SUCCESS &&
                   automatic(fm, NLSCOEF_MAX),
                 "the cap for reltol = 0");
  fails += check(FIRKodeWFtolerances(mem, ewt) == FIRK_SUCCESS &&
                   automatic(fm, NLSCOEF_MAX),
                 "the cap with a user weight function");
  fails += check(marker(mem) == 1, "WriteParameters marks the automatic value");

  /* the setter is an EXPERT option */
  fails +=
    check(FIRKodeSetNonlinConvCoef(mem, SUN_RCONST(0.05)) == FIRK_ILL_INPUT &&
            automatic(fm, NLSCOEF_MAX),
          "SetNonlinConvCoef is rejected in BEGINNER mode");

  if (FIRKodeSetMode(mem, FIRK_EXPERT) ||
      FIRKodeSStolerances(mem, reltols[1], atol))
  {
    printf("  FAIL: EXPERT setup failed\n");
    return 1;
  }
  fails += check(automatic(fm, rule(reltols[1])),
                 "EXPERT mode keeps the automatic value");

  /* a fixed value, which the tolerances no longer change */
  fails += check(FIRKodeSetNonlinConvCoef(mem, SUN_RCONST(0.05)) == FIRK_SUCCESS &&
                   fixed(fm, SUN_RCONST(0.05)),
                 "SetNonlinConvCoef(0.05) fixes the value");
  fails += check(FIRKodeSStolerances(mem, reltols[0], atol) == FIRK_SUCCESS &&
                   fixed(fm, SUN_RCONST(0.05)),
                 "new tolerances keep the fixed value");
  fails += check(marker(mem) == 0, "WriteParameters prints the fixed value");

  /* a non-positive value restores the rule */
  fails += check(FIRKodeSetNonlinConvCoef(mem, ZERO) == FIRK_SUCCESS &&
                   automatic(fm, rule(reltols[0])),
                 "SetNonlinConvCoef(0) restores the rule");
  fails += check(FIRKodeSetNonlinConvCoef(mem, SUN_RCONST(0.05)) == FIRK_SUCCESS &&
                   FIRKodeSetNonlinConvCoef(mem, -ONE) == FIRK_SUCCESS &&
                   automatic(fm, rule(reltols[0])),
                 "SetNonlinConvCoef(-1) restores the rule");

  /* so does a return to BEGINNER mode */
  fails += check(FIRKodeSetNonlinConvCoef(mem, SUN_RCONST(0.05)) == FIRK_SUCCESS &&
                   FIRKodeSetMode(mem, FIRK_BEGINNER) == FIRK_SUCCESS &&
                   automatic(fm, rule(reltols[0])),
                 "SetMode(BEGINNER) restores the rule");

  /* the option key */
  if (FIRKodeSetMode(mem, FIRK_EXPERT))
  {
    printf("  FAIL: SetMode(EXPERT) failed\n");
    return 1;
  }
  fails +=
    check(set_option(mem, "firkode.nonlin_conv_coef", "0.02") == FIRK_SUCCESS &&
            fixed(fm, SUN_RCONST(0.02)),
          "the key nonlin_conv_coef fixes the value");
  fails += check(set_option(mem, "firkode.nonlin_conv_coef", "0") == FIRK_SUCCESS &&
                   automatic(fm, rule(reltols[0])),
                 "the key nonlin_conv_coef 0 restores the rule");
  fails += check(marker(mem) == 1, "WriteParameters marks the rule again");

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(abstol);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("\nAll tests passed\n");
  return 0;
}
