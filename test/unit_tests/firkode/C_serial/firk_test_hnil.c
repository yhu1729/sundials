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
 * Unit test of the t + h = t warnings (design hazard H32): y' = -y from
 * t0 = 1/u with fixed steps h = 1/4, so that t + h = t at every step. The
 * integration ends with FIRK_TOO_MUCH_WORK after 20 steps, and the warnings
 * go through the logger: mxhnil warnings and one notice that they stop, for
 * the default mxhnil = 10 and for FIRKodeSetMaxHnilWarns(3). Built without
 * warning records (SUNDIALS_LOGGING_LEVEL below the warning level), only the
 * return value is checked.
 *
 * Pass criterion: the return value and the warning counts match.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "sundials/priv/sundials_logger_macros.h" /* the logging levels */
#include "sundials/sundials_config.h"
#include "sundials/sundials_logger.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ONE SUN_RCONST(1.0)

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)t;
  (void)user_data;
  N_VScale(-ONE, y, ydot);
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

/* the lines of a file that contain a string */
static int count_lines(const char* path, const char* text)
{
  char line[2048];
  int n    = 0;
  FILE* fp = fopen(path, "r");
  if (fp == NULL) { return -1; }
  while (fgets(line, sizeof(line), fp) != NULL)
  {
    if (strstr(line, text) != NULL) { n++; }
  }
  fclose(fp);
  return n;
}

/* one integration with mxhnil (0: the default 10) and its warnings in path */
static int run(int mxhnil, const char* path)
{
  const sunrealtype t0 = ONE / SUN_UNIT_ROUNDOFF;
  int fails            = 0, flag, nwarn, ndone, expected;
  sunrealtype tret;
  SUNContext sunctx = NULL;
  SUNLogger logger  = NULL;
  N_Vector y;
  SUNMatrix A;
  SUNLinearSolver LS;
  void* mem;
  char what[160];
  FILE* fp;

  fp = fopen(path, "w"); /* a fresh file */
  if (fp == NULL) { return 1; }
  fclose(fp);
  if (SUNLogger_Create(SUN_COMM_NULL, 0, &logger) ||
      SUNLogger_SetWarningFilename(logger, path) ||
      SUNContext_Create(SUN_COMM_NULL, &sunctx) ||
      SUNContext_SetLogger(sunctx, logger))
  {
    return 1;
  }

  y = N_VNew_Serial(1, sunctx);
  N_VConst(ONE, y);
  A   = SUNDenseMatrix(1, 1, sunctx);
  LS  = SUNLinSol_Dense(y, A, sunctx);
  mem = FIRKodeCreate(sunctx);
  if (FIRKodeInit(mem, f, t0, y) ||
      FIRKodeSStolerances(mem, SUN_RCONST(1.0e-4), SUN_RCONST(1.0e-8)) ||
      FIRKodeSetLinearSolver(mem, LS, A) ||
      FIRKodeSetFixedStep(mem, SUN_RCONST(0.25)) ||
      FIRKodeSetMaxNumSteps(mem, 20) ||
      (mxhnil > 0 && FIRKodeSetMaxHnilWarns(mem, mxhnil)))
  {
    return 1;
  }

  flag = FIRKodeEvolve(mem, t0 + SUN_RCONST(1000.0) * t0 * SUN_UNIT_ROUNDOFF, y,
                       &tret, FIRK_NORMAL);
  snprintf(what, sizeof(what), "mxhnil %d: FIRK_TOO_MUCH_WORK (flag %d)",
           mxhnil, flag);
  fails += check(flag == FIRK_TOO_MUCH_WORK, what);

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);
  SUNLogger_Destroy(&logger);

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_WARNING
  expected = (mxhnil > 0) ? mxhnil : 10;
  nwarn    = count_lines(path, "t + h = t on the next step");
  ndone    = count_lines(path, "will not be issued again");
  snprintf(what, sizeof(what),
           "mxhnil %d: %d warnings (expected %d) and %d notice(s) in the log",
           mxhnil, nwarn, expected, ndone);
  fails += check(nwarn == expected && ndone == 1, what);
#else
  (void)nwarn;
  (void)ndone;
  (void)expected;
#endif
  return fails;
}

int main(int argc, char* argv[])
{
  const char* path = (argc > 1) ? argv[1] : "firk_test_hnil.log";
  int fails        = 0;

  fails += run(0, path);
  fails += run(3, path);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("All tests passed\n");
  return 0;
}
