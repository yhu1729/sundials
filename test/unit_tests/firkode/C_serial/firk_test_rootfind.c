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
 * Unit test of the rootfinding (design hazard H22):
 *
 * - y' = -y, y(0) = 1, with fixed steps h = 1/2 and g = y - 1/2: the root
 *   ln 2 is located on the collocation polynomial that FIRKodeGetDky
 *   evaluates, so g at that polynomial changes sign between tret - ttol and
 *   tret, with ttol = 100 u (|t| + |h|) (white box); the linear interpolant
 *   of the step ends would put the root 0.03 away;
 * - the oscillator y1' = y2, y2' = -y1, y(0) = (1, 0), with g = y1 = cos t:
 *   without a direction the roots pi/2 (decreasing, rootsfound -1) and
 *   3 pi/2 (increasing, +1) are returned; FIRKodeSetRootDirection with +1
 *   returns 3 pi/2 only and with -1 pi/2 only;
 * - a root at t0: g = y2 = -sin t vanishes at t0 = 0; it is not returned, and
 *   the first root is pi (increasing);
 * - g identically zero: no root is returned, and one warning that a root
 *   function is still identically zero goes through the logger after the
 *   first step, none after FIRKodeSetNoInactiveRootWarn or for the root at
 *   t0. Built without warning records (SUNDIALS_LOGGING_LEVEL below the
 *   warning level), the warnings are not counted;
 * - a failure met by the rootfinding checks is reported by its own code: with
 *   fixed steps 0.1 and g = y1 where |y1| >= 0.1 and 0 elsewhere, g is
 *   exactly 0 at the step end t = 1.5, which is returned as a root; the next
 *   FIRKodeEvolve then forms y'(1.5) for the check of g beyond that root, and
 *   a recoverable failure of the right-hand side there returns
 *   FIRK_UNREC_RHSFUNC_ERR, an unrecoverable one FIRK_RHSFUNC_FAIL, with the
 *   right-hand side message and not the rootfinding one;
 * - with no root functions (FIRKodeRootInit with nrtfn = 0) the setters and
 *   getter of the rootfinding return FIRK_ILL_INPUT.
 *
 * Pass criterion: every return value, root time (to 1e-5 for ln 2, to 1e-8
 * for the oscillator; 1e-3 in single precision), rootsfound entry, error
 * message and warning count matches, and ttol and the sign change are as
 * stated.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>

#include "sundials/sundials_context.h"

#include "firkode/firkode.h"
#include "firkode/firkode_impl.h"
#include "firkode/firkode_ls.h"
#include "firkode/firkode_root_impl.h"
#include "nvector/nvector_serial.h"
#include "sundials/priv/sundials_logger_macros.h" /* the logging levels */
#include "sundials/sundials_config.h"
#include "sundials/sundials_logger.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ZERO SUN_RCONST(0.0)
#define HALF SUN_RCONST(0.5)
#define ONE  SUN_RCONST(1.0)
#define PI   SUN_RCONST(3.141592653589793238462643383279502884)
#define LN2  SUN_RCONST(0.693147180559945309417232121458176568)

#if defined(SUNDIALS_SINGLE_PRECISION)
#define TOL      SUN_RCONST(1.0e-5)
#define TOL_LN2  SUN_RCONST(1.0e-3)
#define TOL_ROOT SUN_RCONST(1.0e-3)
#else
#define TOL      SUN_RCONST(1.0e-10)
#define TOL_LN2  SUN_RCONST(1.0e-5)
#define TOL_ROOT SUN_RCONST(1.0e-8)
#endif

/* fifteen fixed steps of 0.1 reach t = 1.5 up to accumulated roundoff */
#define TOL_END (SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF)

/* y' = -y */
static int f_decay(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)t;
  (void)user_data;
  N_VScale(-ONE, y, ydot);
  return 0;
}

/* y1' = y2, y2' = -y1 */
/* failures injected into f_osc: 1 recoverable, 2 unrecoverable */
static int fail_mode = 0;

static int f_osc(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)t;
  (void)user_data;
  NV_Ith_S(ydot, 0) = NV_Ith_S(y, 1);
  NV_Ith_S(ydot, 1) = -NV_Ith_S(y, 0);
  if (fail_mode == 1) { return 1; }
  if (fail_mode == 2) { return -1; }
  return 0;
}

static int g_half(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)
{
  (void)t;
  (void)user_data;
  gout[0] = NV_Ith_S(y, 0) - HALF;
  return 0;
}

static int g_y1(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)
{
  (void)t;
  (void)user_data;
  gout[0] = NV_Ith_S(y, 0);
  return 0;
}

static int g_y2(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)
{
  (void)t;
  (void)user_data;
  gout[0] = NV_Ith_S(y, 1);
  return 0;
}

static int g_zero(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)
{
  (void)t;
  (void)y;
  (void)user_data;
  gout[0] = ZERO;
  return 0;
}

/* g = y1 where |y1| >= 1/10 and 0 elsewhere: exactly 0 at a step end */
static int g_dead(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)
{
  (void)t;
  (void)user_data;
  gout[0] = (SUNRabs(NV_Ith_S(y, 0)) < SUN_RCONST(0.1)) ? ZERO : NV_Ith_S(y, 0);
  return 0;
}

/* the last message of the error handler */
static char last_msg[512];

static void err_handler(int line, const char* func, const char* file,
                        const char* msg, SUNErrCode err_code,
                        void* err_user_data, SUNContext sunctx)
{
  (void)line;
  (void)func;
  (void)file;
  (void)err_code;
  (void)err_user_data;
  (void)sunctx;
  strncpy(last_msg, msg, sizeof(last_msg) - 1);
  last_msg[sizeof(last_msg) - 1] = '\0';
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

/* one integrator with its context, logger and solvers */
struct run
{
  SUNContext sunctx;
  SUNLogger logger;
  N_Vector y;
  SUNMatrix A;
  SUNLinearSolver LS;
  void* mem;
};

/* n = 1: the decay problem with fixed steps 1/2; n = 2: the oscillator */
static int setup(struct run* R, int n, FIRKRootFn g, const char* path)
{
  FILE* fp = fopen(path, "w"); /* a fresh file */
  if (fp == NULL) { return 1; }
  fclose(fp);
  memset(R, 0, sizeof(*R));
  if (SUNLogger_Create(SUN_COMM_NULL, 0, &R->logger) ||
      SUNLogger_SetWarningFilename(R->logger, path) ||
      SUNContext_Create(SUN_COMM_NULL, &R->sunctx) ||
      SUNContext_SetLogger(R->sunctx, R->logger))
  {
    return 1;
  }
  R->y = N_VNew_Serial(n, R->sunctx);
  N_VConst(ZERO, R->y);
  NV_Ith_S(R->y, 0) = ONE;
  R->A              = SUNDenseMatrix(n, n, R->sunctx);
  R->LS             = SUNLinSol_Dense(R->y, R->A, R->sunctx);
  R->mem            = FIRKodeCreate(R->sunctx);
  if (R->mem == NULL ||
      FIRKodeInit(R->mem, (n == 1) ? f_decay : f_osc, ZERO, R->y) ||
      FIRKodeSStolerances(R->mem, TOL, TOL) ||
      FIRKodeSetLinearSolver(R->mem, R->LS, R->A) ||
      FIRKodeRootInit(R->mem, 1, g))
  {
    return 1;
  }
  if (n == 1 && FIRKodeSetFixedStep(R->mem, HALF)) { return 1; }
  return 0;
}

static void teardown(struct run* R)
{
  FIRKodeFree(&R->mem);
  SUNLinSolFree(R->LS);
  SUNMatDestroy(R->A);
  N_VDestroy(R->y);
  SUNContext_Free(&R->sunctx);
  SUNLogger_Destroy(&R->logger);
}

/* the warnings about root functions identically zero, if they are logged */
static int inactive_warnings(const char* path, int expected, const char* label)
{
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_WARNING
  char what[160];
  int n = count_lines(path, "root functions identically 0");
  snprintf(what, sizeof(what), "%s: %d inactive-root warning(s), expected %d",
           label, n, expected);
  return check(n == expected, what);
#else
  (void)path;
  (void)expected;
  (void)label;
  return 0;
#endif
}

/* one call to FIRKodeEvolve and the checks of its return */
static int expect(struct run* R, sunrealtype tout, int flag_expected,
                  sunrealtype t_expected, sunrealtype tol, int dir_expected,
                  const char* label, sunrealtype* tret_out)
{
  int flag, rootsfound = 0;
  sunrealtype tret;
  char what[200];

  flag = FIRKodeEvolve(R->mem, tout, R->y, &tret, FIRK_NORMAL);
  if (tret_out != NULL) { *tret_out = tret; }
  if (flag == FIRK_ROOT_RETURN) { FIRKodeGetRootInfo(R->mem, &rootsfound); }
  snprintf(what, sizeof(what),
           "%s: return %d at t = " SUN_FORMAT_G
           " with rootsfound %d, expected %d at " SUN_FORMAT_G " with %d",
           label, flag, tret, rootsfound, flag_expected, t_expected,
           dir_expected);
  return check(flag == flag_expected && SUNRabs(tret - t_expected) <= tol &&
                 rootsfound == dir_expected,
               what);
}

/* y' = -y and g = y - 1/2 with fixed steps 1/2 */
static int test_collocation(const char* path)
{
  int fails = 0;
  sunrealtype tret, ttol, ttol_expected, gsign[2];
  FIRKodeMem fm;
  char what[200];
  struct run R;

  if (setup(&R, 1, g_half, path)) { return 1; }
  fails += expect(&R, SUN_RCONST(2.0), FIRK_ROOT_RETURN, LN2, TOL_LN2, -1,
                  "g = y - 1/2", &tret);
  fm            = (FIRKodeMem)R.mem;
  ttol          = fm->root_mem->ttol;
  ttol_expected = (SUNRabs(fm->tcur) + SUNRabs(fm->h)) * fm->uround *
                  SUN_RCONST(100.0);
  snprintf(what, sizeof(what),
           "ttol = " SUN_FORMAT_E
           ", expected 100 u (|t| + |h|) = " SUN_FORMAT_E,
           ttol, ttol_expected);
  fails += check(ttol == ttol_expected, what);

  /* g on the collocation polynomial changes sign within ttol of tret */
  FIRKodeGetDky(R.mem, tret, 0, R.y);
  gsign[1] = NV_Ith_S(R.y, 0) - HALF;
  FIRKodeGetDky(R.mem, tret - ttol, 0, R.y);
  gsign[0] = NV_Ith_S(R.y, 0) - HALF;
  snprintf(what, sizeof(what),
           "g at the polynomial: " SUN_FORMAT_E " at tret - ttol, " SUN_FORMAT_E
           " at tret",
           gsign[0], gsign[1]);
  fails += check(gsign[0] > ZERO && gsign[1] <= ZERO, what);
  fails += inactive_warnings(path, 0, "g = y - 1/2");
  teardown(&R);
  return fails;
}

/* the oscillator: directions, a root at t0, a root function that is zero */
static int test_oscillator(const char* path)
{
  int fails = 0, dir;
  struct run R;

  /* both directions */
  if (setup(&R, 2, g_y1, path)) { return 1; }
  fails += expect(&R, SUN_RCONST(7.0), FIRK_ROOT_RETURN, HALF * PI, TOL_ROOT,
                  -1, "g = cos t", NULL);
  fails += expect(&R, SUN_RCONST(7.0), FIRK_ROOT_RETURN, SUN_RCONST(1.5) * PI,
                  TOL_ROOT, 1, "g = cos t", NULL);
  fails += expect(&R, SUN_RCONST(7.0), FIRK_SUCCESS, SUN_RCONST(7.0), ZERO, 0,
                  "g = cos t", NULL);
  teardown(&R);

  /* increasing only */
  if (setup(&R, 2, g_y1, path)) { return 1; }
  dir = 1;
  fails += check(FIRKodeSetRootDirection(R.mem, &dir) == FIRK_SUCCESS,
                 "FIRKodeSetRootDirection succeeds");
  fails += expect(&R, SUN_RCONST(7.0), FIRK_ROOT_RETURN, SUN_RCONST(1.5) * PI,
                  TOL_ROOT, 1, "g = cos t, increasing", NULL);
  fails += expect(&R, SUN_RCONST(7.0), FIRK_SUCCESS, SUN_RCONST(7.0), ZERO, 0,
                  "g = cos t, increasing", NULL);
  teardown(&R);

  /* decreasing only */
  if (setup(&R, 2, g_y1, path)) { return 1; }
  dir = -1;
  fails += check(FIRKodeSetRootDirection(R.mem, &dir) == FIRK_SUCCESS,
                 "FIRKodeSetRootDirection succeeds");
  fails += expect(&R, SUN_RCONST(7.0), FIRK_ROOT_RETURN, HALF * PI, TOL_ROOT,
                  -1, "g = cos t, decreasing", NULL);
  fails += expect(&R, SUN_RCONST(7.0), FIRK_SUCCESS, SUN_RCONST(7.0), ZERO, 0,
                  "g = cos t, decreasing", NULL);
  teardown(&R);

  /* a root at t0 is not returned; g moves away from zero in the first step */
  if (setup(&R, 2, g_y2, path)) { return 1; }
  fails += expect(&R, SUN_RCONST(4.0), FIRK_ROOT_RETURN, PI, TOL_ROOT, 1,
                  "g = -sin t", NULL);
  fails += expect(&R, SUN_RCONST(4.0), FIRK_SUCCESS, SUN_RCONST(4.0), ZERO, 0,
                  "g = -sin t", NULL);
  fails += inactive_warnings(path, 0, "g = -sin t");
  teardown(&R);

  /* a root function that stays zero: one warning, unless switched off */
  if (setup(&R, 2, g_zero, path)) { return 1; }
  fails += expect(&R, SUN_RCONST(4.0), FIRK_SUCCESS, SUN_RCONST(4.0), ZERO, 0,
                  "g = 0", NULL);
  teardown(&R);
  fails += inactive_warnings(path, 1, "g = 0");

  if (setup(&R, 2, g_zero, path)) { return 1; }
  fails += check(FIRKodeSetNoInactiveRootWarn(R.mem) == FIRK_SUCCESS,
                 "FIRKodeSetNoInactiveRootWarn succeeds");
  fails += expect(&R, SUN_RCONST(4.0), FIRK_SUCCESS, SUN_RCONST(4.0), ZERO, 0,
                  "g = 0, no warning", NULL);
  teardown(&R);
  fails += inactive_warnings(path, 0, "g = 0, no warning");
  return fails;
}

/* a failure of the right-hand side met by the rootfinding checks */
static int test_ydot_failure(const char* path)
{
  struct run R;
  int fails = 0, flag;
  sunrealtype tret;

  if (setup(&R, 2, g_dead, path) || FIRKodeSetFixedStep(R.mem, SUN_RCONST(0.1)) ||
      SUNContext_PushErrHandler(R.sunctx, err_handler, NULL))
  {
    return 1;
  }
  fails += expect(&R, SUN_RCONST(7.0), FIRK_ROOT_RETURN, SUN_RCONST(1.5),
                  TOL_END, -1, "g = 0 at the step end", NULL);

  fail_mode   = 1;
  last_msg[0] = '\0';
  flag        = FIRKodeEvolve(R.mem, SUN_RCONST(7.0), R.y, &tret, FIRK_NORMAL);
  fails += check(flag == FIRK_UNREC_RHSFUNC_ERR,
                 "a recoverable RHS failure at the root is "
                 "FIRK_UNREC_RHSFUNC_ERR");
  fails += check(strstr(last_msg, "right-hand side") != NULL &&
                   strstr(last_msg, "rootfinding") == NULL,
                 "the message names the right-hand side, not the rootfinding");

  fail_mode   = 2;
  last_msg[0] = '\0';
  flag        = FIRKodeEvolve(R.mem, SUN_RCONST(7.0), R.y, &tret, FIRK_NORMAL);
  fails +=
    check(flag == FIRK_RHSFUNC_FAIL,
          "an unrecoverable RHS failure at the root is FIRK_RHSFUNC_FAIL");
  fails += check(strstr(last_msg, "right-hand side") != NULL,
                 "the message names the right-hand side");
  fail_mode = 0;

  /* no root functions: the setters and getter are refused */
  {
    int dir = 1, found = 0;
    fails += check(FIRKodeRootInit(R.mem, 0, NULL) == FIRK_SUCCESS,
                   "RootInit with no root functions");
    fails += check(FIRKodeSetRootDirection(R.mem, &dir) == FIRK_ILL_INPUT &&
                     FIRKodeSetNoInactiveRootWarn(R.mem) == FIRK_ILL_INPUT &&
                     FIRKodeGetRootInfo(R.mem, &found) == FIRK_ILL_INPUT,
                   "the rootfinding setters and getter are FIRK_ILL_INPUT "
                   "without root functions");
  }
  teardown(&R);
  return fails;
}

int main(int argc, char* argv[])
{
  const char* path = (argc > 1) ? argv[1] : "firk_test_rootfind.log";
  int fails        = 0;

  fails += test_collocation(path);
  fails += test_oscillator(path);
  fails += test_ydot_failure(path);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("All tests passed\n");
  return 0;
}
