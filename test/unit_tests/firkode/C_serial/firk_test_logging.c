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
 * Unit test of the log records (design Table 8.4, hazard H31), built at
 * logging level 4 or more: integrations of M y' = M diag(lambda) y, lambda =
 * -1, -10, -100, -1000, write their info and debug records to the file
 * argv[1], through
 *
 * - dense Jacobian and mass solvers, a root function and a refiltered
 *   estimate after rejections;
 * - SPGMR with jtimes and psolve, and a PCG mass solver (the conversion of
 *   the stacked solve);
 * - recoverable RHS failures in the Newton iteration;
 * - an unrecoverable RHS failure, a psolve failure and a failing mass
 *   solve, each of which ends its integration;
 * - the other exits of FIRKodeEvolve: too much work, an error test failure
 *   with maxnef = 1, a convergence failure with maxncf = 1 after a diverging
 *   Newton iteration (a zero Jacobian) and after the iteration limit
 *   (maxcor = 1), and a failing preconditioner setup.
 *
 * The file is then read back: every begin-<region> has its end-<region>,
 * properly nested, on every exit path. firk_test_logging.py parses the same
 * file with suntools.logs.
 *
 * Pass criterion: each integration ends as expected and the records are
 * balanced and nested.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_logger.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunlinsol/sunlinsol_pcg.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ZERO SUN_RCONST(0.0)
#define HALF SUN_RCONST(0.5)
#define ONE  SUN_RCONST(1.0)
#define NEQ  4

static const sunrealtype lambda[NEQ] = {SUN_RCONST(-1.0), SUN_RCONST(-10.0),
                                        SUN_RCONST(-100.0), SUN_RCONST(-1000.0)};

/* failure injection: the RHS fails with rhs_fail_value on every call
   rhs_fail_every-th call after rhs_fail_after calls; psolve and the mass
   psolve fail from their fail_at-th call */
static long int nrhs = 0, rhs_fail_after = 0, rhs_fail_every = 0;
static int rhs_fail_value = 0;
static long int npsolve = 0, psolve_fail_at = 0;
static long int nmpsolve = 0, mpsolve_fail_at = 0;

/* the other exits: limits (when positive), a zero Jacobian, which makes the
   Newton iteration diverge on the stiff components, and a failing psetup */
static long int opt_mxstep = 0;
static int opt_maxnef = 0, opt_maxncf = 0, opt_maxcor = 0;
static sunrealtype opt_h0 = ZERO;
static int opt_zero_jac = 0, psetup_fail = 0;

static sunrealtype mdiag(int i) { return (sunrealtype)(i + 1); }

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  int i;
  (void)t;
  (void)user_data;
  nrhs++;
  if (rhs_fail_every > 0 && nrhs > rhs_fail_after &&
      (nrhs - rhs_fail_after) % rhs_fail_every == 0)
  {
    return rhs_fail_value;
  }
  for (i = 0; i < NEQ; i++)
  {
    NV_Ith_S(ydot, i) = mdiag(i) * lambda[i] * NV_Ith_S(y, i);
  }
  return 0;
}

static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int i;
  (void)t;
  (void)y;
  (void)fy;
  (void)user_data;
  (void)tmp1;
  (void)tmp2;
  (void)tmp3;
  if (opt_zero_jac) { return SUNMatZero(J) == SUN_SUCCESS ? 0 : -1; }
  for (i = 0; i < NEQ; i++) { SM_ELEMENT_D(J, i, i) = mdiag(i) * lambda[i]; }
  return 0;
}

static int mass(sunrealtype t, SUNMatrix M, void* user_data, N_Vector tmp1,
                N_Vector tmp2, N_Vector tmp3)
{
  int i;
  (void)t;
  (void)user_data;
  (void)tmp1;
  (void)tmp2;
  (void)tmp3;
  SUNMatZero(M);
  for (i = 0; i < NEQ; i++) { SM_ELEMENT_D(M, i, i) = mdiag(i); }
  return 0;
}

static int jtimes(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y,
                  N_Vector fy, void* user_data, N_Vector tmp)
{
  int i;
  (void)t;
  (void)y;
  (void)fy;
  (void)user_data;
  (void)tmp;
  for (i = 0; i < NEQ; i++)
  {
    NV_Ith_S(Jv, i) = mdiag(i) * lambda[i] * NV_Ith_S(v, i);
  }
  return 0;
}

static int psetup(sunrealtype t, N_Vector y, N_Vector fy, sunbooleantype jok,
                  sunbooleantype* jcur, sunrealtype gamma, void* user_data)
{
  (void)t;
  (void)y;
  (void)fy;
  (void)jok;
  (void)gamma;
  (void)user_data;
  *jcur = SUNTRUE;
  return psetup_fail ? -1 : 0;
}

static int psolve(sunrealtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z,
                  sunrealtype gamma, sunrealtype delta, int lr, void* user_data)
{
  int i;
  (void)t;
  (void)y;
  (void)fy;
  (void)delta;
  (void)lr;
  (void)user_data;
  npsolve++;
  if (psolve_fail_at > 0 && npsolve >= psolve_fail_at) { return -1; }
  for (i = 0; i < NEQ; i++)
  {
    NV_Ith_S(z, i) = NV_Ith_S(r, i) / (mdiag(i) * (ONE - gamma * lambda[i]));
  }
  return 0;
}

static int mtimes(N_Vector v, N_Vector Mv, sunrealtype t, void* mtimes_data)
{
  int i;
  (void)t;
  (void)mtimes_data;
  for (i = 0; i < NEQ; i++) { NV_Ith_S(Mv, i) = mdiag(i) * NV_Ith_S(v, i); }
  return 0;
}

static int mpsolve(sunrealtype t, N_Vector r, N_Vector z, sunrealtype delta,
                   int lr, void* user_data)
{
  int i;
  (void)t;
  (void)delta;
  (void)lr;
  (void)user_data;
  nmpsolve++;
  if (mpsolve_fail_at > 0 && nmpsolve >= mpsolve_fail_at) { return -1; }
  for (i = 0; i < NEQ; i++) { NV_Ith_S(z, i) = NV_Ith_S(r, i) / mdiag(i); }
  return 0;
}

static int g(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)
{
  (void)t;
  (void)user_data;
  gout[0] = NV_Ith_S(y, 0) - HALF;
  return 0;
}

/* an integration to tout through roots; returns the last flag */
static int integrate(SUNContext sunctx, int krylov, sunrealtype tout)
{
  int flag;
  sunrealtype tret;
  N_Vector y;
  SUNMatrix A = NULL, M = NULL;
  SUNLinearSolver LS, MLS;
  void* mem;

  nrhs = npsolve = nmpsolve = 0;
  y                         = N_VNew_Serial(NEQ, sunctx);
  N_VConst(ONE, y);
  if (krylov)
  {
    LS  = SUNLinSol_SPGMR(y, SUN_PREC_LEFT, 0, sunctx);
    MLS = SUNLinSol_PCG(y, SUN_PREC_LEFT, 0, sunctx);
  }
  else
  {
    A   = SUNDenseMatrix(NEQ, NEQ, sunctx);
    M   = SUNDenseMatrix(NEQ, NEQ, sunctx);
    LS  = SUNLinSol_Dense(y, A, sunctx);
    MLS = SUNLinSol_Dense(y, M, sunctx);
  }
  mem  = FIRKodeCreate(sunctx);
  flag = FIRKodeInit(mem, f, ZERO, y);
  if (!flag)
  {
    flag = FIRKodeSStolerances(mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-6));
  }
  if (!flag) { flag = FIRKodeSetLinearSolver(mem, LS, A); }
  if (!flag) { flag = FIRKodeSetMassLinearSolver(mem, MLS, M, SUNFALSE); }
  if (!flag) { flag = FIRKodeRootInit(mem, 1, g); }
  if (!flag && krylov)
  {
    flag = FIRKodeSetJacTimes(mem, NULL, jtimes);
    if (!flag)
    {
      flag = FIRKodeSetPreconditioner(mem, psetup_fail ? psetup : NULL, psolve);
    }
    if (!flag) { flag = FIRKodeSetMassTimes(mem, NULL, mtimes, NULL); }
    if (!flag) { flag = FIRKodeSetMassPreconditioner(mem, NULL, mpsolve); }
  }
  else if (!flag)
  {
    flag = FIRKodeSetJacFn(mem, jac);
    if (!flag) { flag = FIRKodeSetMassFn(mem, mass); }
  }
  if (!flag && opt_mxstep > 0)
  {
    flag = FIRKodeSetMaxNumSteps(mem, opt_mxstep);
  }
  if (!flag && opt_h0 > ZERO) { flag = FIRKodeSetInitStep(mem, opt_h0); }
  if (!flag && (opt_maxnef > 0 || opt_maxncf > 0 || opt_maxcor > 0))
  {
    flag = FIRKodeSetMode(mem, FIRK_EXPERT);
    if (!flag && opt_maxnef > 0)
    {
      flag = FIRKodeSetMaxErrTestFails(mem, opt_maxnef);
    }
    if (!flag && opt_maxncf > 0)
    {
      flag = FIRKodeSetMaxConvFails(mem, opt_maxncf);
    }
    if (!flag && opt_maxcor > 0)
    {
      flag = FIRKodeSetMaxNonlinIters(mem, opt_maxcor);
    }
  }
  if (!flag)
  {
    do {
      flag = FIRKodeEvolve(mem, tout, y, &tret, FIRK_NORMAL);
    }
    while (flag == FIRK_ROOT_RETURN);
  }

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNLinSolFree(MLS);
  if (A) { SUNMatDestroy(A); }
  if (M) { SUNMatDestroy(M); }
  N_VDestroy(y);
  return flag;
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

/* every begin-<region> has its end-<region>, properly nested */
static int check_balance(const char* path, long int* nrecords)
{
  char line[4096], stack[64][96];
  int depth = 0, bad = 0;
  long int n = 0;
  FILE* fp   = fopen(path, "r");

  if (fp == NULL) { return 1; }
  while (fgets(line, sizeof(line), fp) != NULL)
  {
    char *p, *q;
    int k;
    /* "[LEVEL][rank r][scope][label] payload" */
    p = line;
    for (k = 0; k < 3 && p != NULL; k++)
    {
      p = strchr(p, ']');
      if (p != NULL) { p++; }
    }
    if (p == NULL || *p != '[') { continue; }
    p++;
    q = strchr(p, ']');
    if (q == NULL) { continue; }
    *q = '\0';
    n++;
    if (strncmp(p, "begin-", 6) == 0)
    {
      if (depth == 64) { bad++; }
      else
      {
        strncpy(stack[depth], p + 6, sizeof(stack[0]) - 1);
        stack[depth][sizeof(stack[0]) - 1] = '\0';
        depth++;
      }
    }
    else if (strncmp(p, "end-", 4) == 0)
    {
      if (depth == 0 || strcmp(stack[depth - 1], p + 4) != 0)
      {
        if (bad == 0)
        {
          printf("  unbalanced record %s after %ld records (open: %s)\n", p, n,
                 depth > 0 ? stack[depth - 1] : "none");
        }
        bad++;
      }
      else { depth--; }
    }
  }
  fclose(fp);
  if (depth != 0)
  {
    printf("  %d region(s) left open, the innermost %s\n", depth,
           stack[depth - 1]);
    bad++;
  }
  *nrecords = n;
  return bad;
}

int main(int argc, char* argv[])
{
  SUNContext sunctx = NULL;
  SUNLogger logger  = NULL;
  const char* path  = (argc > 1) ? argv[1] : "firk_test_logging.log";
  int fails         = 0, flag;
  long int nrecords = 0;
  char what[160];
  FILE* fp;

  /* a fresh file for all records */
  fp = fopen(path, "w");
  if (fp == NULL) { return 1; }
  fclose(fp);

  if (SUNLogger_Create(SUN_COMM_NULL, 0, &logger) ||
      SUNLogger_SetInfoFilename(logger, path) ||
      SUNLogger_SetDebugFilename(logger, path) ||
      SUNContext_Create(SUN_COMM_NULL, &sunctx) ||
      SUNContext_SetLogger(sunctx, logger))
  {
    return 1;
  }

  /* dense solvers, a root, rejections and refiltered estimates */
  flag = integrate(sunctx, 0, SUN_RCONST(2.0));
  fails += check(flag == FIRK_SUCCESS, "the dense integration succeeds");

  /* Krylov solvers, the stacked conversion with a PCG mass solver */
  flag = integrate(sunctx, 1, SUN_RCONST(0.5));
  fails += check(flag == FIRK_SUCCESS, "the Krylov integration succeeds");

  /* recoverable RHS failures in the Newton iteration */
  rhs_fail_after = 20;
  rhs_fail_every = 7;
  rhs_fail_value = 1;
  flag           = integrate(sunctx, 0, SUN_RCONST(0.5));
  snprintf(what, sizeof(what),
           "recoverable RHS failures end with success or a failure flag (%d)",
           flag);
  fails += check(flag == FIRK_SUCCESS || flag < 0, what);

  /* an unrecoverable RHS failure at the 40th call */
  rhs_fail_after = 0;
  rhs_fail_every = 40;
  rhs_fail_value = -1;
  flag           = integrate(sunctx, 0, SUN_RCONST(0.5));
  snprintf(what, sizeof(what), "an unrecoverable RHS failure ends it (%d)", flag);
  fails += check(flag < 0, what);
  rhs_fail_every = 0;

  /* an unrecoverable psolve failure in the stacked solve */
  psolve_fail_at = 20;
  flag           = integrate(sunctx, 1, SUN_RCONST(0.5));
  snprintf(what, sizeof(what), "an unrecoverable psolve failure ends it (%d)",
           flag);
  fails += check(flag < 0, what);
  psolve_fail_at = 0;

  /* an unrecoverable mass psolve failure */
  mpsolve_fail_at = 15;
  flag            = integrate(sunctx, 1, SUN_RCONST(0.5));
  snprintf(what, sizeof(what), "a failing mass solve ends it (%d)", flag);
  fails += check(flag < 0, what);
  mpsolve_fail_at = 0;

  /* too much work */
  opt_mxstep = 3;
  flag       = integrate(sunctx, 0, SUN_RCONST(2.0));
  snprintf(what, sizeof(what), "mxstep = 3 ends with FIRK_TOO_MUCH_WORK (%d)",
           flag);
  fails += check(flag == FIRK_TOO_MUCH_WORK, what);
  opt_mxstep = 0;

  /* an error test failure with maxnef = 1 */
  opt_h0     = ONE;
  opt_maxnef = 1;
  flag       = integrate(sunctx, 0, SUN_RCONST(2.0));
  snprintf(what, sizeof(what),
           "h0 = 1 with maxnef = 1 ends with FIRK_ERR_FAILURE (%d)", flag);
  fails += check(flag == FIRK_ERR_FAILURE, what);
  opt_maxnef = 0;

  /* a diverging Newton iteration with maxncf = 1 */
  opt_zero_jac = 1;
  opt_maxncf   = 1;
  flag         = integrate(sunctx, 0, SUN_RCONST(2.0));
  snprintf(what, sizeof(what),
           "a zero Jacobian with maxncf = 1 ends with FIRK_CONV_FAILURE (%d)",
           flag);
  fails += check(flag == FIRK_CONV_FAILURE, what);
  opt_zero_jac = 0;

  /* the iteration limit with maxncf = 1 */
  opt_maxcor = 1;
  flag       = integrate(sunctx, 0, SUN_RCONST(2.0));
  snprintf(what, sizeof(what),
           "maxcor = 1 with maxncf = 1 ends with FIRK_CONV_FAILURE (%d)", flag);
  fails += check(flag == FIRK_CONV_FAILURE, what);
  opt_maxcor = opt_maxncf = 0;
  opt_h0                  = ZERO;

  /* a failing preconditioner setup */
  psetup_fail = 1;
  flag        = integrate(sunctx, 1, SUN_RCONST(0.5));
  snprintf(what, sizeof(what),
           "a failing psetup ends with FIRK_LSETUP_FAIL (%d)", flag);
  fails += check(flag == FIRK_LSETUP_FAIL, what);
  psetup_fail = 0;

  /* close the file, then read it back */
  SUNContext_Free(&sunctx);
  SUNLogger_Destroy(&logger);
  fails += check(check_balance(path, &nrecords) == 0,
                 "H31: the records are balanced and nested");
  snprintf(what, sizeof(what), "the log has records (%ld)", nrecords);
  fails += check(nrecords > 100, what);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("All tests passed\n");
  return 0;
}
