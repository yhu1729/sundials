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
 * Unit test of the counters across FIRKodeReset and FIRKodeReInit (design
 * hazards H3 and H28) on M y' = M diag(lambda) y with M = diag(1, 2, 3, 4),
 * lambda = -1, -10, -100, -1000, y(0) = 1 and the root of y_1 - 1/2, with
 * dense linear and mass solvers, and with matrix-free SPGMR (user jtimes and
 * exact diagonal preconditioner) and PCG (user mass-times and exact
 * diagonal preconditioner), so that the counters of the iterative solvers
 * are exercised too:
 *
 * - H3: Reset keeps every counter (integrator, FIRKLS, stage, mass, root),
 *   which keep growing from their values; ReInit zeroes them all; every
 *   counter that the configuration exercises is positive before Reset;
 * - H28: the constant mass matrix is set up once per problem: not again
 *   after Reset, anew after ReInit and after the mass routine is attached
 *   again (dense configuration); the integration stays accurate throughout.
 *
 * Pass criterion: every check holds and the solution is within 1e-4 of
 * exp(lambda t) at each output.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
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

static sunrealtype mass_diag(int i) { return (sunrealtype)(i + 1); }

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  int i;
  (void)t;
  (void)user_data;
  for (i = 0; i < NEQ; i++)
  {
    NV_Ith_S(ydot, i) = mass_diag(i) * lambda[i] * NV_Ith_S(y, i);
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
  for (i = 0; i < NEQ; i++)
  {
    SM_ELEMENT_D(J, i, i) = mass_diag(i) * lambda[i];
  }
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
  for (i = 0; i < NEQ; i++) { SM_ELEMENT_D(M, i, i) = mass_diag(i); }
  return 0;
}

/* the Krylov configuration: M - gamma J = diag(m_i (1 - gamma lambda_i)) */
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
  return 0;
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
  for (i = 0; i < NEQ; i++)
  {
    NV_Ith_S(z, i) = NV_Ith_S(r, i) / (mass_diag(i) * (ONE - gamma * lambda[i]));
  }
  return 0;
}

static int jtsetup(sunrealtype t, N_Vector y, N_Vector fy, void* user_data)
{
  (void)t;
  (void)y;
  (void)fy;
  (void)user_data;
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
    NV_Ith_S(Jv, i) = mass_diag(i) * lambda[i] * NV_Ith_S(v, i);
  }
  return 0;
}

static int mtsetup(sunrealtype t, void* mtimes_data)
{
  (void)t;
  (void)mtimes_data;
  return 0;
}

static int mtimes(N_Vector v, N_Vector Mv, sunrealtype t, void* mtimes_data)
{
  int i;
  (void)t;
  (void)mtimes_data;
  for (i = 0; i < NEQ; i++) { NV_Ith_S(Mv, i) = mass_diag(i) * NV_Ith_S(v, i); }
  return 0;
}

static int mpsetup(sunrealtype t, void* user_data)
{
  (void)t;
  (void)user_data;
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
  for (i = 0; i < NEQ; i++) { NV_Ith_S(z, i) = NV_Ith_S(r, i) / mass_diag(i); }
  return 0;
}

static int g(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)
{
  (void)t;
  (void)user_data;
  gout[0] = NV_Ith_S(y, 0) - HALF;
  return 0;
}

/* every counter with a getter */
#define NCOUNT 28
static const char* names[NCOUNT] =
  {"steps",          "step attempts",    "RHS evals",
   "LS setups",      "error test fails", "step solve fails",
   "Newton iters",   "Newton fails",     "root fn evals",
   "Jac evals",      "Jac step",         "prec evals",
   "prec solves",    "block LS iters",   "block LS fails",
   "jtsetup evals",  "jtimes evals",     "LS RHS evals",
   "stage LS iters", "stage LS fails",   "block solves",
   "mass setups",    "mass mult setups", "mass mults",
   "mass solves",    "mass prec evals",  "mass prec solves",
   "mass iters"};

static int read_counters(void* mem, long int c[NCOUNT])
{
  int k = 0, flag = 0;
  flag |= FIRKodeGetNumSteps(mem, &c[k++]);
  flag |= FIRKodeGetNumStepAttempts(mem, &c[k++]);
  flag |= FIRKodeGetNumRhsEvals(mem, &c[k++]);
  flag |= FIRKodeGetNumLinSolvSetups(mem, &c[k++]);
  flag |= FIRKodeGetNumErrTestFails(mem, &c[k++]);
  flag |= FIRKodeGetNumStepSolveFails(mem, &c[k++]);
  flag |= FIRKodeGetNumNonlinSolvIters(mem, &c[k++]);
  flag |= FIRKodeGetNumNonlinSolvConvFails(mem, &c[k++]);
  flag |= FIRKodeGetNumGEvals(mem, &c[k++]);
  flag |= FIRKodeGetNumJacEvals(mem, &c[k++]);
  flag |= FIRKodeGetJacNumSteps(mem, &c[k++]);
  flag |= FIRKodeGetNumPrecEvals(mem, &c[k++]);
  flag |= FIRKodeGetNumPrecSolves(mem, &c[k++]);
  flag |= FIRKodeGetNumLinIters(mem, &c[k++]);
  flag |= FIRKodeGetNumLinConvFails(mem, &c[k++]);
  flag |= FIRKodeGetNumJTSetupEvals(mem, &c[k++]);
  flag |= FIRKodeGetNumJtimesEvals(mem, &c[k++]);
  flag |= FIRKodeGetNumLinRhsEvals(mem, &c[k++]);
  flag |= FIRKodeGetNumStageLinIters(mem, &c[k++]);
  flag |= FIRKodeGetNumStageLinConvFails(mem, &c[k++]);
  flag |= FIRKodeGetNumBlockSolves(mem, &c[k++]);
  flag |= FIRKodeGetNumMassSetups(mem, &c[k++]);
  flag |= FIRKodeGetNumMassMultSetups(mem, &c[k++]);
  flag |= FIRKodeGetNumMassMult(mem, &c[k++]);
  flag |= FIRKodeGetNumMassSolves(mem, &c[k++]);
  flag |= FIRKodeGetNumMassPrecEvals(mem, &c[k++]);
  flag |= FIRKodeGetNumMassPrecSolves(mem, &c[k++]);
  flag |= FIRKodeGetNumMassIters(mem, &c[k++]);
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

/* integrate to tout through the root; the solution within 1e-4 */
static int evolve(void* mem, sunrealtype tout, N_Vector y, const char* what,
                  sunrealtype* tret_ptr)
{
  int i, flag, fails = 0;
  sunrealtype tret = ZERO, err = ZERO;
  char msg[160];

  do {
    flag = FIRKodeEvolve(mem, tout, y, &tret, FIRK_NORMAL);
  }
  while (flag == FIRK_ROOT_RETURN);
  for (i = 0; i < NEQ; i++)
  {
    err = SUNMAX(err, SUNRabs(NV_Ith_S(y, i) - SUNRexp(lambda[i] * tret)));
  }
  snprintf(msg, sizeof(msg), "%s: flag %d, t %g, error %g", what, flag,
           (double)tret, (double)err);
  fails += check(flag == FIRK_SUCCESS && err <= SUN_RCONST(1.0e-4), msg);
  *tret_ptr = tret;
  return fails;
}

/* the counters a configuration exercises: the integrator's, the root
   evaluations, and those of its solvers */
static int exercised(sunbooleantype krylov, int k)
{
  /* steps, attempts, RHS evals, LS setups, Newton iters, root fn evals,
     block solves, mass solves */
  if (k <= 3 || k == 6 || k == 8 || k == 20 || k == 24) { return 1; }
  /* dense: Jac evals, mass setups, mass mults */
  if (!krylov) { return k == 9 || k == 21 || k == 23; }
  /* Krylov: prec evals and solves, block LS iters, jtsetup and jtimes evals,
     stage LS iters, mass mult setups, mass mults, mass prec evals and solves,
     mass iters */
  return (k >= 11 && k <= 13) || k == 15 || k == 16 || k == 18 ||
         (k >= 22 && k <= 27 && k != 24);
}

static int run_config(SUNContext sunctx, sunbooleantype krylov)
{
  void* mem          = NULL;
  N_Vector y         = NULL;
  SUNMatrix A        = NULL;
  SUNMatrix M        = NULL;
  SUNLinearSolver LS = NULL, MLS = NULL;
  long int a[NCOUNT], b[NCOUNT];
  int fails = 0, k, flag;
  char msg[160];
  sunrealtype t;

  printf("%s solvers:\n", krylov ? "Krylov" : "dense");
  y = N_VNew_Serial(NEQ, sunctx);
  N_VConst(ONE, y);
  mem  = FIRKodeCreate(sunctx);
  flag = FIRKodeInit(mem, f, ZERO, y) ||
         FIRKodeSStolerances(mem, SUN_RCONST(1.0e-7), SUN_RCONST(1.0e-7)) ||
         FIRKodeRootInit(mem, 1, g);
  if (krylov)
  {
    LS   = SUNLinSol_SPGMR(y, SUN_PREC_LEFT, 0, sunctx);
    MLS  = SUNLinSol_PCG(y, SUN_PREC_LEFT, 0, sunctx);
    flag = flag || FIRKodeSetLinearSolver(mem, LS, NULL) ||
           FIRKodeSetJacTimes(mem, jtsetup, jtimes) ||
           FIRKodeSetPreconditioner(mem, psetup, psolve) ||
           FIRKodeSetMassLinearSolver(mem, MLS, NULL, SUNFALSE) ||
           FIRKodeSetMassTimes(mem, mtsetup, mtimes, NULL) ||
           FIRKodeSetMassPreconditioner(mem, mpsetup, mpsolve);
  }
  else
  {
    A    = SUNDenseMatrix(NEQ, NEQ, sunctx);
    M    = SUNDenseMatrix(NEQ, NEQ, sunctx);
    LS   = SUNLinSol_Dense(y, A, sunctx);
    MLS  = SUNLinSol_Dense(y, M, sunctx);
    flag = flag || FIRKodeSetLinearSolver(mem, LS, A) ||
           FIRKodeSetJacFn(mem, jac) ||
           FIRKodeSetMassLinearSolver(mem, MLS, M, SUNFALSE) ||
           FIRKodeSetMassFn(mem, mass);
  }
  if (flag) { return 1; }

  fails += evolve(mem, ONE, y, "the integration to t = 1", &t);
  fails += check(read_counters(mem, a) == 0, "the getters succeed");
  for (k = 0; k < NCOUNT; k++)
  {
    if (!exercised(krylov, k)) { continue; }
    snprintf(msg, sizeof(msg), "%s is positive before Reset (%ld)", names[k],
             a[k]);
    fails += check(a[k] > 0, msg);
  }
  if (!krylov) { fails += check(a[21] == 1, "one mass setup before Reset"); }

  /* H3: Reset keeps every counter */
  fails += check(FIRKodeReset(mem, t, y) == FIRK_SUCCESS, "Reset succeeds");
  read_counters(mem, b);
  for (k = 0; k < NCOUNT; k++)
  {
    snprintf(msg, sizeof(msg), "H3: Reset keeps %s (%ld -> %ld)", names[k],
             a[k], b[k]);
    fails += check(b[k] == a[k], msg);
  }

  /* ... and the counters keep growing from their values; the mass matrix
     is not set up again (H28) */
  fails += evolve(mem, SUN_RCONST(1.5), y, "the integration after Reset", &t);
  read_counters(mem, b);
  for (k = 0; k < NCOUNT; k++)
  {
    snprintf(msg, sizeof(msg), "H3: %s grows from its value (%ld -> %ld)",
             names[k], a[k], b[k]);
    fails += check(b[k] >= a[k], msg);
  }
  snprintf(msg, sizeof(msg),
           "H3: steps (%ld -> %ld) and setups (%ld -> %ld) after Reset", a[0],
           b[0], a[3], b[3]);
  fails += check(b[0] > a[0] && b[3] > a[3], msg);
  if (!krylov)
  {
    snprintf(msg, sizeof(msg), "H28: no mass setup after Reset (%ld -> %ld)",
             a[21], b[21]);
    fails += check(b[21] == a[21], msg);
  }

  /* H3: ReInit zeroes every counter; H28: the new problem has its mass
     matrix set up anew */
  fails += check(FIRKodeReInit(mem, t, y) == FIRK_SUCCESS, "ReInit succeeds");
  read_counters(mem, b);
  for (k = 0; k < NCOUNT; k++)
  {
    snprintf(msg, sizeof(msg), "H3: ReInit zeroes %s (%ld)", names[k], b[k]);
    fails += check(b[k] == 0, msg);
  }
  fails += evolve(mem, SUN_RCONST(2.0), y, "the integration after ReInit", &t);
  read_counters(mem, b);
  if (krylov) { fails += check(b[0] > 0, "steps after ReInit"); }
  else
  {
    snprintf(msg, sizeof(msg), "H28: one mass setup after ReInit (%ld)", b[21]);
    fails += check(b[0] > 0 && b[21] == 1, msg);

    /* attaching the mass routine again sets the mass matrix up anew */
    fails += check(FIRKodeSetMassFn(mem, mass) == FIRKLS_SUCCESS,
                   "SetMassFn succeeds");
    fails += check(FIRKodeReset(mem, t, y) == FIRK_SUCCESS, "Reset succeeds");
    fails += evolve(mem, SUN_RCONST(2.5), y,
                    "the integration after SetMassFn and Reset", &t);
    read_counters(mem, b);
    snprintf(msg, sizeof(msg), "H28: one more mass setup after SetMassFn (%ld)",
             b[21]);
    fails += check(b[21] == 2, msg);
  }

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNLinSolFree(MLS);
  SUNMatDestroy(A);
  SUNMatDestroy(M);
  N_VDestroy(y);
  return fails;
}

int main(void)
{
  SUNContext sunctx = NULL;
  int fails         = 0;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }
  fails += run_config(sunctx, SUNFALSE);
  fails += run_config(sunctx, SUNTRUE);
  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("All tests passed\n");
  return 0;
}
