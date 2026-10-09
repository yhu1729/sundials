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
 * Unit test of the operations on the stage stack (design Section 8.3,
 * Table 8.8) on y' = diag(lambda) y, lambda = -1, -10, -100, -1000, with the
 * three-stage method, a matrix-free SPGMR (user jtimes and psolve) and a
 * dense block solver. White box: the operations of the integrator's stacks
 * are replaced right after FIRKodeInit, so that the stacked vectors the
 * internal FGMRES clones from them carry the replacements too:
 *
 * - H9, H10: no forbidden operation of Table 8.8 (and no other operation
 *   outside the whitelist) is applied to a stacked vector; lengths come from
 *   N_VGetLength;
 * - H15: N_VLinearCombination and N_VScaleAddMulti, and N_VDotProdMulti
 *   (Table 8.8), are applied to stacked vectors only inside the stacked
 *   FGMRES solve;
 * - H11: SUNLinSolSetZeroGuess(.., SUNTRUE) precedes every stacked FGMRES
 *   solve;
 * - H16: the weight view is never cloned (and ASan sees its destruction);
 * - H14: no vector is cloned during the steps (a counting clone on y0);
 *   a clone failure at any point of FIRKodeInit returns FIRK_MEM_FAIL and
 *   leaves an integrator that FIRKodeFree releases;
 * - one stage has no stack.
 *
 * Pass criterion: no trapped operation, no clone per step, a zero guess
 * before every stacked solve, every injected clone failure reported, and
 * accurate solutions.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>

#include "firkode/firkode.h"
#include "firkode/firkode_impl.h" /* the stacks, firkSetNumStagesAny */
#include "firkode/firkode_ls.h"
#include "firkode/firkode_ls_impl.h" /* the stacked FGMRES */
#include "nvector/nvector_manyvector.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunlinsol/sunlinsol_spfgmr.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define NEQ  4

static const sunrealtype lambda[NEQ] = {SUN_RCONST(-1.0), SUN_RCONST(-10.0),
                                        SUN_RCONST(-100.0), SUN_RCONST(-1000.0)};

/*---------------------------------------------------------------
  The problem
  ---------------------------------------------------------------*/

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  int i;
  (void)t;
  (void)user_data;
  for (i = 0; i < NEQ; i++) { NV_Ith_S(ydot, i) = lambda[i] * NV_Ith_S(y, i); }
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
  for (i = 0; i < NEQ; i++) { NV_Ith_S(Jv, i) = lambda[i] * NV_Ith_S(v, i); }
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
    NV_Ith_S(z, i) = NV_Ith_S(r, i) / (ONE - gamma * lambda[i]);
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
  for (i = 0; i < NEQ; i++) { SM_ELEMENT_D(J, i, i) = lambda[i]; }
  return 0;
}

/*---------------------------------------------------------------
  A counting clone for y0 and its clones (H14), with fault injection
  ---------------------------------------------------------------*/

static long int nclones = 0;
static long int fail_at = 0; /* clone (1-based) that fails, 0 none */
static long int nfailed = 0;

static N_Vector clone_counting(N_Vector w)
{
  nclones++;
  if (nclones == fail_at)
  {
    nfailed++;
    return NULL;
  }
  return N_VClone_Serial(w); /* copies the operations of w */
}

/*---------------------------------------------------------------
  The operations of the stacks (H9, H10, H15, H16)
  ---------------------------------------------------------------*/

static long int ntrapped              = 0; /* forbidden operations     */
static long int nmulti_outside        = 0; /* H15 outside the FGMRES   */
static long int nview_clones          = 0; /* clones of the weight view */
static int in_fgmr                    = 0;
static long int nfgmr_solves          = 0; /* stacked FGMRES solves     */
static long int nfgmr_no_zero         = 0; /* H11: without a zero guess */
static const char* first_trap         = NULL;
static int (*fgmr_solve)(SUNLinearSolver, SUNMatrix, N_Vector, N_Vector,
                         sunrealtype) = NULL;

static void trap(const char* op)
{
  if (ntrapped == 0) { first_trap = op; }
  ntrapped++;
}

static sunrealtype trap_maxnorm(N_Vector x)
{
  (void)x;
  trap("N_VMaxNorm");
  return ZERO;
}

static sunrealtype trap_min(N_Vector x)
{
  (void)x;
  trap("N_VMin");
  return ZERO;
}

static sunrealtype trap_l1norm(N_Vector x)
{
  (void)x;
  trap("N_VL1Norm");
  return ZERO;
}

static sunrealtype trap_wl2norm(N_Vector x, N_Vector w)
{
  (void)x;
  (void)w;
  trap("N_VWL2Norm");
  return ZERO;
}

static sunrealtype trap_wrmsnormmask(N_Vector x, N_Vector w, N_Vector id)
{
  (void)x;
  (void)w;
  (void)id;
  trap("N_VWrmsNormMask");
  return ZERO;
}

static sunbooleantype trap_invtest(N_Vector x, N_Vector z)
{
  (void)x;
  (void)z;
  trap("N_VInvTest");
  return SUNTRUE;
}

static sunbooleantype trap_constrmask(N_Vector c, N_Vector x, N_Vector m)
{
  (void)c;
  (void)x;
  (void)m;
  trap("N_VConstrMask");
  return SUNTRUE;
}

static sunrealtype trap_minquotient(N_Vector num, N_Vector denom)
{
  (void)num;
  (void)denom;
  trap("N_VMinQuotient");
  return ZERO;
}

static sunindextype trap_getlocallength(N_Vector x)
{
  (void)x;
  trap("N_VGetLocalLength");
  return 0;
}

static sunrealtype* trap_getarraypointer(N_Vector x)
{
  (void)x;
  trap("N_VGetArrayPointer");
  return NULL;
}

static void trap_abs(N_Vector x, N_Vector z)
{
  (void)x;
  (void)z;
  trap("N_VAbs");
}

static void trap_inv(N_Vector x, N_Vector z)
{
  (void)x;
  (void)z;
  trap("N_VInv");
}

static void trap_addconst(N_Vector x, sunrealtype b, N_Vector z)
{
  (void)x;
  (void)b;
  (void)z;
  trap("N_VAddConst");
}

static void trap_compare(sunrealtype c, N_Vector x, N_Vector z)
{
  (void)c;
  (void)x;
  (void)z;
  trap("N_VCompare");
}

static SUNErrCode linearcombination_counting(int nvec, sunrealtype* c,
                                             N_Vector* X, N_Vector z)
{
  if (!in_fgmr) { nmulti_outside++; }
  return N_VLinearCombination_ManyVector(nvec, c, X, z);
}

static SUNErrCode scaleaddmulti_counting(int nvec, sunrealtype* a, N_Vector x,
                                         N_Vector* Y, N_Vector* Z)
{
  if (!in_fgmr) { nmulti_outside++; }
  return N_VScaleAddMulti_ManyVector(nvec, a, x, Y, Z);
}

static SUNErrCode dotprodmulti_counting(int nvec, N_Vector x, N_Vector* Y,
                                        sunrealtype* dotprods)
{
  if (!in_fgmr) { nmulti_outside++; }
  return N_VDotProdMulti_ManyVector(nvec, x, Y, dotprods);
}

static N_Vector clone_view(N_Vector w)
{
  nview_clones++;
  return N_VClone_ManyVector(w);
}

static void replace_stack_ops(N_Vector stk)
{
  stk->ops->nvmaxnorm           = trap_maxnorm;
  stk->ops->nvmin               = trap_min;
  stk->ops->nvl1norm            = trap_l1norm;
  stk->ops->nvwl2norm           = trap_wl2norm;
  stk->ops->nvwrmsnormmask      = trap_wrmsnormmask;
  stk->ops->nvinvtest           = trap_invtest;
  stk->ops->nvconstrmask        = trap_constrmask;
  stk->ops->nvminquotient       = trap_minquotient;
  stk->ops->nvgetlocallength    = trap_getlocallength;
  stk->ops->nvgetarraypointer   = trap_getarraypointer;
  stk->ops->nvabs               = trap_abs;
  stk->ops->nvinv               = trap_inv;
  stk->ops->nvaddconst          = trap_addconst;
  stk->ops->nvcompare           = trap_compare;
  stk->ops->nvlinearcombination = linearcombination_counting;
  stk->ops->nvscaleaddmulti     = scaleaddmulti_counting;
  stk->ops->nvdotprodmulti      = dotprodmulti_counting;
}

/* the stacked FGMRES solve, marked */
static int fgmr_solve_marked(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                             N_Vector b, sunrealtype tol)
{
  int retval;
  nfgmr_solves++;
  if (!((SUNLinearSolverContent_SPFGMR)S->content)->zeroguess)
  {
    nfgmr_no_zero++;
  }
  in_fgmr = 1;
  retval  = fgmr_solve(S, A, x, b, tol);
  in_fgmr = 0;
  return retval;
}

/*---------------------------------------------------------------
  Checks
  ---------------------------------------------------------------*/

static int check(int ok, const char* what)
{
  if (!ok)
  {
    printf("  FAIL: %s\n", what);
    return 1;
  }
  return 0;
}

static sunrealtype exact_error(N_Vector y, sunrealtype t)
{
  int i;
  sunrealtype err = ZERO;
  for (i = 0; i < NEQ; i++)
  {
    err = SUNMAX(err, SUNRabs(NV_Ith_S(y, i) - SUNRexp(lambda[i] * t)));
  }
  return err;
}

/* y0 with the counting clone */
static N_Vector new_y0(SUNContext sunctx)
{
  N_Vector y = N_VNew_Serial(NEQ, sunctx);
  if (y == NULL) { return NULL; }
  y->ops->nvclone = clone_counting;
  N_VConst(ONE, y);
  return y;
}

/* the stack operations over an integration, with SPGMR (dense = 0) or a
   dense block solver (dense = 1) */
static int check_ops(SUNContext sunctx, int dense)
{
  const char* name = dense ? "dense" : "SPGMR";
  int fails        = 0;
  long int nclones_step1;
  sunrealtype tret;
  N_Vector y;
  SUNMatrix A = NULL;
  SUNLinearSolver LS;
  FIRKodeMem fm;
  FIRKLsMem ls;
  void* mem;
  char what[200];

  ntrapped = nmulti_outside = nview_clones = 0;
  nfgmr_solves = nfgmr_no_zero = 0;
  first_trap                   = NULL;
  nclones                      = 0;
  fail_at                      = 0;

  y   = new_y0(sunctx);
  mem = FIRKodeCreate(sunctx);
  if (FIRKodeInit(mem, f, ZERO, y) ||
      FIRKodeSStolerances(mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-6)))
  {
    return 1;
  }
  if (dense)
  {
    A  = SUNDenseMatrix(NEQ, NEQ, sunctx);
    LS = SUNLinSol_Dense(y, A, sunctx);
    if (FIRKodeSetLinearSolver(mem, LS, A) || FIRKodeSetJacFn(mem, jac))
    {
      return 1;
    }
  }
  else
  {
    LS = SUNLinSol_SPGMR(y, SUN_PREC_LEFT, 0, sunctx);
    if (FIRKodeSetLinearSolver(mem, LS, NULL) ||
        FIRKodeSetJacTimes(mem, NULL, jtimes) ||
        FIRKodeSetPreconditioner(mem, NULL, psolve))
    {
      return 1;
    }
  }

  /* the stacks exist after Init: replace their operations before the
     stacked FGMRES clones them, then mark its solve */
  fm = (FIRKodeMem)mem;
  fails += check(fm->g_stk != NULL && fm->delta_stk != NULL && fm->ewt_stk != NULL,
                 "the stacks exist after FIRKodeInit");
  replace_stack_ops(fm->g_stk);
  replace_stack_ops(fm->delta_stk);
  replace_stack_ops(fm->ewt_stk);
  fm->ewt_stk->ops->nvclone = clone_view;
  if (fm->linit(fm) != 0) { return fails + 1; }
  ls                     = (FIRKLsMem)fm->lmem;
  fgmr_solve             = ls->SPFGMR->ops->solve;
  ls->SPFGMR->ops->solve = fgmr_solve_marked;

  fails += check(FIRKodeEvolve(mem, ONE, y, &tret, FIRK_ONE_STEP) == FIRK_SUCCESS,
                 "the first step succeeds");
  nclones_step1 = nclones;
  fails += check(FIRKodeEvolve(mem, SUN_RCONST(0.1), y, &tret, FIRK_NORMAL) ==
                   FIRK_SUCCESS,
                 "the integration succeeds");

  snprintf(what, sizeof(what),
           "%s: H9: %ld forbidden operations on stacked vectors (first: %s)",
           name, ntrapped, first_trap ? first_trap : "none");
  fails += check(ntrapped == 0, what);
  snprintf(what, sizeof(what),
           "%s: H15: %ld multi-vector operations on stacks outside FGMRES",
           name, nmulti_outside);
  fails += check(nmulti_outside == 0, what);
  snprintf(what, sizeof(what),
           "%s: H11: %ld of %ld stacked solves without a zero guess", name,
           nfgmr_no_zero, nfgmr_solves);
  fails += check(nfgmr_solves > 0 && nfgmr_no_zero == 0, what);
  snprintf(what, sizeof(what), "%s: H16: %ld clones of the weight view", name,
           nview_clones);
  fails += check(nview_clones == 0, what);
  snprintf(what, sizeof(what),
           "%s: H14: %ld clones during the steps after the first", name,
           nclones - nclones_step1);
  fails += check(nclones == nclones_step1, what);
  snprintf(what, sizeof(what), "%s: the solution error %g", name,
           (double)exact_error(y, tret));
  fails += check(exact_error(y, tret) <= SUN_RCONST(1.0e-4), what);

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  if (A) { SUNMatDestroy(A); }
  N_VDestroy(y);
  return fails;
}

/* H14: a clone failure at each clone of FIRKodeInit */
static int check_clone_failures(SUNContext sunctx)
{
  int fails = 0, flag;
  long int k, ninit;
  N_Vector y;
  void* mem;
  char what[200];

  /* the clones of a successful Init */
  y       = new_y0(sunctx);
  nclones = 0;
  fail_at = 0;
  mem     = FIRKodeCreate(sunctx);
  if (FIRKodeInit(mem, f, ZERO, y)) { return 1; }
  ninit = nclones;
  FIRKodeFree(&mem);

  for (k = 1; k <= ninit; k++)
  {
    nclones = 0;
    nfailed = 0;
    fail_at = k;
    mem     = FIRKodeCreate(sunctx);
    flag    = FIRKodeInit(mem, f, ZERO, y);
    snprintf(what, sizeof(what),
             "H14: a failure of clone %ld of %ld in Init returns "
             "FIRK_MEM_FAIL (flag %d)",
             k, ninit, flag);
    fails += check(nfailed == 1 && flag == FIRK_MEM_FAIL, what);
    FIRKodeFree(&mem);
  }
  fail_at = 0;
  N_VDestroy(y);
  return fails;
}

/* one stage: no stack, the "stacks" are the stage vectors */
static int check_one_stage(SUNContext sunctx)
{
  int fails = 0;
  sunrealtype tret;
  N_Vector y;
  SUNMatrix A;
  SUNLinearSolver LS;
  FIRKodeMem fm;
  void* mem;
  char what[200];

  fail_at = 0;
  y       = new_y0(sunctx);
  A       = SUNDenseMatrix(NEQ, NEQ, sunctx);
  LS      = SUNLinSol_Dense(y, A, sunctx);
  mem     = FIRKodeCreate(sunctx);
  if (FIRKodeInit(mem, f, ZERO, y) ||
      FIRKodeSStolerances(mem, SUN_RCONST(1.0e-3), SUN_RCONST(1.0e-3)) ||
      FIRKodeSetLinearSolver(mem, LS, A) || FIRKodeSetJacFn(mem, jac) ||
      firkSetNumStagesAny(mem, 1))
  {
    return 1;
  }
  fails += check(FIRKodeEvolve(mem, SUN_RCONST(0.1), y, &tret, FIRK_NORMAL) ==
                   FIRK_SUCCESS,
                 "one stage: the integration succeeds");
  fm = (FIRKodeMem)mem;
  fails += check(fm->g_stk == fm->G[0] && fm->delta_stk == fm->delta[0] &&
                   fm->ewt_stk == fm->ewt,
                 "one stage: no stack");
  snprintf(what, sizeof(what), "one stage: the solution error %g",
           (double)exact_error(y, tret));
  fails += check(exact_error(y, tret) <= SUN_RCONST(1.0e-2), what);

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  return fails;
}

int main(void)
{
  SUNContext sunctx = NULL;
  int fails         = 0;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  fails += check_ops(sunctx, 0);
  fails += check_ops(sunctx, 1);
  fails += check_clone_failures(sunctx);
  fails += check_one_stage(sunctx);

  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("All tests passed\n");
  return 0;
}
