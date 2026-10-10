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
 * Unit test for the FIRKLS protocol (design Listing 8.4) with a recording mock
 * around SUNLinSol_SPGMR, on y' = diag(lambda) y with lambda = -1, -10, -100,
 * -1000, a user jtimes and jtsetup and an exact block preconditioner, for the
 * one-stage method (block Newton solves) and the three-stage method (stacked
 * Newton solves, so the block solver sees only the estimate):
 *
 * - H6: jtsetup runs once per Newton iteration that reaches a solver and once
 *   per estimate solve with this matrix-free block solver: for s = 1 once per
 *   block solver call (the shortcut of design Listing 8.4, which occurs,
 *   takes neither), for s = 3 njtsetup = nni + (estimate solves); with a
 *   dense block solver njtsetup = nni for s = 1 and 3, none for an estimate;
 * - H11: SUNLinSolSetZeroGuess(.., SUNTRUE) precedes every block and mass
 *   solve (firk_test_stack_ops checks the stacked solves);
 * - H17: the estimate solve is never skipped: at least one per step;
 * - H18: SUNLS_RES_REDUCED is accepted at Newton iteration 1 and is a
 *   recoverable failure at iteration 2;
 * - H29: the user's psetup and linsys receive jok = !jbad and whatever they
 *   report in jcur is recorded, under both reuse policies: after a Newton
 *   failure, Jacobian data reported current is reused (jok = SUNTRUE) and
 *   data reported not recomputed is requested again (jok = SUNFALSE); every
 *   psetup counts in npe, every linsys reporting jcur in nje;
 * - linear-solution scaling: with a dense solver the s = 1 Newton correction
 *   from the matrix of the last setup is scaled by 2 / (1 + gamma / gammap),
 *   the estimate solve is not, and FIRKodeSetLinearSolutionScaling turns it
 *   off and rejects a matrix-free solver;
 * - the solve-flag mapping of Listing 8.4: injected recoverable codes fail a
 *   Newton or estimate solve recoverably, unrecoverable ones (including
 *   SUN_ERR_EXT_FAIL, R8) end the integration with FIRK_LSOLVE_FAIL, and a
 *   preconditioner block (s = 3, dense) accepts a recoverable failure and
 *   counts it in nblockfails;
 * - the stacked preconditioner with a matrix-free block solver (s = 3): only
 *   the user's psolve runs, s times per stage iteration, never the user's
 *   Krylov iteration (H37), with the tolerance eps_L ||r_i||_WRMS nrmfac, on
 *   the sides of the user's pretype (BOTH: left, then right; PCG: left;
 *   NONE: no preconditioner); a recoverable psolve failure is accepted and
 *   counted, an unrecoverable one ends the integration;
 * - the stacked operator: J v from the saved Jacobian before a user jtimes,
 *   SUNMatMatvecSetup once per Jacobian evaluation, the difference-quotient
 *   J v of a zero block is 0 without an RHS call (H7), the stacked tolerance
 *   uses the nrmfac of FIRKodeSetLSNormFactor, and a matrix-free estimate
 *   uses gammap, the gamma of the last setup (PERIODIC policy);
 * - R7: a matrix without SUNMatMatvec and no jtimes: the integration runs
 *   with difference-quotient products in the stacked operator;
 * - H2: a recoverable RHS failure in the dense difference-quotient Jacobian,
 *   which perturbs y in place, and in the band one, which perturbs a copy,
 *   leaves the solution unchanged and the step is retried;
 * - H4: stage-solver setters between Evolve calls (SetStageSolverMaxl,
 *   SetStageSolverGSType, SetMode restoring the defaults) take effect at the
 *   next stacked solve without freeing a live solver (the same object stays
 *   attached), SetStageSolverMaxRestarts edits the live solver; the
 *   Gram-Schmidt type is validated and is an EXPERT option;
 * - input validation: a rejected SetPreconditioner changes nothing; a user
 *   linsys set between Evolve calls drops the saved Jacobian (stale for the
 *   stacked operator), and the default linsys set back saves one again,
 *   evaluated before its first use; one SUNLinearSolver object for both the
 *   system and the mass matrix is FIRKLS_ILL_INPUT in either attach order,
 *   leaving the solver already attached in place;
 * - purpose-aware mass solves (M = I, a mock PCG): a reduced residual is
 *   accepted for a derivative of y (the initial step) and is a recoverable
 *   failure of the enclosing solve in the units conversion; the other codes
 *   map as in the block solves.
 *
 * Pass criterion: every count matches.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>

#include "../firk_test_mocks.h"
#include "firkode/firkode.h"
#include "firkode/firkode_impl.h" /* firkSetNumStagesAny */
#include "firkode/firkode_ls.h"
#include "firkode/firkode_ls_impl.h" /* nblockfails */
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_band.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunlinsol/sunlinsol_pcg.h"
#include "sunlinsol/sunlinsol_spfgmr.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunmatrix/sunmatrix_band.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define NEQ  4

static const sunrealtype lambda[NEQ] = {SUN_RCONST(-1.0), SUN_RCONST(-10.0),
                                        SUN_RCONST(-100.0), SUN_RCONST(-1000.0)};
static long int njtsetup_calls       = 0;

/* the record of the user's psetup or linsys calls (H29) */
#define NREC 8
static long int nsetup_calls = 0;
static sunbooleantype setup_jok[NREC];
static sunbooleantype report_jcur = SUNTRUE; /* the jcur reported */
static sunrealtype psetup_gamma   = 0;       /* gamma of the last psetup */

/* the record of the user's psolve calls; a call outside a solve of the
   (mock) user solver comes from the stacked preconditioner (H37) */
static FirkMockLSContent* psolve_mock = NULL; /* the mock of the run */
static void* psolve_mem               = NULL; /* set to check delta */
static N_Vector psolve_w              = NULL; /* error weights      */
static long int npsolve_stacked       = 0;
static long int npsolve_inner         = 0;
static int psolve_lr[2 * NREC];          /* lr of the stacked calls */
static sunrealtype psolve_delta_err = 0; /* relative error of delta */
static long int psolve_fail_at      = 0; /* call (1-based) to fail  */
static int psolve_fail_value        = 0;

/* the psolve calls of estimate solves (inside a solve of the mock user
   solver, s = 3) */
static void* gamma_mem          = NULL; /* set to check gamma         */
static long int nest_gamma_bad  = 0;    /* gamma != psetup's gamma    */
static long int nstk_gamma_bad  = 0;    /* stacked: gamma != psetup's */
static long int nstk_gamma_diff = 0;    /* stacked calls, gamma != gammap */
static long int nest_gamma_diff = 0;    /* calls with gamma != gammap */

static void psolve_reset(void)
{
  psolve_mock      = NULL;
  psolve_mem       = NULL;
  npsolve_stacked  = 0;
  npsolve_inner    = 0;
  psolve_delta_err = ZERO;
  psolve_fail_at   = 0;
  gamma_mem        = NULL;
  nest_gamma_bad   = 0;
  nstk_gamma_bad   = 0;
  nstk_gamma_diff  = 0;
  nest_gamma_diff  = 0;
}

/* H2: the RHS recognizes the calls of the difference-quotient Jacobians
   that perturb one component, against a snapshot of yn taken by the first
   call at that time: the dense one perturbs the integrator's yn in place,
   the band one a copy, which is the temporary vector h2_ytemp */
static void* h2_mem           = NULL; /* set to watch yn       */
static N_Vector h2_ytemp      = NULL; /* the band copy of yn   */
static N_Vector h2_snap       = NULL; /* yn at the time h2_t   */
static sunrealtype h2_t       = 0;
static sunbooleantype h2_have = SUNFALSE;
static long int h2_dq_calls   = 0; /* calls with one perturbed component */
static long int h2_fail_at    = 0; /* DQ call (1-based) to fail          */
static long int h2_failed     = 0;
static long int h2_corrupt    = 0; /* calls with yn changed elsewhere    */

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  int i, ndiff = 0;
  (void)user_data;
  for (i = 0; i < NEQ; i++) { NV_Ith_S(ydot, i) = lambda[i] * NV_Ith_S(y, i); }

  if (h2_mem != NULL && (y == ((FIRKodeMem)h2_mem)->yn || y == h2_ytemp))
  {
    const sunbooleantype is_yn = (y == ((FIRKodeMem)h2_mem)->yn);
    if (h2_have && t == h2_t)
    {
      for (i = 0; i < NEQ; i++)
      {
        if (NV_Ith_S(y, i) != NV_Ith_S(h2_snap, i)) { ndiff++; }
      }
    }
    if (is_yn && (!h2_have || t != h2_t || ndiff == 0))
    {
      N_VScale(ONE, y, h2_snap);
      h2_t    = t;
      h2_have = SUNTRUE;
    }
    else if (h2_have && t == h2_t && ndiff == 1)
    {
      h2_dq_calls++;
      if (h2_dq_calls == h2_fail_at)
      {
        h2_failed++;
        return 1;
      }
    }
    else if (is_yn) { h2_corrupt++; }
  }
  return 0;
}

static int jtsetup(sunrealtype t, N_Vector y, N_Vector fy, void* user_data)
{
  (void)t;
  (void)y;
  (void)fy;
  (void)user_data;
  njtsetup_calls++;
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

static int psetup(sunrealtype t, N_Vector y, N_Vector fy, sunbooleantype jok,
                  sunbooleantype* jcur, sunrealtype gamma, void* user_data)
{
  (void)t;
  (void)y;
  (void)fy;
  (void)user_data;
  if (nsetup_calls < NREC) { setup_jok[nsetup_calls] = jok; }
  nsetup_calls++;
  psetup_gamma = gamma;
  *jcur        = report_jcur;
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

/* A = I - gamma diag(lambda), recorded as psetup */
static int linsys(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix A,
                  SUNMatrix M, sunbooleantype jok, sunbooleantype* jcur,
                  sunrealtype gamma, void* user_data, N_Vector tmp1,
                  N_Vector tmp2, N_Vector tmp3)
{
  int i;
  (void)t;
  (void)y;
  (void)fy;
  (void)M;
  (void)user_data;
  (void)tmp1;
  (void)tmp2;
  (void)tmp3;
  SUNMatZero(A);
  for (i = 0; i < NEQ; i++) { SM_ELEMENT_D(A, i, i) = ONE - gamma * lambda[i]; }
  if (nsetup_calls < NREC) { setup_jok[nsetup_calls] = jok; }
  nsetup_calls++;
  *jcur = report_jcur;
  return 0;
}

/* the exact solve with I - gamma diag(lambda) */
static int psolve(sunrealtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z,
                  sunrealtype gamma, sunrealtype delta, int lr, void* user_data)
{
  int i;
  long int ncall;
  sunrealtype expected;
  (void)t;
  (void)y;
  (void)fy;
  (void)user_data;

  /* a stacked call: its side, its tolerance eps_L ||r||_WRMS nrmfac, and
     the gamma of the last psetup */
  if (psolve_mock == NULL || !psolve_mock->in_solve)
  {
    if (npsolve_stacked < 2 * NREC) { psolve_lr[npsolve_stacked] = lr; }
    npsolve_stacked++;
    if (gamma_mem != NULL)
    {
      if (gamma != psetup_gamma) { nstk_gamma_bad++; }
      if (((FIRKodeMem)gamma_mem)->gamma != ((FIRKodeMem)gamma_mem)->gammap)
      {
        nstk_gamma_diff++;
      }
    }
    if (psolve_mem != NULL)
    {
      FIRKodeGetErrWeights(psolve_mem, psolve_w);
      expected = SUN_RCONST(0.05) * N_VWrmsNorm(r, psolve_w) *
                 SUNRsqrt((sunrealtype)NEQ);
      if (expected > ZERO)
      {
        psolve_delta_err = SUNMAX(psolve_delta_err,
                                  SUNRabs(delta - expected) / expected);
      }
    }
  }
  else
  {
    npsolve_inner++;
    /* an estimate solve uses the gamma of the last setup */
    if (gamma_mem != NULL)
    {
      if (gamma != psetup_gamma) { nest_gamma_bad++; }
      if (((FIRKodeMem)gamma_mem)->gamma != ((FIRKodeMem)gamma_mem)->gammap)
      {
        nest_gamma_diff++;
      }
    }
  }
  ncall = npsolve_stacked + npsolve_inner;

  for (i = 0; i < NEQ; i++)
  {
    NV_Ith_S(z, i) = NV_Ith_S(r, i) / (ONE - gamma * lambda[i]);
  }
  return (ncall == psolve_fail_at) ? psolve_fail_value : 0;
}

static int check(int ok, const char* what, int s)
{
  if (!ok)
  {
    printf("  FAIL (s = %d): %s\n", s, what);
    return 1;
  }
  return 0;
}

struct run
{
  void* mem;
  N_Vector y;
  SUNLinearSolver LS;
  SUNMatrix A;
  FirkMockLSContent* mock;
};

static int setup(SUNContext sunctx, struct run* R, int s, sunrealtype rtol)
{
  int i;
  R->y = N_VNew_Serial(NEQ, sunctx);
  for (i = 0; i < NEQ; i++) { NV_Ith_S(R->y, i) = ONE; }
  R->LS = FirkMockLS_New(SUNLinSol_SPGMR(R->y, SUN_PREC_LEFT, 0, sunctx), sunctx);
  if (R->LS == NULL) { return 1; }
  R->mock        = FirkMockLS_Content(R->LS);
  R->A           = NULL;
  R->mem         = FIRKodeCreate(sunctx);
  njtsetup_calls = 0;
  nsetup_calls   = 0;
  psolve_reset();
  psolve_mock = R->mock;
  return FIRKodeInit(R->mem, f, ZERO, R->y) ||
         FIRKodeSStolerances(R->mem, rtol, rtol) ||
         FIRKodeSetLinearSolver(R->mem, R->LS, NULL) ||
         FIRKodeSetJacTimes(R->mem, jtsetup, jtimes) ||
         FIRKodeSetPreconditioner(R->mem, psetup, psolve) ||
         ((s != FIRK_DEFAULT_STAGES) && firkSetNumStagesAny(R->mem, s));
}

static void teardown(struct run* R)
{
  psolve_reset();
  FIRKodeFree(&R->mem);
  SUNLinSolFree(R->LS);
  SUNMatDestroy(R->A);
  N_VDestroy(R->y);
}

/* H6, H11, H17 over an integration */
static int check_counts(SUNContext sunctx, int s)
{
  int fails    = 0;
  long int nni = 0, nst = 0, njtsetup = 0, nbsolves = 0, nest;
  sunrealtype tret;
  char what[200];
  struct run R;

  /* the one-stage method (order 1) at a looser tolerance */
  if (setup(sunctx, &R, s, (s == 1) ? SUN_RCONST(1.0e-3) : SUN_RCONST(1.0e-6)))
  {
    return 1;
  }
  fails += check(FIRKodeEvolve(R.mem, SUN_RCONST(0.1), R.y, &tret,
                               FIRK_NORMAL) == FIRK_SUCCESS,
                 "the integration succeeds", s);
  FIRKodeGetNumNonlinSolvIters(R.mem, &nni);
  FIRKodeGetNumSteps(R.mem, &nst);
  FIRKodeGetNumJTSetupEvals(R.mem, &njtsetup);
  FIRKodeGetNumBlockSolves(R.mem, &nbsolves);

  if (s == 1)
  {
    /* every Newton iteration and estimate is one block solve, counted with
       the shortcuts; jtsetup precedes exactly the solver calls */
    nest = nbsolves - nni;
    snprintf(what, sizeof(what),
             "H6: njtsetup %ld = block solver calls %ld (nni %ld, estimate "
             "solves %ld, calls %ld)",
             njtsetup, R.mock->nsolves, nni, nest, njtsetup_calls);
    fails += check(njtsetup == R.mock->nsolves && njtsetup_calls == njtsetup,
                   what, s);
    snprintf(what, sizeof(what),
             "H6: the k >= 2 shortcut skips the solver (%ld block solves, %ld "
             "solver calls)",
             nbsolves, R.mock->nsolves);
    fails += check(nbsolves > R.mock->nsolves, what, s);
  }
  else
  {
    /* the block solver sees only the estimates; the preconditioner applies
       the user's psolve directly */
    nest = R.mock->nsolves;
    snprintf(what, sizeof(what),
             "H6: njtsetup %ld = nni %ld + estimate solves %ld (calls %ld)",
             njtsetup, nni, nest, njtsetup_calls);
    fails += check(njtsetup == nni + nest && njtsetup_calls == njtsetup, what, s);
  }
  fails += check(R.mock->nsolves_no_zero == 0,
                 "H11: a zero guess precedes every block solve", s);
  snprintf(what, sizeof(what), "H17: estimate solves %ld >= steps %ld", nest,
           nst);
  fails += check(nest >= nst, what, s);

  teardown(&R);
  return fails;
}

/* H6: with a matrix-based block solver (a dense mock) and the user's jtsetup
   and jtimes, jtsetup runs once per Newton iteration and never for the
   estimate solves */
static int check_jtsetup_matrixbased(SUNContext sunctx, int s)
{
  int i, fails = 0;
  long int nni = 0, njtsetup = 0, nst = 0;
  sunrealtype tret, rtol;
  char what[200];
  struct run R;

  R.y = N_VNew_Serial(NEQ, sunctx);
  for (i = 0; i < NEQ; i++) { NV_Ith_S(R.y, i) = ONE; }
  R.A  = SUNDenseMatrix(NEQ, NEQ, sunctx);
  R.LS = FirkMockLS_New(SUNLinSol_Dense(R.y, R.A, sunctx), sunctx);
  if (R.LS == NULL) { return 1; }
  R.mock         = FirkMockLS_Content(R.LS);
  R.mem          = FIRKodeCreate(sunctx);
  njtsetup_calls = 0;
  rtol           = (s == 1) ? SUN_RCONST(1.0e-3) : SUN_RCONST(1.0e-6);
  if (FIRKodeInit(R.mem, f, ZERO, R.y) || FIRKodeSStolerances(R.mem, rtol, rtol) ||
      FIRKodeSetLinearSolver(R.mem, R.LS, R.A) ||
      FIRKodeSetJacTimes(R.mem, jtsetup, jtimes) ||
      ((s != FIRK_DEFAULT_STAGES) && firkSetNumStagesAny(R.mem, s)))
  {
    return 1;
  }
  fails += check(FIRKodeEvolve(R.mem, SUN_RCONST(0.1), R.y, &tret,
                               FIRK_NORMAL) == FIRK_SUCCESS,
                 "matrix-based: the integration succeeds", s);
  FIRKodeGetNumNonlinSolvIters(R.mem, &nni);
  FIRKodeGetNumJTSetupEvals(R.mem, &njtsetup);
  FIRKodeGetNumSteps(R.mem, &nst);
  /* at least one estimate solve per step (H17) */
  snprintf(what, sizeof(what),
           "H6: matrix-based: njtsetup %ld = nni %ld, none for the estimate "
           "solves of %ld steps (calls %ld)",
           njtsetup, nni, nst, njtsetup_calls);
  fails += check(nst > 0 && njtsetup == nni && njtsetup_calls == njtsetup, what,
                 s);
  teardown(&R);
  return fails;
}

/* H18: an injected SUNLS_RES_REDUCED at the first or second Newton iteration
   of the first step (s = 1: one block solve per iteration) */
static int check_res_reduced(SUNContext sunctx, long int which,
                             sunbooleantype accepted)
{
  int fails     = 0;
  long int ncfn = -1, nni = 0;
  sunrealtype tret;
  struct run R;

  /* tolerances at which the first Newton solve needs a second iteration */
  if (setup(sunctx, &R, 1, SUN_RCONST(1.0e-5))) { return 1; }
  R.mock->inject_first = which;
  R.mock->inject_last  = which;
  R.mock->inject_value = SUNLS_RES_REDUCED;

  fails +=
    check(FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP) == FIRK_SUCCESS,
          "the first step succeeds", 1);
  FIRKodeGetNumStepSolveFails(R.mem, &ncfn);
  FIRKodeGetNumNonlinSolvIters(R.mem, &nni);
  fails += check(nni >= 2, "the first Newton solve needs two iterations", 1);
  if (accepted)
  {
    fails += check(ncfn == 0, "H18: RES_REDUCED accepted at iteration 1", 1);
  }
  else
  {
    fails += check(ncfn == 1, "H18: RES_REDUCED a failure at iteration 2", 1);
  }

  teardown(&R);
  return fails;
}

/* H29: the first step with a Newton failure injected at its first solve
   (s = 1: one block solve per iteration), with the user's psetup (SPGMR) or
   linsys (dense) reporting jcur = report */
static int check_jcur(SUNContext sunctx, int use_linsys, int policy,
                      sunbooleantype report)
{
  int fails    = 0;
  long int npe = -1, nje = -1;
  sunrealtype tret;
  char what[160];
  struct run R;

  if (use_linsys)
  {
    R.y = N_VNew_Serial(NEQ, sunctx);
    N_VConst(ONE, R.y);
    R.A  = SUNDenseMatrix(NEQ, NEQ, sunctx);
    R.LS = FirkMockLS_New(SUNLinSol_Dense(R.y, R.A, sunctx), sunctx);
    if (R.LS == NULL) { return 1; }
    R.mock       = FirkMockLS_Content(R.LS);
    R.mem        = FIRKodeCreate(sunctx);
    nsetup_calls = 0;
    if (FIRKodeInit(R.mem, f, ZERO, R.y) ||
        FIRKodeSStolerances(R.mem, SUN_RCONST(1.0e-5), SUN_RCONST(1.0e-5)) ||
        FIRKodeSetLinearSolver(R.mem, R.LS, R.A) ||
        FIRKodeSetLinSysFn(R.mem, linsys) || firkSetNumStagesAny(R.mem, 1))
    {
      return 1;
    }
  }
  else if (setup(sunctx, &R, 1, SUN_RCONST(1.0e-5))) { return 1; }

  if (FIRKodeSetMode(R.mem, FIRK_EXPERT) || FIRKodeSetReusePolicy(R.mem, policy))
  {
    return 1;
  }
  report_jcur          = report;
  R.mock->inject_first = 1;
  R.mock->inject_last  = 1;
  R.mock->inject_value = SUNLS_CONV_FAIL;

  fails +=
    check(FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP) == FIRK_SUCCESS,
          "H29: the first step succeeds", 1);
  FIRKodeGetNumPrecEvals(R.mem, &npe);
  FIRKodeGetNumJacEvals(R.mem, &nje);

  snprintf(what, sizeof(what),
           "H29 (%s, policy %d, jcur reported %d): two setups, the first "
           "with jok = 0 (calls %ld)",
           use_linsys ? "linsys" : "psetup", policy, (int)report, nsetup_calls);
  fails += check(nsetup_calls >= 2 && !setup_jok[0], what, 1);
  if (nsetup_calls >= 2)
  {
    snprintf(what, sizeof(what),
             "H29 (%s, policy %d, jcur reported %d): the retry has jok = %d",
             use_linsys ? "linsys" : "psetup", policy, (int)report,
             (int)setup_jok[1]);
    fails += check(setup_jok[1] == report, what, 1);
  }
  snprintf(what, sizeof(what),
           "H29 (%s, policy %d, jcur reported %d): npe %ld, nje %ld for %ld "
           "calls",
           use_linsys ? "linsys" : "psetup", policy, (int)report, npe, nje,
           nsetup_calls);
  if (use_linsys)
  {
    fails += check(nje == (report ? nsetup_calls : 0), what, 1);
  }
  else { fails += check(npe == nsetup_calls, what, 1); }

  report_jcur = SUNTRUE;
  teardown(&R);
  return fails;
}

/* Linear-solution scaling: a direct block solve with gamma = 2 gammap after
   one step of the one-stage method (internal interface) */
static int check_scalesol(SUNContext sunctx)
{
  int fails = 0, i, k;
  sunrealtype tret, gammap, expected, err;
  N_Vector b = NULL, x = NULL;
  FIRKodeMem fm;
  struct run R;
  /* scaling on, Newton: 2 / (1 + 2); off, Newton: 1; on, estimate: 1 */
  const sunbooleantype onoff[3] = {SUNTRUE, SUNFALSE, SUNTRUE};
  const int purpose[3]          = {FIRKLS_BLOCK_NEWTON, FIRKLS_BLOCK_NEWTON,
                                   FIRKLS_BLOCK_ESTIMATE};
  const sunrealtype factor[3]   = {SUN_RCONST(2.0) / SUN_RCONST(3.0), ONE, ONE};
  const char* what[3]           = {"scaled Newton correction",
                                   "Newton correction without scaling",
                                   "estimate without scaling"};

  R.y = N_VNew_Serial(NEQ, sunctx);
  N_VConst(ONE, R.y);
  R.A   = SUNDenseMatrix(NEQ, NEQ, sunctx);
  R.LS  = SUNLinSol_Dense(R.y, R.A, sunctx);
  R.mem = FIRKodeCreate(sunctx);
  b     = N_VClone(R.y);
  x     = N_VClone(R.y);
  if (FIRKodeInit(R.mem, f, ZERO, R.y) ||
      FIRKodeSStolerances(R.mem, SUN_RCONST(1.0e-3), SUN_RCONST(1.0e-3)) ||
      FIRKodeSetLinearSolver(R.mem, R.LS, R.A) || FIRKodeSetJacFn(R.mem, jac) ||
      firkSetNumStagesAny(R.mem, 1))
  {
    return 1;
  }
  fails +=
    check(FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP) == FIRK_SUCCESS,
          "scaling: the first step succeeds", 1);

  /* a solve with the matrix I - gammap J of the last setup at gamma = 2
     gammap; b = 1 */
  fm     = (FIRKodeMem)R.mem;
  gammap = fm->gammap;
  for (k = 0; k < 3; k++)
  {
    fails +=
      check(FIRKodeSetLinearSolutionScaling(R.mem, onoff[k]) == FIRKLS_SUCCESS,
            "scaling: the setter accepts a dense solver in BEGINNER mode", 1);
    fm->gamma = SUN_RCONST(2.0) * gammap;
    N_VConst(ONE, b);
    fails += check(fm->lsolve_blk(fm, b, x, purpose[k], 1) == 0,
                   "scaling: the block solve succeeds", 1);
    err = ZERO;
    for (i = 0; i < NEQ; i++)
    {
      expected = factor[k] / (ONE - gammap * lambda[i]);
      err = SUNMAX(err, SUNRabs(NV_Ith_S(x, i) - expected) / SUNRabs(expected));
    }
    fails += check(err <= SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF, what[k], 1);
  }
  fm->gamma = gammap;
  teardown(&R);
  N_VDestroy(b);
  N_VDestroy(x);

  /* a matrix-free solver has no matrix of the last setup */
  if (setup(sunctx, &R, 1, SUN_RCONST(1.0e-3))) { return fails + 1; }
  fails +=
    check(FIRKodeSetLinearSolutionScaling(R.mem, SUNTRUE) == FIRKLS_ILL_INPUT,
          "scaling: a matrix-free solver is rejected", 1);
  teardown(&R);
  return fails;
}

/* The solve-flag mapping: the first step with code value injected at the
   first solve of the (mock) block solver: s = 1 with SPGMR (a Newton solve),
   s = 3 with SPGMR (the estimate), s = 3 with a dense solver (a
   preconditioner block) */
struct flag_case
{
  int s;
  int dense;
  int value;
  int expect_flag;
  long int expect_ncfn;
  long int expect_blockfails;
  const char* what;
};

static int check_flags(SUNContext sunctx)
{
  static const struct flag_case cases[] =
    {{1, 0, SUNLS_PSET_FAIL_REC, FIRK_SUCCESS, 1, 0,
      "PSET_FAIL_REC fails a Newton solve recoverably"},
     {1, 0, SUNLS_RECOV_FAILURE, FIRK_SUCCESS, 1, 0,
      "RECOV_FAILURE fails a Newton solve recoverably"},
     {1, 0, SUN_ERR_EXT_FAIL, FIRK_LSOLVE_FAIL, 1, 0,
      "SUN_ERR_EXT_FAIL in a Newton solve is unrecoverable (R8)"},
     {1, 0, SUNLS_PSET_FAIL_UNREC, FIRK_LSOLVE_FAIL, 1, 0,
      "PSET_FAIL_UNREC in a Newton solve is unrecoverable"},
     {3, 0, SUNLS_RES_REDUCED, FIRK_SUCCESS, 1, 0,
      "RES_REDUCED fails the estimate recoverably"},
     {3, 1, SUNLS_CONV_FAIL, FIRK_SUCCESS, 0, 1,
      "CONV_FAIL in a preconditioner block is accepted and counted"},
     {3, 1, SUN_ERR_EXT_FAIL, FIRK_LSOLVE_FAIL, 1, 0,
      "SUN_ERR_EXT_FAIL in a preconditioner block is unrecoverable"}};
  int fails = 0, c, flag;
  long int ncfn;
  sunrealtype tret;
  char what[200];
  struct run R;

  for (c = 0; c < (int)(sizeof(cases) / sizeof(cases[0])); c++)
  {
    if (cases[c].dense)
    {
      R.y = N_VNew_Serial(NEQ, sunctx);
      N_VConst(ONE, R.y);
      R.A  = SUNDenseMatrix(NEQ, NEQ, sunctx);
      R.LS = FirkMockLS_New(SUNLinSol_Dense(R.y, R.A, sunctx), sunctx);
      if (R.LS == NULL) { return fails + 1; }
      R.mock = FirkMockLS_Content(R.LS);
      R.mem  = FIRKodeCreate(sunctx);
      if (FIRKodeInit(R.mem, f, ZERO, R.y) ||
          FIRKodeSStolerances(R.mem, SUN_RCONST(1.0e-5), SUN_RCONST(1.0e-5)) ||
          FIRKodeSetLinearSolver(R.mem, R.LS, R.A) || FIRKodeSetJacFn(R.mem, jac))
      {
        return fails + 1;
      }
    }
    else if (setup(sunctx, &R, cases[c].s, SUN_RCONST(1.0e-5)))
    {
      return fails + 1;
    }
    R.mock->inject_first = 1;
    R.mock->inject_last  = 1;
    R.mock->inject_value = cases[c].value;

    flag = FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP);
    ncfn = -1;
    FIRKodeGetNumStepSolveFails(R.mem, &ncfn);
    snprintf(what, sizeof(what),
             "%s: flag %d (expected %d), ncfn %ld, nblockfails %ld",
             cases[c].what, flag, cases[c].expect_flag, ncfn,
             ((FIRKLsMem)((FIRKodeMem)R.mem)->lmem)->nblockfails);
    fails += check(flag == cases[c].expect_flag && ncfn == cases[c].expect_ncfn &&
                     ((FIRKLsMem)((FIRKodeMem)R.mem)->lmem)->nblockfails ==
                       cases[c].expect_blockfails,
                   what, cases[c].s);
    teardown(&R);
  }
  return fails;
}

/* H37 and the tolerance of the stacked preconditioner over an integration
   with the (mock) SPGMR, left preconditioned, s = 3 */
static int check_stacked_psolve(SUNContext sunctx)
{
  int fails     = 0;
  long int nsli = 0, nps = 0, nli = 0, k;
  sunrealtype tret;
  char what[200];
  struct run R;
  sunbooleantype all_left = SUNTRUE;

  if (setup(sunctx, &R, 3, SUN_RCONST(1.0e-6))) { return 1; }
  psolve_w   = N_VClone(R.y);
  psolve_mem = R.mem;
  fails += check(FIRKodeEvolve(R.mem, SUN_RCONST(0.1), R.y, &tret,
                               FIRK_NORMAL) == FIRK_SUCCESS,
                 "H37: the integration succeeds", 3);
  FIRKodeGetNumStageLinIters(R.mem, &nsli);
  FIRKodeGetNumPrecSolves(R.mem, &nps);
  FIRKodeGetNumLinIters(R.mem, &nli);
  for (k = 0; k < SUNMIN(npsolve_stacked, 2 * NREC); k++)
  {
    all_left = all_left && (psolve_lr[k] == SUN_PREC_LEFT);
  }

  snprintf(what, sizeof(what),
           "H37: stacked psolve calls %ld = s x stage iterations %ld",
           npsolve_stacked, nsli);
  fails += check(nsli > 0 && npsolve_stacked == 3 * nsli, what, 3);
  snprintf(what, sizeof(what),
           "H37: nps %ld = stacked %ld + user solver %ld psolve calls", nps,
           npsolve_stacked, npsolve_inner);
  fails += check(nps == npsolve_stacked + npsolve_inner, what, 3);
  snprintf(what, sizeof(what),
           "H37: block iterations %ld only from the user solver's solves (%ld)",
           nli, R.mock->niters);
  fails += check(nli == R.mock->niters, what, 3);
  fails += check(all_left, "stacked psolve calls on the user's (left) side", 3);
  snprintf(what, sizeof(what),
           "stacked psolve tolerance eps_L ||r_i|| nrmfac, relative error %g",
           (double)psolve_delta_err);
  fails += check(psolve_delta_err <= SUN_RCONST(10.0) * SUN_UNIT_ROUNDOFF, what,
                 3);

  N_VDestroy(psolve_w);
  psolve_w = NULL;
  teardown(&R);
  return fails;
}

/* The sides of the stacked psolve calls for the user's pretype: the first
   stacked solve of the first step (before any estimate) */
static int check_pretypes(SUNContext sunctx)
{
  int fails = 0, c, i, n;
  sunrealtype tret;
  char what[200];
  struct run R;

  struct
  {
    const char* name;
    int solver; /* 0 SPGMR, 1 SPFGMR, 2 PCG */
    int pretype;
    int nsides; /* psolve calls per block, 0 for none */
    int first;  /* side of the first call of a block  */
  } cases[] = {{"SPGMR right", 0, SUN_PREC_RIGHT, 1, SUN_PREC_RIGHT},
               {"SPGMR both", 0, SUN_PREC_BOTH, 2, SUN_PREC_LEFT},
               {"SPGMR none", 0, SUN_PREC_NONE, 0, 0},
               {"SPFGMR left (right only)", 1, SUN_PREC_LEFT, 1, SUN_PREC_RIGHT},
               {"PCG right (left only)", 2, SUN_PREC_RIGHT, 1, SUN_PREC_LEFT}};

  for (c = 0; c < (int)(sizeof(cases) / sizeof(cases[0])); c++)
  {
    R.y = N_VNew_Serial(NEQ, sunctx);
    N_VConst(ONE, R.y);
    R.A    = NULL;
    R.LS   = (cases[c].solver == 0)
               ? SUNLinSol_SPGMR(R.y, cases[c].pretype, 0, sunctx)
             : (cases[c].solver == 1)
               ? SUNLinSol_SPFGMR(R.y, cases[c].pretype, 0, sunctx)
               : SUNLinSol_PCG(R.y, cases[c].pretype, 0, sunctx);
    R.mock = NULL;
    R.mem  = FIRKodeCreate(sunctx);
    psolve_reset();
    if (FIRKodeInit(R.mem, f, ZERO, R.y) ||
        FIRKodeSStolerances(R.mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-6)) ||
        FIRKodeSetLinearSolver(R.mem, R.LS, NULL) ||
        FIRKodeSetJacTimes(R.mem, jtsetup, jtimes) ||
        FIRKodeSetPreconditioner(R.mem, psetup, psolve))
    {
      return fails + 1;
    }
    fails += check(FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP) ==
                     FIRK_SUCCESS,
                   "pretype: the first step succeeds", 3);

    n = 3 * cases[c].nsides;
    if (n == 0)
    {
      snprintf(what, sizeof(what), "pretype %s: no psolve calls (%ld)",
               cases[c].name, npsolve_stacked);
      fails += check(npsolve_stacked == 0, what, 3);
    }
    else
    {
      int ok = (npsolve_stacked >= n);
      for (i = 0; ok && i < n; i++)
      {
        ok = (psolve_lr[i] == ((cases[c].nsides == 2 && i % 2 == 1)
                                 ? SUN_PREC_RIGHT
                                 : cases[c].first));
      }
      snprintf(what, sizeof(what),
               "pretype %s: the sides of the first stacked psolve calls",
               cases[c].name);
      fails += check(ok, what, 3);
    }
    teardown(&R);
  }
  return fails;
}

/* a recoverable psolve failure in the stacked preconditioner is accepted and
   counted; an unrecoverable one ends the integration */
static int check_psolve_failures(SUNContext sunctx)
{
  int fails = 0, flag;
  long int ncfn;
  sunrealtype tret;
  char what[200];
  struct run R;

  if (setup(sunctx, &R, 3, SUN_RCONST(1.0e-6))) { return 1; }
  psolve_fail_at    = 1;
  psolve_fail_value = 1;
  flag              = FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP);
  ncfn              = -1;
  FIRKodeGetNumStepSolveFails(R.mem, &ncfn);
  snprintf(what, sizeof(what),
           "a recoverable psolve failure is accepted: flag %d, ncfn %ld, "
           "nblockfails %ld",
           flag, ncfn, ((FIRKLsMem)((FIRKodeMem)R.mem)->lmem)->nblockfails);
  fails += check(flag == FIRK_SUCCESS && ncfn == 0 &&
                   ((FIRKLsMem)((FIRKodeMem)R.mem)->lmem)->nblockfails == 1,
                 what, 3);
  teardown(&R);

  if (setup(sunctx, &R, 3, SUN_RCONST(1.0e-6))) { return fails + 1; }
  psolve_fail_at    = 1;
  psolve_fail_value = -1;
  flag              = FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP);
  snprintf(what, sizeof(what),
           "an unrecoverable psolve failure ends the integration: flag %d", flag);
  fails += check(flag == FIRK_LSOLVE_FAIL, what, 3);
  teardown(&R);
  return fails;
}

/* a user jtimes that counts its calls */
static long int njtimes_user = 0;

static int jtimes_counting(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y,
                           N_Vector fy, void* user_data, N_Vector tmp)
{
  njtimes_user++;
  return jtimes(v, Jv, t, y, fy, user_data, tmp);
}

/* a SUNMatMatvecSetup that counts its calls */
static long int nmatvecsetup = 0;

static SUNErrCode matvecsetup_counting(SUNMatrix A)
{
  (void)A;
  nmatvecsetup++;
  return SUN_SUCCESS;
}

/* The stacked operator (s = 3) */
static int check_stacked_operator(SUNContext sunctx)
{
  int fails    = 0, i, j, k, ok;
  long int nje = 0, nje0 = 0, nfeDQ0 = 0, nfeDQ1 = 0, nsli = -1;
  sunbooleantype jcur;
  sunrealtype tret, h, zexp, zmax, err;
  N_Vector v = NULL, z = NULL, Jv = NULL, zero = NULL;
  FIRKodeMem fm;
  char what[200];
  struct run R;

  /* a dense solver with a user jtimes: J v from the saved Jacobian, whose
     products are set up once per Jacobian evaluation */
  R.y = N_VNew_Serial(NEQ, sunctx);
  N_VConst(ONE, R.y);
  R.A   = SUNDenseMatrix(NEQ, NEQ, sunctx);
  R.LS  = SUNLinSol_Dense(R.y, R.A, sunctx);
  R.mem = FIRKodeCreate(sunctx);
  if (FIRKodeInit(R.mem, f, ZERO, R.y) ||
      FIRKodeSStolerances(R.mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-6)) ||
      FIRKodeSetLinearSolver(R.mem, R.LS, R.A) || FIRKodeSetJacFn(R.mem, jac) ||
      FIRKodeSetJacTimes(R.mem, NULL, jtimes_counting))
  {
    return 1;
  }
  njtimes_user = 0;
  fails += check(FIRKodeEvolve(R.mem, SUN_RCONST(0.1), R.y, &tret,
                               FIRK_NORMAL) == FIRK_SUCCESS,
                 "operator: the integration succeeds", 3);
  snprintf(what, sizeof(what),
           "J v from the saved Jacobian, not the user's jtimes (%ld calls)",
           njtimes_user);
  fails += check(njtimes_user == 0, what, 3);

  /* setups with a new Jacobian, twice, and with the saved one */
  fm                                              = (FIRKodeMem)R.mem;
  ((FIRKLsMem)fm->lmem)->savedJ->ops->matvecsetup = matvecsetup_counting;
  nmatvecsetup                                    = 0;
  FIRKodeGetNumJacEvals(R.mem, &nje0);
  fails += check(fm->lsetup(fm, SUNTRUE, &jcur) == 0 &&
                   fm->lsetup(fm, SUNTRUE, &jcur) == 0 &&
                   fm->lsetup(fm, SUNFALSE, &jcur) == 0,
                 "operator: the setups succeed", 3);
  FIRKodeGetNumJacEvals(R.mem, &nje);
  snprintf(what, sizeof(what),
           "SUNMatMatvecSetup once per Jacobian evaluation: %ld for %ld",
           nmatvecsetup, nje - nje0);
  fails += check(nje - nje0 == 2 && nmatvecsetup == 2, what, 3);
  teardown(&R);

  /* the difference-quotient J v of a zero block (H7), with SPGMR and no
     preconditioner, after one step */
  R.y = N_VNew_Serial(NEQ, sunctx);
  N_VConst(ONE, R.y);
  R.A   = NULL;
  R.LS  = SUNLinSol_SPGMR(R.y, SUN_PREC_NONE, 0, sunctx);
  R.mem = FIRKodeCreate(sunctx);
  if (FIRKodeInit(R.mem, f, ZERO, R.y) ||
      FIRKodeSStolerances(R.mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-6)) ||
      FIRKodeSetLinearSolver(R.mem, R.LS, NULL))
  {
    return fails + 1;
  }
  fails +=
    check(FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP) == FIRK_SUCCESS,
          "H7: the first step succeeds", 3);
  fm   = (FIRKodeMem)R.mem;
  Jv   = N_VClone(R.y);
  zero = N_VClone(R.y);
  v    = N_VClone(fm->g_stk);
  z    = N_VClone(fm->g_stk);

  N_VConst(ONE, Jv);
  N_VConst(ZERO, zero);
  FIRKodeGetNumLinRhsEvals(R.mem, &nfeDQ0);
  /* every component 0 (N_VMaxNorm would not see a NaN) */
  ok = (firkLsDQJtimes(zero, Jv, fm->tn, fm->yn, fm->fn, fm, R.y) == 0);
  for (k = 0; k < NEQ; k++) { ok = ok && (NV_Ith_S(Jv, k) == ZERO); }
  fails += check(ok, "H7: the difference-quotient J 0 is 0", 3);
  FIRKodeGetNumLinRhsEvals(R.mem, &nfeDQ1);
  fails += check(nfeDQ1 == nfeDQ0, "H7: J 0 without an RHS call", 3);

  /* z = (I - h A (x) J) v with a zero first stage block, at the current
     linearization point: f(tn, yn) is evaluated lazily, at the start of
     the next step (D12) */
  f(fm->tn, fm->yn, fm->fn, NULL);
  N_VConst(ONE, v);
  N_VConst(ZERO, firkStage(fm, v, 0));
  fails += check(firkLsStackedATimes(fm, v, z) == 0,
                 "H7: the stacked product succeeds", 3);
  h    = fm->h;
  zmax = ZERO;
  err  = ZERO;
  ok   = 1;
  for (i = 0; i < 3; i++)
  {
    for (k = 0; k < NEQ; k++)
    {
      zexp = NV_Ith_S(firkStage(fm, v, i), k);
      for (j = 0; j < 3; j++)
      {
        zexp -= h * fm->tbl->A[i * 3 + j] * lambda[k] *
                NV_Ith_S(firkStage(fm, v, j), k);
      }
      /* finite: not NaN and not infinite */
      ok = ok &&
           (NV_Ith_S(firkStage(fm, z, i), k) == NV_Ith_S(firkStage(fm, z, i), k)) &&
           (SUNRabs(NV_Ith_S(firkStage(fm, z, i), k)) <= SUN_BIG_REAL);
      err  = SUNMAX(err, SUNRabs(NV_Ith_S(firkStage(fm, z, i), k) - zexp));
      zmax = SUNMAX(zmax, SUNRabs(zexp));
    }
  }
  snprintf(what, sizeof(what),
           "H7: the stacked product with a zero block, error %g of %g",
           (double)err, (double)zmax);
  /* the difference quotient is exact for this linear f up to rounding,
     amplified by 1 / sig */
  fails += check(ok &&
                   err <= SUN_RCONST(10.0) * SUNRsqrt(SUN_UNIT_ROUNDOFF) * zmax,
                 what, 3);
  N_VDestroy(v);
  N_VDestroy(z);
  N_VDestroy(Jv);
  N_VDestroy(zero);
  teardown(&R);

  /* the stacked tolerance uses the user's nrmfac: with a huge factor FGMRES
     stops before its first iteration */
  if (setup(sunctx, &R, 3, SUN_RCONST(1.0e-6))) { return fails + 1; }
  fails +=
    check(FIRKodeSetLSNormFactor(R.mem, SUN_RCONST(1.0e12)) == FIRKLS_SUCCESS,
          "nrmfac: the setter succeeds", 3);
  (void)FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP);
  FIRKodeGetNumStageLinIters(R.mem, &nsli);
  snprintf(what, sizeof(what),
           "the stacked tolerance uses SetLSNormFactor: %ld stage iterations",
           nsli);
  fails += check(nsli == 0, what, 3);
  teardown(&R);

  /* a matrix-free estimate, and the stacked block preconditioner, with the
     gamma of the last setup */
  if (setup(sunctx, &R, 3, SUN_RCONST(1.0e-6))) { return fails + 1; }
  if (FIRKodeSetMode(R.mem, FIRK_EXPERT) ||
      FIRKodeSetReusePolicy(R.mem, FIRK_POLICY_PERIODIC))
  {
    return fails + 1;
  }
  gamma_mem = R.mem;
  fails += check(FIRKodeEvolve(R.mem, SUN_RCONST(0.1), R.y, &tret,
                               FIRK_NORMAL) == FIRK_SUCCESS,
                 "gammap: the integration succeeds", 3);
  snprintf(what, sizeof(what),
           "the estimate uses gammap: %ld psolve calls with another gamma, "
           "%ld calls while gamma != gammap",
           nest_gamma_bad, nest_gamma_diff);
  fails += check(nest_gamma_bad == 0 && nest_gamma_diff > 0, what, 3);
  snprintf(what, sizeof(what),
           "the stacked preconditioner uses gammap: %ld psolve calls with "
           "another gamma, %ld calls while gamma != gammap",
           nstk_gamma_bad, nstk_gamma_diff);
  fails += check(nstk_gamma_bad == 0 && nstk_gamma_diff > 0, what, 3);
  teardown(&R);
  return fails;
}

/* R7: a dense matrix stripped of SUNMatMatvec, no jtimes, s = 3 */
static int check_matvec_fallback(SUNContext sunctx)
{
  int fails       = 0, i;
  long int nfeDQ0 = 0, nfeDQ1 = 0;
  sunrealtype tret, err;
  char what[200];
  struct run R;

  R.y = N_VNew_Serial(NEQ, sunctx);
  N_VConst(ONE, R.y);
  R.A   = SUNDenseMatrix(NEQ, NEQ, sunctx);
  R.LS  = SUNLinSol_Dense(R.y, R.A, sunctx);
  R.mem = FIRKodeCreate(sunctx);
  if (FIRKodeInit(R.mem, f, ZERO, R.y) ||
      FIRKodeSStolerances(R.mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-6)) ||
      FIRKodeSetLinearSolver(R.mem, R.LS, R.A) || FIRKodeSetJacFn(R.mem, jac))
  {
    return 1;
  }
  R.A->ops->matvec = NULL;
  if (check(FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP) == FIRK_SUCCESS,
            "R7: a matrix without SUNMatMatvec is accepted", 3))
  {
    teardown(&R);
    return 1;
  }

  /* the saved Jacobian, a clone of A, without the product as well */
  ((FIRKLsMem)((FIRKodeMem)R.mem)->lmem)->savedJ->ops->matvec = NULL;
  FIRKodeGetNumLinRhsEvals(R.mem, &nfeDQ0);
  fails += check(FIRKodeEvolve(R.mem, SUN_RCONST(0.1), R.y, &tret,
                               FIRK_NORMAL) == FIRK_SUCCESS,
                 "R7: the integration succeeds", 3);
  FIRKodeGetNumLinRhsEvals(R.mem, &nfeDQ1);
  err = ZERO;
  for (i = 0; i < NEQ; i++)
  {
    err = SUNMAX(err, SUNRabs(NV_Ith_S(R.y, i) - SUNRexp(lambda[i] * tret)));
  }
  snprintf(what, sizeof(what),
           "R7: difference-quotient products (%ld RHS evaluations), error %g",
           nfeDQ1 - nfeDQ0, (double)err);
  fails += check(nfeDQ1 > nfeDQ0 && err <= SUN_RCONST(1.0e-4), what, 3);
  teardown(&R);
  return fails;
}

/* H2: the second difference-quotient column of the first Jacobian fails
   recoverably (dense or band solver without a Jacobian function, s = 3);
   with the full bandwidth each group of the band one is one column */
static int check_dq_restore(SUNContext sunctx, sunbooleantype band)
{
  int fails = 0, flag;
  sunrealtype tret;
  char what[200];
  struct run R;

  R.y = N_VNew_Serial(NEQ, sunctx);
  N_VConst(ONE, R.y);
  if (band)
  {
    R.A  = SUNBandMatrix(NEQ, NEQ - 1, NEQ - 1, sunctx);
    R.LS = SUNLinSol_Band(R.y, R.A, sunctx);
  }
  else
  {
    R.A  = SUNDenseMatrix(NEQ, NEQ, sunctx);
    R.LS = SUNLinSol_Dense(R.y, R.A, sunctx);
  }
  R.mem = FIRKodeCreate(sunctx);
  if (FIRKodeInit(R.mem, f, ZERO, R.y) ||
      FIRKodeSStolerances(R.mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-6)) ||
      FIRKodeSetLinearSolver(R.mem, R.LS, R.A))
  {
    return 1;
  }
  h2_snap     = N_VClone(R.y);
  h2_mem      = R.mem;
  h2_ytemp    = band ? ((FIRKodeMem)R.mem)->tempv2 : NULL;
  h2_have     = SUNFALSE;
  h2_dq_calls = 0;
  h2_failed   = 0;
  h2_corrupt  = 0;
  h2_fail_at  = 2;

  flag = FIRKodeEvolve(R.mem, SUN_RCONST(0.1), R.y, &tret, FIRK_NORMAL);
  snprintf(what, sizeof(what),
           "H2 (%s): flag %d, %ld failed DQ calls, %ld calls with yn changed "
           "in another component",
           band ? "band" : "dense", flag, h2_failed, h2_corrupt);
  fails += check(flag == FIRK_SUCCESS && h2_failed == 1 && h2_corrupt == 0,
                 what, 3);

  h2_mem   = NULL;
  h2_ytemp = NULL;
  N_VDestroy(h2_snap);
  h2_snap = NULL;
  teardown(&R);
  return fails;
}

/* H4: the stacked SPFGMR of a run (internal) */
static SUNLinearSolverContent_SPFGMR stacked_content(void* mem)
{
  SUNLinearSolver S = ((FIRKLsMem)((FIRKodeMem)mem)->lmem)->SPFGMR;
  return (S == NULL) ? NULL : (SUNLinearSolverContent_SPFGMR)S->content;
}

static int check_stage_setters(SUNContext sunctx)
{
  int fails = 0;
  sunrealtype tret;
  SUNLinearSolverContent_SPFGMR c;
  SUNLinearSolver S0;
  char what[200];
  struct run R;

  if (setup(sunctx, &R, 3, SUN_RCONST(1.0e-6))) { return 1; }
  fails += check(FIRKodeSetStageSolverGSType(R.mem, SUN_CLASSICAL_GS) ==
                   FIRK_ILL_INPUT,
                 "H4: SetStageSolverGSType is an EXPERT option", 3);
  fails += check(FIRKodeSetMode(R.mem, FIRK_EXPERT) == FIRK_SUCCESS,
                 "H4: SetMode(EXPERT) succeeds", 3);
  fails += check(FIRKodeEvolve(R.mem, SUN_RCONST(0.01), R.y, &tret,
                               FIRK_NORMAL) == FIRK_SUCCESS,
                 "H4: the first Evolve succeeds", 3);

  /* the number of restarts edits the live solver */
  S0 = ((FIRKLsMem)((FIRKodeMem)R.mem)->lmem)->SPFGMR;
  fails += check(FIRKodeSetStageSolverMaxRestarts(R.mem, 3) == FIRKLS_SUCCESS,
                 "H4: SetStageSolverMaxRestarts succeeds", 3);
  c = stacked_content(R.mem);
  fails += check(S0 != NULL &&
                   ((FIRKLsMem)((FIRKodeMem)R.mem)->lmem)->SPFGMR == S0 &&
                   c != NULL && c->max_restarts == 3,
                 "H4: SetStageSolverMaxRestarts edits the live solver", 3);

  /* a new dimension and Gram-Schmidt type between Evolve calls */
  fails += check(FIRKodeSetStageSolverMaxl(R.mem, 2) == FIRKLS_SUCCESS &&
                   FIRKodeSetStageSolverGSType(R.mem, SUN_CLASSICAL_GS) ==
                     FIRKLS_SUCCESS,
                 "H4: the setters succeed", 3);
  fails += check(((FIRKLsMem)((FIRKodeMem)R.mem)->lmem)->SPFGMR == S0 &&
                   stacked_content(R.mem) != NULL,
                 "H4: the stacked solver is not freed by a setter", 3);
  fails += check(FIRKodeSetStageSolverGSType(R.mem, 99) == FIRKLS_ILL_INPUT,
                 "H4: an invalid Gram-Schmidt type is rejected", 3);
  fails += check(FIRKodeEvolve(R.mem, SUN_RCONST(0.05), R.y, &tret,
                               FIRK_NORMAL) == FIRK_SUCCESS,
                 "H4: Evolve after the setters succeeds", 3);
  c = stacked_content(R.mem);
  snprintf(what, sizeof(what),
           "H4: the new solver has maxl %d, gstype %d, restarts %d",
           c ? c->maxl : -1, c ? c->gstype : -1, c ? c->max_restarts : -1);
  fails += check(c != NULL && c->maxl == 2 && c->gstype == SUN_CLASSICAL_GS &&
                   c->max_restarts == 3,
                 what, 3);

  /* after ReInit, BEGINNER mode restores the defaults: dimension
     max(5, 3s), modified Gram-Schmidt */
  fails += check(FIRKodeReInit(R.mem, tret, R.y) == FIRK_SUCCESS &&
                   FIRKodeSetMode(R.mem, FIRK_BEGINNER) == FIRK_SUCCESS,
                 "H4: ReInit and SetMode(BEGINNER) succeed", 3);
  fails += check(FIRKodeEvolve(R.mem, SUN_RCONST(0.1), R.y, &tret,
                               FIRK_NORMAL) == FIRK_SUCCESS,
                 "H4: Evolve after SetMode succeeds", 3);
  c = stacked_content(R.mem);
  snprintf(what, sizeof(what),
           "H4: the default solver has maxl %d, gstype %d again",
           c ? c->maxl : -1, c ? c->gstype : -1);
  fails += check(c != NULL && c->maxl == 9 && c->gstype == SUN_MODIFIED_GS,
                 what, 3);
  teardown(&R);
  return fails;
}

/* the largest error from exp(lambda t) */
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

static int check_input_validation(SUNContext sunctx)
{
  int fails     = 0;
  long int nje0 = 0, nje = 0, nfeDQ0 = 0, nfeDQ = 0;
  sunrealtype tret;
  SUNMatrix J = NULL;
  FIRKLsMem ls;
  char what[200];
  struct run R;

  R.y = N_VNew_Serial(NEQ, sunctx);
  N_VConst(ONE, R.y);
  R.A   = SUNDenseMatrix(NEQ, NEQ, sunctx);
  R.LS  = SUNLinSol_Dense(R.y, R.A, sunctx);
  R.mem = FIRKodeCreate(sunctx);
  if (FIRKodeInit(R.mem, f, ZERO, R.y) ||
      FIRKodeSStolerances(R.mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-6)) ||
      FIRKodeSetLinearSolver(R.mem, R.LS, R.A) || FIRKodeSetJacFn(R.mem, jac))
  {
    return 1;
  }
  ls = (FIRKLsMem)((FIRKodeMem)R.mem)->lmem;

  /* a dense solver takes no preconditioner, and the rejection changes
     nothing */
  fails +=
    check(FIRKodeSetPreconditioner(R.mem, psetup, psolve) == FIRKLS_ILL_INPUT,
          "SetPreconditioner rejects a dense solver", 3);
  fails += check(ls->pset == NULL && ls->psolve == NULL,
                 "a rejected SetPreconditioner changes nothing", 3);

  fails += check(FIRKodeEvolve(R.mem, SUN_RCONST(0.01), R.y, &tret,
                               FIRK_NORMAL) == FIRK_SUCCESS,
                 "linsys: the first Evolve succeeds", 3);

  /* a user linsys between Evolve calls: no saved Jacobian, so the stacked
     operator uses difference quotients */
  report_jcur = SUNTRUE;
  fails += check(FIRKodeSetLinSysFn(R.mem, linsys) == FIRKLS_SUCCESS,
                 "SetLinSysFn(linsys) succeeds", 3);
  FIRKodeGetJac(R.mem, &J);
  fails += check(J == NULL, "a user linsys drops the saved Jacobian", 3);
  FIRKodeGetNumLinRhsEvals(R.mem, &nfeDQ0);
  fails += check(FIRKodeEvolve(R.mem, SUN_RCONST(0.05), R.y, &tret,
                               FIRK_NORMAL) == FIRK_SUCCESS,
                 "linsys: Evolve with the user linsys succeeds", 3);
  FIRKodeGetNumLinRhsEvals(R.mem, &nfeDQ);
  snprintf(what, sizeof(what),
           "linsys: difference-quotient products (%ld), error %g",
           nfeDQ - nfeDQ0, (double)exact_error(R.y, tret));
  fails += check(nfeDQ > nfeDQ0 && exact_error(R.y, tret) <= SUN_RCONST(1.0e-4),
                 what, 3);

  /* the default linsys back: a saved Jacobian again, evaluated before its
     first use */
  fails += check(FIRKodeSetLinSysFn(R.mem, NULL) == FIRKLS_SUCCESS,
                 "SetLinSysFn(NULL) succeeds", 3);
  FIRKodeGetJac(R.mem, &J);
  fails += check(J != NULL, "the default linsys saves a Jacobian again", 3);
  FIRKodeGetNumJacEvals(R.mem, &nje0);
  fails += check(FIRKodeEvolve(R.mem, SUN_RCONST(0.1), R.y, &tret,
                               FIRK_NORMAL) == FIRK_SUCCESS,
                 "linsys: Evolve with the default linsys succeeds", 3);
  FIRKodeGetNumJacEvals(R.mem, &nje);
  snprintf(what, sizeof(what),
           "linsys: the Jacobian evaluated again (%ld), error %g", nje - nje0,
           (double)exact_error(R.y, tret));
  fails += check(nje > nje0 && exact_error(R.y, tret) <= SUN_RCONST(1.0e-4),
                 what, 3);
  teardown(&R);
  return fails;
}

/* one SUNLinearSolver object for the system and the mass matrix */
static int check_shared_solver(SUNContext sunctx)
{
  int fails = 0;
  SUNMatrix M;
  SUNLinearSolver MLS;
  FIRKodeMem mem;
  struct run R;

  R.y = N_VNew_Serial(NEQ, sunctx);
  N_VConst(ONE, R.y);
  R.A   = SUNDenseMatrix(NEQ, NEQ, sunctx);
  R.LS  = SUNLinSol_Dense(R.y, R.A, sunctx);
  M     = SUNDenseMatrix(NEQ, NEQ, sunctx);
  MLS   = SUNLinSol_Dense(R.y, M, sunctx);
  R.mem = FIRKodeCreate(sunctx);
  if (FIRKodeInit(R.mem, f, ZERO, R.y) ||
      FIRKodeSStolerances(R.mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-6)) ||
      FIRKodeSetLinearSolver(R.mem, R.LS, R.A))
  {
    return 1;
  }
  mem = (FIRKodeMem)R.mem;

  /* the system solver offered as the mass solver */
  fails += check(FIRKodeSetMassLinearSolver(R.mem, R.LS, M, SUNFALSE) ==
                   FIRKLS_ILL_INPUT,
                 "the system solver is rejected as the mass solver", 3);
  fails += check(mem->mass_mem == NULL && !mem->mass_set,
                 "a rejected mass solver attaches nothing", 3);

  /* a distinct mass solver, then offered as the system solver */
  fails +=
    check(FIRKodeSetMassLinearSolver(R.mem, MLS, M, SUNFALSE) == FIRKLS_SUCCESS,
          "a distinct mass solver is accepted", 3);
  fails += check(FIRKodeSetMassEpsLin(R.mem, -ONE) == FIRKLS_ILL_INPUT &&
                   FIRKodeSetMassEpsLin(R.mem, SUN_RCONST(0.1)) == FIRKLS_SUCCESS &&
                   ((FIRKLsMassMem)mem->mass_mem)->eplifac == SUN_RCONST(0.1) &&
                   FIRKodeSetMassEpsLin(R.mem, ZERO) == FIRKLS_SUCCESS &&
                   ((FIRKLsMassMem)mem->mass_mem)->eplifac == FIRKLS_EPLIN,
                 "SetMassEpsLin rejects a negative factor, 0 restores the "
                 "default",
                 3);
  fails += check(FIRKodeSetLinearSolver(R.mem, MLS, R.A) == FIRKLS_ILL_INPUT,
                 "the mass solver is rejected as the system solver", 3);
  fails += check(mem->lmem != NULL && ((FIRKLsMem)mem->lmem)->LS == R.LS,
                 "a rejected system solver leaves the attached one", 3);

  teardown(&R);
  SUNLinSolFree(MLS);
  SUNMatDestroy(M);
  return fails;
}

/* M = I, matrix-free */
static int mtimes_identity(N_Vector v, N_Vector Mv, sunrealtype t,
                           void* mtimes_data)
{
  (void)t;
  (void)mtimes_data;
  N_VScale(ONE, v, Mv);
  return 0;
}

/* a run with SPGMR, M = I and a mock PCG mass solver, s = 3 */
static int setup_mass(SUNContext sunctx, struct run* R, SUNLinearSolver* MLS,
                      FirkMockLSContent** mmock)
{
  if (setup(sunctx, R, 3, SUN_RCONST(1.0e-6))) { return 1; }
  *MLS = FirkMockLS_New(SUNLinSol_PCG(R->y, SUN_PREC_NONE, 0, sunctx), sunctx);
  if (*MLS == NULL) { return 1; }
  *mmock = FirkMockLS_Content(*MLS);
  return FIRKodeSetMassLinearSolver(R->mem, *MLS, NULL, SUNFALSE) ||
         FIRKodeSetMassTimes(R->mem, NULL, mtimes_identity, NULL);
}

static int check_mass_flags(SUNContext sunctx)
{
  int fails = 0, flag, c;
  long int ncfn;
  sunrealtype tret;
  N_Vector b;
  FIRKodeMem fm;
  SUNLinearSolver MLS;
  FirkMockLSContent* mmock;
  char what[200];
  struct run R;

  static const struct
  {
    int value;
    int deriv; /* expected return of a derivative solve */
    int conv;  /* and of a conversion solve             */
  } cases[] = {{SUNLS_RES_REDUCED, 0, 1},
               {SUNLS_RECOV_FAILURE, 1, 1},
               {SUNLS_PSET_FAIL_REC, 1, 1},
               {SUN_ERR_EXT_FAIL, -1, -1},
               {SUNLS_PSET_FAIL_UNREC, -1, -1}};

  /* direct mass solves after one step */
  if (setup_mass(sunctx, &R, &MLS, &mmock)) { return 1; }
  fails +=
    check(FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP) == FIRK_SUCCESS,
          "mass: the first step succeeds", 3);
  fm = (FIRKodeMem)R.mem;
  b  = N_VClone(R.y);
  for (c = 0; c < (int)(sizeof(cases) / sizeof(cases[0])); c++)
  {
    int rd, rc;
    mmock->inject_first = mmock->nsolves + 1;
    mmock->inject_last  = mmock->nsolves + 2;
    mmock->inject_value = cases[c].value;
    N_VConst(ONE, b);
    rd = fm->msolve(fm, b, SUN_RCONST(1.0));
    N_VConst(ONE, b);
    rc = firkLsMassSolve(fm, b, SUN_RCONST(1.0), SUNFALSE);
    snprintf(what, sizeof(what),
             "mass flag %d: derivative solve %d (expected %d), conversion %d "
             "(expected %d)",
             cases[c].value, rd, cases[c].deriv, rc, cases[c].conv);
    fails += check(rd == cases[c].deriv && rc == cases[c].conv, what, 3);
  }
  mmock->inject_first = 0;
  snprintf(what, sizeof(what), "H11: %ld of %ld mass solves without a zero guess",
           mmock->nsolves_no_zero, mmock->nsolves);
  fails += check(mmock->nsolves > 0 && mmock->nsolves_no_zero == 0, what, 3);
  N_VDestroy(b);
  teardown(&R);
  SUNLinSolFree(MLS);

  /* the first mass solve of a run is a derivative of y in the initial step:
     a reduced residual is accepted */
  if (setup_mass(sunctx, &R, &MLS, &mmock)) { return fails + 1; }
  mmock->inject_first = 1;
  mmock->inject_last  = 1;
  mmock->inject_value = SUNLS_RES_REDUCED;
  flag                = FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP);
  ncfn                = -1;
  FIRKodeGetNumStepSolveFails(R.mem, &ncfn);
  snprintf(what, sizeof(what),
           "mass: RES_REDUCED accepted in the initial step: flag %d, ncfn %ld",
           flag, ncfn);
  fails += check(flag == FIRK_SUCCESS && ncfn == 0, what, 3);
  teardown(&R);
  SUNLinSolFree(MLS);

  /* with a given initial step the first mass solve converts the first stacked
     solve: a reduced residual fails it recoverably */
  if (setup_mass(sunctx, &R, &MLS, &mmock) ||
      FIRKodeSetInitStep(R.mem, SUN_RCONST(1.0e-3)))
  {
    return fails + 1;
  }
  mmock->inject_first = 1;
  mmock->inject_last  = 1;
  mmock->inject_value = SUNLS_RES_REDUCED;
  flag                = FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP);
  ncfn                = -1;
  FIRKodeGetNumStepSolveFails(R.mem, &ncfn);
  snprintf(what, sizeof(what),
           "mass: RES_REDUCED fails the conversion recoverably: flag %d, ncfn "
           "%ld",
           flag, ncfn);
  fails += check(flag == FIRK_SUCCESS && ncfn == 1, what, 3);
  teardown(&R);
  SUNLinSolFree(MLS);
  return fails;
}

int main(void)
{
  SUNContext sunctx = NULL;
  int fails         = 0;
  int policy, use_linsys;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  fails += check_counts(sunctx, 1);
  fails += check_counts(sunctx, 3);
  fails += check_jtsetup_matrixbased(sunctx, 1);
  fails += check_jtsetup_matrixbased(sunctx, 3);
  fails += check_res_reduced(sunctx, 1, SUNTRUE);
  fails += check_res_reduced(sunctx, 2, SUNFALSE);
  for (policy = FIRK_POLICY_RADAU5; policy <= FIRK_POLICY_PERIODIC; policy++)
  {
    for (use_linsys = 0; use_linsys <= 1; use_linsys++)
    {
      fails += check_jcur(sunctx, use_linsys, policy, SUNTRUE);
      fails += check_jcur(sunctx, use_linsys, policy, SUNFALSE);
    }
  }
  fails += check_scalesol(sunctx);
  fails += check_flags(sunctx);
  fails += check_stacked_psolve(sunctx);
  fails += check_pretypes(sunctx);
  fails += check_psolve_failures(sunctx);
  fails += check_stacked_operator(sunctx);
  fails += check_matvec_fallback(sunctx);
  fails += check_dq_restore(sunctx, SUNFALSE);
  fails += check_dq_restore(sunctx, SUNTRUE);
  fails += check_stage_setters(sunctx);
  fails += check_input_validation(sunctx);
  fails += check_shared_solver(sunctx);
  fails += check_mass_flags(sunctx);

  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("All tests passed\n");
  return 0;
}
