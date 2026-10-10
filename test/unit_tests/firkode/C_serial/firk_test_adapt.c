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
 * Unit test for the selection of the step size controller (design Table 9.1):
 *
 * - the default controller is of type SUN_ADAPTCONTROLLER_H, owned by
 *   FIRKODE, and has the parameters of ImpGus with k1 = k2 = 1 (RADAU5's
 *   predictive controller, design Table 5.9);
 * - FIRKodeSetAdaptController installs a controller of type
 *   SUN_ADAPTCONTROLLER_H without taking ownership, rejects a controller of
 *   any other type with FIRK_ILL_INPUT (keeping the current one), and
 *   reinstalls the owned default for C = NULL;
 * - the default stores accepted errors with a floor of 1e-2, while the
 *   current error estimate and explicitly selected controllers retain the
 *   actual error; forward/backward steps and controller replacement obey
 *   the same rule;
 * - FIRKodeSetAdaptControllerByName installs, owned by FIRKODE, the controller
 *   of each name of the Soderlind family (its parameters equal those of the
 *   constructor of that name), and rejects ImExGus, unknown names and NULL
 *   with FIRK_ILL_INPUT (keeping the current controller); a user controller it
 *   replaces is not freed;
 * - firkAdapt (design Listing 8.6) gives the closed forms of the next step for
 *   a mock controller proposing h_C = factor h: with the Newton count factor
 *   kappa = 0.9 min(1, 15 / (14 + k_nls)) and h_I = kappa |h| dsm^(-1/4),
 *   min(0.9 h_C, h_I) (BEGINNER, accepted), h_I without calling the controller
 *   (BEGINNER, rejected), kappa h_C (EXPERT, accepted), min(kappa h_C, h_I)
 *   (EXPERT, rejected); the caps etamax and etamin; the dead band only with a
 *   Newton rate of at most theta_jac, then with keep_setup; hmax last;
 * - H13: dsm = 0 reaches the controller as 1e-10 and gives a finite step;
 *   on y' = 3 t^2, whose solution t^3 the method and its estimate reproduce,
 *   the step after the first grows by kappa (1e-10)^(-1/4), the floor's
 *   I-controller proposal, below the first-step cap 1e4, and the later steps
 *   by the growth cap 8 (not in single precision, where the roundoff of the
 *   estimate is above the floor);
 * - H36: during an integration the controller is asked for a step at most once
 *   per attempt and always before its history update: in BEGINNER mode once per
 *   accepted step, in EXPERT mode also once per failed error test; on the
 *   van der Pol oscillator (mu = 10), whose error varies, the predictive
 *   proposal 0.9 h_C of the default controller is more than 1% below the
 *   I-controller's kappa |h| dsm^(-1/4) on some accepted steps of BEGINNER
 *   mode, so that it sets their next step, and the I-controller itself, in
 *   its place, never is;
 * - R18: at t0 = 0.01 / u the solution sin(1e3 (t - t0)), which no step of at
 *   least RADAU5's minimum 10 u |tn| = 0.1 resolves: no attempt is shorter than
 *   that minimum and the run ends with FIRK_ERR_FAILURE before maxnef error
 *   test failures.
 *
 * Pass criterion: every return value, message, installed controller,
 * ownership flag, step factor and call count matches.
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "firkode/firkode.h"
#include "firkode/firkode_impl.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "sunadaptcontroller/sunadaptcontroller_soderlind.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ONE SUN_RCONST(1.0)
#include "sundials/sundials_context.h"

static char last_msg[256];

static void capture(int line, const char* func, const char* file, const char* msg,
                    SUNErrCode err_code, void* err_user_data, SUNContext sunctx)
{
  (void)line;
  (void)func;
  (void)file;
  (void)err_code;
  (void)err_user_data;
  (void)sunctx;
  snprintf(last_msg, sizeof(last_msg), "%s", msg);
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

/* the installed controller and its ownership flag */
static SUNAdaptController controller(void* mem)
{
  return ((FIRKodeMem)mem)->hadapt_mem->hcontroller;
}

static sunbooleantype owned(void* mem)
{
  return ((FIRKodeMem)mem)->hadapt_mem->owncontroller;
}

static int same_parameters(SUNAdaptController A, SUNAdaptController B);

static int check_types(SUNContext sunctx, void* mem)
{
  int fails = 0;
  SUNAdaptController C0, C_none, C_pi, C_ref;

  /* the default */
  C0 = controller(mem);
  fails += check(C0 != NULL &&
                   SUNAdaptController_GetType(C0) == SUN_ADAPTCONTROLLER_H &&
                   owned(mem),
                 "the default controller is of type H and owned");
  fails += check(((FIRKodeMem)mem)->hadapt_mem->floor_history,
                 "the default controller floors stored errors");
  C_ref = SUNAdaptController_ImpGus(sunctx);
  if (C_ref == NULL ||
      SUNAdaptController_SetParams_ImpGus(C_ref, SUN_RCONST(1.0),
                                          SUN_RCONST(1.0)) != SUN_SUCCESS)
  {
    return fails + 1;
  }
  fails += check(same_parameters(C0, C_ref),
                 "the default controller is ImpGus with k1 = k2 = 1");
  SUNAdaptController_Destroy(C_ref);

  /* a controller without a type is rejected and the default kept */
  C_none = SUNAdaptController_NewEmpty(sunctx);
  if (C_none == NULL) { return fails + 1; }
  last_msg[0] = '\0';
  fails += check(FIRKodeSetAdaptController(mem, C_none) == FIRK_ILL_INPUT,
                 "a controller of type NONE is rejected");
  fails += check(strcmp(last_msg, "The controller must be of type "
                                  "SUN_ADAPTCONTROLLER_H.") == 0,
                 "the type message");
  fails += check(controller(mem) == C0 && owned(mem),
                 "a rejected controller keeps the current one");
  SUNAdaptController_Destroy(C_none);

  /* a user controller of type H is installed, not owned */
  C_pi = SUNAdaptController_PI(sunctx);
  if (C_pi == NULL) { return fails + 1; }
  fails += check(FIRKodeSetAdaptController(mem, C_pi) == FIRK_SUCCESS &&
                   controller(mem) == C_pi && !owned(mem) &&
                   !((FIRKodeMem)mem)->hadapt_mem->floor_history,
                 "a user controller is installed and not owned");

  /* NULL reinstalls an owned default; the user's controller stays alive */
  fails += check(FIRKodeSetAdaptController(mem, NULL) == FIRK_SUCCESS &&
                   controller(mem) != NULL && controller(mem) != C_pi &&
                   owned(mem) && ((FIRKodeMem)mem)->hadapt_mem->floor_history,
                 "NULL reinstalls the owned default");
  fails += check(SUNAdaptController_GetType(C_pi) == SUN_ADAPTCONTROLLER_H,
                 "the user controller is not freed");
  SUNAdaptController_Destroy(C_pi);

  return fails;
}

/* the names accepted by FIRKodeSetAdaptControllerByName */
#define NUM_NAMES 10
static const char* names[NUM_NAMES] = {"Soderlind", "PID",    "PI",    "I",
                                       "ExpGus",    "ImpGus", "H0211", "H0321",
                                       "H211",      "H312"};

static SUNAdaptController create(SUNContext sunctx, int i)
{
  switch (i)
  {
  case 0: return SUNAdaptController_Soderlind(sunctx);
  case 1: return SUNAdaptController_PID(sunctx);
  case 2: return SUNAdaptController_PI(sunctx);
  case 3: return SUNAdaptController_I(sunctx);
  case 4: return SUNAdaptController_ExpGus(sunctx);
  case 5: return SUNAdaptController_ImpGus(sunctx);
  case 6: return SUNAdaptController_H0211(sunctx);
  case 7: return SUNAdaptController_H0321(sunctx);
  case 8: return SUNAdaptController_H211(sunctx);
  case 9: return SUNAdaptController_H312(sunctx);
  default: return NULL;
  }
}

/* the parameters of two controllers of the Soderlind family are equal */
static int same_parameters(SUNAdaptController A, SUNAdaptController B)
{
  SUNAdaptControllerContent_Soderlind a, b;

  if (A == NULL || B == NULL) { return 0; }
  a = (SUNAdaptControllerContent_Soderlind)A->content;
  b = (SUNAdaptControllerContent_Soderlind)B->content;
  return a->k1 == b->k1 && a->k2 == b->k2 && a->k3 == b->k3 && a->k4 == b->k4 &&
         a->k5 == b->k5 && a->historysize == b->historysize;
}

static int check_names(SUNContext sunctx, void* mem)
{
  int i, fails = 0;
  char what[128];
  SUNAdaptController C, C_ref, C_user;
  const char* bad_names[] = {"ImExGus", "pid", "", "Unknown"};
  const char* expected =
    "Unknown controller name (must be Soderlind, PID, PI, I, ExpGus, ImpGus, "
    "H0211, H0321, H211 or H312).";

  for (i = 0; i < NUM_NAMES; i++)
  {
    C_ref = create(sunctx, i);
    snprintf(what, sizeof(what), "the controller %s is installed and owned",
             names[i]);
    fails +=
      check(FIRKodeSetAdaptControllerByName(mem, names[i]) == FIRK_SUCCESS &&
              owned(mem) && !((FIRKodeMem)mem)->hadapt_mem->floor_history &&
              SUNAdaptController_GetType(controller(mem)) ==
                SUN_ADAPTCONTROLLER_H &&
              same_parameters(controller(mem), C_ref),
            what);
    SUNAdaptController_Destroy(C_ref);
  }

  /* rejected names keep the current controller */
  C = controller(mem);
  for (i = 0; i < (int)(sizeof(bad_names) / sizeof(bad_names[0])); i++)
  {
    last_msg[0] = '\0';
    snprintf(what, sizeof(what), "the name '%s' is rejected", bad_names[i]);
    fails += check(FIRKodeSetAdaptControllerByName(mem, bad_names[i]) ==
                       FIRK_ILL_INPUT &&
                     strcmp(last_msg, expected) == 0 && controller(mem) == C &&
                     owned(mem),
                   what);
  }
  last_msg[0] = '\0';
  fails += check(FIRKodeSetAdaptControllerByName(mem, NULL) == FIRK_ILL_INPUT &&
                   strcmp(last_msg, expected) == 0 && controller(mem) == C,
                 "a NULL name is rejected");

  /* a user controller replaced by name is not freed */
  C_user = SUNAdaptController_PI(sunctx);
  if (C_user == NULL) { return fails + 1; }
  fails += check(FIRKodeSetAdaptController(mem, C_user) == FIRK_SUCCESS &&
                   FIRKodeSetAdaptControllerByName(mem, "I") == FIRK_SUCCESS &&
                   controller(mem) != C_user && owned(mem) &&
                   SUNAdaptController_GetType(C_user) == SUN_ADAPTCONTROLLER_H,
                 "a user controller replaced by name is not freed");
  SUNAdaptController_Destroy(C_user);

  return fails;
}

/* ---------------------------------------------------------------------------
 * A recording mock controller: proposes factor * h, or wraps an inner
 * controller, and records its calls as 'E' (EstimateStep) and 'U' (UpdateH)
 * ---------------------------------------------------------------------------*/

static sunrealtype mock_factor = SUN_RCONST(2.0);
static sunrealtype mock_dsm    = SUN_RCONST(-1.0);
static char mock_calls[8192];
static int mock_ncalls = 0;
static void* mock_mem  = NULL; /* set to compare with the I-controller */
static int mock_nbind  = 0;    /* proposals with 0.9 h_C < 0.99 h_I    */

static void mock_record(char c)
{
  if (mock_ncalls < (int)sizeof(mock_calls) - 1)
  {
    mock_calls[mock_ncalls++] = c;
    mock_calls[mock_ncalls]   = '\0';
  }
}

static SUNAdaptController_Type mock_gettype(SUNAdaptController C)
{
  (void)C;
  return SUN_ADAPTCONTROLLER_H;
}

static SUNErrCode mock_estimate(SUNAdaptController C, sunrealtype h, int p,
                                sunrealtype dsm, sunrealtype* hnew)
{
  mock_record('E');
  mock_dsm = dsm;
  if (C->content != NULL)
  {
    SUNErrCode retval =
      SUNAdaptController_EstimateStep((SUNAdaptController)C->content, h, p, dsm,
                                      hnew);
    if (retval == SUN_SUCCESS && mock_mem != NULL)
    {
      /* the I-controller proposal of design Listing 8.6 */
      FIRKodeMem fm     = (FIRKodeMem)mock_mem;
      sunrealtype kappa = fm->hadapt_mem->safety;
      sunrealtype h_I;
      if (fm->hadapt_mem->newt_safety)
      {
        kappa *=
          SUNMIN(SUN_RCONST(1.0),
                 (SUN_RCONST(2.0) * (sunrealtype)fm->maxcor + SUN_RCONST(1.0)) /
                   (SUN_RCONST(2.0) * (sunrealtype)fm->maxcor +
                    (sunrealtype)fm->nni_last));
      }
      h_I = kappa * SUNRabs(h) *
            SUNRpowerR(dsm,
                       -SUN_RCONST(1.0) / ((sunrealtype)p + SUN_RCONST(1.0)));
      if (fm->hadapt_mem->safety * SUNRabs(*hnew) < SUN_RCONST(0.99) * h_I)
      {
        mock_nbind++;
      }
    }
    return retval;
  }
  *hnew = mock_factor * h;
  return SUN_SUCCESS;
}

static SUNErrCode mock_updateh(SUNAdaptController C, sunrealtype h,
                               sunrealtype dsm)
{
  mock_record('U');
  if (C->content != NULL)
  {
    return SUNAdaptController_UpdateH((SUNAdaptController)C->content, h, dsm);
  }
  return SUN_SUCCESS;
}

static SUNErrCode mock_reset(SUNAdaptController C)
{
  if (C->content != NULL)
  {
    return SUNAdaptController_Reset((SUNAdaptController)C->content);
  }
  return SUN_SUCCESS;
}

static SUNAdaptController mock_new(SUNContext sunctx, SUNAdaptController inner)
{
  SUNAdaptController C = SUNAdaptController_NewEmpty(sunctx);
  if (C == NULL) { return NULL; }
  C->ops->gettype      = mock_gettype;
  C->ops->estimatestep = mock_estimate;
  C->ops->updateh      = mock_updateh;
  C->ops->reset        = mock_reset;
  C->content           = inner;
  return C;
}

static void mock_free(SUNAdaptController C)
{
  if (C->content != NULL)
  {
    SUNAdaptController_Destroy((SUNAdaptController)C->content);
  }
  C->content = NULL;
  SUNAdaptController_DestroyEmpty(C);
}

static int close_to(sunrealtype a, sunrealtype b)
{
  return SUNRabs(a - b) <= SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF * SUNRabs(b);
}

/* one call of firkAdapt with h = 0.1, p = 3, k_nls = 2 (kappa = 0.84375) */
static sunrealtype adapt(void* mem, int mode, sunrealtype dsm,
                         sunbooleantype accepted, sunrealtype factor,
                         sunrealtype theta)
{
  FIRKodeMem fm = (FIRKodeMem)mem;

  fm->mode               = mode;
  fm->h                  = SUN_RCONST(0.1);
  fm->hadapt_mem->p      = 3;
  fm->hadapt_mem->etamax = SUN_RCONST(8.0);
  fm->nni_last           = 2;
  fm->theta              = theta;
  mock_factor            = factor;
  mock_ncalls            = 0;
  if (firkAdapt(fm, dsm, accepted) != FIRK_SUCCESS) { return -SUN_RCONST(1.0); }
  return fm->eta;
}

static int check_closed_forms(SUNContext sunctx)
{
  int fails = 0;
  void* mem = FIRKodeCreate(sunctx);
  FIRKodeMem fm;
  SUNAdaptController C   = mock_new(sunctx, NULL);
  const sunrealtype k    = SUN_RCONST(0.84375); /* 0.9 * 15 / 16 */
  const sunrealtype big  = SUN_RCONST(1.0e-1);  /* a rate above theta_jac */
  const sunrealtype fac2 = SUN_RCONST(2.0);
  const sunrealtype dsm1 = SUN_RCONST(0.01); /* dsm^(-1/4) = sqrt(10) */
  const sunrealtype dsm2 = SUN_RCONST(16.0); /* dsm^(-1/4) = 1/2 */
  sunrealtype eta, hI1;

  if (mem == NULL || C == NULL ||
      FIRKodeSetMode(mem, FIRK_EXPERT) != FIRK_SUCCESS ||
      FIRKodeSetAdaptController(mem, C) != FIRK_SUCCESS)
  {
    return 1;
  }
  fm  = (FIRKodeMem)mem;
  hI1 = k * SUNRsqrt(SUN_RCONST(10.0)); /* h_I / h for dsm1 */

  /* BEGINNER, accepted: min(0.9 h_C, h_I) */
  eta = adapt(mem, FIRK_BEGINNER, dsm1, SUNTRUE, fac2, big);
  fails += check(close_to(eta, SUN_RCONST(1.8)) && mock_ncalls == 1,
                 "BEGINNER accepted: 0.9 h_C below h_I");
  eta = adapt(mem, FIRK_BEGINNER, dsm1, SUNTRUE, SUN_RCONST(10.0), big);
  fails += check(close_to(eta, hI1), "BEGINNER accepted: the h_I minimum");

  /* BEGINNER, rejected: h_I, without the controller */
  eta = adapt(mem, FIRK_BEGINNER, dsm2, SUNFALSE, fac2, big);
  fails += check(close_to(eta, SUN_RCONST(0.5) * k) && mock_ncalls == 0,
                 "BEGINNER rejected: h_I, the controller not called");

  /* EXPERT, accepted: kappa h_C */
  eta = adapt(mem, FIRK_EXPERT, dsm1, SUNTRUE, fac2, big);
  fails += check(close_to(eta, fac2 * k) && mock_ncalls == 1,
                 "EXPERT accepted: kappa h_C");

  /* EXPERT, rejected: min(kappa h_C, h_I), then etamin */
  eta = adapt(mem, FIRK_EXPERT, dsm2, SUNFALSE, fac2, big);
  fails += check(close_to(eta, SUN_RCONST(0.5) * k),
                 "EXPERT rejected: h_I below kappa h_C");
  eta = adapt(mem, FIRK_EXPERT, dsm2, SUNFALSE, SUN_RCONST(0.3), big);
  fails += check(close_to(eta, SUN_RCONST(0.3) * k),
                 "EXPERT rejected: kappa h_C below h_I");
  eta = adapt(mem, FIRK_EXPERT, dsm2, SUNFALSE, SUN_RCONST(0.1), big);
  fails += check(close_to(eta, fm->hadapt_mem->etamin),
                 "EXPERT rejected: the reduction bound etamin");

  /* the growth cap etamax */
  eta = adapt(mem, FIRK_EXPERT, dsm1, SUNTRUE, SUN_RCONST(100.0), big);
  fails += check(close_to(eta, SUN_RCONST(8.0)), "the growth cap etamax");

  /* the dead band [1, 1.2], only with a rate of at most theta_jac */
  eta = adapt(mem, FIRK_EXPERT, dsm1, SUNTRUE, SUN_RCONST(1.1) / k,
              SUN_RCONST(1.0e-4));
  fails += check(eta == SUN_RCONST(1.0) && fm->keep_setup,
                 "dead band: h and the setup kept");
  eta = adapt(mem, FIRK_EXPERT, dsm1, SUNTRUE, SUN_RCONST(1.1) / k,
              SUN_RCONST(1.0e-2));
  fails += check(close_to(eta, SUN_RCONST(1.1)) && !fm->keep_setup,
                 "no dead band with a rate above theta_jac");

  /* hmax applies last */
  fm->hmax_inv = ONE / SUN_RCONST(0.15);
  eta          = adapt(mem, FIRK_EXPERT, dsm1, SUNTRUE, fac2, big);
  fails += check(close_to(eta, SUN_RCONST(1.5)), "hmax applies last");
  fm->hmax_inv = SUN_RCONST(0.0);

  /* H13: dsm = 0 reaches the controller as the floor 1e-10 */
  eta = adapt(mem, FIRK_BEGINNER, SUN_RCONST(0.0), SUNTRUE, fac2, big);
  fails += check(mock_dsm == SUN_RCONST(1.0e-10) && close_to(eta, SUN_RCONST(1.8)),
                 "H13: dsm = 0 is floored, the step is finite");

  FIRKodeFree(&mem);
  mock_free(C);
  return fails;
}

/* H36 on y' = -y from y(0) = 1 to t = 10 with a first step of 1, which fails
   the error test at tolerances 1e-5 (usable in every precision); the mock
   wraps the default controller */
static int f_decay(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)t;
  (void)user_data;
  NV_Ith_S(ydot, 0) = -NV_Ith_S(y, 0);
  return 0;
}

static int check_call_order(SUNContext sunctx, int mode)
{
  int i, fails = 0, nE = 0, nU = 0, order_ok = 1;
  long int nst = 0, netf = 0;
  sunrealtype tret;
  void* mem                = FIRKodeCreate(sunctx);
  N_Vector y               = N_VNew_Serial(1, sunctx);
  SUNMatrix A              = SUNDenseMatrix(1, 1, sunctx);
  SUNLinearSolver LS       = SUNLinSol_Dense(y, A, sunctx);
  SUNAdaptController inner = SUNAdaptController_ImpGus(sunctx);
  SUNAdaptController C;
  char what[128];

  /* the default controller, ImpGus with k1 = k2 = 1 */
  if (inner == NULL ||
      SUNAdaptController_SetParams_ImpGus(inner, ONE, ONE) != SUN_SUCCESS)
  {
    return 1;
  }
  C = mock_new(sunctx, inner);

  NV_Ith_S(y, 0) = ONE;
  if (FIRKodeInit(mem, f_decay, SUN_RCONST(0.0), y) ||
      FIRKodeSStolerances(mem, SUN_RCONST(1.0e-5), SUN_RCONST(1.0e-5)) ||
      FIRKodeSetLinearSolver(mem, LS, A) || FIRKodeSetMode(mem, FIRK_EXPERT) ||
      FIRKodeSetAdaptController(mem, C) || FIRKodeSetInitStep(mem, ONE))
  {
    return 1;
  }
  ((FIRKodeMem)mem)->mode = mode; /* keeps the mock, unlike SetMode */
  mock_ncalls             = 0;

  fails += check(FIRKodeEvolve(mem, SUN_RCONST(10.0), y, &tret, FIRK_NORMAL) ==
                   FIRK_SUCCESS,
                 "the integration succeeds");
  FIRKodeGetNumSteps(mem, &nst);
  FIRKodeGetNumErrTestFails(mem, &netf);

  /* every U follows an E of its attempt */
  for (i = 0; i < mock_ncalls; i++)
  {
    if (mock_calls[i] == 'E') { nE++; }
    else
    {
      nU++;
      if (i == 0 || mock_calls[i - 1] != 'E') { order_ok = 0; }
    }
  }
  snprintf(what, sizeof(what),
           "%s: E calls %d, U calls %d, steps %ld, error test fails %ld",
           (mode == FIRK_BEGINNER) ? "BEGINNER" : "EXPERT", nE, nU, nst, netf);
  fails += check(netf > 0, "the first step fails the error test");
  fails += check(order_ok && nU == nst &&
                   nE == ((mode == FIRK_BEGINNER) ? nst : nst + netf),
                 what);

  FIRKodeFree(&mem);
  mock_free(C);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  return fails;
}

#if !defined(SUNDIALS_SINGLE_PRECISION)
/* H13: a polynomial exact solution, so that dsm is at roundoff; the estimate
   (as RADAU5's) and firkAdapt both floor it at 1e-10 */
static int f_cubic(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)y;
  (void)user_data;
  NV_Ith_S(ydot, 0) = SUN_RCONST(3.0) * t * t;
  return 0;
}

static int check_polynomial(SUNContext sunctx)
{
  int k, fails = 0;
  long int nni1 = 0;
  sunrealtype tret, h[6], expected, kappa;
  const sunrealtype tol = SUN_RCONST(1.0e-6);
  void* mem             = FIRKodeCreate(sunctx);
  N_Vector y            = N_VNew_Serial(1, sunctx);
  SUNMatrix A           = SUNDenseMatrix(1, 1, sunctx);
  SUNLinearSolver LS    = SUNLinSol_Dense(y, A, sunctx);
  char what[160];

  NV_Ith_S(y, 0) = SUN_RCONST(1.0);
  if (FIRKodeInit(mem, f_cubic, SUN_RCONST(1.0), y) ||
      FIRKodeSStolerances(mem, tol, tol) || FIRKodeSetLinearSolver(mem, LS, A) ||
      FIRKodeSetInitStep(mem, SUN_RCONST(1.0e-3)))
  {
    return 1;
  }
  for (k = 0; k < 6; k++)
  {
    if (FIRKodeEvolve(mem, SUN_RCONST(1.0e30), y, &tret, FIRK_ONE_STEP) ||
        FIRKodeGetLastStep(mem, &h[k]))
    {
      return fails + 1;
    }
    if (k == 0) { FIRKodeGetNumNonlinSolvIters(mem, &nni1); }
  }
  for (k = 1; k < 6; k++)
  {
    /* the Newton count safety factor of the first step (maxcor = 7) */
    kappa = SUN_RCONST(0.9) *
            SUNMIN(SUN_RCONST(1.0),
                   SUN_RCONST(15.0) / (SUN_RCONST(14.0) + (sunrealtype)nni1));
    expected = (k == 1)
                 ? kappa * SUNRpowerR(SUN_RCONST(1.0e-10), -SUN_RCONST(0.25))
                 : SUN_RCONST(8.0);
    snprintf(what, sizeof(what), "H13: step %d grows by %g, expected %g", k + 1,
             (double)(h[k] / h[k - 1]), (double)expected);
    fails += check(isfinite((double)h[k]) && SUNRabs(h[k] / h[k - 1] - expected) <=
                                               SUN_RCONST(1.0e-5) * expected,
                   what);
  }

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  return fails;
}
#endif

/* H36: the predictive proposal of ImpGus with k1 = k2 = 1 (the default) sets
   some accepted steps of BEGINNER mode on the van der Pol oscillator with
   mu = 10, and the I-controller in its place none */
static int f_vdp(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  const sunrealtype y1 = NV_Ith_S(y, 0), y2 = NV_Ith_S(y, 1);
  (void)t;
  (void)user_data;
  NV_Ith_S(ydot, 0) = y2;
  NV_Ith_S(ydot, 1) = SUN_RCONST(10.0) * (ONE - y1 * y1) * y2 - y1;
  return 0;
}

static int count_binding(SUNContext sunctx, sunbooleantype predictive,
                         int* nbind, long int* nst)
{
  int flag;
  sunrealtype tret;
  void* mem                = FIRKodeCreate(sunctx);
  N_Vector y               = N_VNew_Serial(2, sunctx);
  SUNMatrix A              = SUNDenseMatrix(2, 2, sunctx);
  SUNLinearSolver LS       = SUNLinSol_Dense(y, A, sunctx);
  SUNAdaptController inner = predictive ? SUNAdaptController_ImpGus(sunctx)
                                        : SUNAdaptController_I(sunctx);
  SUNAdaptController C;

  if (inner == NULL ||
      (predictive &&
       SUNAdaptController_SetParams_ImpGus(inner, ONE, ONE) != SUN_SUCCESS))
  {
    return 1;
  }
  C              = mock_new(sunctx, inner);
  NV_Ith_S(y, 0) = SUN_RCONST(2.0);
  NV_Ith_S(y, 1) = SUN_RCONST(0.0);
  if (FIRKodeInit(mem, f_vdp, SUN_RCONST(0.0), y) ||
      FIRKodeSStolerances(mem, SUN_RCONST(1.0e-5), SUN_RCONST(1.0e-5)) ||
      FIRKodeSetLinearSolver(mem, LS, A) || FIRKodeSetMode(mem, FIRK_EXPERT) ||
      FIRKodeSetAdaptController(mem, C))
  {
    return 1;
  }
  ((FIRKodeMem)mem)->mode = FIRK_BEGINNER; /* keeps the mock */
  mock_mem                = mem;
  mock_nbind              = 0;
  flag = FIRKodeEvolve(mem, SUN_RCONST(20.0), y, &tret, FIRK_NORMAL);
  FIRKodeGetNumSteps(mem, nst);
  *nbind   = mock_nbind;
  mock_mem = NULL;

  FIRKodeFree(&mem);
  mock_free(C);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  return (flag == FIRK_SUCCESS) ? 0 : 1;
}

static int check_predictive(SUNContext sunctx)
{
  int fails = 0, nbind_gus = 0, nbind_i = 0;
  long int nst_gus = 0, nst_i = 0;
  char what[160];

  fails += check(count_binding(sunctx, SUNTRUE, &nbind_gus, &nst_gus) == 0,
                 "van der Pol with ImpGus succeeds");
  fails += check(count_binding(sunctx, SUNFALSE, &nbind_i, &nst_i) == 0,
                 "van der Pol with the I-controller succeeds");
  snprintf(what, sizeof(what),
           "H36: the predictive proposal sets %d of %ld accepted steps (the "
           "I-controller %d of %ld)",
           nbind_gus, nst_gus, nbind_i, nst_i);
  fails += check(nbind_gus > 0 && nbind_i == 0, what);
  return fails;
}

/* The stored error floor belongs to the default controller, not to the
   integrator's error test or to ownership of a controller. A constant RHS
   gives a near-zero error in every precision and in either time direction. */
static int f_constant(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)t;
  (void)y;
  (void)user_data;
  N_VConst(ONE, ydot);
  return 0;
}

static int check_history_floor(SUNContext sunctx, int selection,
                               sunrealtype direction)
{
  int k, fails = 0;
  sunrealtype tret, expected;
  void* mem                          = NULL;
  N_Vector y                         = NULL;
  SUNMatrix A                        = NULL;
  SUNLinearSolver LS                 = NULL;
  SUNAdaptController user_controller = NULL;
  SUNAdaptControllerContent_Soderlind content;
  FIRKodeMem fm;
  char what[160];

  mem = FIRKodeCreate(sunctx);
  y   = N_VNew_Serial(1, sunctx);
  A   = SUNDenseMatrix(1, 1, sunctx);
  if (mem == NULL || y == NULL || A == NULL)
  {
    fails++;
    goto cleanup;
  }
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (LS == NULL)
  {
    fails++;
    goto cleanup;
  }
  fm = (FIRKodeMem)mem;
  N_VConst(ONE, y);
  if (FIRKodeInit(mem, f_constant, SUN_RCONST(0.0), y) ||
      FIRKodeSStolerances(mem, SUN_RCONST(1.0e-3), SUN_RCONST(1.0e-3)) ||
      FIRKodeSetLinearSolver(mem, LS, A) ||
      FIRKodeSetInitStep(mem, direction * SUN_RCONST(1.0e-3)))
  {
    fails++;
    goto cleanup;
  }
  if (selection > 0 && FIRKodeSetMode(mem, FIRK_EXPERT))
  {
    fails++;
    goto cleanup;
  }
  if (selection == 2)
  {
    user_controller = SUNAdaptController_ImpGus(sunctx);
    if (user_controller == NULL ||
        FIRKodeSetAdaptController(mem, user_controller))
    {
      fails++;
      goto cleanup;
    }
  }
  else if (selection == 3 && FIRKodeSetAdaptControllerByName(mem, "ImpGus"))
  {
    fails++;
    goto cleanup;
  }

  for (k = 0; k < 3; k++)
  {
    if (FIRKodeEvolve(mem, direction * SUN_RCONST(10.0), y, &tret, FIRK_ONE_STEP))
    {
      fails++;
      goto cleanup;
    }
    content  = (SUNAdaptControllerContent_Soderlind)controller(mem)->content;
    expected = (selection < 2) ? SUN_RCONST(1.0e-2) : fm->dsm;
    fails += check(fm->dsm < SUN_RCONST(1.0e-2),
                   "the actual accepted error remains below the history floor");
    fails += check(close_to(content->ep, expected),
                   "only the default controller receives floored history");
    fails += check(close_to(content->hp, fm->h),
                   "history retains the signed accepted step");
    snprintf(what, sizeof(what),
             "the accepted solution matches the analytic solution: "
             "y=%g, t=%g, error=%g",
             (double)NV_Ith_S(y, 0), (double)tret,
             (double)(NV_Ith_S(y, 0) - (ONE + tret)));
    /* Inexact stage solves need not reproduce the constant solution to
       roundoff; compare against the requested mixed tolerance. */
    fails += check(SUNRabs(NV_Ith_S(y, 0) - (ONE + tret)) <=
                     SUN_RCONST(1.0e-3) * (ONE + SUNRabs(ONE + tret)),
                   what);
  }

  /* Reset clears history without changing which controller uses the floor. */
  fails += check(FIRKodeReset(mem, tret, y) == FIRK_SUCCESS,
                 "Reset succeeds with controller history");
  if (FIRKodeSetInitStep(mem, direction * SUN_RCONST(1.0e-3)) ||
      FIRKodeEvolve(mem, tret + direction, y, &tret, FIRK_ONE_STEP))
  {
    fails++;
    goto cleanup;
  }
  content  = (SUNAdaptControllerContent_Soderlind)controller(mem)->content;
  expected = (selection < 2) ? SUN_RCONST(1.0e-2) : fm->dsm;
  fails += check(close_to(content->ep, expected),
                 "the history policy survives Reset");

  /* A rejected large step must not become the stored previous error. The
     eventual accepted error above the floor is stored without alteration. */
  N_VConst(ONE, y);
  if (FIRKodeInit(mem, f_decay, SUN_RCONST(0.0), y) ||
      FIRKodeSStolerances(mem, SUN_RCONST(1.0e-5), SUN_RCONST(1.0e-5)) ||
      FIRKodeSetLinearSolver(mem, LS, A) ||
      FIRKodeSetInitStep(mem, direction) ||
      FIRKodeEvolve(mem, direction * SUN_RCONST(10.0), y, &tret, FIRK_ONE_STEP))
  {
    fails++;
    goto cleanup;
  }
  content = (SUNAdaptControllerContent_Soderlind)controller(mem)->content;
  fails +=
    check(fm->netf > 0 && fm->dsm >= SUN_RCONST(1.0e-2) && fm->dsm <= ONE,
          "a large step is rejected before an accepted error above the floor");
  fails += check(close_to(content->ep, fm->dsm),
                 "history stores the accepted error above the floor unchanged");

cleanup:
  FIRKodeFree(&mem);
  if (user_controller != NULL) { SUNAdaptController_Destroy(user_controller); }
  if (LS != NULL) { SUNLinSolFree(LS); }
  if (A != NULL) { SUNMatDestroy(A); }
  if (y != NULL) { N_VDestroy(y); }
  return fails;
}

/* R18: the smallest step is 10 u |tn|. At t0 = 0.01 / u it is 0.1, much
   longer than the period of y = sin(1e3 (t - t0)), in every precision */
static sunrealtype t0_fast = SUN_RCONST(0.0);

static int f_fast(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)y;
  (void)user_data;
  NV_Ith_S(ydot,
           0) = SUN_RCONST(1.0e3) *
                (sunrealtype)cos((double)(SUN_RCONST(1.0e3) * (t - t0_fast)));
  return 0;
}

static int check_min_step(SUNContext sunctx)
{
  int fails     = 0, flag;
  long int netf = 0;
  sunrealtype tret, hmin_rel;
  void* mem          = FIRKodeCreate(sunctx);
  N_Vector y         = N_VNew_Serial(1, sunctx);
  SUNMatrix A        = SUNDenseMatrix(1, 1, sunctx);
  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);

  t0_fast        = SUN_RCONST(0.01) / SUN_UNIT_ROUNDOFF;
  NV_Ith_S(y, 0) = SUN_RCONST(0.0);
  if (FIRKodeInit(mem, f_fast, t0_fast, y) ||
      FIRKodeSStolerances(mem, SUN_RCONST(1.0e-4), SUN_RCONST(1.0e-4)) ||
      FIRKodeSetLinearSolver(mem, LS, A))
  {
    return 1;
  }

  flag = FIRKodeEvolve(mem, t0_fast + ONE, y, &tret, FIRK_NORMAL);
  FIRKodeGetNumErrTestFails(mem, &netf);
  hmin_rel = SUN_RCONST(10.0) * SUN_UNIT_ROUNDOFF * SUNRabs(tret);

  fails += check(flag == FIRK_ERR_FAILURE, "R18: the run fails its error test");
  fails += check(SUNRabs(((FIRKodeMem)mem)->h) >= SUN_RCONST(0.999) * hmin_rel,
                 "R18: no attempt below 10 u |tn|");
  fails += check(netf < MXNEF, "R18: the guard ends the run before maxnef");

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  return fails;
}

int main(void)
{
  SUNContext sunctx = NULL;
  void* mem         = NULL;
  int fails         = 0, selection;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }
  if (SUNContext_PushErrHandler(sunctx, capture, NULL)) { return 1; }

  mem = FIRKodeCreate(sunctx);
  if (mem == NULL) { return 1; }
  if (FIRKodeSetMode(mem, FIRK_EXPERT) != FIRK_SUCCESS) { return 1; }

  fails += check_types(sunctx, mem);
  fails += check_names(sunctx, mem);
  fails += check_closed_forms(sunctx);
  fails += check_call_order(sunctx, FIRK_BEGINNER);
  fails += check_call_order(sunctx, FIRK_EXPERT);
  fails += check_predictive(sunctx);
  for (selection = 0; selection < 4; selection++)
  {
    fails += check_history_floor(sunctx, selection, ONE);
    fails += check_history_floor(sunctx, selection, -ONE);
  }
#if !defined(SUNDIALS_SINGLE_PRECISION)
  fails += check_polynomial(sunctx);
#endif
  fails += check_min_step(sunctx);

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
