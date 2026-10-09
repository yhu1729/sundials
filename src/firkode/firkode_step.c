/* -----------------------------------------------------------------
 * Programmer(s): Yifan Hu @ UMBC
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * This is the implementation file for one FIRKODE step attempt:
 * the stage solve, the handling of its failures, the error estimate
 * and test, and the completion of an accepted step.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <firkode/firkode.h>
#include <sundials/sundials_types.h>

#include "firkode_impl.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/*===============================================================
  Private function prototypes
  ===============================================================*/

static int firkStepAttempt(FIRKodeMem firk_mem, sunbooleantype attempt_rejected,
                           int* nflagPtr, int* ncfPtr, int* nefPtr,
                           sunrealtype* dsmPtr);
static int firkHandleNFlag(FIRKodeMem firk_mem, int nls, int* nflagPtr,
                           int* ncfPtr);
static int firkErrorEstimate(FIRKodeMem firk_mem, sunbooleantype attempt_rejected,
                             sunrealtype* dsmPtr);
static int firkDoErrorTest(FIRKodeMem firk_mem, int* nflagPtr, int* nefPtr,
                           sunrealtype dsm);
static int firkCompleteStep(FIRKodeMem firk_mem, sunrealtype dsm);

/*===============================================================
  Step attempt
  ===============================================================*/

/*---------------------------------------------------------------
  firkStep

  Performs one internal step from tn to tn + h, with retries on
  nonlinear solver or error test failures (design Listing 8.2).
  Every attempt is counted in nst_attempts and logged by one
  begin-step-attempt and one end-step-attempt record.
  ---------------------------------------------------------------*/
int firkStep(FIRKodeMem firk_mem)
{
  int nflag, kflag, ncf, nef;
  sunrealtype dsm;
  sunbooleantype attempt_rejected;

  SUNDIALS_MARK_FUNCTION_BEGIN(FIRK_PROFILER);

  nflag            = FIRST_CALL;
  ncf              = 0;
  nef              = 0;
  attempt_rejected = SUNFALSE;

  for (;;)
  {
    firk_mem->nst_attempts++;

    firk_mem->h     = firk_mem->hprime;
    firk_mem->gamma = firk_mem->h * firk_mem->tbl->gamma0;
    firk_mem->gamrat =
      (firk_mem->gammap != ZERO) ? firk_mem->gamma / firk_mem->gammap : ONE;

    SUNLogInfo(FIRK_LOGGER, "begin-step-attempt",
               "step = %li, tn = " SUN_FORMAT_G ", h = " SUN_FORMAT_G
               ", s = %i",
               firk_mem->nst + 1, firk_mem->tn, firk_mem->h, firk_mem->s);

    kflag = firkStepAttempt(firk_mem, attempt_rejected, &nflag, &ncf, &nef, &dsm);

    /* a failed attempt ends its record here; an accepted one below, once
       completed (as CVODE, which reports a step after completing it) */
    SUNLogInfoIf(kflag == TRY_AGAIN || kflag == FIRK_ERR_FAILURE, FIRK_LOGGER,
                 "end-step-attempt",
                 "status = failed error test, dsm = " SUN_FORMAT_G, dsm);
    SUNLogInfoIf(kflag != FIRK_SUCCESS && kflag != TRY_AGAIN &&
                   kflag != FIRK_ERR_FAILURE,
                 FIRK_LOGGER, "end-step-attempt",
                 "status = failed solve, nflag = %i", nflag);

    if (kflag != PREDICT_AGAIN && kflag != TRY_AGAIN) { break; }

    /* RADAU5's REJECT: this step was rejected, by its Newton iteration or
       by its error test */
    attempt_rejected = SUNTRUE;
  }

  if (kflag == FIRK_SUCCESS)
  {
    kflag = firkCompleteStep(firk_mem, dsm);
    SUNLogInfoIf(kflag == FIRK_SUCCESS, FIRK_LOGGER, "end-step-attempt",
                 "status = success, dsm = " SUN_FORMAT_G, dsm);
    SUNLogInfoIf(kflag != FIRK_SUCCESS, FIRK_LOGGER, "end-step-attempt",
                 "status = failed completion, retval = %i", kflag);
  }

  SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
  return kflag;
}

/*---------------------------------------------------------------
  firkStepAttempt

  One attempt of the step: predictor, stage solve, error estimate
  and error test. Returns FIRK_SUCCESS, PREDICT_AGAIN after a
  recoverable failure of the solve or the estimate, TRY_AGAIN after
  a failed error test, or a negative FIRKODE return code; nflag,
  ncf and nef carry the outcome to the next attempt, dsm is the
  norm of the error estimate (0 in fixed-step mode).
  ---------------------------------------------------------------*/
static int firkStepAttempt(FIRKodeMem firk_mem, sunbooleantype attempt_rejected,
                           int* nflagPtr, int* ncfPtr, int* nefPtr,
                           sunrealtype* dsmPtr)
{
  int kflag, eflag, nls, retval;
  sunbooleantype jbad, need_setup;

  *dsmPtr = ZERO;

  /* the reuse policy decides the setup (design decision D9) */
  firkReusePolicy(firk_mem, *nflagPtr, &jbad, &need_setup);

  /* predict the stage increments */
  kflag = firkPredict(firk_mem);
  if (kflag != FIRK_SUCCESS) { return kflag; }

  /* set up the linear solver, then solve the stage system */
  nls = FIRK_NLS_SUCCESS;
  if (need_setup)
  {
    retval           = firk_mem->lsetup(firk_mem, jbad, &(firk_mem->jcur));
    firk_mem->gamrat = ONE;
    if (retval < 0) { return FIRK_LSETUP_FAIL; }
    if (retval > 0)
    {
      nls                = FIRK_NLS_LSETUP_RECVR;
      firk_mem->nls_hfac = firk_mem->hadapt_mem->etacf;
    }
    else if (firk_mem->jcur) { firk_mem->jac_current = SUNTRUE; }
  }
  if (nls == FIRK_NLS_SUCCESS) { nls = firkNewton(firk_mem); }
  kflag = firkHandleNFlag(firk_mem, nls, nflagPtr, ncfPtr);
  if (kflag != FIRK_SUCCESS) { return kflag; }

  if (firk_mem->fixedstep) { return FIRK_SUCCESS; }

  /* error estimate; a recoverable failure is handled as one of the solve,
     with the step factor of a recoverable solve failure */
  eflag = firkErrorEstimate(firk_mem, attempt_rejected, dsmPtr);
  if (eflag > 0)
  {
    firk_mem->nls_hfac = firk_mem->hadapt_mem->etacf;
    return firkHandleNFlag(firk_mem, FIRK_NLS_LSOLVE_RECVR, nflagPtr, ncfPtr);
  }
  if (eflag < 0) { return eflag; }

  /* error test */
  return firkDoErrorTest(firk_mem, nflagPtr, nefPtr, *dsmPtr);
}

/*---------------------------------------------------------------
  firkHandleNFlag

  Handles the status of the stage solve, or of a recoverable failure
  of the setup or the error estimate (design Listing 8.2). Returns
  FIRK_SUCCESS on success, PREDICT_AGAIN to retry the step, or a
  negative FIRKODE return code. An adaptive step is retried with the
  step reduced by the factor nls_hfac of the failure; in fixed-step
  mode the step is retried once, unreduced, when the Jacobian data
  was not current, so that the next attempt evaluates it (R3).
  ---------------------------------------------------------------*/
static int firkHandleNFlag(FIRKodeMem firk_mem, int nls, int* nflagPtr,
                           int* ncfPtr)
{
  sunbooleantype too_many;

  if (nls == FIRK_NLS_SUCCESS) { return FIRK_SUCCESS; }

  firk_mem->ncfn++;

  /* unrecoverable failures */
  if (nls < 0) { return nls; }

  /* recoverable failure; the policy of the next attempt sees whether J was
     current and how far gamma was from its setup */
  (*ncfPtr)++;
  firk_mem->hadapt_mem->etamax = ONE;
  firk_mem->jac_stale_fail     = !firk_mem->jac_current;
  firk_mem->dgamma_fail        = (firk_mem->gammap != ZERO)
                                   ? SUNRabs(firk_mem->gamma / firk_mem->gammap - ONE)
                                   : ONE;

  if (firk_mem->fixedstep)
  {
    too_many = (*ncfPtr > 1) || firk_mem->jac_current;
  }
  else
  {
    too_many = (SUNRabs(firk_mem->h) <= firkHminEff(firk_mem) * ONEPSM) ||
               (*ncfPtr == firk_mem->maxncf);
  }
  if (too_many)
  {
    if (nls == FIRK_NLS_RHS_RECVR) { return FIRK_REPTD_RHSFUNC_ERR; }
    return FIRK_CONV_FAILURE;
  }

  /* reduce the step size (not in fixed-step mode) and retry */
  if (!firk_mem->fixedstep)
  {
    firk_mem->eta    = SUNMAX(firk_mem->nls_hfac,
                              firkHminEff(firk_mem) / SUNRabs(firk_mem->h));
    firk_mem->hprime = firk_mem->h * firk_mem->eta;
  }
  *nflagPtr = PREV_CONV_FAIL;

  return PREDICT_AGAIN;
}

/*---------------------------------------------------------------
  firkErrorEstimate

  Computes the RADAU5-style filtered local error estimate
    ele = (M - gamma J)^{-1} [ gamma f(tn,yn) + sum_i e_i M Z_i ]
  and its WRMS norm dsm (design Listing 8.5). On the first step, and
  after any rejected attempt of this step, an estimate with dsm >= 1
  is filtered again with f evaluated at yn + ele, as in RADAU5; a
  recoverable failure of that RHS evaluation or solve keeps the first
  estimate (H20). Returns FIRK_SUCCESS, +1 after a recoverable failure
  of the first solve, or a negative failure code.
  ---------------------------------------------------------------*/
static int firkErrorEstimateBody(FIRKodeMem firk_mem,
                                 sunbooleantype attempt_rejected,
                                 sunrealtype* dsmPtr)
{
  int i, s, retval;
  N_Vector w;
  sunrealtype dsm;
  SUNDIALS_MAYBE_UNUSED sunbooleantype refiltered;

  s = firk_mem->s;

  /* w = sum_i e_i Z_i, then M w if a mass matrix is present */
  for (i = 0; i < s; i++)
  {
    firk_mem->cvals[i] = firk_mem->tbl->e[i];
    firk_mem->Xvecs[i] = firk_mem->Z[i];
  }
  retval = N_VLinearCombination(s, firk_mem->cvals, firk_mem->Xvecs,
                                firk_mem->tempv2);
  if (retval != 0) { return FIRK_VECTOROP_ERR; }
  w = firk_mem->tempv2;
  if (firk_mem->mass_set)
  {
    retval = firk_mem->mmult(firk_mem, firk_mem->tempv2, firk_mem->tempv3);
    if (retval != 0) { return FIRK_MASSMULT_FAIL; }
    w = firk_mem->tempv3;
  }

  /* filtered estimate */
  N_VLinearSum(firk_mem->gamma, firk_mem->fn, ONE, w, firk_mem->ele);
  retval = firk_mem->lsolve_blk(firk_mem, firk_mem->ele, firk_mem->ele,
                                FIRKLS_BLOCK_ESTIMATE, 0);
  if (retval < 0) { return FIRK_LSOLVE_FAIL; }
  if (retval > 0) { return 1; }
  dsm = N_VWrmsNorm(firk_mem->ele, firk_mem->ewt);

  /* the refilter, into tempv1 so that a failure keeps the first estimate */
  refiltered = SUNFALSE;
  if (dsm >= ONE && firk_mem->refilter &&
      (firk_mem->firststage || attempt_rejected))
  {
    N_VLinearSum(ONE, firk_mem->yn, ONE, firk_mem->ele, firk_mem->tempv1);
    retval = firkRhs(firk_mem, firk_mem->tn, firk_mem->tempv1, firk_mem->ftemp);
    if (retval < 0) { return retval; }
    if (retval == FIRK_SUCCESS)
    {
      N_VLinearSum(firk_mem->gamma, firk_mem->ftemp, ONE, w, firk_mem->tempv1);
      retval = firk_mem->lsolve_blk(firk_mem, firk_mem->tempv1,
                                    firk_mem->tempv1, FIRKLS_BLOCK_ESTIMATE, 0);
      if (retval < 0) { return FIRK_LSOLVE_FAIL; }
      if (retval == 0)
      {
        N_VScale(ONE, firk_mem->tempv1, firk_mem->ele);
        dsm        = N_VWrmsNorm(firk_mem->ele, firk_mem->ewt);
        refiltered = SUNTRUE;
      }
    }
  }

  SUNLogDebug(FIRK_LOGGER, "error-estimate",
              "est-norm = " SUN_FORMAT_G ", refiltered = %i", dsm,
              (int)refiltered);

  *dsmPtr = SUNMAX(dsm, DSM_FLOOR);

  return FIRK_SUCCESS;
}

/* the profiling region of the estimate (design Table 8.5) */
static int firkErrorEstimate(FIRKodeMem firk_mem,
                             sunbooleantype attempt_rejected, sunrealtype* dsmPtr)
{
  int retval;
  SUNDIALS_MARK_FUNCTION_BEGIN(FIRK_PROFILER);
  retval = firkErrorEstimateBody(firk_mem, attempt_rejected, dsmPtr);
  SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
  return retval;
}

/*---------------------------------------------------------------
  firkDoErrorTest

  Selects the next step size and performs the local error test
  (design Listing 8.6). Returns FIRK_SUCCESS if the test passes,
  TRY_AGAIN to retry the step, or FIRK_ERR_FAILURE /
  FIRK_CONTROLLER_ERR.
  ---------------------------------------------------------------*/
static int firkDoErrorTest(FIRKodeMem firk_mem, int* nflagPtr, int* nefPtr,
                           sunrealtype dsm)
{
  int retval;
  FIRKodeHAdaptMem hadapt_mem = firk_mem->hadapt_mem;

  /* the next step size, once per attempt and before the controller's
     history update in firkCompleteStep (H36) */
  retval = firkAdapt(firk_mem, dsm, dsm <= ONE);
  if (retval != FIRK_SUCCESS) { return retval; }

  if (dsm <= ONE) { return FIRK_SUCCESS; }

  /* test failed */
  (*nefPtr)++;
  firk_mem->netf++;
  *nflagPtr = PREV_ERR_FAIL;

  if ((*nefPtr == firk_mem->maxnef) ||
      (SUNRabs(firk_mem->h) <= firkHminEff(firk_mem) * ONEPSM))
  {
    return FIRK_ERR_FAILURE;
  }

  /* RADAU5 reduces a rejected first step by 10; repeated failures are
     bounded by etamxf (an ARKODE safeguard) */
  if (firk_mem->firststage)
  {
    firk_mem->eta = SUNMIN(firk_mem->eta, SUN_RCONST(0.1));
  }
  else if (*nefPtr >= hadapt_mem->small_nef)
  {
    firk_mem->eta = SUNMIN(firk_mem->eta, hadapt_mem->etamxf);
  }
  firk_mem->eta = SUNMAX(firk_mem->eta,
                         firkHminEff(firk_mem) / SUNRabs(firk_mem->h));

  /* no growth on the next accepted attempt of this step */
  hadapt_mem->etamax = ONE;

  firk_mem->hprime = firk_mem->h * firk_mem->eta;

  return TRY_AGAIN;
}

/*---------------------------------------------------------------
  firkCompleteStep

  Updates the integrator state after a successful step: notifies
  the step size controller, advances the solution and the time, and
  stores the stage increments for dense output and the next
  predictor. A failure leaves the state of the last step intact.
  ---------------------------------------------------------------*/
static int firkCompleteStep(FIRKodeMem firk_mem, sunrealtype dsm)
{
  int i, s, retval;
  sunrealtype troundoff;
  N_Vector* atmp;

  s = firk_mem->s;

  /* notify the step size controller */
  if (!firk_mem->fixedstep && firk_mem->hadapt_mem->hcontroller != NULL)
  {
    retval = SUNAdaptController_UpdateH(firk_mem->hadapt_mem->hcontroller,
                                        firk_mem->h, dsm);
    if (retval != SUN_SUCCESS)
    {
      firkProcessError(firk_mem, FIRK_CONTROLLER_ERR, __LINE__, __func__,
                       __FILE__, "Failure updating controller object");
      return FIRK_CONTROLLER_ERR;
    }
  }

  /* the new solution yn + sum_i d_i Z_i, in place */
  firk_mem->cvals[0] = ONE;
  firk_mem->Xvecs[0] = firk_mem->yn;
  for (i = 0; i < s; i++)
  {
    firk_mem->cvals[i + 1] = firk_mem->tbl->d[i];
    firk_mem->Xvecs[i + 1] = firk_mem->Z[i];
  }
  retval = N_VLinearCombination(s + 1, firk_mem->cvals, firk_mem->Xvecs,
                                firk_mem->yn);
  if (retval != 0) { return FIRK_VECTOROP_ERR; }
  firk_mem->fn_is_current  = SUNFALSE;
  firk_mem->ydn_is_current = SUNFALSE;
  firk_mem->jac_current    = SUNFALSE;

  firk_mem->tcur = firk_mem->tn + firk_mem->h;

  if (firk_mem->tstopset)
  {
    troundoff = FUZZ_FACTOR * firk_mem->uround *
                (SUNRabs(firk_mem->tcur) + SUNRabs(firk_mem->h));
    if (SUNRabs(firk_mem->tcur - firk_mem->tstop) <= troundoff)
    {
      firk_mem->tcur = firk_mem->tstop;
    }
  }

  /* keep the accepted stage increments (Z <-> Zprev) for dense output and
     the next predictor; no stack holds Z or Zprev */
  atmp                = firk_mem->Zprev;
  firk_mem->Zprev     = firk_mem->Z;
  firk_mem->Z         = atmp;
  firk_mem->have_prev = SUNTRUE;
  firk_mem->hold      = firk_mem->h;

  firk_mem->nst++;
  firk_mem->tn = firk_mem->tcur;
  if (firk_mem->fixedstep)
  {
    /* return to the fixed step size, which the step to tstop may have
       shortened */
    firk_mem->eta    = ONE;
    firk_mem->hprime = firk_mem->hin;
  }
  else { firk_mem->hprime = firk_mem->h * firk_mem->eta; }
  firk_mem->dsm = dsm;

  firk_mem->hadapt_mem->etamax = firk_mem->hadapt_mem->growth;

  firk_mem->initsetup  = SUNFALSE;
  firk_mem->firststage = SUNFALSE;

  return FIRK_SUCCESS;
}
