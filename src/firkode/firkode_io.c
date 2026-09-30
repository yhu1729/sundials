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
 * This is the implementation file for the optional input and output
 * functions of the FIRKODE integrator, and for its time step
 * adaptivity utilities.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include <sundials/sundials_types.h>

#include "firkode_impl.h"
#include "firkode_ls_impl.h"
#include "firkode_root_impl.h"
#include "sundials_utils.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)

/* access the integrator memory or return an error */
#define FIRK_ACCESS_MEM(firkode_mem, firk_mem)                          \
  if ((firkode_mem) == NULL)                                            \
  {                                                                     \
    firkProcessError(NULL, FIRK_MEM_NULL, __LINE__, __func__, __FILE__, \
                     MSGFIRK_NO_MEM);                                   \
    return FIRK_MEM_NULL;                                               \
  }                                                                     \
  (firk_mem) = (FIRKodeMem)(firkode_mem)

/*===============================================================
  Method selection
  ===============================================================*/

int FIRKodeSetNumStages(void* firkode_mem, int s)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);

  if (s < 1 || s > FIRK_MAX_STAGES)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_STAGES);
    return FIRK_ILL_INPUT;
  }

  if (firk_mem->MallocDone && !firk_mem->initsetup && s != firk_mem->s)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_STAGES_LOCKED);
    return FIRK_ILL_INPUT;
  }

  firk_mem->s = s;
  return FIRK_SUCCESS;
}

int FIRKodeSetOrder(void* firkode_mem, int ord)
{
  int s;

  if (ord <= 0) { s = FIRK_DEFAULT_STAGES; }
  else { s = ord / 2 + 1; } /* smallest s with 2s-1 >= ord; no overflow */

  return FIRKodeSetNumStages(firkode_mem, s);
}

int FIRKodeGetNumStages(void* firkode_mem, int* s)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *s = firk_mem->s;
  return FIRK_SUCCESS;
}

int FIRKodeGetCurrentTable(void* firkode_mem, FIRKodeTable* T)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);

  /* construct the table on demand so it can be inspected before the
     first step */
  if (firk_mem->T == NULL || firk_mem->T->s != firk_mem->s)
  {
    if (firk_mem->T != NULL) { FIRKodeTable_Free(firk_mem->T); }
    firk_mem->T = FIRKodeTable_RadauIIA(firk_mem->s);
    if (firk_mem->T == NULL)
    {
      firkProcessError(firk_mem, FIRK_TABLE_FAIL, __LINE__, __func__, __FILE__,
                       MSGFIRK_TABLE_FAIL);
      return FIRK_TABLE_FAIL;
    }
  }

  *T = firk_mem->T;
  return FIRK_SUCCESS;
}

/*===============================================================
  Optional inputs
  ===============================================================*/

int FIRKodeSetUserData(void* firkode_mem, void* user_data)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->user_data = user_data;
  return FIRK_SUCCESS;
}

int FIRKodeSetMaxNumSteps(void* firkode_mem, long int mxsteps)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  /* passing mxsteps = 0 sets the default, mxsteps < 0 disables the test */
  if (mxsteps == 0) { firk_mem->mxstep = MXSTEP_DEFAULT; }
  else { firk_mem->mxstep = mxsteps; }
  return FIRK_SUCCESS;
}

int FIRKodeSetMaxHnilWarns(void* firkode_mem, int mxhnil)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->mxhnil = mxhnil;
  return FIRK_SUCCESS;
}

int FIRKodeSetInitStep(void* firkode_mem, sunrealtype hin)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (firk_mem->fixedstep)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "The initial step cannot be set in fixed-step mode.");
    return FIRK_ILL_INPUT;
  }
  firk_mem->hin = hin;
  return FIRK_SUCCESS;
}

int FIRKodeSetMinStep(void* firkode_mem, sunrealtype hmin)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);

  if (hmin < ZERO)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_NEG_HMIN);
    return FIRK_ILL_INPUT;
  }

  if (hmin == ZERO)
  {
    firk_mem->hmin = HMIN_DEFAULT;
    return FIRK_SUCCESS;
  }

  if (hmin * firk_mem->hmax_inv > ONE)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_HMIN_HMAX);
    return FIRK_ILL_INPUT;
  }

  firk_mem->hmin = hmin;
  return FIRK_SUCCESS;
}

int FIRKodeSetMaxStep(void* firkode_mem, sunrealtype hmax)
{
  sunrealtype hmax_inv;
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);

  if (hmax < ZERO)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_NEG_HMAX);
    return FIRK_ILL_INPUT;
  }

  if (hmax == ZERO)
  {
    firk_mem->hmax_inv = HMAX_INV_DEFAULT;
    return FIRK_SUCCESS;
  }

  hmax_inv = ONE / hmax;
  if (hmax_inv * firk_mem->hmin > ONE)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_HMIN_HMAX);
    return FIRK_ILL_INPUT;
  }

  firk_mem->hmax_inv = hmax_inv;
  return FIRK_SUCCESS;
}

int FIRKodeSetFixedStep(void* firkode_mem, sunrealtype hfixed)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);

  if (hfixed == ZERO)
  {
    /* return to adaptive mode with default step bounds */
    firk_mem->fixedstep = SUNFALSE;
    firk_mem->hin       = ZERO;
    firk_mem->hmin      = HMIN_DEFAULT;
    firk_mem->hmax_inv  = HMAX_INV_DEFAULT;
    return FIRK_SUCCESS;
  }

  firk_mem->fixedstep = SUNTRUE;
  firk_mem->hin       = hfixed;
  firk_mem->hmin      = SUNRabs(hfixed);
  firk_mem->hmax_inv  = ONE / SUNRabs(hfixed);
  return FIRK_SUCCESS;
}

int FIRKodeSetStopTime(void* firkode_mem, sunrealtype tstop)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);

  /* if the integration has started, check that tstop is ahead of tcur */
  if (firk_mem->nst > 0 && !firk_mem->initsetup)
  {
    if ((tstop - firk_mem->tcur) * firk_mem->h < ZERO)
    {
      firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                       MSGFIRK_BAD_TSTOP, tstop, firk_mem->tcur);
      return FIRK_ILL_INPUT;
    }
  }

  firk_mem->tstop    = tstop;
  firk_mem->tstopset = SUNTRUE;
  return FIRK_SUCCESS;
}

int FIRKodeSetInterpolateStopTime(void* firkode_mem, sunbooleantype interp)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->tstopinterp = interp;
  return FIRK_SUCCESS;
}

int FIRKodeClearStopTime(void* firkode_mem)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->tstopset = SUNFALSE;
  return FIRK_SUCCESS;
}

int FIRKodeSetMaxErrTestFails(void* firkode_mem, int maxnef)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->maxnef = (maxnef <= 0) ? MXNEF : maxnef;
  return FIRK_SUCCESS;
}

int FIRKodeSetMaxConvFails(void* firkode_mem, int maxncf)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->maxncf = (maxncf <= 0) ? MXNCF : maxncf;
  return FIRK_SUCCESS;
}

int FIRKodeSetLSetupFrequency(void* firkode_mem, long int msbp)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (msbp < 0)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "A negative setup frequency is illegal.");
    return FIRK_ILL_INPUT;
  }
  firk_mem->msbp = (msbp == 0) ? MSBP_DEFAULT : msbp;
  return FIRK_SUCCESS;
}

int FIRKodeSetDeltaGammaMaxLSetup(void* firkode_mem, sunrealtype dgmax_lsetup)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->dgmax_lsetup = (dgmax_lsetup < ZERO) ? DGMAX_LSETUP_DEFAULT
                                                 : dgmax_lsetup;
  return FIRK_SUCCESS;
}

int FIRKodeSetPredictorMethod(void* firkode_mem, int method)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (method != FIRK_PREDICT_EXTRAPOLATE && method != FIRK_PREDICT_TRIVIAL)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Illegal predictor method.");
    return FIRK_ILL_INPUT;
  }
  firk_mem->predictor = method;
  return FIRK_SUCCESS;
}

int FIRKodeSetErrorRefilter(void* firkode_mem, sunbooleantype onoff)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->refilter = onoff;
  return FIRK_SUCCESS;
}

/*===============================================================
  Optional step adaptivity inputs
  ===============================================================*/

int FIRKodeSetAdaptController(void* firkode_mem, SUNAdaptController C)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  return firkReplaceAdaptController(firk_mem, C, SUNFALSE);
}

int FIRKodeSetSafetyFactor(void* firkode_mem, sunrealtype safety)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (safety > ONE)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Illegal safety factor");
    return FIRK_ILL_INPUT;
  }
  firk_mem->hadapt_mem->safety = (safety <= ZERO) ? SAFETY : safety;
  return FIRK_SUCCESS;
}

int FIRKodeSetErrorBias(void* firkode_mem, sunrealtype bias)
{
  int retval;
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);

  if (firk_mem->hadapt_mem->hcontroller == NULL)
  {
    firkProcessError(firk_mem, FIRK_MEM_NULL, __LINE__, __func__, __FILE__,
                     "SUNAdaptController NULL -- must be set before setting "
                     "the error bias");
    return FIRK_MEM_NULL;
  }

  retval = SUNAdaptController_SetErrorBias(firk_mem->hadapt_mem->hcontroller,
                                           (bias < ONE) ? -ONE : bias);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRK_CONTROLLER_ERR, __LINE__, __func__,
                     __FILE__, "SUNAdaptController_SetErrorBias failure");
    return FIRK_CONTROLLER_ERR;
  }
  return FIRK_SUCCESS;
}

int FIRKodeSetMaxGrowth(void* firkode_mem, sunrealtype mx_growth)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->hadapt_mem->growth = (mx_growth <= ONE) ? GROWTH : mx_growth;
  return FIRK_SUCCESS;
}

int FIRKodeSetMinReduction(void* firkode_mem, sunrealtype eta_min)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->hadapt_mem->etamin = (eta_min >= ONE || eta_min <= ZERO) ? ETAMIN
                                                                     : eta_min;
  return FIRK_SUCCESS;
}

int FIRKodeSetFixedStepBounds(void* firkode_mem, sunrealtype lb, sunrealtype ub)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if ((lb <= ONE) && (ub >= ONE))
  {
    firk_mem->hadapt_mem->lbound = lb;
    firk_mem->hadapt_mem->ubound = ub;
  }
  else
  {
    firk_mem->hadapt_mem->lbound = HFIXED_LB;
    firk_mem->hadapt_mem->ubound = HFIXED_UB;
  }
  return FIRK_SUCCESS;
}

int FIRKodeSetMaxFirstGrowth(void* firkode_mem, sunrealtype etamx1)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->hadapt_mem->etamx1 = (etamx1 <= ONE) ? ETAMX1 : etamx1;
  return FIRK_SUCCESS;
}

int FIRKodeSetMaxEFailGrowth(void* firkode_mem, sunrealtype etamxf)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->hadapt_mem->etamxf = ((etamxf <= ZERO) || (etamxf > ONE)) ? ETAMXF
                                                                      : etamxf;
  return FIRK_SUCCESS;
}

int FIRKodeSetSmallNumEFails(void* firkode_mem, int small_nef)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->hadapt_mem->small_nef = (small_nef <= 0) ? SMALL_NEF : small_nef;
  return FIRK_SUCCESS;
}

int FIRKodeSetMaxCFailGrowth(void* firkode_mem, sunrealtype etacf)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->hadapt_mem->etacf = ((etacf <= ZERO) || (etacf > ONE)) ? ETACF
                                                                   : etacf;
  return FIRK_SUCCESS;
}

int FIRKodeSetNewtonCountSafety(void* firkode_mem, sunbooleantype onoff)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->hadapt_mem->newt_safety = onoff;
  return FIRK_SUCCESS;
}

/*===============================================================
  Optional nonlinear solver inputs
  ===============================================================*/

int FIRKodeSetMaxNonlinIters(void* firkode_mem, int maxcor)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->maxcor = (maxcor <= 0) ? MAXCOR_DEFAULT : maxcor;
  if (firk_mem->NLS != NULL)
  {
    if (SUNNonlinSolSetMaxIters(firk_mem->NLS, firk_mem->maxcor) != SUN_SUCCESS)
    {
      firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                       "Setting maximum number of nonlinear iterations failed");
      return FIRK_ILL_INPUT;
    }
  }
  return FIRK_SUCCESS;
}

int FIRKodeSetNonlinConvCoef(void* firkode_mem, sunrealtype nlscoef)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->nlscoef = (nlscoef <= ZERO) ? NLSCOEF_DEFAULT : nlscoef;
  return FIRK_SUCCESS;
}

int FIRKodeSetNonlinDivergenceRate(void* firkode_mem, sunrealtype theta_max)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->theta_max = (theta_max <= ZERO) ? THETA_MAX_DEFAULT : theta_max;
  return FIRK_SUCCESS;
}

int FIRKodeSetJacBadConvRate(void* firkode_mem, sunrealtype theta_jbad)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  firk_mem->theta_jbad = (theta_jbad <= ZERO) ? THETA_JBAD_DEFAULT : theta_jbad;
  return FIRK_SUCCESS;
}

/*===============================================================
  Optional outputs
  ===============================================================*/

int FIRKodeGetNumSteps(void* firkode_mem, long int* nsteps)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *nsteps = firk_mem->nst;
  return FIRK_SUCCESS;
}

int FIRKodeGetNumStepAttempts(void* firkode_mem, long int* step_attempts)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *step_attempts = firk_mem->nst_attempts;
  return FIRK_SUCCESS;
}

int FIRKodeGetNumRhsEvals(void* firkode_mem, long int* nfevals)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *nfevals = firk_mem->nfe;
  return FIRK_SUCCESS;
}

int FIRKodeGetNumLinSolvSetups(void* firkode_mem, long int* nlinsetups)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *nlinsetups = firk_mem->nsetups;
  return FIRK_SUCCESS;
}

int FIRKodeGetNumErrTestFails(void* firkode_mem, long int* netfails)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *netfails = firk_mem->netf;
  return FIRK_SUCCESS;
}

int FIRKodeGetNumStepSolveFails(void* firkode_mem, long int* nncfails)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *nncfails = firk_mem->ncfn;
  return FIRK_SUCCESS;
}

int FIRKodeGetNumNonlinSolvIters(void* firkode_mem, long int* nniters)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *nniters = firk_mem->nni;
  return FIRK_SUCCESS;
}

int FIRKodeGetNumNonlinSolvConvFails(void* firkode_mem, long int* nnfails)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *nnfails = firk_mem->nnf;
  return FIRK_SUCCESS;
}

int FIRKodeGetNonlinSolvStats(void* firkode_mem, long int* nniters,
                              long int* nnfails)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *nniters = firk_mem->nni;
  *nnfails = firk_mem->nnf;
  return FIRK_SUCCESS;
}

int FIRKodeGetActualInitStep(void* firkode_mem, sunrealtype* hinused)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *hinused = firk_mem->h0u;
  return FIRK_SUCCESS;
}

int FIRKodeGetLastStep(void* firkode_mem, sunrealtype* hlast)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *hlast = firk_mem->hold;
  return FIRK_SUCCESS;
}

int FIRKodeGetCurrentStep(void* firkode_mem, sunrealtype* hcur)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *hcur = firk_mem->next_h;
  return FIRK_SUCCESS;
}

int FIRKodeGetCurrentTime(void* firkode_mem, sunrealtype* tcur)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *tcur = firk_mem->tcur;
  return FIRK_SUCCESS;
}

int FIRKodeGetCurrentState(void* firkode_mem, N_Vector* y)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *y = firk_mem->yn;
  return FIRK_SUCCESS;
}

int FIRKodeGetCurrentGamma(void* firkode_mem, sunrealtype* gamma)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *gamma = firk_mem->gamma;
  return FIRK_SUCCESS;
}

int FIRKodeGetTolScaleFactor(void* firkode_mem, sunrealtype* tolsfac)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *tolsfac = firk_mem->tolsf;
  return FIRK_SUCCESS;
}

int FIRKodeGetErrWeights(void* firkode_mem, N_Vector eweight)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  N_VScale(ONE, firk_mem->ewt, eweight);
  return FIRK_SUCCESS;
}

int FIRKodeGetEstLocalErrors(void* firkode_mem, N_Vector ele)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  N_VScale(ONE, firk_mem->ele, ele);
  return FIRK_SUCCESS;
}

int FIRKodeGetIntegratorStats(void* firkode_mem, long int* nsteps,
                              long int* nfevals, long int* nlinsetups,
                              long int* netfails, sunrealtype* hinused,
                              sunrealtype* hlast, sunrealtype* hcur,
                              sunrealtype* tcur)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *nsteps     = firk_mem->nst;
  *nfevals    = firk_mem->nfe;
  *nlinsetups = firk_mem->nsetups;
  *netfails   = firk_mem->netf;
  *hinused    = firk_mem->h0u;
  *hlast      = firk_mem->hold;
  *hcur       = firk_mem->next_h;
  *tcur       = firk_mem->tcur;
  return FIRK_SUCCESS;
}

int FIRKodeGetUserData(void* firkode_mem, void** user_data)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *user_data = firk_mem->user_data;
  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  FIRKodePrintAllStats

  Prints all integrator statistics.
  ---------------------------------------------------------------*/
int FIRKodePrintAllStats(void* firkode_mem, FILE* outfile, SUNOutputFormat fmt)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  FIRKLsMassMem firkls_massmem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);

  if (fmt != SUN_OUTPUTFORMAT_TABLE && fmt != SUN_OUTPUTFORMAT_CSV)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Invalid formatting option.");
    return FIRK_ILL_INPUT;
  }

  /* step and method stats */
  sunfprintf_real(outfile, fmt, SUNTRUE, "Current time", firk_mem->tcur);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Steps", firk_mem->nst);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Step attempts",
                  firk_mem->nst_attempts);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Stability limited steps", 0);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Accuracy limited steps",
                  firk_mem->hadapt_mem->nst_acc);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Error test fails", firk_mem->netf);
  sunfprintf_long(outfile, fmt, SUNFALSE, "NLS step fails", firk_mem->ncfn);
  sunfprintf_real(outfile, fmt, SUNFALSE, "Initial step size", firk_mem->h0u);
  sunfprintf_real(outfile, fmt, SUNFALSE, "Last step size", firk_mem->hold);
  sunfprintf_real(outfile, fmt, SUNFALSE, "Current step size", firk_mem->next_h);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Stages", firk_mem->s);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Method order",
                  (firk_mem->T != NULL) ? firk_mem->T->q : 2 * firk_mem->s - 1);

  /* function evaluations */
  sunfprintf_long(outfile, fmt, SUNFALSE, "RHS fn evals", firk_mem->nfe);

  /* nonlinear solver stats */
  sunfprintf_long(outfile, fmt, SUNFALSE, "NLS iters", firk_mem->nni);
  sunfprintf_long(outfile, fmt, SUNFALSE, "NLS fails", firk_mem->nnf);
  if (firk_mem->nst > 0)
  {
    sunfprintf_real(outfile, fmt, SUNFALSE, "NLS iters per step",
                    (sunrealtype)firk_mem->nni / (sunrealtype)firk_mem->nst);
  }

  /* linear solver stats */
  sunfprintf_long(outfile, fmt, SUNFALSE, "LS setups", firk_mem->nsetups);
  if (firk_mem->lmem)
  {
    firkls_mem = (FIRKLsMem)(firk_mem->lmem);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Jac fn evals", firkls_mem->nje);
    sunfprintf_long(outfile, fmt, SUNFALSE, "LS RHS fn evals", firkls_mem->nfeDQ);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Prec setup evals", firkls_mem->npe);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Prec solves", firkls_mem->nps);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Block LS solves", firkls_mem->nbs);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Block LS iters", firkls_mem->nli);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Block LS fails", firkls_mem->ncfl);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Stage LS iters", firkls_mem->nsli);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Stage LS fails", firkls_mem->nscf);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Jac-times setups",
                    firkls_mem->njtsetup);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Jac-times evals",
                    firkls_mem->njtimes);
    if (firk_mem->nni > 0)
    {
      sunfprintf_real(outfile, fmt, SUNFALSE, "Stage LS iters per NLS iter",
                      (sunrealtype)firkls_mem->nsli / (sunrealtype)firk_mem->nni);
      sunfprintf_real(outfile, fmt, SUNFALSE, "Jac evals per NLS iter",
                      (sunrealtype)firkls_mem->nje / (sunrealtype)firk_mem->nni);
      sunfprintf_real(outfile, fmt, SUNFALSE, "Prec evals per NLS iter",
                      (sunrealtype)firkls_mem->npe / (sunrealtype)firk_mem->nni);
    }
  }

  /* mass matrix stats */
  if (firk_mem->mass_mem)
  {
    firkls_massmem = (FIRKLsMassMem)(firk_mem->mass_mem);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass setups",
                    firkls_massmem->nmsetups);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass solves",
                    firkls_massmem->nmsolves);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass mult evals",
                    firkls_massmem->nmvevals);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass LS iters",
                    firkls_massmem->nmli);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass LS fails",
                    firkls_massmem->nmcfails);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass prec setups",
                    firkls_massmem->nmpe);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass prec solves",
                    firkls_massmem->nmps);
  }

  /* rootfinding stats */
  if (firk_mem->root_mem)
  {
    sunfprintf_long(outfile, fmt, SUNFALSE, "Root fn evals",
                    firk_mem->root_mem->nge);
  }

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  FIRKodeWriteParameters

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int FIRKodeWriteParameters(void* firkode_mem, FILE* fp)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);

  fprintf(fp, "FIRKODE solver parameters:\n");
  fprintf(fp, "  Number of stages = %i\n", firk_mem->s);
  fprintf(fp, "  Method order = %i\n", 2 * firk_mem->s - 1);
  if (firk_mem->itol == FIRK_SS)
  {
    fprintf(fp, "  Solver relative tolerance = " SUN_FORMAT_G "\n",
            firk_mem->reltol);
    fprintf(fp, "  Solver absolute tolerance = " SUN_FORMAT_G "\n",
            firk_mem->Sabstol);
  }
  else if (firk_mem->itol == FIRK_SV)
  {
    fprintf(fp, "  Solver relative tolerance = " SUN_FORMAT_G "\n",
            firk_mem->reltol);
    fprintf(fp, "  Vector-valued solver absolute tolerance\n");
  }
  else if (firk_mem->itol == FIRK_WF)
  {
    fprintf(fp, "  User provided error weight function\n");
  }
  if (firk_mem->fixedstep)
  {
    fprintf(fp, "  Fixed step size = " SUN_FORMAT_G "\n", firk_mem->hin);
  }
  else
  {
    if (firk_mem->hin != ZERO)
    {
      fprintf(fp, "  Initial step size = " SUN_FORMAT_G "\n", firk_mem->hin);
    }
    if (firk_mem->hmin != ZERO)
    {
      fprintf(fp, "  Minimum step size = " SUN_FORMAT_G "\n", firk_mem->hmin);
    }
    if (firk_mem->hmax_inv != ZERO)
    {
      fprintf(fp, "  Maximum step size = " SUN_FORMAT_G "\n",
              ONE / firk_mem->hmax_inv);
    }
    fprintf(fp, "  Predictor method = %i\n", firk_mem->predictor);
    fprintf(fp, "  Refilter error estimate = %i\n", firk_mem->refilter);
    firkPrintAdaptMem(firk_mem->hadapt_mem, fp);
  }
  fprintf(fp, "  Maximum number of steps = %li\n", firk_mem->mxstep);
  fprintf(fp, "  Maximum error test failures = %i\n", firk_mem->maxnef);
  fprintf(fp, "  Maximum convergence failures = %i\n", firk_mem->maxncf);
  fprintf(fp, "  Maximum nonlinear iterations = %i\n", firk_mem->maxcor);
  fprintf(fp, "  Nonlinear convergence coefficient = " SUN_FORMAT_G "\n",
          firk_mem->nlscoef);
  fprintf(fp, "  Nonlinear divergence rate = " SUN_FORMAT_G "\n",
          firk_mem->theta_max);
  fprintf(fp, "  Linear solver setup frequency = %li\n", firk_mem->msbp);
  fprintf(fp, "  Linear solver setup threshold = " SUN_FORMAT_G "\n",
          firk_mem->dgmax_lsetup);
  fprintf(fp, "\n");

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  FIRKodeGetReturnFlagName

  Returns the name of the FIRKODE return flag.
  ---------------------------------------------------------------*/
char* FIRKodeGetReturnFlagName(long int flag)
{
  char* name;

  name = (char*)malloc(26 * sizeof(char));
  if (name == NULL) { return NULL; }

  switch (flag)
  {
  case FIRK_SUCCESS: sprintf(name, "FIRK_SUCCESS"); break;
  case FIRK_TSTOP_RETURN: sprintf(name, "FIRK_TSTOP_RETURN"); break;
  case FIRK_ROOT_RETURN: sprintf(name, "FIRK_ROOT_RETURN"); break;
  case FIRK_WARNING: sprintf(name, "FIRK_WARNING"); break;
  case FIRK_TOO_MUCH_WORK: sprintf(name, "FIRK_TOO_MUCH_WORK"); break;
  case FIRK_TOO_MUCH_ACC: sprintf(name, "FIRK_TOO_MUCH_ACC"); break;
  case FIRK_ERR_FAILURE: sprintf(name, "FIRK_ERR_FAILURE"); break;
  case FIRK_CONV_FAILURE: sprintf(name, "FIRK_CONV_FAILURE"); break;
  case FIRK_LINIT_FAIL: sprintf(name, "FIRK_LINIT_FAIL"); break;
  case FIRK_LSETUP_FAIL: sprintf(name, "FIRK_LSETUP_FAIL"); break;
  case FIRK_LSOLVE_FAIL: sprintf(name, "FIRK_LSOLVE_FAIL"); break;
  case FIRK_RHSFUNC_FAIL: sprintf(name, "FIRK_RHSFUNC_FAIL"); break;
  case FIRK_FIRST_RHSFUNC_ERR: sprintf(name, "FIRK_FIRST_RHSFUNC_ERR"); break;
  case FIRK_REPTD_RHSFUNC_ERR: sprintf(name, "FIRK_REPTD_RHSFUNC_ERR"); break;
  case FIRK_UNREC_RHSFUNC_ERR: sprintf(name, "FIRK_UNREC_RHSFUNC_ERR"); break;
  case FIRK_RTFUNC_FAIL: sprintf(name, "FIRK_RTFUNC_FAIL"); break;
  case FIRK_NLS_INIT_FAIL: sprintf(name, "FIRK_NLS_INIT_FAIL"); break;
  case FIRK_NLS_SETUP_FAIL: sprintf(name, "FIRK_NLS_SETUP_FAIL"); break;
  case FIRK_NLS_FAIL: sprintf(name, "FIRK_NLS_FAIL"); break;
  case FIRK_MASSINIT_FAIL: sprintf(name, "FIRK_MASSINIT_FAIL"); break;
  case FIRK_MASSSETUP_FAIL: sprintf(name, "FIRK_MASSSETUP_FAIL"); break;
  case FIRK_MASSSOLVE_FAIL: sprintf(name, "FIRK_MASSSOLVE_FAIL"); break;
  case FIRK_MASSMULT_FAIL: sprintf(name, "FIRK_MASSMULT_FAIL"); break;
  case FIRK_MEM_FAIL: sprintf(name, "FIRK_MEM_FAIL"); break;
  case FIRK_MEM_NULL: sprintf(name, "FIRK_MEM_NULL"); break;
  case FIRK_ILL_INPUT: sprintf(name, "FIRK_ILL_INPUT"); break;
  case FIRK_NO_MALLOC: sprintf(name, "FIRK_NO_MALLOC"); break;
  case FIRK_BAD_K: sprintf(name, "FIRK_BAD_K"); break;
  case FIRK_BAD_T: sprintf(name, "FIRK_BAD_T"); break;
  case FIRK_BAD_DKY: sprintf(name, "FIRK_BAD_DKY"); break;
  case FIRK_TOO_CLOSE: sprintf(name, "FIRK_TOO_CLOSE"); break;
  case FIRK_VECTOROP_ERR: sprintf(name, "FIRK_VECTOROP_ERR"); break;
  case FIRK_CONTROLLER_ERR: sprintf(name, "FIRK_CONTROLLER_ERR"); break;
  case FIRK_TABLE_FAIL: sprintf(name, "FIRK_TABLE_FAIL"); break;
  case FIRK_CONTEXT_ERR: sprintf(name, "FIRK_CONTEXT_ERR"); break;
  case FIRK_UNRECOGNIZED_ERR: sprintf(name, "FIRK_UNRECOGNIZED_ERR"); break;
  default: sprintf(name, "NONE");
  }

  return name;
}

/*===============================================================
  Time step adaptivity utilities
  ===============================================================*/

/*---------------------------------------------------------------
  firkAdaptInit

  Creates a zeroed adaptivity memory structure (defaults are set
  by FIRKodeCreate).
  ---------------------------------------------------------------*/
FIRKodeHAdaptMem firkAdaptInit(void)
{
  FIRKodeHAdaptMem hadapt_mem;

  hadapt_mem = (FIRKodeHAdaptMem)malloc(sizeof(struct FIRKodeHAdaptMemRec));
  if (hadapt_mem == NULL) { return NULL; }
  memset(hadapt_mem, 0, sizeof(struct FIRKodeHAdaptMemRec));
  return hadapt_mem;
}

/*---------------------------------------------------------------
  firkPrintAdaptMem

  Outputs the time step adaptivity parameters.
  ---------------------------------------------------------------*/
void firkPrintAdaptMem(FIRKodeHAdaptMem hadapt_mem, FILE* outfile)
{
  if (hadapt_mem == NULL) { return; }

  fprintf(outfile, "  Step adaptivity: etamx1 = " SUN_FORMAT_G "\n",
          hadapt_mem->etamx1);
  fprintf(outfile, "  Step adaptivity: etamxf = " SUN_FORMAT_G "\n",
          hadapt_mem->etamxf);
  fprintf(outfile, "  Step adaptivity: etamin = " SUN_FORMAT_G "\n",
          hadapt_mem->etamin);
  fprintf(outfile, "  Step adaptivity: small_nef = %i\n", hadapt_mem->small_nef);
  fprintf(outfile, "  Step adaptivity: etacf = " SUN_FORMAT_G "\n",
          hadapt_mem->etacf);
  fprintf(outfile, "  Step adaptivity: safety = " SUN_FORMAT_G "\n",
          hadapt_mem->safety);
  fprintf(outfile, "  Step adaptivity: growth = " SUN_FORMAT_G "\n",
          hadapt_mem->growth);
  fprintf(outfile, "  Step adaptivity: lbound = " SUN_FORMAT_G "\n",
          hadapt_mem->lbound);
  fprintf(outfile, "  Step adaptivity: ubound = " SUN_FORMAT_G "\n",
          hadapt_mem->ubound);
  fprintf(outfile, "  Step adaptivity: Newton count safety = %i\n",
          hadapt_mem->newt_safety);
  if (hadapt_mem->hcontroller != NULL)
  {
    (void)SUNAdaptController_Write(hadapt_mem->hcontroller, outfile);
  }
}

/*---------------------------------------------------------------
  firkAdapt

  Time step adaptivity wrapper: computes firk_mem->eta from the
  controller estimate, the safety factors, the step growth bounds,
  and the no-change deadband that avoids refactorizations.
  ---------------------------------------------------------------*/
int firkAdapt(FIRKodeMem firk_mem, sunrealtype dsm)
{
  int retval;
  sunrealtype h_acc, hcur, fac;
  FIRKodeHAdaptMem hadapt_mem = firk_mem->hadapt_mem;

  hcur = firk_mem->h;

  if (firk_mem->fixedstep || hadapt_mem->hcontroller == NULL)
  {
    firk_mem->eta = ONE;
    return FIRK_SUCCESS;
  }

  retval = SUNAdaptController_EstimateStep(hadapt_mem->hcontroller, hcur,
                                           hadapt_mem->p, dsm, &h_acc);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRK_CONTROLLER_ERR, __LINE__, __func__,
                     __FILE__, "SUNAdaptController_EstimateStep failure.");
    return FIRK_CONTROLLER_ERR;
  }

  SUNLogDebug(FIRK_LOGGER, "new-step-before-bounds",
              "h_acc_controller = " SUN_FORMAT_G, h_acc);

  /* safety factor, scaled by the Newton iteration count (RADAU5) */
  fac = hadapt_mem->safety;
  if (hadapt_mem->newt_safety)
  {
    fac = SUNMIN(fac, (TWO * (sunrealtype)firk_mem->maxcor + ONE) /
                        (TWO * (sunrealtype)firk_mem->maxcor +
                         (sunrealtype)firk_mem->nni_last));
  }
  h_acc *= fac;

  /* enforce bounds on the step size growth and reduction */
  h_acc = SUNMIN(SUNRabs(h_acc), SUNRabs(hadapt_mem->etamax * hcur));
  h_acc = SUNMAX(h_acc, SUNRabs(hadapt_mem->etamin * hcur));
  hadapt_mem->nst_acc++;

  /* deadband: keep h unchanged (and hence the factorization) */
  if (dsm <= ONE)
  {
    if ((h_acc > SUNRabs(hcur * hadapt_mem->lbound * ONEMSM)) &&
        (h_acc < SUNRabs(hcur * hadapt_mem->ubound * ONEPSM)))
    {
      h_acc = SUNRabs(hcur);
    }
  }
  h_acc = SUNRcopysign(h_acc, hcur);

  firk_mem->eta = h_acc / hcur;
  firk_mem->eta = SUNMAX(firk_mem->eta, firk_mem->hmin / SUNRabs(hcur));
  firk_mem->eta /= SUNMAX(ONE,
                          SUNRabs(hcur) * firk_mem->hmax_inv * firk_mem->eta);

  SUNLogDebug(FIRK_LOGGER, "new-step-eta", "eta = " SUN_FORMAT_G, firk_mem->eta);

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkReplaceAdaptController

  Replaces the current SUNAdaptController; a NULL input installs
  the default implicit Gustafsson controller.
  ---------------------------------------------------------------*/
int firkReplaceAdaptController(FIRKodeMem firk_mem, SUNAdaptController C,
                               sunbooleantype take_ownership)
{
  int retval;

  if (firk_mem->hadapt_mem->owncontroller &&
      (firk_mem->hadapt_mem->hcontroller != NULL))
  {
    retval = SUNAdaptController_Destroy(firk_mem->hadapt_mem->hcontroller);
    firk_mem->hadapt_mem->owncontroller = SUNFALSE;
    if (retval != SUN_SUCCESS)
    {
      firkProcessError(firk_mem, FIRK_MEM_FAIL, __LINE__, __func__, __FILE__,
                       "SUNAdaptController_Destroy failure");
      return FIRK_MEM_FAIL;
    }
  }
  firk_mem->hadapt_mem->hcontroller = NULL;

  if (C == NULL)
  {
    C = SUNAdaptController_ImpGus(firk_mem->sunctx);
    if (C == NULL)
    {
      firkProcessError(firk_mem, FIRK_MEM_FAIL, __LINE__, __func__, __FILE__,
                       "SUNAdaptController_ImpGus allocation failure");
      return FIRK_MEM_FAIL;
    }
    firk_mem->hadapt_mem->owncontroller = SUNTRUE;
  }
  else { firk_mem->hadapt_mem->owncontroller = take_ownership; }

  firk_mem->hadapt_mem->hcontroller = C;

  return FIRK_SUCCESS;
}
