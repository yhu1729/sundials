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
 * functions of the FIRKODE integrator.
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
  Mode selection
  ===============================================================*/

/* Returns SUNTRUE if the EXPERT option set by the function func may be set,
   and reports an error otherwise (design decision D8) */
sunbooleantype firkExpertAllowed(FIRKodeMem firk_mem, const char* func)
{
  if (firk_mem->mode == FIRK_EXPERT) { return SUNTRUE; }
  firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, func, __FILE__,
                   MSGFIRK_REQUIRES_EXPERT, func);
  return SUNFALSE;
}

int FIRKodeSetMode(void* firkode_mem, int mode)
{
  int retval;
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);

  if (mode != FIRK_BEGINNER && mode != FIRK_EXPERT)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_MODE);
    return FIRK_ILL_INPUT;
  }

  /* the mode is fixed once a step has been attempted (ReInit clears this) */
  if (firk_mem->nst_attempts > 0)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_MODE_LOCKED);
    return FIRK_ILL_INPUT;
  }

  /* BEGINNER runs use the defaults of every EXPERT option: those of the
     linear solver interface first, whose restore can only fail in the option
     setters of the stacked solver, before anything of the integrator is
     touched, and whose code is mapped to the integrator's; on any failure
     the mode is unchanged */
  if (mode == FIRK_BEGINNER && firk_mem->mode == FIRK_EXPERT)
  {
    if (firk_mem->lmem != NULL &&
        firkLsSetExpertDefaults(firk_mem) != FIRKLS_SUCCESS)
    {
      return FIRK_LINIT_FAIL;
    }
    retval = firkSetExpertDefaults(firk_mem);
    if (retval != FIRK_SUCCESS) { return retval; }
  }

  firk_mem->mode = mode;
  return FIRK_SUCCESS;
}

int FIRKodeGetMode(void* firkode_mem, int* mode_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *mode_ptr = firk_mem->mode;
  return FIRK_SUCCESS;
}

/*===============================================================
  Method selection
  ===============================================================*/

/* The stage count for the order ord: the smallest s with 2s - 1 >= ord,
   FIRK_DEFAULT_STAGES for ord <= 0, and -1 above the largest order */
int firkOrderToStages(int ord)
{
  if (ord <= 0) { return FIRK_DEFAULT_STAGES; }
  if (ord > 2 * FIRK_MAX_STAGES - 1) { return -1; }
  return ord / 2 + 1;
}

/* Selects any of the built-in Radau IIA methods, 1 <= s <= FIRK_MAX_STAGES.
   Version 0.1.0 exposes only s = 3 (design decision D7); the general-s code
   path is kept for the later versions and tested through this function. */
int firkSetNumStagesAny(void* firkode_mem, int s)
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

  if (firk_mem->tbl == NULL || firk_mem->tbl->s != s)
  {
    FIRKodeTable_Free(firk_mem->tbl);
    firk_mem->tbl = FIRKodeTable_RadauIIA(s);
    if (firk_mem->tbl == NULL)
    {
      firkProcessError(firk_mem, FIRK_MEM_FAIL, __LINE__, __func__, __FILE__,
                       MSGFIRK_MEM_FAIL);
      return FIRK_MEM_FAIL;
    }
  }
  firk_mem->s = s;
  return FIRK_SUCCESS;
}

int FIRKodeSetNumStages(void* firkode_mem, int s)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }

  if (s != FIRK_SUPPORTED_STAGES)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_ONLY_STAGES);
    return FIRK_ILL_INPUT;
  }
  return firkSetNumStagesAny(firkode_mem, s);
}

int FIRKodeSetOrder(void* firkode_mem, int ord)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }

  if (firkOrderToStages(ord) != FIRK_SUPPORTED_STAGES)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_ONLY_ORDER);
    return FIRK_ILL_INPUT;
  }
  return firkSetNumStagesAny(firkode_mem, FIRK_SUPPORTED_STAGES);
}

int FIRKodeGetNumStages(void* firkode_mem, int* s_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *s_ptr = firk_mem->s;
  return FIRK_SUCCESS;
}

int FIRKodeGetCurrentTable(void* firkode_mem, FIRKodeTable* tbl_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);

  /* the integrator's table, not owned by the caller */
  *tbl_ptr = firk_mem->tbl;
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

  /* every callback that receives user_data, whenever this is called (H5) */
  if (firk_mem->user_efun) { firk_mem->e_data = user_data; }
  if (firk_mem->root_mem != NULL) { firk_mem->root_mem->root_data = user_data; }
  firkLsSetUserData(firk_mem, user_data);
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
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
  firk_mem->maxnef = (maxnef <= 0) ? MXNEF : maxnef;
  return FIRK_SUCCESS;
}

int FIRKodeSetMaxConvFails(void* firkode_mem, int maxncf)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
  firk_mem->maxncf = (maxncf <= 0) ? MXNCF : maxncf;
  return FIRK_SUCCESS;
}

int FIRKodeSetReusePolicy(void* firkode_mem, int policy)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
  if (policy != FIRK_POLICY_RADAU5 && policy != FIRK_POLICY_PERIODIC)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_POLICY);
    return FIRK_ILL_INPUT;
  }
  firk_mem->reuse_policy = policy;
  return FIRK_SUCCESS;
}

int FIRKodeSetLSetupFrequency(void* firkode_mem, long int msbp)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
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
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
  firk_mem->dgmax_lsetup = (dgmax_lsetup < ZERO) ? DGMAX_LSETUP_DEFAULT
                                                 : dgmax_lsetup;
  return FIRK_SUCCESS;
}

int FIRKodeSetPredictorMethod(void* firkode_mem, int method)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
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
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
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
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }

  /* the step size is computed with SUNAdaptController_EstimateStep, which
     leaves the step unchanged for a controller of any other type */
  if (C != NULL && SUNAdaptController_GetType(C) != SUN_ADAPTCONTROLLER_H)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_CONTROLLER_TYPE);
    return FIRK_ILL_INPUT;
  }

  return firkReplaceAdaptController(firk_mem, C, SUNFALSE);
}

/* Creates a controller of the Soderlind family by its name, the names of
   ARKodeSetAdaptControllerByName without ImExGus (a controller for IMEX
   methods, whose object library FIRKODE does not include) */
int FIRKodeSetAdaptControllerByName(void* firkode_mem, const char* cname)
{
  int retval;
  SUNAdaptController C = NULL;
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }

  if (cname == NULL)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_CONTROLLER_NAME);
    return FIRK_ILL_INPUT;
  }

  if (strcmp(cname, "Soderlind") == 0)
  {
    C = SUNAdaptController_Soderlind(firk_mem->sunctx);
  }
  else if (strcmp(cname, "PID") == 0)
  {
    C = SUNAdaptController_PID(firk_mem->sunctx);
  }
  else if (strcmp(cname, "PI") == 0)
  {
    C = SUNAdaptController_PI(firk_mem->sunctx);
  }
  else if (strcmp(cname, "I") == 0)
  {
    C = SUNAdaptController_I(firk_mem->sunctx);
  }
  else if (strcmp(cname, "ExpGus") == 0)
  {
    C = SUNAdaptController_ExpGus(firk_mem->sunctx);
  }
  else if (strcmp(cname, "ImpGus") == 0)
  {
    C = SUNAdaptController_ImpGus(firk_mem->sunctx);
  }
  else if (strcmp(cname, "H0211") == 0)
  {
    C = SUNAdaptController_H0211(firk_mem->sunctx);
  }
  else if (strcmp(cname, "H0321") == 0)
  {
    C = SUNAdaptController_H0321(firk_mem->sunctx);
  }
  else if (strcmp(cname, "H211") == 0)
  {
    C = SUNAdaptController_H211(firk_mem->sunctx);
  }
  else if (strcmp(cname, "H312") == 0)
  {
    C = SUNAdaptController_H312(firk_mem->sunctx);
  }
  else
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_CONTROLLER_NAME);
    return FIRK_ILL_INPUT;
  }

  if (C == NULL)
  {
    firkProcessError(firk_mem, FIRK_MEM_FAIL, __LINE__, __func__, __FILE__,
                     "SUNAdaptController allocation failure");
    return FIRK_MEM_FAIL;
  }

  /* FIRKODE owns the new controller; it is destroyed if not installed */
  retval = firkReplaceAdaptController(firk_mem, C, SUNTRUE);
  if (retval != FIRK_SUCCESS) { (void)SUNAdaptController_Destroy(C); }
  return retval;
}

int FIRKodeSetSafetyFactor(void* firkode_mem, sunrealtype safety)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
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
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }

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
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
  firk_mem->hadapt_mem->growth = (mx_growth <= ONE) ? GROWTH : mx_growth;
  return FIRK_SUCCESS;
}

int FIRKodeSetMinReduction(void* firkode_mem, sunrealtype eta_min)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
  firk_mem->hadapt_mem->etamin = (eta_min >= ONE || eta_min <= ZERO) ? ETAMIN
                                                                     : eta_min;
  return FIRK_SUCCESS;
}

int FIRKodeSetFixedStepBounds(void* firkode_mem, sunrealtype lb, sunrealtype ub)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
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
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
  firk_mem->hadapt_mem->etamx1 = (etamx1 <= ONE) ? ETAMX1 : etamx1;
  return FIRK_SUCCESS;
}

int FIRKodeSetMaxEFailGrowth(void* firkode_mem, sunrealtype etamxf)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
  firk_mem->hadapt_mem->etamxf = ((etamxf <= ZERO) || (etamxf > ONE)) ? ETAMXF
                                                                      : etamxf;
  return FIRK_SUCCESS;
}

int FIRKodeSetSmallNumEFails(void* firkode_mem, int small_nef)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
  firk_mem->hadapt_mem->small_nef = (small_nef <= 0) ? SMALL_NEF : small_nef;
  return FIRK_SUCCESS;
}

int FIRKodeSetMaxCFailGrowth(void* firkode_mem, sunrealtype etacf)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
  firk_mem->hadapt_mem->etacf = ((etacf <= ZERO) || (etacf > ONE)) ? ETACF
                                                                   : etacf;
  return FIRK_SUCCESS;
}

int FIRKodeSetNewtonCountSafety(void* firkode_mem, sunbooleantype onoff)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
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
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
  firk_mem->maxcor = (maxcor <= 0) ? MAXCOR_DEFAULT : maxcor;
  return FIRK_SUCCESS;
}

int FIRKodeSetNonlinConvCoef(void* firkode_mem, sunrealtype nlscoef)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
  /* a non-positive value restores the default, RADAU5's rule */
  firk_mem->nlscoef_user = (nlscoef <= ZERO) ? NLSCOEF_AUTO : nlscoef;
  firkSetNlsCoef(firk_mem);
  return FIRK_SUCCESS;
}

int FIRKodeSetNonlinDivergenceRate(void* firkode_mem, sunrealtype theta_max)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
  /* the stop test eta dyno <= eps_nls, eta = theta / (1 - theta), needs
     theta < 1 for every rate that is not a divergence */
  if (theta_max >= ONE)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_THETA_MAX);
    return FIRK_ILL_INPUT;
  }
  firk_mem->theta_max = (theta_max <= ZERO) ? THETA_MAX_DEFAULT : theta_max;
  return FIRK_SUCCESS;
}

int FIRKodeSetJacBadConvRate(void* firkode_mem, sunrealtype theta_jbad)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
  firk_mem->theta_jbad = (theta_jbad <= ZERO) ? THETA_JBAD_DEFAULT : theta_jbad;
  return FIRK_SUCCESS;
}

int FIRKodeSetNonlinRateExponent(void* firkode_mem, sunrealtype rate_exp)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
  firk_mem->rate_exp = (rate_exp <= ZERO) ? RATE_EXP_DEFAULT : rate_exp;
  return FIRK_SUCCESS;
}

int FIRKodeSetJacReuseConvRate(void* firkode_mem, sunrealtype theta_jac)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  if (!firkExpertAllowed(firk_mem, __func__)) { return FIRK_ILL_INPUT; }
  if (theta_jac >= ONE)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_THETA_JAC);
    return FIRK_ILL_INPUT;
  }
  firk_mem->theta_jac = (theta_jac <= ZERO) ? THETA_JAC_DEFAULT : theta_jac;
  return FIRK_SUCCESS;
}

/*===============================================================
  Optional outputs
  ===============================================================*/

int FIRKodeGetNumSteps(void* firkode_mem, long int* nsteps_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *nsteps_ptr = firk_mem->nst;
  return FIRK_SUCCESS;
}

int FIRKodeGetNumStepAttempts(void* firkode_mem, long int* step_attempts_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *step_attempts_ptr = firk_mem->nst_attempts;
  return FIRK_SUCCESS;
}

int FIRKodeGetNumRhsEvals(void* firkode_mem, long int* nfevals_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *nfevals_ptr = firk_mem->nfe;
  return FIRK_SUCCESS;
}

int FIRKodeGetNumLinSolvSetups(void* firkode_mem, long int* nlinsetups_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *nlinsetups_ptr = firk_mem->nsetups;
  return FIRK_SUCCESS;
}

int FIRKodeGetNumErrTestFails(void* firkode_mem, long int* netfails_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *netfails_ptr = firk_mem->netf;
  return FIRK_SUCCESS;
}

int FIRKodeGetNumStepSolveFails(void* firkode_mem, long int* nncfails_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *nncfails_ptr = firk_mem->ncfn;
  return FIRK_SUCCESS;
}

int FIRKodeGetNumNonlinSolvIters(void* firkode_mem, long int* nniters_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *nniters_ptr = firk_mem->nni;
  return FIRK_SUCCESS;
}

int FIRKodeGetNumNonlinSolvConvFails(void* firkode_mem, long int* nnfails_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *nnfails_ptr = firk_mem->nnf;
  return FIRK_SUCCESS;
}

int FIRKodeGetNonlinSolvStats(void* firkode_mem, long int* nniters_ptr,
                              long int* nnfails_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *nniters_ptr = firk_mem->nni;
  *nnfails_ptr = firk_mem->nnf;
  return FIRK_SUCCESS;
}

int FIRKodeGetLastNonlinConvRate(void* firkode_mem, sunrealtype* theta_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *theta_ptr = firk_mem->theta;
  return FIRK_SUCCESS;
}

int FIRKodeGetActualInitStep(void* firkode_mem, sunrealtype* hinused_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *hinused_ptr = firk_mem->h0u;
  return FIRK_SUCCESS;
}

int FIRKodeGetLastStep(void* firkode_mem, sunrealtype* hlast_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *hlast_ptr = firk_mem->hold;
  return FIRK_SUCCESS;
}

int FIRKodeGetCurrentStep(void* firkode_mem, sunrealtype* hcur_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *hcur_ptr = firk_mem->next_h;
  return FIRK_SUCCESS;
}

int FIRKodeGetCurrentTime(void* firkode_mem, sunrealtype* tcur_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *tcur_ptr = firk_mem->tcur;
  return FIRK_SUCCESS;
}

int FIRKodeGetCurrentState(void* firkode_mem, N_Vector* y_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *y_ptr = firk_mem->yn;
  return FIRK_SUCCESS;
}

int FIRKodeGetCurrentGamma(void* firkode_mem, sunrealtype* gamma_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *gamma_ptr = firk_mem->gamma;
  return FIRK_SUCCESS;
}

int FIRKodeGetTolScaleFactor(void* firkode_mem, sunrealtype* tolsfac_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *tolsfac_ptr = firk_mem->tolsf;
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

int FIRKodeGetIntegratorStats(void* firkode_mem, long int* nsteps_ptr,
                              long int* nfevals_ptr, long int* nlinsetups_ptr,
                              long int* netfails_ptr, sunrealtype* hinused_ptr,
                              sunrealtype* hlast_ptr, sunrealtype* hcur_ptr,
                              sunrealtype* tcur_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *nsteps_ptr     = firk_mem->nst;
  *nfevals_ptr    = firk_mem->nfe;
  *nlinsetups_ptr = firk_mem->nsetups;
  *netfails_ptr   = firk_mem->netf;
  *hinused_ptr    = firk_mem->h0u;
  *hlast_ptr      = firk_mem->hold;
  *hcur_ptr       = firk_mem->next_h;
  *tcur_ptr       = firk_mem->tcur;
  return FIRK_SUCCESS;
}

int FIRKodeGetUserData(void* firkode_mem, void** user_data_ptr)
{
  FIRKodeMem firk_mem;
  FIRK_ACCESS_MEM(firkode_mem, firk_mem);
  *user_data_ptr = firk_mem->user_data;
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
  sunfprintf_long(outfile, fmt, SUNFALSE, "Error test fails", firk_mem->netf);
  sunfprintf_long(outfile, fmt, SUNFALSE, "NLS step fails", firk_mem->ncfn);
  sunfprintf_real(outfile, fmt, SUNFALSE, "Initial step size", firk_mem->h0u);
  sunfprintf_real(outfile, fmt, SUNFALSE, "Last step size", firk_mem->hold);
  sunfprintf_real(outfile, fmt, SUNFALSE, "Current step size", firk_mem->next_h);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Stages", firk_mem->s);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Method order", firk_mem->tbl->q);

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
    sunfprintf_long(outfile, fmt, SUNFALSE, "Block LS solves",
                    firkls_mem->nblocksolves);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Block LS iters", firkls_mem->nli);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Block LS fails", firkls_mem->ncfl);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Block prec fails",
                    firkls_mem->nblockfails);
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
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass-times evals",
                    firkls_massmem->nmmults);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass LS iters",
                    firkls_massmem->nmli);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass LS fails",
                    firkls_massmem->nmcfl);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass Prec setup evals",
                    firkls_massmem->nmpe);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass Prec solves",
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
  fprintf(fp, "  Mode = %s\n",
          (firk_mem->mode == FIRK_EXPERT) ? "EXPERT" : "BEGINNER");
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
    else { fprintf(fp, "  Initial step size = estimated\n"); }
    fprintf(fp, "  Minimum step size = " SUN_FORMAT_G "\n", firk_mem->hmin);
    if (firk_mem->hmax_inv != ZERO)
    {
      fprintf(fp, "  Maximum step size = " SUN_FORMAT_G "\n",
              ONE / firk_mem->hmax_inv);
    }
    else { fprintf(fp, "  Maximum step size = unbounded\n"); }
    fprintf(fp, "  Predictor method = %i\n", firk_mem->predictor);
    fprintf(fp, "  Refilter error estimate = %i\n", firk_mem->refilter);
    firkPrintAdaptMem(firk_mem->hadapt_mem, fp);
  }
  if (firk_mem->tstopset)
  {
    fprintf(fp, "  Stop time = " SUN_FORMAT_G ", interpolated = %i\n",
            firk_mem->tstop, (int)firk_mem->tstopinterp);
  }
  else { fprintf(fp, "  Stop time = none\n"); }
  fprintf(fp, "  Maximum number of steps = %li\n", firk_mem->mxstep);
  fprintf(fp, "  Maximum t + h = t warnings = %i\n", firk_mem->mxhnil);
  fprintf(fp, "  Maximum error test failures = %i\n", firk_mem->maxnef);
  fprintf(fp, "  Maximum convergence failures = %i\n", firk_mem->maxncf);
  fprintf(fp, "  Maximum nonlinear iterations = %i\n", firk_mem->maxcor);
  fprintf(fp, "  Nonlinear convergence coefficient = " SUN_FORMAT_G "%s\n",
          firk_mem->nlscoef,
          (firk_mem->nlscoef_user <= ZERO) ? " (automatic)" : "");
  fprintf(fp, "  Nonlinear divergence rate = " SUN_FORMAT_G "\n",
          firk_mem->theta_max);
  fprintf(fp, "  Nonlinear rate exponent = " SUN_FORMAT_G "\n",
          firk_mem->rate_exp);
  fprintf(fp, "  Jacobian reuse rate = " SUN_FORMAT_G "\n", firk_mem->theta_jac);
  fprintf(fp, "  Jacobian bad rate = " SUN_FORMAT_G "\n", firk_mem->theta_jbad);
  fprintf(fp, "  Reuse policy = %s\n",
          (firk_mem->reuse_policy == FIRK_POLICY_RADAU5) ? "RADAU5" : "PERIODIC");
  fprintf(fp, "  Linear solver setup frequency = %li\n", firk_mem->msbp);
  fprintf(fp, "  Linear solver setup threshold = " SUN_FORMAT_G "\n",
          firk_mem->dgmax_lsetup);
  if (firk_mem->root_mem != NULL && firk_mem->root_mem->nrtfn > 0)
  {
    fprintf(fp, "  Root functions = %i\n", firk_mem->root_mem->nrtfn);
  }
  firkLsWriteParameters(firk_mem, fp);
  fprintf(fp, "\n");

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  FIRKodeGetReturnFlagName

  Returns the name of the FIRKODE return flag.
  ---------------------------------------------------------------*/
char* FIRKodeGetReturnFlagName(long int flag)
{
  const char* str;
  char* name;

  switch (flag)
  {
  case FIRK_SUCCESS: str = "FIRK_SUCCESS"; break;
  case FIRK_TSTOP_RETURN: str = "FIRK_TSTOP_RETURN"; break;
  case FIRK_ROOT_RETURN: str = "FIRK_ROOT_RETURN"; break;
  case FIRK_WARNING: str = "FIRK_WARNING"; break;
  case FIRK_TOO_MUCH_WORK: str = "FIRK_TOO_MUCH_WORK"; break;
  case FIRK_TOO_MUCH_ACC: str = "FIRK_TOO_MUCH_ACC"; break;
  case FIRK_ERR_FAILURE: str = "FIRK_ERR_FAILURE"; break;
  case FIRK_CONV_FAILURE: str = "FIRK_CONV_FAILURE"; break;
  case FIRK_LINIT_FAIL: str = "FIRK_LINIT_FAIL"; break;
  case FIRK_LSETUP_FAIL: str = "FIRK_LSETUP_FAIL"; break;
  case FIRK_LSOLVE_FAIL: str = "FIRK_LSOLVE_FAIL"; break;
  case FIRK_RHSFUNC_FAIL: str = "FIRK_RHSFUNC_FAIL"; break;
  case FIRK_FIRST_RHSFUNC_ERR: str = "FIRK_FIRST_RHSFUNC_ERR"; break;
  case FIRK_REPTD_RHSFUNC_ERR: str = "FIRK_REPTD_RHSFUNC_ERR"; break;
  case FIRK_UNREC_RHSFUNC_ERR: str = "FIRK_UNREC_RHSFUNC_ERR"; break;
  case FIRK_RTFUNC_FAIL: str = "FIRK_RTFUNC_FAIL"; break;
  case FIRK_MEM_FAIL: str = "FIRK_MEM_FAIL"; break;
  case FIRK_MEM_NULL: str = "FIRK_MEM_NULL"; break;
  case FIRK_ILL_INPUT: str = "FIRK_ILL_INPUT"; break;
  case FIRK_NO_MALLOC: str = "FIRK_NO_MALLOC"; break;
  case FIRK_BAD_K: str = "FIRK_BAD_K"; break;
  case FIRK_BAD_T: str = "FIRK_BAD_T"; break;
  case FIRK_BAD_DKY: str = "FIRK_BAD_DKY"; break;
  case FIRK_TOO_CLOSE: str = "FIRK_TOO_CLOSE"; break;
  case FIRK_VECTOROP_ERR: str = "FIRK_VECTOROP_ERR"; break;
  case FIRK_MASSINIT_FAIL: str = "FIRK_MASSINIT_FAIL"; break;
  case FIRK_MASSSETUP_FAIL: str = "FIRK_MASSSETUP_FAIL"; break;
  case FIRK_MASSSOLVE_FAIL: str = "FIRK_MASSSOLVE_FAIL"; break;
  case FIRK_MASSMULT_FAIL: str = "FIRK_MASSMULT_FAIL"; break;
  case FIRK_TABLE_FAIL: str = "FIRK_TABLE_FAIL"; break;
  case FIRK_CONTROLLER_ERR: str = "FIRK_CONTROLLER_ERR"; break;
  case FIRK_UNRECOGNIZED_ERR: str = "FIRK_UNRECOGNIZED_ERR"; break;
  default: str = "NONE";
  }

  name = (char*)malloc(strlen(str) + 1);
  if (name == NULL) { return NULL; }
  memcpy(name, str, strlen(str) + 1);
  return name;
}
