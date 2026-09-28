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
 * This is the implementation file for the main FIRKODE integrator.
 *
 * FIRKODE integrates M y' = f(t,y) with the s-stage Radau IIA
 * collocation methods. Each step solves the coupled stage system
 *   M z_i = h sum_j a_ij f(t_n + c_j h, y_n + z_j),  i = 1..s,
 * for the stage increments z_i with a simplified Newton iteration.
 * The Newton systems are solved by an internal flexible GMRES
 * iteration on the stacked (s x N) vector of increments, using the
 * user's N x N linear solver on (M - h gamma0 J) as a block-diagonal
 * preconditioner; the same solver provides the RADAU5-style filtered
 * error estimate. The collocation polynomial of the last step gives
 * the stage predictor and dense output.
 * -----------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <firkode/firkode.h>
#include <nvector/nvector_manyvector.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_types.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#include "firkode_impl.h"
#include "firkode_root_impl.h"
#include "sundials_utils.h"

#define ZERO SUN_RCONST(0.0)
#define HALF SUN_RCONST(0.5)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)
#define FOUR SUN_RCONST(4.0)

/*===============================================================
  Private function prototypes
  ===============================================================*/

static sunbooleantype firkCheckNvector(N_Vector tmpl);
static sunbooleantype firkAllocVectors(FIRKodeMem firk_mem, N_Vector tmpl);
static void firkFreeVectors(FIRKodeMem firk_mem);
static sunbooleantype firkAllocStageVectors(FIRKodeMem firk_mem);
static void firkFreeStageVectors(FIRKodeMem firk_mem);

static int firkReInitCommon(FIRKodeMem firk_mem, sunrealtype t0, N_Vector y0,
                            sunbooleantype reset_counters);
static int firkInitialSetup(FIRKodeMem firk_mem, sunrealtype tout);

static int firkHin(FIRKodeMem firk_mem, sunrealtype tout);
static sunrealtype firkUpperBoundH0(FIRKodeMem firk_mem, sunrealtype tdist);
static int firkYddNorm(FIRKodeMem firk_mem, sunrealtype hg, sunrealtype* yddnrm);

static int firkStep(FIRKodeMem firk_mem);
static int firkPredict(FIRKodeMem firk_mem);
static int firkNls(FIRKodeMem firk_mem, int nflag);
static int firkHandleNFlag(FIRKodeMem firk_mem, int* nflagPtr, int* ncfPtr);
static int firkErrorEstimate(FIRKodeMem firk_mem, sunbooleantype after_reject,
                             sunrealtype* dsmPtr);
static int firkDoErrorTest(FIRKodeMem firk_mem, int* nflagPtr, int* nefPtr,
                           sunrealtype dsm);
static int firkCompleteStep(FIRKodeMem firk_mem, sunrealtype dsm);
static int firkHandleFailure(FIRKodeMem firk_mem, int flag);

static int firkEwtSetSS(FIRKodeMem firk_mem, N_Vector ycur, N_Vector weight);
static int firkEwtSetSV(FIRKodeMem firk_mem, N_Vector ycur, N_Vector weight);

/*===============================================================
  Exported functions
  ===============================================================*/

/*---------------------------------------------------------------
  FIRKodeCreate

  Creates an internal memory block for a problem to be solved by
  FIRKODE and sets all optional inputs to their default values.
  ---------------------------------------------------------------*/
void* FIRKodeCreate(SUNContext sunctx)
{
  FIRKodeMem firk_mem;
  int retval;

  if (sunctx == NULL)
  {
    firkProcessError(NULL, 0, __LINE__, __func__, __FILE__, MSGFIRK_NULL_SUNCTX);
    return NULL;
  }

  firk_mem = (FIRKodeMem)malloc(sizeof(struct FIRKodeMemRec));
  if (firk_mem == NULL)
  {
    firkProcessError(NULL, 0, __LINE__, __func__, __FILE__, MSGFIRK_MEM_FAIL);
    return NULL;
  }
  memset(firk_mem, 0, sizeof(struct FIRKodeMemRec));

  firk_mem->sunctx = sunctx;
  firk_mem->uround = SUN_UNIT_ROUNDOFF;

  /* problem specification defaults */
  firk_mem->f         = NULL;
  firk_mem->user_data = NULL;
  firk_mem->itol      = FIRK_NN;
  firk_mem->atolmin0  = SUNTRUE;
  firk_mem->user_efun = SUNFALSE;
  firk_mem->efun      = NULL;
  firk_mem->e_data    = NULL;

  /* method defaults */
  firk_mem->s         = FIRK_DEFAULT_STAGES;
  firk_mem->s_alloc   = 0;
  firk_mem->T         = NULL;
  firk_mem->predictor = FIRK_PREDICT_EXTRAPOLATE;
  firk_mem->refilter  = SUNTRUE;

  /* integrator defaults */
  firk_mem->mxstep       = MXSTEP_DEFAULT;
  firk_mem->mxhnil       = MXHNIL_DEFAULT;
  firk_mem->hin          = ZERO;
  firk_mem->hmin         = HMIN_DEFAULT;
  firk_mem->hmax_inv     = HMAX_INV_DEFAULT;
  firk_mem->fixedstep    = SUNFALSE;
  firk_mem->tstopset     = SUNFALSE;
  firk_mem->tstopinterp  = SUNFALSE;
  firk_mem->maxnef       = MXNEF;
  firk_mem->maxncf       = MXNCF;
  firk_mem->msbp         = MSBP_DEFAULT;
  firk_mem->dgmax_lsetup = DGMAX_LSETUP_DEFAULT;
  firk_mem->convfail     = FIRK_NO_FAILURES;

  /* nonlinear solver defaults */
  firk_mem->NLS        = NULL;
  firk_mem->ownNLS     = SUNFALSE;
  firk_mem->maxcor     = MAXCOR_DEFAULT;
  firk_mem->nlscoef    = NLSCOEF_DEFAULT;
  firk_mem->theta_max  = THETA_MAX_DEFAULT;
  firk_mem->theta_jbad = THETA_JBAD_DEFAULT;
  firk_mem->eta_nls    = ONE;
  firk_mem->theta      = ZERO;

  /* adaptivity defaults */
  firk_mem->hadapt_mem = firkAdaptInit();
  if (firk_mem->hadapt_mem == NULL)
  {
    firkProcessError(NULL, 0, __LINE__, __func__, __FILE__, MSGFIRK_MEM_FAIL);
    free(firk_mem);
    return NULL;
  }
  firk_mem->hadapt_mem->etamx1      = ETAMX1;
  firk_mem->hadapt_mem->etamxf      = ETAMXF;
  firk_mem->hadapt_mem->etamin      = ETAMIN;
  firk_mem->hadapt_mem->etacf       = ETACF;
  firk_mem->hadapt_mem->small_nef   = SMALL_NEF;
  firk_mem->hadapt_mem->safety      = SAFETY;
  firk_mem->hadapt_mem->growth      = GROWTH;
  firk_mem->hadapt_mem->lbound      = HFIXED_LB;
  firk_mem->hadapt_mem->ubound      = HFIXED_UB;
  firk_mem->hadapt_mem->newt_safety = SUNTRUE;
  firk_mem->hadapt_mem->etamax      = ETAMX1;
  retval = firkReplaceAdaptController(firk_mem, NULL, SUNTRUE);
  if (retval != FIRK_SUCCESS)
  {
    free(firk_mem->hadapt_mem);
    free(firk_mem);
    return NULL;
  }

  /* linear solver and mass matrix hooks */
  firk_mem->linit      = NULL;
  firk_mem->lreinit    = NULL;
  firk_mem->lsetup     = NULL;
  firk_mem->lsolve_stk = NULL;
  firk_mem->lsolve_blk = NULL;
  firk_mem->lfree      = NULL;
  firk_mem->lmem       = NULL;
  firk_mem->minit      = NULL;
  firk_mem->msetup     = NULL;
  firk_mem->mmult      = NULL;
  firk_mem->msolve     = NULL;
  firk_mem->mfree      = NULL;
  firk_mem->mass_mem   = NULL;
  firk_mem->mass_set   = SUNFALSE;

  /* rootfinding */
  firk_mem->root_mem = NULL;

  firk_mem->VabstolMallocDone = SUNFALSE;
  firk_mem->MallocDone        = SUNFALSE;

  return (void*)firk_mem;
}

/*---------------------------------------------------------------
  FIRKodeInit

  Allocates and initializes memory for a problem.
  ---------------------------------------------------------------*/
int FIRKodeInit(void* firkode_mem, FIRKRhsFn f, sunrealtype t0, N_Vector y0)
{
  FIRKodeMem firk_mem;
  int retval;

  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRK_MEM_NULL, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MEM);
    return FIRK_MEM_NULL;
  }
  firk_mem = (FIRKodeMem)firkode_mem;

  SUNDIALS_MARK_FUNCTION_BEGIN(FIRK_PROFILER);

  if (y0 == NULL)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_NULL_Y0);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_ILL_INPUT;
  }

  if (f == NULL)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_NULL_F);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_ILL_INPUT;
  }

  if (!firkCheckNvector(y0))
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_NVECTOR);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_ILL_INPUT;
  }

  /* allocate the N-length work vectors (using y0 as a template) */
  if (firk_mem->MallocDone) { firkFreeVectors(firk_mem); }
  if (!firkAllocVectors(firk_mem, y0))
  {
    firkProcessError(firk_mem, FIRK_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSGFIRK_MEM_FAIL);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_MEM_FAIL;
  }

  /* stage vectors and the nonlinear solver are allocated in the initial
     setup once the number of stages is final */
  firkFreeStageVectors(firk_mem);

  firk_mem->f = f;

  /* the linear solver hooks are checked in the initial setup */
  firk_mem->MallocDone = SUNTRUE;

  retval = firkReInitCommon(firk_mem, t0, y0, SUNTRUE);

  SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
  return retval;
}

/*---------------------------------------------------------------
  FIRKodeReInit

  Re-initializes FIRKODE for a new problem with the same sizes;
  all counters are reset.
  ---------------------------------------------------------------*/
int FIRKodeReInit(void* firkode_mem, sunrealtype t0, N_Vector y0)
{
  FIRKodeMem firk_mem;
  int retval;

  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRK_MEM_NULL, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MEM);
    return FIRK_MEM_NULL;
  }
  firk_mem = (FIRKodeMem)firkode_mem;

  SUNDIALS_MARK_FUNCTION_BEGIN(FIRK_PROFILER);

  if (firk_mem->MallocDone == SUNFALSE)
  {
    firkProcessError(firk_mem, FIRK_NO_MALLOC, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MALLOC);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_NO_MALLOC;
  }

  if (y0 == NULL)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_NULL_Y0);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_ILL_INPUT;
  }

  retval = firkReInitCommon(firk_mem, t0, y0, SUNTRUE);

  SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
  return retval;
}

/*---------------------------------------------------------------
  FIRKodeReset

  Resets the integrator state to (tR, yR) while keeping all
  counters and optional inputs.
  ---------------------------------------------------------------*/
int FIRKodeReset(void* firkode_mem, sunrealtype tR, N_Vector yR)
{
  FIRKodeMem firk_mem;
  int retval;

  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRK_MEM_NULL, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MEM);
    return FIRK_MEM_NULL;
  }
  firk_mem = (FIRKodeMem)firkode_mem;

  SUNDIALS_MARK_FUNCTION_BEGIN(FIRK_PROFILER);

  if (firk_mem->MallocDone == SUNFALSE)
  {
    firkProcessError(firk_mem, FIRK_NO_MALLOC, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MALLOC);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_NO_MALLOC;
  }

  if (yR == NULL)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_NULL_Y0);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_ILL_INPUT;
  }

  retval = firkReInitCommon(firk_mem, tR, yR, SUNFALSE);

  SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
  return retval;
}

/*---------------------------------------------------------------
  Tolerance specification functions
  ---------------------------------------------------------------*/

int FIRKodeSStolerances(void* firkode_mem, sunrealtype reltol, sunrealtype abstol)
{
  FIRKodeMem firk_mem;

  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRK_MEM_NULL, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MEM);
    return FIRK_MEM_NULL;
  }
  firk_mem = (FIRKodeMem)firkode_mem;

  if (firk_mem->MallocDone == SUNFALSE)
  {
    firkProcessError(firk_mem, FIRK_NO_MALLOC, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MALLOC);
    return FIRK_NO_MALLOC;
  }

  if (reltol < ZERO)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_RELTOL);
    return FIRK_ILL_INPUT;
  }

  if (abstol < ZERO)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_ABSTOL);
    return FIRK_ILL_INPUT;
  }

  firk_mem->reltol   = reltol;
  firk_mem->Sabstol  = abstol;
  firk_mem->atolmin0 = (abstol == ZERO);
  firk_mem->itol     = FIRK_SS;

  firk_mem->user_efun = SUNFALSE;
  firk_mem->efun      = firkEwtSet;
  firk_mem->e_data    = NULL; /* set to firk_mem in the initial setup */

  return FIRK_SUCCESS;
}

int FIRKodeSVtolerances(void* firkode_mem, sunrealtype reltol, N_Vector abstol)
{
  FIRKodeMem firk_mem;
  sunrealtype atolmin;

  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRK_MEM_NULL, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MEM);
    return FIRK_MEM_NULL;
  }
  firk_mem = (FIRKodeMem)firkode_mem;

  if (firk_mem->MallocDone == SUNFALSE)
  {
    firkProcessError(firk_mem, FIRK_NO_MALLOC, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MALLOC);
    return FIRK_NO_MALLOC;
  }

  if (reltol < ZERO)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_RELTOL);
    return FIRK_ILL_INPUT;
  }

  if (abstol == NULL)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_NULL_ABSTOL);
    return FIRK_ILL_INPUT;
  }

  if (abstol->ops->nvmin == NULL)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Missing N_VMin routine from N_Vector");
    return FIRK_ILL_INPUT;
  }
  atolmin = N_VMin(abstol);
  if (atolmin < ZERO)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_ABSTOL);
    return FIRK_ILL_INPUT;
  }

  if (!(firk_mem->VabstolMallocDone))
  {
    firk_mem->Vabstol = N_VClone(firk_mem->ewt);
    if (firk_mem->Vabstol == NULL)
    {
      firkProcessError(firk_mem, FIRK_MEM_FAIL, __LINE__, __func__, __FILE__,
                       MSGFIRK_MEM_FAIL);
      return FIRK_MEM_FAIL;
    }
    firk_mem->VabstolMallocDone = SUNTRUE;
  }

  firk_mem->reltol = reltol;
  N_VScale(ONE, abstol, firk_mem->Vabstol);
  firk_mem->atolmin0 = (atolmin == ZERO);
  firk_mem->itol     = FIRK_SV;

  firk_mem->user_efun = SUNFALSE;
  firk_mem->efun      = firkEwtSet;
  firk_mem->e_data    = NULL;

  return FIRK_SUCCESS;
}

int FIRKodeWFtolerances(void* firkode_mem, FIRKEwtFn efun)
{
  FIRKodeMem firk_mem;

  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRK_MEM_NULL, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MEM);
    return FIRK_MEM_NULL;
  }
  firk_mem = (FIRKodeMem)firkode_mem;

  if (firk_mem->MallocDone == SUNFALSE)
  {
    firkProcessError(firk_mem, FIRK_NO_MALLOC, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MALLOC);
    return FIRK_NO_MALLOC;
  }

  if (efun == NULL)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "efun = NULL illegal.");
    return FIRK_ILL_INPUT;
  }

  firk_mem->itol      = FIRK_WF;
  firk_mem->user_efun = SUNTRUE;
  firk_mem->efun      = efun;
  firk_mem->e_data    = NULL; /* set to user_data in the initial setup */

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  FIRKodeEvolve

  Main solver function. In FIRK_NORMAL mode the solver steps until
  it reaches or passes tout and then interpolates to obtain y(tout).
  In FIRK_ONE_STEP mode it takes one internal step and returns.
  ---------------------------------------------------------------*/
int FIRKodeEvolve(void* firkode_mem, sunrealtype tout, N_Vector yout,
                  sunrealtype* tret, int itask)
{
  FIRKodeMem firk_mem;
  long int nstloc;
  int retval, hflag, kflag, istate, ier, ewtsetOK, ir, irfndp;
  sunrealtype troundoff, tout_hin, rh, nrm;
  sunbooleantype inactive_roots;
  FIRKodeRootMem rootmem;

  /*
   * -------------------------------------
   * 1. Check and process inputs
   * -------------------------------------
   */

  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRK_MEM_NULL, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MEM);
    return FIRK_MEM_NULL;
  }
  firk_mem = (FIRKodeMem)firkode_mem;

  SUNDIALS_MARK_FUNCTION_BEGIN(FIRK_PROFILER);

  if (firk_mem->MallocDone == SUNFALSE)
  {
    firkProcessError(firk_mem, FIRK_NO_MALLOC, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MALLOC);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_NO_MALLOC;
  }

  if (yout == NULL)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_YOUT_NULL);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_ILL_INPUT;
  }

  if (tret == NULL)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_TRET_NULL);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_ILL_INPUT;
  }

  if ((itask != FIRK_NORMAL) && (itask != FIRK_ONE_STEP))
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_ITASK);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_ILL_INPUT;
  }

  /*
   * ----------------------------------------
   * 2. Initializations performed only at
   *    the first step after (Re)Init/Reset:
   *    - initial setup
   *    - evaluate f(t0, y0)
   *    - compute initial step size
   *    - check for approach to tstop
   * ----------------------------------------
   */

  if (firk_mem->initsetup)
  {
    firk_mem->tretlast = *tret = firk_mem->tcur;

    ier = firkInitialSetup(firk_mem, tout);
    if (ier != FIRK_SUCCESS)
    {
      SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
      return ier;
    }

    /* evaluate f at (t0, y0) */
    retval = firk_mem->f(firk_mem->tcur, firk_mem->yn, firk_mem->fn,
                         firk_mem->user_data);
    firk_mem->nfe++;
    if (retval < 0)
    {
      firkProcessError(firk_mem, FIRK_RHSFUNC_FAIL, __LINE__, __func__,
                       __FILE__, MSGFIRK_RHSFUNC_FAILED, firk_mem->tcur);
      SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
      return FIRK_RHSFUNC_FAIL;
    }
    if (retval > 0)
    {
      firkProcessError(firk_mem, FIRK_FIRST_RHSFUNC_ERR, __LINE__, __func__,
                       __FILE__, MSGFIRK_RHSFUNC_FIRST);
      SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
      return FIRK_FIRST_RHSFUNC_ERR;
    }
    firk_mem->fn_is_current  = SUNTRUE;
    firk_mem->ydn_is_current = !firk_mem->mass_set;

    /* test input tstop for legality */
    if (firk_mem->tstopset)
    {
      if ((firk_mem->tstop - firk_mem->tcur) * (tout - firk_mem->tcur) <= ZERO)
      {
        firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                         MSGFIRK_BAD_TSTOP, firk_mem->tstop, firk_mem->tcur);
        SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
        return FIRK_ILL_INPUT;
      }
    }

    /* set initial h (from hin or firkHin) */
    firk_mem->h = firk_mem->hin;
    if ((firk_mem->h != ZERO) && ((tout - firk_mem->tcur) * firk_mem->h < ZERO))
    {
      firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                       firk_mem->fixedstep ? MSGFIRK_BAD_HFIXED : MSGFIRK_BAD_H0);
      SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
      return FIRK_ILL_INPUT;
    }
    if (firk_mem->h == ZERO)
    {
      tout_hin = tout;
      if (firk_mem->tstopset &&
          (tout - firk_mem->tcur) * (tout - firk_mem->tstop) > ZERO)
      {
        tout_hin = firk_mem->tstop;
      }
      hflag = firkHin(firk_mem, tout_hin);
      if (hflag != FIRK_SUCCESS)
      {
        istate = firkHandleFailure(firk_mem, hflag);
        SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
        return istate;
      }
    }

    /* enforce hmax and hmin */
    rh = SUNRabs(firk_mem->h) * firk_mem->hmax_inv;
    if (rh > ONE) { firk_mem->h /= rh; }
    if (SUNRabs(firk_mem->h) < firk_mem->hmin)
    {
      firk_mem->h *= firk_mem->hmin / SUNRabs(firk_mem->h);
    }

    /* check for approach to tstop */
    if (firk_mem->tstopset)
    {
      if ((firk_mem->tcur + firk_mem->h - firk_mem->tstop) * firk_mem->h > ZERO)
      {
        firk_mem->h = (firk_mem->tstop - firk_mem->tcur) *
                      (ONE - FOUR * firk_mem->uround);
      }
    }

    firk_mem->h0u    = firk_mem->h;
    firk_mem->hprime = firk_mem->h;
    firk_mem->eta    = ONE;

    /* check for zeros of the root function g at and near t0 */
    if (firk_mem->root_mem != NULL && firk_mem->root_mem->nrtfn > 0)
    {
      retval = firkRootCheck1(firk_mem);
      if (retval == FIRK_RTFUNC_FAIL)
      {
        firkProcessError(firk_mem, FIRK_RTFUNC_FAIL, __LINE__, __func__,
                         __FILE__, MSGFIRK_RTFUNC_FAILED, firk_mem->tcur);
        SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
        return FIRK_RTFUNC_FAIL;
      }
      if (retval < 0)
      {
        SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
        return retval;
      }
    }
  }

  /*
   * ------------------------------------------------------
   * 3. At following steps, perform stop tests:
   *    - check if we passed tstop
   *    - check if we passed tout (NORMAL mode)
   *    - check if current time was returned (ONE_STEP mode)
   *    - check if we are close to tstop
   * -------------------------------------------------------
   */

  if (!firk_mem->initsetup)
  {
    troundoff = FUZZ_FACTOR * firk_mem->uround *
                (SUNRabs(firk_mem->tcur) + SUNRabs(firk_mem->h));

    /* check for a root in the last step taken, other than the last root
       found, if any; in ONE_STEP mode return y(tcur) if it was not returned
       because of an intervening root */
    if (firk_mem->root_mem != NULL && firk_mem->root_mem->nrtfn > 0)
    {
      rootmem = firk_mem->root_mem;
      irfndp  = rootmem->irfnd;

      retval = firkRootCheck2(firk_mem);
      if (retval == CLOSERT)
      {
        firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                         MSGFIRK_CLOSE_ROOTS, rootmem->tlo);
        SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
        return FIRK_ILL_INPUT;
      }
      else if (retval == FIRK_RTFUNC_FAIL)
      {
        firkProcessError(firk_mem, FIRK_RTFUNC_FAIL, __LINE__, __func__,
                         __FILE__, MSGFIRK_RTFUNC_FAILED, rootmem->tlo);
        SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
        return FIRK_RTFUNC_FAIL;
      }
      else if (retval == RTFOUND)
      {
        /* firkRootCheck2 leaves y(tlo) in ycur */
        firk_mem->tretlast = *tret = rootmem->tlo;
        N_VScale(ONE, firk_mem->ycur, yout);
        SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
        return FIRK_ROOT_RETURN;
      }

      /* if tcur is distinct from tretlast (within roundoff), check the
         remaining interval for roots */
      if (SUNRabs(firk_mem->tcur - firk_mem->tretlast) > troundoff)
      {
        retval = firkRootCheck3(firk_mem, tout, itask);
        if (retval == FIRK_SUCCESS)
        {
          rootmem->irfnd = 0;
          if ((irfndp == 1) && (itask == FIRK_ONE_STEP))
          {
            firk_mem->tretlast = *tret = firk_mem->tcur;
            N_VScale(ONE, firk_mem->yn, yout);
            SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
            return FIRK_SUCCESS;
          }
        }
        else if (retval == RTFOUND)
        {
          /* firkRootCheck3 leaves y(tlo) in ycur */
          rootmem->irfnd     = 1;
          firk_mem->tretlast = *tret = rootmem->tlo;
          N_VScale(ONE, firk_mem->ycur, yout);
          SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
          return FIRK_ROOT_RETURN;
        }
        else if (retval == FIRK_RTFUNC_FAIL)
        {
          firkProcessError(firk_mem, FIRK_RTFUNC_FAIL, __LINE__, __func__,
                           __FILE__, MSGFIRK_RTFUNC_FAILED, rootmem->tlo);
          SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
          return FIRK_RTFUNC_FAIL;
        }
      }
    }

    /* test for tcur at tstop or near tstop */
    if (firk_mem->tstopset)
    {
      if (SUNRabs(firk_mem->tcur - firk_mem->tstop) <= troundoff)
      {
        if ((tout - firk_mem->tstop) * firk_mem->h >= ZERO ||
            SUNRabs(tout - firk_mem->tstop) <= troundoff)
        {
          if (firk_mem->tstopinterp)
          {
            ier = FIRKodeGetDky(firk_mem, firk_mem->tstop, 0, yout);
            if (ier != FIRK_SUCCESS)
            {
              firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__,
                               __FILE__, MSGFIRK_BAD_TSTOP, firk_mem->tstop,
                               firk_mem->tcur);
              SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
              return FIRK_ILL_INPUT;
            }
          }
          else { N_VScale(ONE, firk_mem->yn, yout); }
          firk_mem->tretlast = *tret = firk_mem->tstop;
          firk_mem->tstopset         = SUNFALSE;
          SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
          return FIRK_TSTOP_RETURN;
        }
      }
      /* if the next step would overtake tstop, adjust the step size */
      else if ((firk_mem->tcur + firk_mem->hprime - firk_mem->tstop) *
                 firk_mem->h >
               ZERO)
      {
        firk_mem->hprime = (firk_mem->tstop - firk_mem->tcur) *
                           (ONE - FOUR * firk_mem->uround);
        firk_mem->eta = firk_mem->hprime / firk_mem->h;
      }
    }

    /* in FIRK_NORMAL mode, test if tout was reached */
    if ((itask == FIRK_NORMAL) && ((firk_mem->tcur - tout) * firk_mem->h >= ZERO))
    {
      firk_mem->tretlast = *tret = tout;
      ier                        = FIRKodeGetDky(firk_mem, tout, 0, yout);
      if (ier != FIRK_SUCCESS)
      {
        firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                         MSGFIRK_BAD_TOUT, tout);
        SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
        return FIRK_ILL_INPUT;
      }
      SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
      return FIRK_SUCCESS;
    }

    /* in FIRK_ONE_STEP mode, test if tcur was returned */
    if (itask == FIRK_ONE_STEP &&
        SUNRabs(firk_mem->tcur - firk_mem->tretlast) > troundoff)
    {
      firk_mem->tretlast = *tret = firk_mem->tcur;
      N_VScale(ONE, firk_mem->yn, yout);
      SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
      return FIRK_SUCCESS;
    }
  }

  /*
   * --------------------------------------------------
   * 4. Looping point for internal steps
   * --------------------------------------------------
   */

  nstloc = 0;
  istate = FIRK_SUCCESS;
  for (;;)
  {
    firk_mem->next_h = firk_mem->hprime;

    /* reset and check ewt (the initial setup loaded it for the first step) */
    if (!firk_mem->initsetup)
    {
      ewtsetOK = firk_mem->efun(firk_mem->yn, firk_mem->ewt, firk_mem->e_data);
      if (ewtsetOK != 0)
      {
        if (firk_mem->itol == FIRK_WF)
        {
          firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__,
                           __FILE__, MSGFIRK_EWT_NOW_FAIL, firk_mem->tcur);
        }
        else
        {
          firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__,
                           __FILE__, MSGFIRK_EWT_NOW_BAD, firk_mem->tcur);
        }
        istate             = FIRK_ILL_INPUT;
        firk_mem->tretlast = *tret = firk_mem->tcur;
        N_VScale(ONE, firk_mem->yn, yout);
        break;
      }
    }

    /* check for too many steps */
    if ((firk_mem->mxstep > 0) && (nstloc >= firk_mem->mxstep))
    {
      firkProcessError(firk_mem, FIRK_TOO_MUCH_WORK, __LINE__, __func__,
                       __FILE__, MSGFIRK_MAX_STEPS, firk_mem->tcur);
      istate             = FIRK_TOO_MUCH_WORK;
      firk_mem->tretlast = *tret = firk_mem->tcur;
      N_VScale(ONE, firk_mem->yn, yout);
      break;
    }

    /* check for too much accuracy requested */
    nrm             = N_VWrmsNorm(firk_mem->yn, firk_mem->ewt);
    firk_mem->tolsf = firk_mem->uround * nrm;
    if (firk_mem->tolsf > ONE)
    {
      firkProcessError(firk_mem, FIRK_TOO_MUCH_ACC, __LINE__, __func__,
                       __FILE__, MSGFIRK_TOO_MUCH_ACC, firk_mem->tcur);
      istate             = FIRK_TOO_MUCH_ACC;
      firk_mem->tretlast = *tret = firk_mem->tcur;
      N_VScale(ONE, firk_mem->yn, yout);
      firk_mem->tolsf *= TWO;
      break;
    }
    else { firk_mem->tolsf = ONE; }

    /* check for h below roundoff level in tcur */
    if (firk_mem->tcur + firk_mem->hprime == firk_mem->tcur)
    {
      firk_mem->nhnil++;
      if (firk_mem->nhnil <= firk_mem->mxhnil)
      {
        firkProcessError(firk_mem, FIRK_WARNING, __LINE__, __func__, __FILE__,
                         MSGFIRK_HNIL, firk_mem->tcur, firk_mem->hprime);
      }
      if (firk_mem->nhnil == firk_mem->mxhnil)
      {
        firkProcessError(firk_mem, FIRK_WARNING, __LINE__, __func__, __FILE__,
                         MSGFIRK_HNIL_DONE);
      }
    }

    /* make sure f(tn, yn) is available for this step */
    if (!firk_mem->fn_is_current)
    {
      retval = firk_mem->f(firk_mem->tn, firk_mem->yn, firk_mem->fn,
                           firk_mem->user_data);
      firk_mem->nfe++;
      if (retval != 0)
      {
        istate = (retval < 0) ? FIRK_RHSFUNC_FAIL : FIRK_UNREC_RHSFUNC_ERR;
        istate = firkHandleFailure(firk_mem, istate);
        firk_mem->tretlast = *tret = firk_mem->tcur;
        N_VScale(ONE, firk_mem->yn, yout);
        break;
      }
      firk_mem->fn_is_current  = SUNTRUE;
      firk_mem->ydn_is_current = !firk_mem->mass_set;
    }

    /* take a step */
    kflag = firkStep(firk_mem);

    /* process failed step cases, and exit loop */
    if (kflag != FIRK_SUCCESS)
    {
      istate             = firkHandleFailure(firk_mem, kflag);
      firk_mem->tretlast = *tret = firk_mem->tcur;
      N_VScale(ONE, firk_mem->yn, yout);
      break;
    }

    nstloc++;

    /* check for a root in the last step taken */
    if (firk_mem->root_mem != NULL && firk_mem->root_mem->nrtfn > 0)
    {
      rootmem = firk_mem->root_mem;
      retval  = firkRootCheck3(firk_mem, tout, itask);
      if (retval == RTFOUND)
      {
        /* firkRootCheck3 leaves y(tlo) in ycur */
        rootmem->irfnd     = 1;
        istate             = FIRK_ROOT_RETURN;
        firk_mem->tretlast = *tret = rootmem->tlo;
        N_VScale(ONE, firk_mem->ycur, yout);
        break;
      }
      else if (retval == FIRK_RTFUNC_FAIL)
      {
        firkProcessError(firk_mem, FIRK_RTFUNC_FAIL, __LINE__, __func__,
                         __FILE__, MSGFIRK_RTFUNC_FAILED, rootmem->tlo);
        istate = FIRK_RTFUNC_FAIL;
        break;
      }

      /* at the end of the first step, warn about root functions that are
         still inactive (identically zero) */
      if (firk_mem->nst == 1)
      {
        inactive_roots = SUNFALSE;
        for (ir = 0; ir < rootmem->nrtfn; ir++)
        {
          if (!rootmem->gactive[ir])
          {
            inactive_roots = SUNTRUE;
            break;
          }
        }
        if ((rootmem->mxgnull > 0) && inactive_roots)
        {
          firkProcessError(firk_mem, FIRK_WARNING, __LINE__, __func__, __FILE__,
                           MSGFIRK_INACTIVE_ROOTS);
        }
      }
    }

    /* check if tcur is at tstop or near tstop */
    if (firk_mem->tstopset)
    {
      troundoff = FUZZ_FACTOR * firk_mem->uround *
                  (SUNRabs(firk_mem->tcur) + SUNRabs(firk_mem->h));

      if (SUNRabs(firk_mem->tcur - firk_mem->tstop) <= troundoff)
      {
        if ((tout - firk_mem->tstop) * firk_mem->h >= ZERO ||
            SUNRabs(tout - firk_mem->tstop) <= troundoff)
        {
          if (firk_mem->tstopinterp)
          {
            (void)FIRKodeGetDky(firk_mem, firk_mem->tstop, 0, yout);
          }
          else { N_VScale(ONE, firk_mem->yn, yout); }
          firk_mem->tretlast = *tret = firk_mem->tstop;
          firk_mem->tstopset         = SUNFALSE;
          istate                     = FIRK_TSTOP_RETURN;
          break;
        }
      }
      /* if the next step would overtake tstop, adjust the step size */
      else if ((firk_mem->tcur + firk_mem->hprime - firk_mem->tstop) *
                 firk_mem->h >
               ZERO)
      {
        firk_mem->hprime = (firk_mem->tstop - firk_mem->tcur) *
                           (ONE - FOUR * firk_mem->uround);
        firk_mem->eta = firk_mem->hprime / firk_mem->h;
      }
    }

    /* in NORMAL mode, check if tout was reached */
    if ((itask == FIRK_NORMAL) && (firk_mem->tcur - tout) * firk_mem->h >= ZERO)
    {
      istate             = FIRK_SUCCESS;
      firk_mem->tretlast = *tret = tout;
      (void)FIRKodeGetDky(firk_mem, tout, 0, yout);
      firk_mem->next_h = firk_mem->hprime;
      break;
    }

    /* in ONE_STEP mode, copy y and exit loop */
    if (itask == FIRK_ONE_STEP)
    {
      istate             = FIRK_SUCCESS;
      firk_mem->tretlast = *tret = firk_mem->tcur;
      N_VScale(ONE, firk_mem->yn, yout);
      firk_mem->next_h = firk_mem->hprime;
      break;
    }
  }

  SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
  return istate;
}

/*---------------------------------------------------------------
  FIRKodeGetDky

  Computes the k-th derivative of the collocation polynomial of
  the last completed step at time t and stores it in dky:
    u(tcur + (theta - 1) hold) = yn + sum_j (L_j(theta) - d_j) Zprev_j,
  with L_j the Lagrange basis on the nodes {0, c_1, ..., c_s}.
  ---------------------------------------------------------------*/
int FIRKodeGetDky(void* firkode_mem, sunrealtype t, int k, N_Vector dky)
{
  FIRKodeMem firk_mem;
  FIRKodeTable T;
  sunrealtype tfuzz, tp, tn1, theta, tpow, fac;
  int i, j, m, s, ier;

  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRK_MEM_NULL, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MEM);
    return FIRK_MEM_NULL;
  }
  firk_mem = (FIRKodeMem)firkode_mem;

  SUNDIALS_MARK_FUNCTION_BEGIN(FIRK_PROFILER);

  if (dky == NULL)
  {
    firkProcessError(firk_mem, FIRK_BAD_DKY, __LINE__, __func__, __FILE__,
                     MSGFIRK_NULL_DKY);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_BAD_DKY;
  }

  if ((k < 0) || (k > firk_mem->s))
  {
    firkProcessError(firk_mem, FIRK_BAD_K, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_K);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_BAD_K;
  }

  /* without a completed step only the current state is available */
  if (!firk_mem->have_prev)
  {
    tfuzz = FUZZ_FACTOR * firk_mem->uround *
            (SUNRabs(firk_mem->tcur) + SUNRabs(firk_mem->h));
    if (k == 0 && SUNRabs(t - firk_mem->tcur) <= tfuzz)
    {
      N_VScale(ONE, firk_mem->yn, dky);
      SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
      return FIRK_SUCCESS;
    }
    firkProcessError(firk_mem, FIRK_BAD_T, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_T, t, firk_mem->tcur, firk_mem->tcur);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_BAD_T;
  }

  /* allow for some slack */
  tfuzz = FUZZ_FACTOR * firk_mem->uround *
          (SUNRabs(firk_mem->tcur) + SUNRabs(firk_mem->hold));
  if (firk_mem->hold < ZERO) { tfuzz = -tfuzz; }
  tp  = firk_mem->tcur - firk_mem->hold - tfuzz;
  tn1 = firk_mem->tcur + tfuzz;
  if ((t - tp) * (t - tn1) > ZERO)
  {
    firkProcessError(firk_mem, FIRK_BAD_T, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_T, t, firk_mem->tcur - firk_mem->hold,
                     firk_mem->tcur);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_BAD_T;
  }

  T     = firk_mem->T;
  s     = firk_mem->s;
  theta = ONE + (t - firk_mem->tcur) / firk_mem->hold;

  if (k == 0)
  {
    firk_mem->cvals[0] = ONE;
    firk_mem->Xvecs[0] = firk_mem->yn;
    for (j = 0; j < s; j++)
    {
      firk_mem->cvals[j + 1] = -T->d[j];
      firk_mem->Xvecs[j + 1] = firk_mem->Zprev[j];
    }
    tpow = ONE;
    for (m = 0; m < s; m++)
    {
      tpow *= theta;
      for (j = 0; j < s; j++) { firk_mem->cvals[j + 1] += tpow * T->P[m][j]; }
    }
    ier = N_VLinearCombination(s + 1, firk_mem->cvals, firk_mem->Xvecs, dky);
    if (ier != 0)
    {
      SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
      return FIRK_VECTOROP_ERR;
    }
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_SUCCESS;
  }

  /* k-th derivative: d^k/dtheta^k theta^(m+1) = (m+1)!/(m+1-k)! theta^(m+1-k) */
  for (j = 0; j < s; j++)
  {
    firk_mem->cvals[j] = ZERO;
    firk_mem->Xvecs[j] = firk_mem->Zprev[j];
  }
  for (m = k - 1; m < s; m++)
  {
    fac = ONE;
    for (i = m + 1; i > m + 1 - k; i--) { fac *= (sunrealtype)i; }
    tpow = SUNRpowerI(theta, m + 1 - k);
    for (j = 0; j < s; j++) { firk_mem->cvals[j] += fac * tpow * T->P[m][j]; }
  }
  ier = N_VLinearCombination(s, firk_mem->cvals, firk_mem->Xvecs, dky);
  if (ier != 0)
  {
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_VECTOROP_ERR;
  }
  N_VScale(SUNRpowerI(firk_mem->hold, -k), dky, dky);

  SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  FIRKodeFree

  Frees all memory allocated by FIRKodeCreate and FIRKodeInit.
  ---------------------------------------------------------------*/
void FIRKodeFree(void** firkode_mem)
{
  FIRKodeMem firk_mem;

  if (firkode_mem == NULL || *firkode_mem == NULL) { return; }
  firk_mem = (FIRKodeMem)(*firkode_mem);

  firkFreeStageVectors(firk_mem);
  if (firk_mem->MallocDone) { firkFreeVectors(firk_mem); }

  if (firk_mem->ownNLS && firk_mem->NLS != NULL)
  {
    SUNNonlinSolFree(firk_mem->NLS);
    firk_mem->NLS    = NULL;
    firk_mem->ownNLS = SUNFALSE;
  }

  if (firk_mem->lfree != NULL) { firk_mem->lfree(firk_mem); }
  if (firk_mem->mfree != NULL) { firk_mem->mfree(firk_mem); }

  (void)firkRootFree(firk_mem);

  if (firk_mem->T != NULL)
  {
    FIRKodeTable_Free(firk_mem->T);
    firk_mem->T = NULL;
  }

  if (firk_mem->hadapt_mem != NULL)
  {
    if (firk_mem->hadapt_mem->owncontroller &&
        firk_mem->hadapt_mem->hcontroller != NULL)
    {
      (void)SUNAdaptController_Destroy(firk_mem->hadapt_mem->hcontroller);
    }
    free(firk_mem->hadapt_mem);
    firk_mem->hadapt_mem = NULL;
  }

  free(*firkode_mem);
  *firkode_mem = NULL;
}

/*===============================================================
  Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  firkCheckNvector

  Checks that all required vector operations are present. The
  stacked stage vector (NVECTOR_MANYVECTOR) additionally requires
  N_VGetLength, and the Krylov iteration on it N_VDotProd.
  ---------------------------------------------------------------*/
static sunbooleantype firkCheckNvector(N_Vector tmpl)
{
  if ((tmpl->ops->nvclone == NULL) || (tmpl->ops->nvdestroy == NULL) ||
      (tmpl->ops->nvlinearsum == NULL) || (tmpl->ops->nvconst == NULL) ||
      (tmpl->ops->nvprod == NULL) || (tmpl->ops->nvdiv == NULL) ||
      (tmpl->ops->nvscale == NULL) || (tmpl->ops->nvabs == NULL) ||
      (tmpl->ops->nvinv == NULL) || (tmpl->ops->nvaddconst == NULL) ||
      (tmpl->ops->nvmaxnorm == NULL) || (tmpl->ops->nvwrmsnorm == NULL) ||
      (tmpl->ops->nvdotprod == NULL) || (tmpl->ops->nvgetlength == NULL))
  {
    return SUNFALSE;
  }
  return SUNTRUE;
}

/*---------------------------------------------------------------
  firkAllocVectors / firkFreeVectors

  Allocate/free the N-length work vectors of the integrator.
  ---------------------------------------------------------------*/
static sunbooleantype firkAllocVectors(FIRKodeMem firk_mem, N_Vector tmpl)
{
  N_Vector* vecs[9];
  int i, j;

  vecs[0] = &firk_mem->yn;
  vecs[1] = &firk_mem->ycur;
  vecs[2] = &firk_mem->fn;
  vecs[3] = &firk_mem->ewt;
  vecs[4] = &firk_mem->ele;
  vecs[5] = &firk_mem->ftemp;
  vecs[6] = &firk_mem->tempv1;
  vecs[7] = &firk_mem->tempv2;
  vecs[8] = &firk_mem->tempv3;

  for (i = 0; i < 9; i++)
  {
    *vecs[i] = N_VClone(tmpl);
    if (*vecs[i] == NULL)
    {
      for (j = 0; j < i; j++)
      {
        N_VDestroy(*vecs[j]);
        *vecs[j] = NULL;
      }
      return SUNFALSE;
    }
  }

  /* without a mass matrix ydn aliases fn; a separate vector is allocated
     in the initial setup when a mass matrix is attached */
  firk_mem->ydn = firk_mem->fn;

  return SUNTRUE;
}

static void firkFreeVectors(FIRKodeMem firk_mem)
{
  if (firk_mem->ydn != NULL && firk_mem->ydn != firk_mem->fn)
  {
    N_VDestroy(firk_mem->ydn);
  }
  firk_mem->ydn = NULL;

  if (firk_mem->yn) { N_VDestroy(firk_mem->yn); }
  if (firk_mem->ycur) { N_VDestroy(firk_mem->ycur); }
  if (firk_mem->fn) { N_VDestroy(firk_mem->fn); }
  if (firk_mem->ewt) { N_VDestroy(firk_mem->ewt); }
  if (firk_mem->ele) { N_VDestroy(firk_mem->ele); }
  if (firk_mem->ftemp) { N_VDestroy(firk_mem->ftemp); }
  if (firk_mem->tempv1) { N_VDestroy(firk_mem->tempv1); }
  if (firk_mem->tempv2) { N_VDestroy(firk_mem->tempv2); }
  if (firk_mem->tempv3) { N_VDestroy(firk_mem->tempv3); }
  firk_mem->yn = firk_mem->ycur = firk_mem->fn = firk_mem->ewt = NULL;
  firk_mem->ele = firk_mem->ftemp = firk_mem->tempv1 = NULL;
  firk_mem->tempv2 = firk_mem->tempv3 = NULL;

  if (firk_mem->VabstolMallocDone)
  {
    N_VDestroy(firk_mem->Vabstol);
    firk_mem->Vabstol           = NULL;
    firk_mem->VabstolMallocDone = SUNFALSE;
  }
}

/*---------------------------------------------------------------
  firkAllocStageVectors / firkFreeStageVectors

  Allocate/free the stage vector arrays, the stacked ManyVector
  wrappers, and the fused-operation scratch arrays for s stages.
  ---------------------------------------------------------------*/
static sunbooleantype firkAllocStageVectors(FIRKodeMem firk_mem)
{
  int i, s;
  N_Vector* ewt_arr;

  s = firk_mem->s;

  firk_mem->Zpred = N_VCloneVectorArray(s, firk_mem->yn);
  firk_mem->Zcor  = N_VCloneVectorArray(s, firk_mem->yn);
  firk_mem->Z     = N_VCloneVectorArray(s, firk_mem->yn);
  firk_mem->Zprev = N_VCloneVectorArray(s, firk_mem->yn);
  firk_mem->F     = N_VCloneVectorArray(s, firk_mem->yn);
  firk_mem->cvals = (sunrealtype*)calloc((size_t)(s + 2), sizeof(sunrealtype));
  firk_mem->Xvecs = (N_Vector*)calloc((size_t)(s + 2), sizeof(N_Vector));
  firk_mem->wtmp  = (sunrealtype*)calloc((size_t)(s + 2), sizeof(sunrealtype));
  if (firk_mem->Zpred == NULL || firk_mem->Zcor == NULL || firk_mem->Z == NULL ||
      firk_mem->Zprev == NULL || firk_mem->F == NULL || firk_mem->cvals == NULL ||
      firk_mem->Xvecs == NULL || firk_mem->wtmp == NULL)
  {
    firkFreeStageVectors(firk_mem);
    return SUNFALSE;
  }

  firk_mem->zpred_stk = N_VNew_ManyVector((sunindextype)s, firk_mem->Zpred,
                                          firk_mem->sunctx);
  firk_mem->zcor_stk  = N_VNew_ManyVector((sunindextype)s, firk_mem->Zcor,
                                          firk_mem->sunctx);

  ewt_arr = (N_Vector*)malloc((size_t)s * sizeof(N_Vector));
  if (ewt_arr == NULL)
  {
    firkFreeStageVectors(firk_mem);
    return SUNFALSE;
  }
  for (i = 0; i < s; i++) { ewt_arr[i] = firk_mem->ewt; }
  firk_mem->ewt_stk = N_VNew_ManyVector((sunindextype)s, ewt_arr,
                                        firk_mem->sunctx);
  free(ewt_arr);

  if (firk_mem->zpred_stk == NULL || firk_mem->zcor_stk == NULL ||
      firk_mem->ewt_stk == NULL)
  {
    firkFreeStageVectors(firk_mem);
    return SUNFALSE;
  }

  firk_mem->s_alloc = s;

  return SUNTRUE;
}

static void firkFreeStageVectors(FIRKodeMem firk_mem)
{
  int s = firk_mem->s_alloc;

  if (firk_mem->zpred_stk) { N_VDestroy(firk_mem->zpred_stk); }
  if (firk_mem->zcor_stk) { N_VDestroy(firk_mem->zcor_stk); }
  if (firk_mem->ewt_stk) { N_VDestroy(firk_mem->ewt_stk); }
  firk_mem->zpred_stk = firk_mem->zcor_stk = firk_mem->ewt_stk = NULL;

  if (s > 0)
  {
    if (firk_mem->Zpred) { N_VDestroyVectorArray(firk_mem->Zpred, s); }
    if (firk_mem->Zcor) { N_VDestroyVectorArray(firk_mem->Zcor, s); }
    if (firk_mem->Z) { N_VDestroyVectorArray(firk_mem->Z, s); }
    if (firk_mem->Zprev) { N_VDestroyVectorArray(firk_mem->Zprev, s); }
    if (firk_mem->F) { N_VDestroyVectorArray(firk_mem->F, s); }
  }
  firk_mem->Zpred = firk_mem->Zcor = firk_mem->Z = NULL;
  firk_mem->Zprev = firk_mem->F = NULL;

  free(firk_mem->cvals);
  free(firk_mem->Xvecs);
  free(firk_mem->wtmp);
  firk_mem->cvals = NULL;
  firk_mem->Xvecs = NULL;
  firk_mem->wtmp  = NULL;

  firk_mem->s_alloc   = 0;
  firk_mem->have_prev = SUNFALSE;
}

/*---------------------------------------------------------------
  firkReInitCommon

  Shared body of FIRKodeInit, FIRKodeReInit and FIRKodeReset.
  ---------------------------------------------------------------*/
static int firkReInitCommon(FIRKodeMem firk_mem, sunrealtype t0, N_Vector y0,
                            sunbooleantype reset_counters)
{
  N_VScale(ONE, y0, firk_mem->yn);
  firk_mem->tn   = t0;
  firk_mem->tcur = t0;

  firk_mem->initsetup      = SUNTRUE;
  firk_mem->firststage     = SUNTRUE;
  firk_mem->have_prev      = SUNFALSE;
  firk_mem->fn_is_current  = SUNFALSE;
  firk_mem->ydn_is_current = SUNFALSE;

  firk_mem->h      = ZERO;
  firk_mem->hprime = ZERO;
  firk_mem->next_h = ZERO;
  firk_mem->eta    = ONE;
  firk_mem->hold   = ZERO;
  firk_mem->h0u    = ZERO;
  firk_mem->gamma  = ZERO;
  firk_mem->gammap = ZERO;
  firk_mem->gamrat = ONE;
  firk_mem->dsm    = ZERO;
  firk_mem->tolsf  = ONE;

  firk_mem->eta_nls    = ONE;
  firk_mem->theta      = ZERO;
  firk_mem->delp       = ZERO;
  firk_mem->delnrm     = ZERO;
  firk_mem->nls_reason = FIRK_NLS_REASON_NONE;
  firk_mem->nni_last   = 0;
  firk_mem->convfail   = FIRK_NO_FAILURES;
  firk_mem->jcur       = SUNFALSE;
  firk_mem->nstlp      = 0;

  if (reset_counters)
  {
    firk_mem->nst          = 0;
    firk_mem->nst_attempts = 0;
    firk_mem->nfe          = 0;
    firk_mem->ncfn         = 0;
    firk_mem->nni          = 0;
    firk_mem->nnf          = 0;
    firk_mem->netf         = 0;
    firk_mem->nsetups      = 0;
    firk_mem->nhnil        = 0;
    if (firk_mem->hadapt_mem) { firk_mem->hadapt_mem->nst_acc = 0; }
    if (firk_mem->lreinit) { firk_mem->lreinit(firk_mem); }
    if (firk_mem->root_mem) { firk_mem->root_mem->nge = 0; }
  }
  if (firk_mem->root_mem) { firk_mem->root_mem->irfnd = 0; }

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkInitialSetup

  Performs input consistency checks and completes the allocation
  of stage-dependent memory before the first step.
  ---------------------------------------------------------------*/
static int firkInitialSetup(FIRKodeMem firk_mem, sunrealtype tout)
{
  int ier;
  sunrealtype tdist, tround;
  SUNNonlinearSolver NLS;

  /* is tout too close to tn? */
  tdist  = SUNRabs(tout - firk_mem->tcur);
  tround = firk_mem->uround * SUNMAX(SUNRabs(firk_mem->tcur), SUNRabs(tout));
  if (tdist == ZERO || tdist < TWO * tround)
  {
    firkProcessError(firk_mem, FIRK_TOO_CLOSE, __LINE__, __func__, __FILE__,
                     MSGFIRK_TOO_CLOSE);
    return FIRK_TOO_CLOSE;
  }

  /* did the user specify tolerances? */
  if (firk_mem->itol == FIRK_NN)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_TOL);
    return FIRK_ILL_INPUT;
  }

  /* built-in error weights with abstol == 0 need N_VMin */
  if ((!firk_mem->user_efun) && (firk_mem->atolmin0) &&
      (!firk_mem->tempv1->ops->nvmin))
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Missing N_VMin routine from N_Vector");
    return FIRK_ILL_INPUT;
  }

  /* set data for efun */
  if (firk_mem->user_efun) { firk_mem->e_data = firk_mem->user_data; }
  else { firk_mem->e_data = firk_mem; }

  /* the fixed step size must point towards tout */
  if (firk_mem->fixedstep && (firk_mem->hin * (tout - firk_mem->tcur) < ZERO))
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_BAD_HFIXED);
    return FIRK_ILL_INPUT;
  }

  /* construct the table for the requested number of stages */
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
  firk_mem->hadapt_mem->p = firk_mem->T->p;
  firk_mem->hadapt_mem->q = firk_mem->T->q;

  /* allocate the stage memory and the Newton solver for s stages */
  if (firk_mem->s_alloc != firk_mem->s)
  {
    firkFreeStageVectors(firk_mem);
    if (firk_mem->ownNLS && firk_mem->NLS != NULL)
    {
      SUNNonlinSolFree(firk_mem->NLS);
      firk_mem->NLS    = NULL;
      firk_mem->ownNLS = SUNFALSE;
    }

    if (!firkAllocStageVectors(firk_mem))
    {
      firkProcessError(firk_mem, FIRK_MEM_FAIL, __LINE__, __func__, __FILE__,
                       MSGFIRK_MEM_FAIL);
      return FIRK_MEM_FAIL;
    }

    NLS = SUNNonlinSol_Newton(firk_mem->zpred_stk, firk_mem->sunctx);
    if (NLS == NULL)
    {
      firkProcessError(firk_mem, FIRK_MEM_FAIL, __LINE__, __func__, __FILE__,
                       MSGFIRK_MEM_FAIL);
      return FIRK_MEM_FAIL;
    }
    ier = firkNlsAttach(firk_mem, NLS);
    if (ier != FIRK_SUCCESS)
    {
      SUNNonlinSolFree(NLS);
      return ier;
    }
    firk_mem->ownNLS = SUNTRUE;
  }
  firk_mem->have_prev = SUNFALSE;

  /* load initial error weights */
  ier = firk_mem->efun(firk_mem->yn, firk_mem->ewt, firk_mem->e_data);
  if (ier != 0)
  {
    if (firk_mem->itol == FIRK_WF)
    {
      firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                       MSGFIRK_EWT_FAIL);
    }
    else
    {
      firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                       MSGFIRK_BAD_EWT);
    }
    return FIRK_ILL_INPUT;
  }

  /* a linear solver is required */
  if (firk_mem->lmem == NULL || firk_mem->linit == NULL)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_LS);
    return FIRK_ILL_INPUT;
  }
  ier = firk_mem->linit(firk_mem);
  if (ier != 0)
  {
    firkProcessError(firk_mem, FIRK_LINIT_FAIL, __LINE__, __func__, __FILE__,
                     MSGFIRK_LINIT_FAIL);
    return FIRK_LINIT_FAIL;
  }

  /* mass matrix initialization and (one-time) setup */
  if (firk_mem->mass_set)
  {
    if (firk_mem->ydn == firk_mem->fn)
    {
      firk_mem->ydn = N_VClone(firk_mem->fn);
      if (firk_mem->ydn == NULL)
      {
        firk_mem->ydn = firk_mem->fn;
        firkProcessError(firk_mem, FIRK_MEM_FAIL, __LINE__, __func__, __FILE__,
                         MSGFIRK_MEM_FAIL);
        return FIRK_MEM_FAIL;
      }
    }
    if (firk_mem->minit != NULL)
    {
      ier = firk_mem->minit(firk_mem);
      if (ier != 0)
      {
        firkProcessError(firk_mem, FIRK_MASSINIT_FAIL, __LINE__, __func__,
                         __FILE__, MSGFIRK_MASSINIT_FAIL);
        return FIRK_MASSINIT_FAIL;
      }
    }
    if (firk_mem->msetup != NULL)
    {
      ier = firk_mem->msetup(firk_mem, firk_mem->tcur, firk_mem->tempv1,
                             firk_mem->tempv2, firk_mem->tempv3);
      if (ier != 0)
      {
        firkProcessError(firk_mem, FIRK_MASSSETUP_FAIL, __LINE__, __func__,
                         __FILE__, MSGFIRK_MASSSETUP_FAIL);
        return FIRK_MASSSETUP_FAIL;
      }
    }
  }

  /* initialize the nonlinear solver (after the linear solver so that the
     lsetup and lsolve hooks are set) */
  ier = firkNlsInit(firk_mem);
  if (ier != 0)
  {
    firkProcessError(firk_mem, FIRK_NLS_INIT_FAIL, __LINE__, __func__, __FILE__,
                     MSGFIRK_NLS_INIT_FAIL);
    return FIRK_NLS_INIT_FAIL;
  }

  /* reset the step size controller and step-size heuristics */
  if (firk_mem->hadapt_mem->hcontroller != NULL)
  {
    ier = SUNAdaptController_Reset(firk_mem->hadapt_mem->hcontroller);
    if (ier != SUN_SUCCESS)
    {
      firkProcessError(firk_mem, FIRK_CONTROLLER_ERR, __LINE__, __func__,
                       __FILE__, "SUNAdaptController_Reset failure");
      return FIRK_CONTROLLER_ERR;
    }
  }
  firk_mem->hadapt_mem->etamax = firk_mem->hadapt_mem->etamx1;
  firk_mem->eta_nls            = ONE;
  firk_mem->gammap             = ZERO;
  firk_mem->nstlp              = 0;

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkGetYdot

  Makes ydn = M^{-1} f(tn, yn) current (ydn aliases fn without a
  mass matrix). Returns FIRK_SUCCESS or FIRK_MASSSOLVE_FAIL.
  ---------------------------------------------------------------*/
int firkGetYdot(FIRKodeMem firk_mem)
{
  int retval;

  if (!firk_mem->fn_is_current)
  {
    retval = firk_mem->f(firk_mem->tcur, firk_mem->yn, firk_mem->fn,
                         firk_mem->user_data);
    firk_mem->nfe++;
    if (retval < 0) { return FIRK_RHSFUNC_FAIL; }
    if (retval > 0) { return FIRK_UNREC_RHSFUNC_ERR; }
    firk_mem->fn_is_current  = SUNTRUE;
    firk_mem->ydn_is_current = SUNFALSE;
  }

  if (!firk_mem->mass_set)
  {
    firk_mem->ydn_is_current = SUNTRUE;
    return FIRK_SUCCESS;
  }
  if (firk_mem->ydn_is_current) { return FIRK_SUCCESS; }

  N_VScale(ONE, firk_mem->fn, firk_mem->ydn);
  retval = firk_mem->msolve(firk_mem, firk_mem->ydn, firk_mem->nlscoef);
  if (retval != 0) { return FIRK_MASSSOLVE_FAIL; }
  firk_mem->ydn_is_current = SUNTRUE;
  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkHin

  Computes a tentative initial step size h0 as the solution of
    (WRMS norm of (h0^2 ydd / 2)) = 1,
  where ydd is an estimated second derivative of y (same algorithm
  as CVODE's cvHin).
  ---------------------------------------------------------------*/
static int firkHin(FIRKodeMem firk_mem, sunrealtype tout)
{
  int retval, sign, count1, count2;
  sunrealtype tdiff, tdist, tround, hlb, hub;
  sunrealtype hg, hgs, hs, hnew, hrat, h0, yddnrm;
  sunbooleantype hgOK;

  retval = firkGetYdot(firk_mem);
  if (retval != FIRK_SUCCESS) { return retval; }

  tdiff  = tout - firk_mem->tcur;
  sign   = (tdiff > ZERO) ? 1 : -1;
  tdist  = SUNRabs(tdiff);
  tround = firk_mem->uround * SUNMAX(SUNRabs(firk_mem->tcur), SUNRabs(tout));

  hlb = HLB_FACTOR * tround;
  hub = firkUpperBoundH0(firk_mem, tdist);

  hg = SUNRsqrt(hlb * hub);

  if (hub < hlb)
  {
    firk_mem->h = (sign == -1) ? -hg : hg;
    return FIRK_SUCCESS;
  }

  hs   = hg;
  hnew = hg;

  for (count1 = 1; count1 <= MAX_ITERS; count1++)
  {
    hgOK = SUNFALSE;

    for (count2 = 1; count2 <= MAX_ITERS; count2++)
    {
      hgs    = hg * sign;
      retval = firkYddNorm(firk_mem, hgs, &yddnrm);
      if (retval < 0) { return retval; }
      if (retval == FIRK_SUCCESS)
      {
        hgOK = SUNTRUE;
        break;
      }
      hg *= SUN_RCONST(0.2);
    }

    if (!hgOK)
    {
      if (count1 <= 2) { return FIRK_REPTD_RHSFUNC_ERR; }
      hnew = hs;
      break;
    }

    hs = hg;

    hnew = (yddnrm * hub * hub > TWO) ? SUNRsqrt(TWO / yddnrm)
                                      : SUNRsqrt(hg * hub);

    if (count1 == MAX_ITERS) { break; }

    hrat = hnew / hg;

    if ((hrat > HALF) && (hrat < TWO)) { break; }

    if ((count1 > 1) && (hrat > TWO))
    {
      hnew = hg;
      break;
    }

    hg = hnew;
  }

  h0 = H_BIAS * hnew;
  if (h0 < hlb) { h0 = hlb; }
  if (h0 > hub) { h0 = hub; }
  if (sign == -1) { h0 = -h0; }
  firk_mem->h = h0;

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkUpperBoundH0

  Sets an upper bound on |h0| based on tdist and |y0|/|y0'|.
  ---------------------------------------------------------------*/
static sunrealtype firkUpperBoundH0(FIRKodeMem firk_mem, sunrealtype tdist)
{
  sunrealtype hub_inv, hub;
  N_Vector temp1, temp2;

  temp1 = firk_mem->tempv1;
  temp2 = firk_mem->tempv2;

  N_VAbs(firk_mem->yn, temp2);
  firk_mem->efun(firk_mem->yn, temp1, firk_mem->e_data);
  N_VInv(temp1, temp1);
  N_VLinearSum(HUB_FACTOR, temp2, ONE, temp1, temp1);

  N_VAbs(firk_mem->ydn, temp2);

  N_VDiv(temp2, temp1, temp1);
  hub_inv = N_VMaxNorm(temp1);

  hub = HUB_FACTOR * tdist;

  if (hub * hub_inv > ONE) { hub = ONE / hub_inv; }

  return hub;
}

/*---------------------------------------------------------------
  firkYddNorm

  Computes a difference-quotient estimate of the second derivative
  of y and returns its WRMS norm.
  ---------------------------------------------------------------*/
static int firkYddNorm(FIRKodeMem firk_mem, sunrealtype hg, sunrealtype* yddnrm)
{
  int retval;

  N_VLinearSum(hg, firk_mem->ydn, ONE, firk_mem->yn, firk_mem->ycur);
  retval = firk_mem->f(firk_mem->tcur + hg, firk_mem->ycur, firk_mem->tempv3,
                       firk_mem->user_data);
  firk_mem->nfe++;
  if (retval < 0) { return FIRK_RHSFUNC_FAIL; }
  if (retval > 0) { return RHSFUNC_RECVR; }

  if (firk_mem->mass_set)
  {
    retval = firk_mem->msolve(firk_mem, firk_mem->tempv3, firk_mem->nlscoef);
    if (retval != 0) { return FIRK_MASSSOLVE_FAIL; }
  }

  N_VLinearSum(ONE / hg, firk_mem->tempv3, -ONE / hg, firk_mem->ydn,
               firk_mem->tempv3);

  *yddnrm = N_VWrmsNorm(firk_mem->tempv3, firk_mem->ewt);

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkStep

  Performs one internal step from tn to tn + h, with retries on
  nonlinear solver or error test failures.
  ---------------------------------------------------------------*/
static int firkStep(FIRKodeMem firk_mem)
{
  int i, s, nflag, kflag, eflag, ncf, nef, retval;
  sunrealtype dsm;
  sunbooleantype after_reject;

  s     = firk_mem->s;
  nflag = FIRST_CALL;
  ncf   = 0;
  nef   = 0;
  dsm   = ZERO;

  for (;;)
  {
    firk_mem->nst_attempts++;

    firk_mem->h     = firk_mem->hprime;
    firk_mem->gamma = firk_mem->h * firk_mem->T->gamma0;
    firk_mem->gamrat =
      (firk_mem->gammap != ZERO) ? firk_mem->gamma / firk_mem->gammap : ONE;
    after_reject = (nflag == PREV_ERR_FAIL);

    SUNLogInfo(FIRK_LOGGER, "begin-step-attempt",
               "step = %li, tn = " SUN_FORMAT_G ", h = " SUN_FORMAT_G
               ", stages = %i",
               firk_mem->nst + 1, firk_mem->tn, firk_mem->h, s);

    /* predict the stage increments */
    retval = firkPredict(firk_mem);
    if (retval != FIRK_SUCCESS) { return retval; }

    /* solve the stage system */
    nflag = firkNls(firk_mem, nflag);

    kflag = firkHandleNFlag(firk_mem, &nflag, &ncf);
    if (kflag == PREDICT_AGAIN)
    {
      SUNLogInfo(FIRK_LOGGER, "end-step-attempt",
                 "status = failed solve, dsm = " SUN_FORMAT_G, dsm);
      continue;
    }
    if (kflag != DO_ERROR_TEST)
    {
      SUNLogInfo(FIRK_LOGGER, "end-step-attempt",
                 "status = failed solve, kflag = %i", kflag);
      return kflag;
    }

    /* stage completion: ycur = yn + sum_i d_i Z_i */
    firk_mem->cvals[0] = ONE;
    firk_mem->Xvecs[0] = firk_mem->yn;
    for (i = 0; i < s; i++)
    {
      firk_mem->cvals[i + 1] = firk_mem->T->d[i];
      firk_mem->Xvecs[i + 1] = firk_mem->Z[i];
    }
    retval = N_VLinearCombination(s + 1, firk_mem->cvals, firk_mem->Xvecs,
                                  firk_mem->ycur);
    if (retval != 0) { return FIRK_VECTOROP_ERR; }

    if (firk_mem->fixedstep)
    {
      dsm = ZERO;
      break;
    }

    /* error estimate */
    eflag = firkErrorEstimate(firk_mem, after_reject, &dsm);
    if (eflag > 0)
    {
      nflag = eflag;
      kflag = firkHandleNFlag(firk_mem, &nflag, &ncf);
      if (kflag == PREDICT_AGAIN)
      {
        SUNLogInfo(FIRK_LOGGER, "end-step-attempt",
                   "status = failed error estimate");
        continue;
      }
      return kflag;
    }
    if (eflag < 0)
    {
      SUNLogInfo(FIRK_LOGGER, "end-step-attempt",
                 "status = failed error estimate, eflag = %i", eflag);
      return eflag;
    }

    /* error test */
    eflag = firkDoErrorTest(firk_mem, &nflag, &nef, dsm);
    if (eflag == TRY_AGAIN)
    {
      SUNLogInfo(FIRK_LOGGER, "end-step-attempt",
                 "status = failed error test, dsm = " SUN_FORMAT_G, dsm);
      continue;
    }
    if (eflag != FIRK_SUCCESS)
    {
      SUNLogInfo(FIRK_LOGGER, "end-step-attempt",
                 "status = failed error test, eflag = %i", eflag);
      return eflag;
    }

    break;
  }

  SUNLogInfo(FIRK_LOGGER, "end-step-attempt",
             "status = success, dsm = " SUN_FORMAT_G, dsm);

  return firkCompleteStep(firk_mem, dsm);
}

/*---------------------------------------------------------------
  firkPredict

  Predicts the stage increments, either as zero or by extrapolating
  the collocation polynomial of the previous step:
    Zpred_i = u_prev(tn + c_i h) - yn
            = sum_j (L_j(theta_i) - d_j) Zprev_j,  theta_i = 1 + c_i h/hold.
  ---------------------------------------------------------------*/
static int firkPredict(FIRKodeMem firk_mem)
{
  int i, j, m, s, retval;
  sunrealtype theta, tpow;
  FIRKodeTable T;

  s = firk_mem->s;
  T = firk_mem->T;

  if (firk_mem->predictor == FIRK_PREDICT_TRIVIAL || !firk_mem->have_prev ||
      firk_mem->hold == ZERO)
  {
    for (i = 0; i < s; i++) { N_VConst(ZERO, firk_mem->Zpred[i]); }
    return FIRK_SUCCESS;
  }

  for (i = 0; i < s; i++)
  {
    theta = ONE + T->c[i] * firk_mem->h / firk_mem->hold;
    for (j = 0; j < s; j++) { firk_mem->wtmp[j] = -T->d[j]; }
    tpow = ONE;
    for (m = 0; m < s; m++)
    {
      tpow *= theta;
      for (j = 0; j < s; j++) { firk_mem->wtmp[j] += tpow * T->P[m][j]; }
    }
    retval = N_VLinearCombination(s, firk_mem->wtmp, firk_mem->Zprev,
                                  firk_mem->Zpred[i]);
    if (retval != 0) { return FIRK_VECTOROP_ERR; }
  }

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkNls

  Solves the stacked stage system with the nonlinear solver.
  ---------------------------------------------------------------*/
static int firkNls(FIRKodeMem firk_mem, int nflag)
{
  int i, flag;
  sunbooleantype callSetup;
  long int nni_inc = 0;
  long int nnf_inc = 0;

  /* decide whether to call the linear solver setup routine */
  if (firk_mem->lsetup)
  {
    firk_mem->convfail = ((nflag == FIRST_CALL) || (nflag == PREV_ERR_FAIL))
                           ? FIRK_NO_FAILURES
                           : FIRK_FAIL_OTHER;

    callSetup = (nflag == PREV_CONV_FAIL) || (nflag == PREV_ERR_FAIL) ||
                (firk_mem->firststage) || (firk_mem->gammap == ZERO) ||
                (firk_mem->nst >= firk_mem->nstlp + firk_mem->msbp) ||
                (SUNRabs(firk_mem->gamrat - ONE) > firk_mem->dgmax_lsetup) ||
                ((firk_mem->theta_jbad < ONE) &&
                 (firk_mem->theta > firk_mem->theta_jbad));
  }
  else { callSetup = SUNFALSE; }

  /* initial guess for the correction to the predictor */
  N_VConst(ZERO, firk_mem->zcor_stk);
  firk_mem->nls_reason = FIRK_NLS_REASON_NONE;

  /* call nonlinear solver setup if it exists */
  if ((firk_mem->NLS)->ops->setup)
  {
    flag = SUNNonlinSolSetup(firk_mem->NLS, firk_mem->zcor_stk, firk_mem);
    if (flag < 0) { return FIRK_NLS_SETUP_FAIL; }
    if (flag > 0) { return SUN_NLS_CONV_RECVR; }
  }

  SUNLogInfo(FIRK_LOGGER, "begin-nonlinear-solve", "tol = " SUN_FORMAT_G,
             firk_mem->nlscoef);

  /* solve the nonlinear system */
  flag = SUNNonlinSolSolve(firk_mem->NLS, firk_mem->zpred_stk,
                           firk_mem->zcor_stk, firk_mem->ewt_stk,
                           firk_mem->nlscoef, callSetup, firk_mem);

  /* increment counters */
  (void)SUNNonlinSolGetNumIters(firk_mem->NLS, &nni_inc);
  firk_mem->nni += nni_inc;
  firk_mem->nni_last = nni_inc;

  (void)SUNNonlinSolGetNumConvFails(firk_mem->NLS, &nnf_inc);
  firk_mem->nnf += nnf_inc;

  if (flag != SUN_SUCCESS)
  {
    SUNLogInfo(FIRK_LOGGER, "end-nonlinear-solve",
               "status = failed, flag = %i, iters = %li", flag, nni_inc);
    return flag;
  }

  /* final stage increments */
  for (i = 0; i < firk_mem->s; i++)
  {
    N_VLinearSum(ONE, firk_mem->Zpred[i], ONE, firk_mem->Zcor[i], firk_mem->Z[i]);
  }

  SUNLogInfo(FIRK_LOGGER, "end-nonlinear-solve",
             "status = success, iters = %li", nni_inc);

  firk_mem->jcur = SUNFALSE;

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkHandleNFlag

  Handles the return value of the nonlinear solve (or a recoverable
  failure in the error estimate). Returns DO_ERROR_TEST on success,
  PREDICT_AGAIN to retry the step with a smaller step size, or a
  negative FIRKODE return code.
  ---------------------------------------------------------------*/
static int firkHandleNFlag(FIRKodeMem firk_mem, int* nflagPtr, int* ncfPtr)
{
  int nflag;

  nflag = *nflagPtr;

  if (nflag == FIRK_SUCCESS) { return DO_ERROR_TEST; }

  firk_mem->ncfn++;

  /* unrecoverable failures */
  if (nflag < 0)
  {
    if (nflag == FIRK_LSETUP_FAIL) { return FIRK_LSETUP_FAIL; }
    else if (nflag == FIRK_LSOLVE_FAIL) { return FIRK_LSOLVE_FAIL; }
    else if (nflag == FIRK_RHSFUNC_FAIL) { return FIRK_RHSFUNC_FAIL; }
    else if (nflag == FIRK_MASSMULT_FAIL) { return FIRK_MASSMULT_FAIL; }
    else if (nflag == FIRK_NLS_SETUP_FAIL) { return FIRK_NLS_SETUP_FAIL; }
    else { return FIRK_NLS_FAIL; }
  }

  /* recoverable failure */
  (*ncfPtr)++;
  firk_mem->hadapt_mem->etamax = ONE;

  if ((SUNRabs(firk_mem->h) <= firk_mem->hmin * ONEPSM) ||
      (*ncfPtr == firk_mem->maxncf))
  {
    if (nflag == RHSFUNC_RECVR) { return FIRK_REPTD_RHSFUNC_ERR; }
    return FIRK_CONV_FAILURE;
  }

  /* reduce the step size and retry */
  if (firk_mem->nls_reason == FIRK_NLS_REASON_PREDICT)
  {
    firk_mem->eta = firk_mem->eta_nlsfail;
  }
  else { firk_mem->eta = firk_mem->hadapt_mem->etacf; }
  firk_mem->eta = SUNMAX(firk_mem->eta, firk_mem->hmin / SUNRabs(firk_mem->h));
  firk_mem->hprime = firk_mem->h * firk_mem->eta;
  *nflagPtr        = PREV_CONV_FAIL;

  return PREDICT_AGAIN;
}

/*---------------------------------------------------------------
  firkErrorEstimate

  Computes the RADAU5-style filtered local error estimate
    ele = (M - gamma J)^{-1} [ gamma f(tn,yn) + sum_i e_i M Z_i ]
  and its WRMS norm dsm. On the first step, or after a rejected
  step, a large estimate is recomputed with f evaluated at yn + ele.
  Returns FIRK_SUCCESS, a positive recoverable flag (RHSFUNC_RECVR
  or SUN_NLS_CONV_RECVR), or a negative failure code.
  ---------------------------------------------------------------*/
static int firkErrorEstimate(FIRKodeMem firk_mem, sunbooleantype after_reject,
                             sunrealtype* dsmPtr)
{
  int i, s, retval;
  N_Vector w;
  sunrealtype dsm;

  s = firk_mem->s;

  /* w = sum_i e_i Z_i, then M w if a mass matrix is present */
  for (i = 0; i < s; i++)
  {
    firk_mem->cvals[i] = firk_mem->T->e[i];
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
  retval = firk_mem->lsolve_blk(firk_mem, firk_mem->ele, firk_mem->ewt);
  if (retval < 0) { return FIRK_LSOLVE_FAIL; }
  if (retval > 0) { return SUN_NLS_CONV_RECVR; }
  dsm = N_VWrmsNorm(firk_mem->ele, firk_mem->ewt);

  SUNLogDebug(FIRK_LOGGER, "error-estimate", "dsm = " SUN_FORMAT_G, dsm);

  /* on the first step or after a rejection, filter a large estimate again
     using f evaluated at the perturbed state */
  if (dsm >= ONE && firk_mem->refilter && (firk_mem->firststage || after_reject))
  {
    N_VLinearSum(ONE, firk_mem->yn, ONE, firk_mem->ele, firk_mem->tempv1);
    retval = firk_mem->f(firk_mem->tn, firk_mem->tempv1, firk_mem->ftemp,
                         firk_mem->user_data);
    firk_mem->nfe++;
    if (retval < 0) { return FIRK_RHSFUNC_FAIL; }
    if (retval > 0) { return RHSFUNC_RECVR; }

    N_VLinearSum(firk_mem->gamma, firk_mem->ftemp, ONE, w, firk_mem->ele);
    retval = firk_mem->lsolve_blk(firk_mem, firk_mem->ele, firk_mem->ewt);
    if (retval < 0) { return FIRK_LSOLVE_FAIL; }
    if (retval > 0) { return SUN_NLS_CONV_RECVR; }
    dsm = N_VWrmsNorm(firk_mem->ele, firk_mem->ewt);

    SUNLogDebug(FIRK_LOGGER, "error-estimate-filtered", "dsm = " SUN_FORMAT_G,
                dsm);
  }

  *dsmPtr = SUNMAX(dsm, DSM_FLOOR);

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkDoErrorTest

  Performs the local error test and selects the next step size.
  Returns FIRK_SUCCESS if the test passes, TRY_AGAIN to retry the
  step, or FIRK_ERR_FAILURE / FIRK_CONTROLLER_ERR.
  ---------------------------------------------------------------*/
static int firkDoErrorTest(FIRKodeMem firk_mem, int* nflagPtr, int* nefPtr,
                           sunrealtype dsm)
{
  int retval;
  FIRKodeHAdaptMem hadapt_mem = firk_mem->hadapt_mem;

  /* consider a change of step size for the next attempt */
  retval = firkAdapt(firk_mem, dsm);
  if (retval != FIRK_SUCCESS) { return retval; }

  firk_mem->eta = SUNMIN(firk_mem->eta, hadapt_mem->etamax);
  firk_mem->eta = SUNMAX(firk_mem->eta, firk_mem->hmin / SUNRabs(firk_mem->h));
  firk_mem->eta /=
    SUNMAX(ONE, SUNRabs(firk_mem->h) * firk_mem->hmax_inv * firk_mem->eta);

  if (dsm <= ONE) { return FIRK_SUCCESS; }

  /* test failed */
  (*nefPtr)++;
  firk_mem->netf++;
  *nflagPtr = PREV_ERR_FAIL;

  if ((*nefPtr == firk_mem->maxnef) ||
      (SUNRabs(firk_mem->h) <= firk_mem->hmin * ONEPSM))
  {
    return FIRK_ERR_FAILURE;
  }

  /* On a rejection the step must shrink. A predictive controller may still
     propose growth (its history refers to the last accepted step), so the
     step is also bounded by the basic error-based reduction (RADAU5). */
  hadapt_mem->etamax = ONE;
  firk_mem->eta =
    SUNMIN(firk_mem->eta,
           hadapt_mem->safety *
             SUNRpowerR(dsm, -ONE / ((sunrealtype)hadapt_mem->p + ONE)));
  firk_mem->eta = SUNMAX(firk_mem->eta, hadapt_mem->etamin);

  if (*nefPtr >= hadapt_mem->small_nef)
  {
    firk_mem->eta = SUNMIN(firk_mem->eta, hadapt_mem->etamxf);
  }

  firk_mem->eta = SUNMIN(firk_mem->eta, hadapt_mem->etamax);
  firk_mem->eta = SUNMAX(firk_mem->eta, firk_mem->hmin / SUNRabs(firk_mem->h));
  firk_mem->eta /=
    SUNMAX(ONE, SUNRabs(firk_mem->h) * firk_mem->hmax_inv * firk_mem->eta);

  firk_mem->hprime = firk_mem->h * firk_mem->eta;

  return TRY_AGAIN;
}

/*---------------------------------------------------------------
  firkCompleteStep

  Updates the integrator state after a successful step: advances
  the time, stores the stage increments for dense output and the
  next predictor, and notifies the step size controller.
  ---------------------------------------------------------------*/
static int firkCompleteStep(FIRKodeMem firk_mem, sunrealtype dsm)
{
  int retval;
  sunrealtype troundoff;
  N_Vector vtmp;
  N_Vector* atmp;

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
     the next predictor; the wrappers zpred_stk/zcor_stk are unaffected */
  atmp                = firk_mem->Zprev;
  firk_mem->Zprev     = firk_mem->Z;
  firk_mem->Z         = atmp;
  firk_mem->have_prev = SUNTRUE;
  firk_mem->hold      = firk_mem->h;

  /* update yn to the new solution (yn <-> ycur); ewt and fn are not
     swapped so the stacked weight wrapper stays valid */
  vtmp                     = firk_mem->yn;
  firk_mem->yn             = firk_mem->ycur;
  firk_mem->ycur           = vtmp;
  firk_mem->fn_is_current  = SUNFALSE;
  firk_mem->ydn_is_current = SUNFALSE;

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

  firk_mem->nst++;
  firk_mem->tn     = firk_mem->tcur;
  firk_mem->hprime = firk_mem->h * firk_mem->eta;
  firk_mem->dsm    = dsm;

  firk_mem->hadapt_mem->etamax = firk_mem->hadapt_mem->growth;

  firk_mem->initsetup  = SUNFALSE;
  firk_mem->firststage = SUNFALSE;

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkHandleFailure

  Prints error messages for all cases of failure by firkHin and
  firkStep and returns the value FIRKodeEvolve returns to the user.
  ---------------------------------------------------------------*/
static int firkHandleFailure(FIRKodeMem firk_mem, int flag)
{
  switch (flag)
  {
  case FIRK_ERR_FAILURE:
    firkProcessError(firk_mem, FIRK_ERR_FAILURE, __LINE__, __func__, __FILE__,
                     MSGFIRK_ERR_FAILS, firk_mem->tcur, firk_mem->h);
    break;
  case FIRK_CONV_FAILURE:
    firkProcessError(firk_mem, FIRK_CONV_FAILURE, __LINE__, __func__, __FILE__,
                     MSGFIRK_CONV_FAILS, firk_mem->tcur, firk_mem->h);
    break;
  case FIRK_LSETUP_FAIL:
    firkProcessError(firk_mem, FIRK_LSETUP_FAIL, __LINE__, __func__, __FILE__,
                     MSGFIRK_SETUP_FAILED, firk_mem->tcur);
    break;
  case FIRK_LSOLVE_FAIL:
    firkProcessError(firk_mem, FIRK_LSOLVE_FAIL, __LINE__, __func__, __FILE__,
                     MSGFIRK_SOLVE_FAILED, firk_mem->tcur);
    break;
  case FIRK_RHSFUNC_FAIL:
    firkProcessError(firk_mem, FIRK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                     MSGFIRK_RHSFUNC_FAILED, firk_mem->tcur);
    break;
  case FIRK_UNREC_RHSFUNC_ERR:
    firkProcessError(firk_mem, FIRK_UNREC_RHSFUNC_ERR, __LINE__, __func__,
                     __FILE__, MSGFIRK_RHSFUNC_UNREC, firk_mem->tcur);
    break;
  case FIRK_REPTD_RHSFUNC_ERR:
    firkProcessError(firk_mem, FIRK_REPTD_RHSFUNC_ERR, __LINE__, __func__,
                     __FILE__, MSGFIRK_RHSFUNC_REPTD, firk_mem->tcur);
    break;
  case FIRK_RTFUNC_FAIL:
    firkProcessError(firk_mem, FIRK_RTFUNC_FAIL, __LINE__, __func__, __FILE__,
                     MSGFIRK_RTFUNC_FAILED, firk_mem->tcur);
    break;
  case FIRK_TOO_CLOSE:
    firkProcessError(firk_mem, FIRK_TOO_CLOSE, __LINE__, __func__, __FILE__,
                     MSGFIRK_TOO_CLOSE);
    break;
  case FIRK_MEM_NULL:
    firkProcessError(NULL, FIRK_MEM_NULL, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MEM);
    break;
  case SUN_ERR_ARG_CORRUPT:
    firkProcessError(firk_mem, FIRK_MEM_NULL, __LINE__, __func__, __FILE__,
                     MSGFIRK_NLS_INPUT_NULL, firk_mem->tcur);
    break;
  case FIRK_NLS_SETUP_FAIL:
    firkProcessError(firk_mem, FIRK_NLS_SETUP_FAIL, __LINE__, __func__,
                     __FILE__, MSGFIRK_NLS_SETUP_FAILED, firk_mem->tcur);
    break;
  case FIRK_NLS_FAIL:
    firkProcessError(firk_mem, FIRK_NLS_FAIL, __LINE__, __func__, __FILE__,
                     MSGFIRK_NLS_FAIL, firk_mem->tcur);
    break;
  case FIRK_MASSSOLVE_FAIL:
    firkProcessError(firk_mem, FIRK_MASSSOLVE_FAIL, __LINE__, __func__,
                     __FILE__, MSGFIRK_MASSSOLVE_FAIL, firk_mem->tcur);
    break;
  case FIRK_MASSMULT_FAIL:
    firkProcessError(firk_mem, FIRK_MASSMULT_FAIL, __LINE__, __func__, __FILE__,
                     MSGFIRK_MASSMULT_FAIL, firk_mem->tcur);
    break;
  case FIRK_CONTROLLER_ERR:
    firkProcessError(firk_mem, FIRK_CONTROLLER_ERR, __LINE__, __func__,
                     __FILE__, MSGFIRK_CONTROLLER_ERR, firk_mem->tcur);
    break;
  case FIRK_VECTOROP_ERR:
    firkProcessError(firk_mem, FIRK_VECTOROP_ERR, __LINE__, __func__, __FILE__,
                     "A vector operation failed.");
    break;
  default:
    firkProcessError(firk_mem, FIRK_UNRECOGNIZED_ERR, __LINE__, __func__,
                     __FILE__,
                     "FIRKODE encountered an unrecognized error. Please "
                     "report this to the SUNDIALS developers at "
                     "sundials-users@llnl.gov");
    return FIRK_UNRECOGNIZED_ERR;
  }

  return flag;
}

/*---------------------------------------------------------------
  Error weight functions
  ---------------------------------------------------------------*/

int firkEwtSet(N_Vector ycur, N_Vector weight, void* data)
{
  FIRKodeMem firk_mem;
  int flag = 0;

  firk_mem = (FIRKodeMem)data;

  switch (firk_mem->itol)
  {
  case FIRK_SS: flag = firkEwtSetSS(firk_mem, ycur, weight); break;
  case FIRK_SV: flag = firkEwtSetSV(firk_mem, ycur, weight); break;
  }

  return flag;
}

static int firkEwtSetSS(FIRKodeMem firk_mem, N_Vector ycur, N_Vector weight)
{
  N_VAbs(ycur, firk_mem->tempv1);
  N_VScale(firk_mem->reltol, firk_mem->tempv1, firk_mem->tempv1);
  N_VAddConst(firk_mem->tempv1, firk_mem->Sabstol, firk_mem->tempv1);
  if (firk_mem->atolmin0)
  {
    if (N_VMin(firk_mem->tempv1) <= ZERO) { return -1; }
  }
  N_VInv(firk_mem->tempv1, weight);
  return 0;
}

static int firkEwtSetSV(FIRKodeMem firk_mem, N_Vector ycur, N_Vector weight)
{
  N_VAbs(ycur, firk_mem->tempv1);
  N_VLinearSum(firk_mem->reltol, firk_mem->tempv1, ONE, firk_mem->Vabstol,
               firk_mem->tempv1);
  if (firk_mem->atolmin0)
  {
    if (N_VMin(firk_mem->tempv1) <= ZERO) { return -1; }
  }
  N_VInv(firk_mem->tempv1, weight);
  return 0;
}

/*---------------------------------------------------------------
  Error message handling
  ---------------------------------------------------------------*/

void firkProcessError(FIRKodeMem firk_mem, int error_code, int line,
                      const char* func, const char* file, const char* msgfmt, ...)
{
  va_list ap;
  size_t msglen;
  char* msg;

  va_start(ap, msgfmt);
  msglen = 1;
  if (msgfmt) { msglen += (size_t)vsnprintf(NULL, 0, msgfmt, ap); }
  va_end(ap);

  msg = (char*)malloc(msglen);
  if (msg == NULL) { return; }

  va_start(ap, msgfmt);
  vsnprintf(msg, msglen, msgfmt, ap);
  va_end(ap);

  do {
    if (firk_mem == NULL)
    {
      SUNGlobalFallbackErrHandler(line, func, file, msg, error_code);
      break;
    }

    if (error_code == FIRK_WARNING)
    {
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_WARNING
      char* file_and_line = sunCombineFileAndLine(line, file);
      SUNLogger_QueueMsg(FIRK_LOGGER, SUN_LOGLEVEL_WARNING, file_and_line, func,
                         msg);
      free(file_and_line);
#endif
      break;
    }

    SUNHandleErrWithMsg(line, func, file, msg, error_code, firk_mem->sunctx);

    (void)SUNContext_GetLastError(firk_mem->sunctx);
  }
  while (0);

  free(msg);
}
