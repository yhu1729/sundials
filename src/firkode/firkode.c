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
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_types.h>

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

static int firkMassSolveDeriv(FIRKodeMem firk_mem, N_Vector b);
static int firkHin(FIRKodeMem firk_mem, sunrealtype tout);
static sunrealtype firkUpperBoundH0(FIRKodeMem firk_mem, sunrealtype tdist);
static int firkYddNorm(FIRKodeMem firk_mem, sunrealtype hg, sunrealtype* yddnrm);

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
  firk_mem->mode    = FIRK_BEGINNER;
  firk_mem->s_alloc = 0;

  /* integrator defaults */
  firk_mem->mxstep      = MXSTEP_DEFAULT;
  firk_mem->mxhnil      = MXHNIL_DEFAULT;
  firk_mem->hin         = ZERO;
  firk_mem->hmin        = HMIN_DEFAULT;
  firk_mem->hmax_inv    = HMAX_INV_DEFAULT;
  firk_mem->fixedstep   = SUNFALSE;
  firk_mem->tstopset    = SUNFALSE;
  firk_mem->tstopinterp = SUNFALSE;

  /* Newton iteration state */
  firk_mem->eta_nls  = ONE;
  firk_mem->theta    = ZERO;
  firk_mem->nls_hfac = ONE;

  /* adaptivity defaults */
  firk_mem->hadapt_mem = firkAdaptInit();
  if (firk_mem->hadapt_mem == NULL)
  {
    firkProcessError(NULL, 0, __LINE__, __func__, __FILE__, MSGFIRK_MEM_FAIL);
    free(firk_mem);
    return NULL;
  }

  /* defaults of the EXPERT options, including the controller */
  retval = firkSetExpertDefaults(firk_mem);
  if (retval != FIRK_SUCCESS)
  {
    FIRKodeTable_Free(firk_mem->tbl);
    free(firk_mem->hadapt_mem);
    free(firk_mem);
    return NULL;
  }

  /* linear solver and mass matrix hooks */
  firk_mem->linit          = NULL;
  firk_mem->lreinit        = NULL;
  firk_mem->lsetup         = NULL;
  firk_mem->lsolve_stk     = NULL;
  firk_mem->lsolve_blk     = NULL;
  firk_mem->lfree          = NULL;
  firk_mem->lmem           = NULL;
  firk_mem->minit          = NULL;
  firk_mem->mreinit        = NULL;
  firk_mem->msetup         = NULL;
  firk_mem->mmult          = NULL;
  firk_mem->msolve         = NULL;
  firk_mem->mfree          = NULL;
  firk_mem->mass_mem       = NULL;
  firk_mem->mass_set       = SUNFALSE;
  firk_mem->mass_iterative = SUNFALSE;

  /* rootfinding */
  firk_mem->root_mem = NULL;

  firk_mem->VabstolMallocDone = SUNFALSE;
  firk_mem->MallocDone        = SUNFALSE;

  return (void*)firk_mem;
}

/*---------------------------------------------------------------
  firkSetExpertDefaults

  Sets every option that only EXPERT mode may change to its default
  (design decision D8): the method, the predictor and refilter, the
  failure limits, the Newton parameters, the reuse frequencies, and
  the step-size controller and heuristics. Used by FIRKodeCreate and
  by FIRKodeSetMode when it returns to BEGINNER mode.
  ---------------------------------------------------------------*/
int firkSetExpertDefaults(FIRKodeMem firk_mem)
{
  firk_mem->s = FIRK_DEFAULT_STAGES;
  if (firk_mem->tbl == NULL || firk_mem->tbl->s != FIRK_DEFAULT_STAGES)
  {
    FIRKodeTable_Free(firk_mem->tbl);
    firk_mem->tbl = FIRKodeTable_RadauIIA(FIRK_DEFAULT_STAGES);
    if (firk_mem->tbl == NULL)
    {
      firkProcessError(firk_mem, FIRK_MEM_FAIL, __LINE__, __func__, __FILE__,
                       MSGFIRK_MEM_FAIL);
      return FIRK_MEM_FAIL;
    }
  }
  firk_mem->predictor = FIRK_PREDICT_EXTRAPOLATE;
  firk_mem->refilter  = SUNTRUE;

  firk_mem->maxnef       = MXNEF;
  firk_mem->maxncf       = MXNCF;
  firk_mem->msbp         = MSBP_DEFAULT;
  firk_mem->dgmax_lsetup = DGMAX_LSETUP_DEFAULT;

  firk_mem->maxcor       = MAXCOR_DEFAULT;
  firk_mem->nlscoef_user = NLSCOEF_AUTO;
  firkSetNlsCoef(firk_mem);
  firk_mem->theta_max    = THETA_MAX_DEFAULT;
  firk_mem->theta_jbad   = THETA_JBAD_DEFAULT;
  firk_mem->theta_jac    = THETA_JAC_DEFAULT;
  firk_mem->rate_exp     = RATE_EXP_DEFAULT;
  firk_mem->reuse_policy = FIRK_POLICY_RADAU5;

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

  /* the default controller (this also restores its error bias) */
  return firkReplaceAdaptController(firk_mem, NULL, SUNTRUE);
}

/*---------------------------------------------------------------
  firkSetNlsCoef

  Sets the Newton stopping tolerance in use (design Section 4.4):
  the value fixed by FIRKodeSetNonlinConvCoef if positive, and
  otherwise RADAU5's rule
    eps_nls = max(10 u / reltol, min(0.03, sqrt(reltol)))
  for a scalar relative tolerance, or its cap 0.03 before the
  tolerances are set, for reltol = 0 and for a user weight function,
  which give no relative accuracy to derive it from. Called whenever
  its inputs change: the tolerance setters, the EXPERT setter and
  the return to the EXPERT defaults.
  ---------------------------------------------------------------*/
void firkSetNlsCoef(FIRKodeMem firk_mem)
{
  sunrealtype eps;

  if (firk_mem->nlscoef_user > ZERO)
  {
    firk_mem->nlscoef = firk_mem->nlscoef_user;
    return;
  }

  eps = NLSCOEF_MAX;
  if ((firk_mem->itol == FIRK_SS || firk_mem->itol == FIRK_SV) &&
      firk_mem->reltol > ZERO)
  {
    eps = SUNMAX(NLSCOEF_UROUND_FAC * firk_mem->uround / firk_mem->reltol,
                 SUNMIN(NLSCOEF_MAX, SUNRsqrt(firk_mem->reltol)));
  }
  firk_mem->nlscoef = eps;
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

  /* the kind of the stages' stack */
  retval = firkStackSetKind(firk_mem, y0);
  if (retval != FIRK_SUCCESS)
  {
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return retval;
  }

  /* allocate the N-length work vectors (using y0 as a template); a second
     call frees the vectors of the first, among them the vector absolute
     tolerances, so the tolerances must be set again */
  if (firk_mem->MallocDone)
  {
    firkFreeVectors(firk_mem);
    firk_mem->MallocDone = SUNFALSE;
    firk_mem->itol       = FIRK_NN;
    firk_mem->user_efun  = SUNFALSE;
    firk_mem->efun       = NULL;
    firk_mem->e_data     = NULL;
    firkSetNlsCoef(firk_mem);
  }
  if (!firkAllocVectors(firk_mem, y0))
  {
    firkProcessError(firk_mem, FIRK_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSGFIRK_MEM_FAIL);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_MEM_FAIL;
  }

  /* the stages and their stacks for the current number of stages: vectors
     are cloned here and at an order change only (design Section 8.3, H14) */
  firkFreeStageVectors(firk_mem);
  if (!firkAllocStageVectors(firk_mem))
  {
    firkFreeVectors(firk_mem);
    firkProcessError(firk_mem, FIRK_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSGFIRK_MEM_FAIL);
    SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
    return FIRK_MEM_FAIL;
  }

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
  firk_mem->e_data    = firk_mem;

  firkSetNlsCoef(firk_mem);

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
  firk_mem->e_data    = firk_mem;

  firkSetNlsCoef(firk_mem);

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
  firk_mem->e_data    = firk_mem->user_data;

  firkSetNlsCoef(firk_mem);

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  FIRKodeEvolve

  Main solver function. In FIRK_NORMAL mode the solver steps until
  it reaches or passes tout and then interpolates to obtain y(tout).
  In FIRK_ONE_STEP mode it takes one internal step and returns.
  ---------------------------------------------------------------*/
int FIRKodeEvolve(void* firkode_mem, sunrealtype tout, N_Vector yout,
                  sunrealtype* tret_ptr, int itask)
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

  if (tret_ptr == NULL)
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
    firk_mem->tretlast = *tret_ptr = firk_mem->tcur;

    ier = firkInitialSetup(firk_mem, tout);
    if (ier != FIRK_SUCCESS)
    {
      SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
      return ier;
    }

    /* evaluate f at (t0, y0) */
    retval = firkRhs(firk_mem, firk_mem->tcur, firk_mem->yn, firk_mem->fn);
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
      if (retval < 0)
      {
        /* the root function, the right-hand side or the mass solve failed:
           reported by its code */
        istate = firkHandleFailure(firk_mem, retval);
        SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
        return istate;
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
      else if (retval < 0)
      {
        /* the root function, the right-hand side or the mass solve failed:
           reported by its code */
        istate = firkHandleFailure(firk_mem, retval);
        SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
        return istate;
      }
      else if (retval == RTFOUND)
      {
        /* firkRootCheck2 leaves y(tlo) in ycur */
        firk_mem->tretlast = *tret_ptr = rootmem->tlo;
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
            firk_mem->tretlast = *tret_ptr = firk_mem->tcur;
            N_VScale(ONE, firk_mem->yn, yout);
            SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
            return FIRK_SUCCESS;
          }
        }
        else if (retval == RTFOUND)
        {
          /* firkRootCheck3 leaves y(tlo) in ycur */
          rootmem->irfnd     = 1;
          firk_mem->tretlast = *tret_ptr = rootmem->tlo;
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
          firk_mem->tretlast = *tret_ptr = firk_mem->tstop;
          firk_mem->tstopset             = SUNFALSE;
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
      firk_mem->tretlast = *tret_ptr = tout;
      ier                            = FIRKodeGetDky(firk_mem, tout, 0, yout);
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
      firk_mem->tretlast = *tret_ptr = firk_mem->tcur;
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
        firk_mem->tretlast = *tret_ptr = firk_mem->tcur;
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
      firk_mem->tretlast = *tret_ptr = firk_mem->tcur;
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
      firk_mem->tretlast = *tret_ptr = firk_mem->tcur;
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
      retval = firkRhs(firk_mem, firk_mem->tn, firk_mem->yn, firk_mem->fn);
      if (retval != FIRK_SUCCESS)
      {
        istate = (retval < 0) ? FIRK_RHSFUNC_FAIL : FIRK_UNREC_RHSFUNC_ERR;
        istate = firkHandleFailure(firk_mem, istate);
        firk_mem->tretlast = *tret_ptr = firk_mem->tcur;
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
      firk_mem->tretlast = *tret_ptr = firk_mem->tcur;
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
        firk_mem->tretlast = *tret_ptr = rootmem->tlo;
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
          firk_mem->tretlast = *tret_ptr = firk_mem->tstop;
          firk_mem->tstopset             = SUNFALSE;
          istate                         = FIRK_TSTOP_RETURN;
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
      firk_mem->tretlast = *tret_ptr = tout;
      (void)FIRKodeGetDky(firk_mem, tout, 0, yout);
      firk_mem->next_h = firk_mem->hprime;
      break;
    }

    /* in ONE_STEP mode, copy y and exit loop */
    if (itask == FIRK_ONE_STEP)
    {
      istate             = FIRK_SUCCESS;
      firk_mem->tretlast = *tret_ptr = firk_mem->tcur;
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
  sunrealtype tfuzz, tp, tn1;
  int ier;

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

  ier = firkDenseEval(firk_mem, t, k, dky);

  SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
  return ier;
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

  if (firk_mem->lfree != NULL) { firk_mem->lfree(firk_mem); }
  if (firk_mem->mfree != NULL) { firk_mem->mfree(firk_mem); }

  (void)firkRootFree(firk_mem);

  FIRKodeTable_Free(firk_mem->tbl);
  firk_mem->tbl = NULL;

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
  int s;

  s = firk_mem->s;

  firk_mem->Z     = firkCloneVectorArray(s, firk_mem->yn);
  firk_mem->Zprev = firkCloneVectorArray(s, firk_mem->yn);
  firk_mem->F     = firkCloneVectorArray(s, firk_mem->yn);
  firk_mem->G     = firkCloneVectorArray(s, firk_mem->yn);
  firk_mem->delta = firkCloneVectorArray(s, firk_mem->yn);
  firk_mem->cvals = (sunrealtype*)calloc((size_t)(s + 2), sizeof(sunrealtype));
  firk_mem->Xvecs = (N_Vector*)calloc((size_t)(s + 2), sizeof(N_Vector));
  firk_mem->wtmp  = (sunrealtype*)calloc((size_t)(s + 2), sizeof(sunrealtype));
  if (firk_mem->Z == NULL || firk_mem->Zprev == NULL || firk_mem->F == NULL ||
      firk_mem->G == NULL || firk_mem->delta == NULL || firk_mem->cvals == NULL ||
      firk_mem->Xvecs == NULL || firk_mem->wtmp == NULL)
  {
    firkFreeStageVectors(firk_mem);
    return SUNFALSE;
  }

  /* the one-stage method has no stack: its "stacks" are the vectors */
  firk_mem->s_alloc = s;
  if (s == 1)
  {
    firk_mem->g_stk     = firk_mem->G[0];
    firk_mem->delta_stk = firk_mem->delta[0];
    firk_mem->ewt_stk   = firk_mem->ewt;
    return SUNTRUE;
  }

  firk_mem->g_stk     = firkStackNew(firk_mem, s, firk_mem->G);
  firk_mem->delta_stk = firkStackNew(firk_mem, s, firk_mem->delta);
  firk_mem->ewt_stk   = firkStackWeightView(firk_mem, s, firk_mem->ewt);

  if (firk_mem->g_stk == NULL || firk_mem->delta_stk == NULL ||
      firk_mem->ewt_stk == NULL)
  {
    firkFreeStageVectors(firk_mem);
    return SUNFALSE;
  }

  return SUNTRUE;
}

static void firkFreeStageVectors(FIRKodeMem firk_mem)
{
  int s = firk_mem->s_alloc;

  /* for one stage the stacks are the vectors themselves */
  if (s >= 2)
  {
    if (firk_mem->g_stk) { N_VDestroy(firk_mem->g_stk); }
    if (firk_mem->delta_stk) { N_VDestroy(firk_mem->delta_stk); }
    if (firk_mem->ewt_stk) { N_VDestroy(firk_mem->ewt_stk); }
  }
  firk_mem->g_stk = firk_mem->delta_stk = firk_mem->ewt_stk = NULL;

  if (s > 0)
  {
    if (firk_mem->Z) { N_VDestroyVectorArray(firk_mem->Z, s); }
    if (firk_mem->Zprev) { N_VDestroyVectorArray(firk_mem->Zprev, s); }
    if (firk_mem->F) { N_VDestroyVectorArray(firk_mem->F, s); }
    if (firk_mem->G) { N_VDestroyVectorArray(firk_mem->G, s); }
    if (firk_mem->delta) { N_VDestroyVectorArray(firk_mem->delta, s); }
  }
  firk_mem->Z = firk_mem->Zprev = firk_mem->F = NULL;
  firk_mem->G = firk_mem->delta = NULL;

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

  firk_mem->eta_nls        = ONE;
  firk_mem->theta          = ZERO;
  firk_mem->nls_hfac       = ONE;
  firk_mem->nni_last       = 0;
  firk_mem->jcur           = SUNFALSE;
  firk_mem->keep_setup     = SUNFALSE;
  firk_mem->jac_current    = SUNFALSE;
  firk_mem->jac_stale_fail = SUNFALSE;
  firk_mem->dgamma_fail    = ONE;
  firk_mem->nstlp          = 0;

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
    if (firk_mem->lreinit) { firk_mem->lreinit(firk_mem); }
    if (firk_mem->mreinit) { firk_mem->mreinit(firk_mem); }
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

  /* the built-in table for the requested number of stages */
  if (firk_mem->tbl == NULL || firk_mem->tbl->s != firk_mem->s)
  {
    FIRKodeTable_Free(firk_mem->tbl);
    firk_mem->tbl = FIRKodeTable_RadauIIA(firk_mem->s);
  }
  if (firk_mem->tbl == NULL)
  {
    firkProcessError(firk_mem, FIRK_TABLE_FAIL, __LINE__, __func__, __FILE__,
                     MSGFIRK_TABLE_FAIL);
    return FIRK_TABLE_FAIL;
  }
  firk_mem->hadapt_mem->p = firk_mem->tbl->p;
  firk_mem->hadapt_mem->q = firk_mem->tbl->q;

  /* allocate the stage memory for s stages */
  if (firk_mem->s_alloc != firk_mem->s)
  {
    firkFreeStageVectors(firk_mem);
    if (!firkAllocStageVectors(firk_mem))
    {
      firkProcessError(firk_mem, FIRK_MEM_FAIL, __LINE__, __func__, __FILE__,
                       MSGFIRK_MEM_FAIL);
      return FIRK_MEM_FAIL;
    }
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
  firkRhs

  Evaluates ydot = f(t, y) (design decision D12): every evaluation
  of the user's f by the integrator goes through this wrapper, which
  counts it in nfe and maps the return value to FIRK_SUCCESS,
  FIRK_RHSFUNC_FAIL (unrecoverable) or RHSFUNC_RECVR (recoverable).
  The difference quotients of FIRKLS count in nfeDQ instead and may
  use the jtimes RHS function, so they call the functions directly.
  ---------------------------------------------------------------*/
int firkRhs(FIRKodeMem firk_mem, sunrealtype t, N_Vector y, N_Vector ydot)
{
  int retval;

  SUNDIALS_MARK_FUNCTION_BEGIN(FIRK_PROFILER);
  retval = firk_mem->f(t, y, ydot, firk_mem->user_data);
  firk_mem->nfe++;
  SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
  if (retval < 0) { return FIRK_RHSFUNC_FAIL; }
  if (retval > 0) { return RHSFUNC_RECVR; }
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
    retval = firkRhs(firk_mem, firk_mem->tcur, firk_mem->yn, firk_mem->fn);
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
  retval = firkMassSolveDeriv(firk_mem, firk_mem->ydn);
  if (retval != FIRK_SUCCESS) { return retval; }
  firk_mem->ydn_is_current = SUNTRUE;
  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkMassSolveDeriv

  Overwrites b, a value or a difference quotient of f, with M^{-1} b,
  a derivative of y that only feeds the heuristics of the initial
  step size estimate and of rootfinding. The residual tolerance of an
  iterative mass solver is relative to the WRMS norm of b, so that the
  result has the relative accuracy set by FIRKodeSetMassEpsLin
  whatever the scale of M and f; the tolerance of the nonlinear
  solver, in the units of y, does not apply to a residual in the
  units of f and may be unattainable for large f, as for nonsmooth
  initial data. Since the result only feeds heuristics, the mass
  solve of the integrator accepts the approximation of an iterative
  solver that reduced the residual without reaching the tolerance;
  any other failure is fatal, as there is no step to retry (design
  Table 8.6). Returns FIRK_SUCCESS or FIRK_MASSSOLVE_FAIL.
  ---------------------------------------------------------------*/
static int firkMassSolveDeriv(FIRKodeMem firk_mem, N_Vector b)
{
  int retval;

  retval = firk_mem->msolve(firk_mem, b, N_VWrmsNorm(b, firk_mem->ewt));
  return (retval == 0) ? FIRK_SUCCESS : FIRK_MASSSOLVE_FAIL;
}

/*---------------------------------------------------------------
  firkHin

  Computes a tentative initial step size h0 as the solution of
    (WRMS norm of (h0^2 ydd / 2)) = 1,
  where ydd is an estimated second derivative of y (same algorithm
  as CVODE's cvHin).
  ---------------------------------------------------------------*/
static int firkHinBody(FIRKodeMem firk_mem, sunrealtype tout)
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

/* the profiling region of the initial step (design Table 8.5) */
static int firkHin(FIRKodeMem firk_mem, sunrealtype tout)
{
  int retval;
  SUNDIALS_MARK_FUNCTION_BEGIN(FIRK_PROFILER);
  retval = firkHinBody(firk_mem, tout);
  SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
  return retval;
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
  retval = firkRhs(firk_mem, firk_mem->tcur + hg, firk_mem->ycur,
                   firk_mem->tempv3);
  if (retval != FIRK_SUCCESS) { return retval; }

  /* ydd = M^{-1} (f(t + hg, y + hg y') - f(t, y)) / hg. An iterative
     mass solver is applied once to the difference quotient of f, so that
     its error is relative to ydd: the errors of separate inexact solves
     for the two values of y' need not cancel, and their difference would
     be divided by hg. A direct solver is exact up to roundoff and solves
     for y'(t + hg) instead. */
  if (firk_mem->mass_iterative)
  {
    N_VLinearSum(ONE / hg, firk_mem->tempv3, -ONE / hg, firk_mem->fn,
                 firk_mem->tempv3);
    retval = firkMassSolveDeriv(firk_mem, firk_mem->tempv3);
    if (retval != FIRK_SUCCESS) { return retval; }
  }
  else
  {
    if (firk_mem->mass_set)
    {
      retval = firkMassSolveDeriv(firk_mem, firk_mem->tempv3);
      if (retval != FIRK_SUCCESS) { return retval; }
    }
    N_VLinearSum(ONE / hg, firk_mem->tempv3, -ONE / hg, firk_mem->ydn,
                 firk_mem->tempv3);
  }

  *yddnrm = N_VWrmsNorm(firk_mem->tempv3, firk_mem->ewt);

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
