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
 * This is the implementation file for FIRKODE's linear solver
 * interface (FIRKLS).
 *
 * The user attaches an N x N SUNLinearSolver (and, for matrix-based
 * solvers, a SUNMatrix template) for the shifted block system
 *   (M - gamma J) x = b,  gamma = h * gamma0.
 * The coupled stage system of the simplified Newton iteration,
 *   (I_s (x) M - h A (x) J) X = B,
 * is solved by an internal flexible GMRES iteration on the stacked
 * ManyVector, with the block solve applied to each stage as a right
 * preconditioner. For single-stage methods (A = [1], gamma0 = 1)
 * the block solve is exact and the Krylov iteration is skipped.
 * The block solve is also used for the filtered error estimate.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_pcg.h>
#include <sunlinsol/sunlinsol_spbcgs.h>
#include <sunlinsol/sunlinsol_spfgmr.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunlinsol/sunlinsol_sptfqmr.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>

#include "firkode_impl.h"
#include "firkode_ls_impl.h"

/* constants */
#define MIN_INC_MULT SUN_RCONST(1000.0)
#define MAX_DQITERS  3 /* max number of stages tried in DQ Jtimes */
#define ZERO         SUN_RCONST(0.0)
#define PT25         SUN_RCONST(0.25)
#define ONE          SUN_RCONST(1.0)
#define TWO          SUN_RCONST(2.0)

/*=================================================================
  PRIVATE FUNCTION PROTOTYPES
  =================================================================*/

static int firkLsLinSys(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix A,
                        SUNMatrix M, sunbooleantype jok, sunbooleantype* jcur,
                        sunrealtype gamma, void* user_data, N_Vector tmp1,
                        N_Vector tmp2, N_Vector tmp3);
static int firkLsSetupStacked(FIRKodeMem firk_mem, FIRKLsMem firkls_mem);
static void firkLsFreeStacked(FIRKLsMem firkls_mem);
static int firkLsMapSolveFlag(FIRKodeMem firk_mem, FIRKLsMem firkls_mem,
                              int retval, sunbooleantype accept_reduced);
static int firkLsSetupMatrix(FIRKodeMem firk_mem, FIRKLsMem firkls_mem,
                             sunbooleantype jbad, sunbooleantype* jcur);
static int firkLsJtSetup(FIRKodeMem firk_mem, FIRKLsMem firkls_mem);
static int firkLsUserPretype(SUNLinearSolver LS);
static int firkLsUserPSolveBlock(FIRKodeMem firk_mem, FIRKLsMem firkls_mem,
                                 N_Vector r, N_Vector z, int pretype);

/*===============================================================
  FIRKLS Exported functions -- Required
  ===============================================================*/

/*---------------------------------------------------------------
  FIRKodeSetLinearSolver specifies the block linear solver
  ---------------------------------------------------------------*/
int FIRKodeSetLinearSolver(void* firkode_mem, SUNLinearSolver LS, SUNMatrix A)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval, LSType;
  sunbooleantype iterative;   /* is the solver iterative?    */
  sunbooleantype matrixbased; /* is a matrix structure used? */

  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRKLS_MEM_NULL, __LINE__, __func__, __FILE__,
                     MSG_LS_FIRKMEM_NULL);
    return FIRKLS_MEM_NULL;
  }
  if (LS == NULL)
  {
    firkProcessError(NULL, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "LS must be non-NULL");
    return FIRKLS_ILL_INPUT;
  }
  firk_mem = (FIRKodeMem)firkode_mem;

  if (firk_mem->MallocDone == SUNFALSE)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MALLOC);
    return FIRKLS_ILL_INPUT;
  }

  /* test if solver is compatible with the LS interface */
  if ((LS->ops->gettype == NULL) || (LS->ops->solve == NULL))
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "LS object is missing a required operation");
    return FIRKLS_ILL_INPUT;
  }

  LSType = SUNLinSolGetType(LS);

  iterative   = (LSType != SUNLINEARSOLVER_DIRECT);
  matrixbased = ((LSType != SUNLINEARSOLVER_ITERATIVE) &&
                 (LSType != SUNLINEARSOLVER_MATRIX_EMBEDDED));

  /* test if vector is compatible with the LS interface */
  if ((firk_mem->tempv1->ops->nvconst == NULL) ||
      (firk_mem->tempv1->ops->nvwrmsnorm == NULL) ||
      (firk_mem->tempv1->ops->nvgetlength == NULL))
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSG_LS_BAD_NVECTOR);
    return FIRKLS_ILL_INPUT;
  }

  /* ensure that A is NULL when LS is matrix-embedded */
  if ((LSType == SUNLINEARSOLVER_MATRIX_EMBEDDED) && (A != NULL))
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Incompatible inputs: matrix-embedded LS requires NULL "
                     "matrix");
    return FIRKLS_ILL_INPUT;
  }

  /* check for compatible LS type, matrix and "atimes" support */
  if (iterative)
  {
    if (!matrixbased && (LSType != SUNLINEARSOLVER_MATRIX_EMBEDDED) &&
        (LS->ops->setatimes == NULL))
    {
      firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                       "Incompatible inputs: iterative LS must support ATimes "
                       "routine");
      return FIRKLS_ILL_INPUT;
    }

    if (matrixbased && (A == NULL))
    {
      firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                       "Incompatible inputs: matrix-iterative LS requires "
                       "non-NULL matrix");
      return FIRKLS_ILL_INPUT;
    }
  }
  else if (A == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Incompatible inputs: direct LS requires non-NULL matrix");
    return FIRKLS_ILL_INPUT;
  }

  /* free any existing system solver attached to FIRKODE */
  if (firk_mem->lfree) { firk_mem->lfree(firk_mem); }

  /* set the linear solver function fields in firk_mem */
  firk_mem->linit      = firkLsInitialize;
  firk_mem->lreinit    = firkLsReInitialize;
  firk_mem->lsetup     = firkLsSetup;
  firk_mem->lsolve_stk = firkLsSolveStacked;
  firk_mem->lsolve_blk = firkLsSolveBlock;
  firk_mem->lfree      = firkLsFree;

  /* allocate memory for FIRKLsMemRec */
  firkls_mem = (FIRKLsMem)malloc(sizeof(struct FIRKLsMemRec));
  if (firkls_mem == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSG_LS_MEM_FAIL);
    return FIRKLS_MEM_FAIL;
  }
  memset(firkls_mem, 0, sizeof(struct FIRKLsMemRec));

  firkls_mem->LS = LS;

  firkls_mem->iterative   = iterative;
  firkls_mem->matrixbased = matrixbased;

  /* the Newton correction of a matrix-based solver is scaled for a change
     of gamma since the setup, as in CVODE */
  firkls_mem->scalesol = matrixbased;

  /* set defaults for Jacobian-related fields */
  if (A != NULL)
  {
    firkls_mem->user_jac = SUNFALSE;
    firkls_mem->jac      = firkLsDQJac;
    firkls_mem->J_data   = firk_mem;
  }
  else
  {
    firkls_mem->user_jac = SUNFALSE;
    firkls_mem->jac      = NULL;
    firkls_mem->J_data   = NULL;
  }

  firkls_mem->jtimesDQ  = SUNTRUE;
  firkls_mem->jtsetup   = NULL;
  firkls_mem->jtimes    = firkLsDQJtimes;
  firkls_mem->jtimesRhs = firk_mem->f;
  firkls_mem->jt_data   = firk_mem;

  firkls_mem->user_linsys = SUNFALSE;
  firkls_mem->linsys      = firkLsLinSys;
  firkls_mem->A_data      = firk_mem;

  /* set defaults for preconditioner-related fields */
  firkls_mem->pset   = NULL;
  firkls_mem->psolve = NULL;
  firkls_mem->P_data = firk_mem->user_data;

  /* initialize counters */
  firkLsInitializeCounters(firkls_mem);

  /* set default values for the rest of the LS parameters */
  firkls_mem->jbad      = SUNTRUE;
  firkls_mem->eplifac   = FIRKLS_EPLIN;
  firkls_mem->sysfac    = ONE;
  firkls_mem->maxl      = 0; /* set from s in firkLsInitialize */
  firkls_mem->maxl_set  = SUNFALSE;
  firkls_mem->last_flag = FIRKLS_SUCCESS;

  /* if LS supports ATimes, attach the FIRKLs routine */
  if (LS->ops->setatimes)
  {
    retval = SUNLinSolSetATimes(LS, firk_mem, firkLsATimes);
    if (retval != SUN_SUCCESS)
    {
      firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__,
                       __FILE__, "Error in calling SUNLinSolSetATimes");
      free(firkls_mem);
      return FIRKLS_SUNLS_FAIL;
    }
  }

  /* if LS supports preconditioning, initialize pset/psol to NULL */
  if (LS->ops->setpreconditioner)
  {
    retval = SUNLinSolSetPreconditioner(LS, firk_mem, NULL, NULL);
    if (retval != SUN_SUCCESS)
    {
      firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__,
                       __FILE__, "Error in calling SUNLinSolSetPreconditioner");
      free(firkls_mem);
      return FIRKLS_SUNLS_FAIL;
    }
  }

  /* when using a SUNMatrix object, store pointer to A and initialize savedJ */
  if (A != NULL)
  {
    firkls_mem->A      = A;
    firkls_mem->savedJ = NULL; /* allocated in firkLsInitialize */
  }

  /* allocate memory for ytemp and x */
  firkls_mem->ytemp = N_VClone(firk_mem->tempv1);
  if (firkls_mem->ytemp == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSG_LS_MEM_FAIL);
    free(firkls_mem);
    return FIRKLS_MEM_FAIL;
  }

  firkls_mem->x = N_VClone(firk_mem->tempv1);
  if (firkls_mem->x == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSG_LS_MEM_FAIL);
    N_VDestroy(firkls_mem->ytemp);
    free(firkls_mem);
    return FIRKLS_MEM_FAIL;
  }

  /* norm conversion factor for the block solver */
  firkls_mem->nrmfac = SUNRsqrt((sunrealtype)N_VGetLength(firkls_mem->ytemp));

  /* stacked solver objects are created in firkLsInitialize */
  firkls_mem->SPFGMR  = NULL;
  firkls_mem->Jv      = NULL;
  firkls_mem->s_alloc = 0;

  /* attach linear solver memory to integrator memory */
  firk_mem->lmem = firkls_mem;

  /* defaults of the options only EXPERT mode may change */
  return firkLsSetExpertDefaults(firk_mem);
}

/*===============================================================
  Optional Set routines
  ===============================================================*/

int FIRKodeSetJacFn(void* firkode_mem, FIRKLsJacFn jac)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if ((jac != NULL) && (firkls_mem->A == NULL))
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Jacobian routine cannot be supplied for NULL SUNMatrix");
    return FIRKLS_ILL_INPUT;
  }

  if (jac != NULL)
  {
    firkls_mem->user_jac = SUNTRUE;
    firkls_mem->jac      = jac;
    firkls_mem->J_data   = firk_mem->user_data;
  }
  else
  {
    firkls_mem->user_jac = SUNFALSE;
    firkls_mem->jac      = firkLsDQJac;
    firkls_mem->J_data   = firk_mem;
  }

  /* reset the internal linear system function (if applicable) */
  if (!firkls_mem->user_linsys)
  {
    firkls_mem->linsys = firkLsLinSys;
    firkls_mem->A_data = firk_mem;
  }

  return FIRKLS_SUCCESS;
}

int FIRKodeSetJacEvalFrequency(void* firkode_mem, long int msbj)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMemExpert(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (msbj < 0)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "A negative evaluation frequency is illegal.");
    return FIRKLS_ILL_INPUT;
  }

  firkls_mem->msbj = (msbj == 0) ? FIRKLS_MSBJ : msbj;

  return FIRKLS_SUCCESS;
}

int FIRKodeSetDeltaGammaMaxBadJac(void* firkode_mem, sunrealtype dgmax_jbad)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMemExpert(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  firkls_mem->dgmax_jbad = (dgmax_jbad <= ZERO) ? FIRKLS_DGMAX : dgmax_jbad;

  return FIRKLS_SUCCESS;
}

int FIRKodeSetEpsLin(void* firkode_mem, sunrealtype eplifac)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (eplifac < ZERO)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "eplifac < 0 illegal.");
    return FIRKLS_ILL_INPUT;
  }

  firkls_mem->eplifac = (eplifac == ZERO) ? FIRKLS_EPLIN : eplifac;

  return FIRKLS_SUCCESS;
}

int FIRKodeSetLinearSolutionScaling(void* firkode_mem, sunbooleantype onoff)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  /* only a matrix-based solver solves with the matrix of the last setup */
  if (!firkls_mem->matrixbased)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Linear solution scaling requires a matrix-based linear "
                     "solver.");
    return FIRKLS_ILL_INPUT;
  }

  firkls_mem->scalesol = onoff;

  return FIRKLS_SUCCESS;
}

int FIRKodeSetLSNormFactor(void* firkode_mem, sunrealtype nrmfac)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (nrmfac > ZERO) { firkls_mem->nrmfac = nrmfac; }
  else if (nrmfac < ZERO)
  {
    N_VConst(ONE, firkls_mem->ytemp);
    firkls_mem->nrmfac =
      SUNRsqrt(N_VDotProd(firkls_mem->ytemp, firkls_mem->ytemp));
  }
  else
  {
    firkls_mem->nrmfac = SUNRsqrt((sunrealtype)N_VGetLength(firkls_mem->ytemp));
  }

  return FIRKLS_SUCCESS;
}

int FIRKodeSetPreconditioner(void* firkode_mem, FIRKLsPrecSetupFn pset,
                             FIRKLsPrecSolveFn psolve)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  SUNPSetupFn firkls_psetup;
  SUNPSolveFn firkls_psolve;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  /* a rejected call changes nothing */
  if (firkls_mem->LS->ops->setpreconditioner == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "SUNLinearSolver object does not support user-supplied "
                     "preconditioning");
    return FIRKLS_ILL_INPUT;
  }

  firkls_mem->pset   = pset;
  firkls_mem->psolve = psolve;

  firkls_psetup = (pset == NULL) ? NULL : firkLsPSetup;
  firkls_psolve = (psolve == NULL) ? NULL : firkLsPSolve;
  retval = SUNLinSolSetPreconditioner(firkls_mem->LS, firk_mem, firkls_psetup,
                                      firkls_psolve);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__, __FILE__,
                     "Error in calling SUNLinSolSetPreconditioner");
    return FIRKLS_SUNLS_FAIL;
  }

  return FIRKLS_SUCCESS;
}

int FIRKodeSetJacTimes(void* firkode_mem, FIRKLsJacTimesSetupFn jtsetup,
                       FIRKLsJacTimesVecFn jtimes)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (jtimes != NULL)
  {
    firkls_mem->jtimesDQ = SUNFALSE;
    firkls_mem->jtsetup  = jtsetup;
    firkls_mem->jtimes   = jtimes;
    firkls_mem->jt_data  = firk_mem->user_data;
  }
  else
  {
    firkls_mem->jtimesDQ  = SUNTRUE;
    firkls_mem->jtsetup   = NULL;
    firkls_mem->jtimes    = firkLsDQJtimes;
    firkls_mem->jtimesRhs = firk_mem->f;
    firkls_mem->jt_data   = firk_mem;
  }

  return FIRKLS_SUCCESS;
}

int FIRKodeSetJacTimesRhsFn(void* firkode_mem, FIRKRhsFn jtimesRhsFn)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (!firkls_mem->jtimesDQ)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Internal finite-difference Jacobian-vector product is "
                     "disabled.");
    return FIRKLS_ILL_INPUT;
  }

  firkls_mem->jtimesRhs = (jtimesRhsFn != NULL) ? jtimesRhsFn : firk_mem->f;

  return FIRKLS_SUCCESS;
}

int FIRKodeSetLinSysFn(void* firkode_mem, FIRKLsLinSysFn linsys)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if ((linsys != NULL) && (firkls_mem->A == NULL))
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Linear system setup routine cannot be supplied for NULL "
                     "SUNMatrix");
    return FIRKLS_ILL_INPUT;
  }

  if (linsys != NULL)
  {
    firkls_mem->user_linsys = SUNTRUE;
    firkls_mem->linsys      = linsys;
    firkls_mem->A_data      = firk_mem->user_data;

    /* no Jacobian is saved with a user linsys: a saved one would be stale
       for the stacked operator, which then uses jtimes or difference
       quotients */
    if (firkls_mem->savedJ != NULL)
    {
      SUNMatDestroy(firkls_mem->savedJ);
      firkls_mem->savedJ = NULL;
    }
  }
  else
  {
    firkls_mem->user_linsys = SUNFALSE;
    firkls_mem->linsys      = firkLsLinSys;
    firkls_mem->A_data      = firk_mem;

    /* the default linsys saves the Jacobian, from its next evaluation on */
    if (firkls_mem->A != NULL && firkls_mem->savedJ == NULL)
    {
      firkls_mem->savedJ = SUNMatClone(firkls_mem->A);
      if (firkls_mem->savedJ == NULL)
      {
        firkProcessError(firk_mem, FIRKLS_MEM_FAIL, __LINE__, __func__,
                         __FILE__, MSG_LS_MEM_FAIL);
        return FIRKLS_MEM_FAIL;
      }
      firkls_mem->jac_saved = SUNFALSE;
    }
  }

  return FIRKLS_SUCCESS;
}

int FIRKodeSetStageSolverMaxl(void* firkode_mem, int maxl)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMemExpert(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (maxl <= 0)
  {
    firkls_mem->maxl_set = SUNFALSE;
    firkls_mem->maxl     = 0;
  }
  else
  {
    firkls_mem->maxl_set = SUNTRUE;
    firkls_mem->maxl     = maxl;
  }

  /* the stacked solver is created again, with the new dimension, at the
     next stacked solve: it may be live now (H4) */
  firkls_mem->stk_stale = SUNTRUE;

  return FIRKLS_SUCCESS;
}

int FIRKodeSetStageSolverMaxRestarts(void* firkode_mem, int maxrs)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMemExpert(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  firkls_mem->maxrs = (maxrs < 0) ? FIRKLS_MAXRS : maxrs;
  if (firkls_mem->SPFGMR != NULL)
  {
    if (SUNLinSol_SPFGMRSetMaxRestarts(firkls_mem->SPFGMR, firkls_mem->maxrs) !=
        SUN_SUCCESS)
    {
      firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__, __FILE__,
                       "Error in calling SUNLinSol_SPFGMRSetMaxRestarts");
      return FIRKLS_SUNLS_FAIL;
    }
  }

  return FIRKLS_SUCCESS;
}

int FIRKodeSetStageSolverGSType(void* firkode_mem, int gstype)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMemExpert(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (gstype != SUN_MODIFIED_GS && gstype != SUN_CLASSICAL_GS)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "gstype must be SUN_MODIFIED_GS or SUN_CLASSICAL_GS.");
    return FIRKLS_ILL_INPUT;
  }

  firkls_mem->gstype = gstype;
  if (firkls_mem->SPFGMR != NULL)
  {
    if (SUNLinSol_SPFGMRSetGSType(firkls_mem->SPFGMR, gstype) != SUN_SUCCESS)
    {
      firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__,
                       __FILE__, "Error in calling SUNLinSol_SPFGMRSetGSType");
      return FIRKLS_SUNLS_FAIL;
    }
  }

  return FIRKLS_SUCCESS;
}

int FIRKodeSetStageEpsLin(void* firkode_mem, sunrealtype eplifac)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMemExpert(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (eplifac < ZERO)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "eplifac < 0 illegal.");
    return FIRKLS_ILL_INPUT;
  }

  firkls_mem->eps_stk = (eplifac == ZERO) ? FIRKLS_EPLIN : eplifac;

  return FIRKLS_SUCCESS;
}

/*---------------------------------------------------------------
  firkLsSetUserData

  Passes a new user_data pointer to the callbacks of the linear and
  mass solver interfaces that receive it (FIRKodeSetUserData, H5):
  the user's Jacobian, linsys, jtimes, preconditioner, mass function
  and mass preconditioner. The mass-times routine keeps the data given
  to FIRKodeSetMassTimes.
  ---------------------------------------------------------------*/
void firkLsSetUserData(FIRKodeMem firk_mem, void* user_data)
{
  FIRKLsMem firkls_mem;
  FIRKLsMassMem firkls_massmem;

  if (firk_mem->lmem != NULL)
  {
    firkls_mem = (FIRKLsMem)firk_mem->lmem;
    if (firkls_mem->user_jac) { firkls_mem->J_data = user_data; }
    if (firkls_mem->user_linsys) { firkls_mem->A_data = user_data; }
    if (!firkls_mem->jtimesDQ) { firkls_mem->jt_data = user_data; }
    firkls_mem->P_data = user_data;
  }
  if (firk_mem->mass_mem != NULL)
  {
    firkls_massmem         = (FIRKLsMassMem)firk_mem->mass_mem;
    firkls_massmem->M_data = user_data;
    firkls_massmem->P_data = user_data;
  }
}

/*---------------------------------------------------------------
  firkLsWriteParameters

  The effective parameters of the linear solver interface and of
  the mass matrix interface, for FIRKodeWriteParameters.
  ---------------------------------------------------------------*/
void firkLsWriteParameters(FIRKodeMem firk_mem, FILE* fp)
{
  FIRKLsMem ls;
  FIRKLsMassMem ms;
  int maxl;

  if (firk_mem->lmem == NULL) { fprintf(fp, "  Linear solver = none\n"); }
  else
  {
    ls = (FIRKLsMem)firk_mem->lmem;
    fprintf(fp, "  Linear solver: iterative = %i, matrix-based = %i\n",
            (int)ls->iterative, (int)ls->matrixbased);
    fprintf(fp, "  Linear solver: Jacobian = %s\n",
            ls->user_linsys   ? "user linear system"
            : (ls->A == NULL) ? "none (matrix-free)"
            : ls->user_jac    ? "user"
                              : "difference quotient");
    fprintf(fp, "  Linear solver: Jacobian-vector product = %s\n",
            ls->jtimesDQ ? "difference quotient" : "user");
    fprintf(fp, "  Linear solver: preconditioner setup = %s, solve = %s\n",
            (ls->pset != NULL) ? "user" : "none",
            (ls->psolve != NULL) ? "user" : "none");
    fprintf(fp, "  Linear solver: tolerance factor = " SUN_FORMAT_G "\n",
            ls->eplifac);
    fprintf(fp, "  Linear solver: norm factor = " SUN_FORMAT_G "\n", ls->nrmfac);
    fprintf(fp, "  Linear solver: solution scaling = %i\n", (int)ls->scalesol);
    fprintf(fp, "  Jacobian evaluation frequency = %li\n", ls->msbj);
    fprintf(fp, "  Bad Jacobian threshold = " SUN_FORMAT_G "\n", ls->dgmax_jbad);
    maxl = ls->maxl_set ? ls->maxl : SUNMAX(5, 3 * firk_mem->s);
    fprintf(fp,
            "  Stage solver: maxl = %i, restarts = %i, Gram-Schmidt = %s, "
            "tolerance factor = " SUN_FORMAT_G "\n",
            maxl, ls->maxrs,
            (ls->gstype == SUN_CLASSICAL_GS) ? "classical" : "modified",
            ls->eps_stk);
  }

  if (firk_mem->mass_mem != NULL)
  {
    ms = (FIRKLsMassMem)firk_mem->mass_mem;
    fprintf(fp, "  Mass solver: iterative = %i, matrix-based = %i\n",
            (int)ms->iterative, (int)(ms->M != NULL));
    fprintf(fp, "  Mass solver: products = %s\n",
            (ms->mtimes != NULL) ? "user" : "matrix");
    fprintf(fp, "  Mass solver: preconditioner setup = %s, solve = %s\n",
            (ms->pset != NULL) ? "user" : "none",
            (ms->psolve != NULL) ? "user" : "none");
    fprintf(fp, "  Mass solver: tolerance factor = " SUN_FORMAT_G "\n",
            ms->eplifac);
    fprintf(fp, "  Mass solver: norm factor = " SUN_FORMAT_G "\n", ms->nrmfac);
  }
}

/*===============================================================
  Optional Get routines
  ===============================================================*/

int FIRKodeGetJac(void* firkode_mem, SUNMatrix* J_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *J_ptr = firkls_mem->savedJ;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetJacTime(void* firkode_mem, sunrealtype* t_J_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *t_J_ptr = firkls_mem->tnlj;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetJacNumSteps(void* firkode_mem, long int* nst_J_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nst_J_ptr = firkls_mem->nstlj;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumJacEvals(void* firkode_mem, long int* njevals_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *njevals_ptr = firkls_mem->nje;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumLinRhsEvals(void* firkode_mem, long int* nfevalsLS_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nfevalsLS_ptr = firkls_mem->nfeDQ;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumPrecEvals(void* firkode_mem, long int* npevals_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *npevals_ptr = firkls_mem->npe;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumPrecSolves(void* firkode_mem, long int* npsolves_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *npsolves_ptr = firkls_mem->nps;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumLinIters(void* firkode_mem, long int* nliters_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nliters_ptr = firkls_mem->nli;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumLinConvFails(void* firkode_mem, long int* nlcfails_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nlcfails_ptr = firkls_mem->ncfl;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumJTSetupEvals(void* firkode_mem, long int* njtsetups_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *njtsetups_ptr = firkls_mem->njtsetup;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumJtimesEvals(void* firkode_mem, long int* njvevals_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *njvevals_ptr = firkls_mem->njtimes;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumBlockSolves(void* firkode_mem, long int* nbsolves_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nbsolves_ptr = firkls_mem->nblocksolves;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumStageLinIters(void* firkode_mem, long int* nsliters_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nsliters_ptr = firkls_mem->nsli;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumStageLinConvFails(void* firkode_mem, long int* nslcfails_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nslcfails_ptr = firkls_mem->nscf;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetLastLinFlag(void* firkode_mem, long int* flag_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *flag_ptr = firkls_mem->last_flag;
  return FIRKLS_SUCCESS;
}

char* FIRKodeGetLinReturnFlagName(long int flag)
{
  const char* str;
  char* name;

  switch (flag)
  {
  case FIRKLS_SUCCESS: str = "FIRKLS_SUCCESS"; break;
  case FIRKLS_MEM_NULL: str = "FIRKLS_MEM_NULL"; break;
  case FIRKLS_LMEM_NULL: str = "FIRKLS_LMEM_NULL"; break;
  case FIRKLS_ILL_INPUT: str = "FIRKLS_ILL_INPUT"; break;
  case FIRKLS_MEM_FAIL: str = "FIRKLS_MEM_FAIL"; break;
  case FIRKLS_MASSMEM_NULL: str = "FIRKLS_MASSMEM_NULL"; break;
  case FIRKLS_JACFUNC_UNRECVR: str = "FIRKLS_JACFUNC_UNRECVR"; break;
  case FIRKLS_JACFUNC_RECVR: str = "FIRKLS_JACFUNC_RECVR"; break;
  case FIRKLS_MASSFUNC_UNRECVR: str = "FIRKLS_MASSFUNC_UNRECVR"; break;
  case FIRKLS_MASSFUNC_RECVR: str = "FIRKLS_MASSFUNC_RECVR"; break;
  case FIRKLS_SUNMAT_FAIL: str = "FIRKLS_SUNMAT_FAIL"; break;
  case FIRKLS_SUNLS_FAIL: str = "FIRKLS_SUNLS_FAIL"; break;
  default:
    /* a negative integrator code that a FIRKLS function can return: the
       setters of EXPERT options return FIRK_ILL_INPUT in BEGINNER mode; a
       positive flag is an LU factorization failure */
    if (flag < 0) { return FIRKodeGetReturnFlagName(flag); }
    str = "NONE";
  }

  name = (char*)malloc(strlen(str) + 1);
  if (name == NULL) { return NULL; }
  memcpy(name, str, strlen(str) + 1);
  return name;
}

/*=================================================================
  FIRKLS private functions
  =================================================================*/

/*-----------------------------------------------------------------
  firkLsATimes

  Generates the block matrix-vector product z = (M - gamma J) v for
  the user's iterative block solver, gamma = gamma_op of the solve.
  J v is obtained from the jtimes routine. During a block solve with
  a mass matrix the product is multiplied by sysfac (see
  firkLsSolveBlock).
  -----------------------------------------------------------------*/
int firkLsATimes(void* firkode_mem, N_Vector v, N_Vector z)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  retval = firkls_mem->jtimes(v, z, firk_mem->tn, firkls_mem->ycur,
                              firkls_mem->fcur, firkls_mem->jt_data,
                              firkls_mem->ytemp);
  firkls_mem->njtimes++;
  if (retval != 0) { return retval; }

  if (firk_mem->mass_set)
  {
    retval = firk_mem->mmult(firk_mem, v, firkls_mem->ytemp);
    if (retval != 0) { return retval; }
    N_VLinearSum(firkls_mem->sysfac, firkls_mem->ytemp,
                 -firkls_mem->sysfac * firkls_mem->gamma_op, z, z);
  }
  else { N_VLinearSum(ONE, v, -firkls_mem->gamma_op, z, z); }

  return 0;
}

/*---------------------------------------------------------------
  firkLsPSetup

  Interfaces the user's psetup routine with the block solver: it
  receives jok = !jbad and reports into firk_mem->jcur, which
  firkLsSetupMatrix records (H29). Every call counts in npe.
  ---------------------------------------------------------------*/
int firkLsPSetup(void* firkode_mem)
{
  int retval;
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  retval = firkls_mem->pset(firk_mem->tn, firkls_mem->ycur, firkls_mem->fcur,
                            !(firkls_mem->jbad), &firk_mem->jcur,
                            firk_mem->gamma, firkls_mem->P_data);
  firkls_mem->npe++;
  return retval;
}

/*-----------------------------------------------------------------
  firkLsPSolve

  Interfaces the user's psolve routine with the block solver, with
  gamma = gamma_op of the solve. During a block solve with a mass
  matrix the preconditioner is that of the system multiplied by
  sysfac, so z is divided by sysfac (see firkLsSolveBlock).
  -----------------------------------------------------------------*/
int firkLsPSolve(void* firkode_mem, N_Vector r, N_Vector z, sunrealtype tol,
                 int lr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  retval = firkls_mem->psolve(firk_mem->tn, firkls_mem->ycur, firkls_mem->fcur,
                              r, z, firkls_mem->gamma_op, tol, lr,
                              firkls_mem->P_data);
  firkls_mem->nps++;
  if (retval == 0 && firkls_mem->sysfac != ONE)
  {
    N_VScale(ONE / firkls_mem->sysfac, z, z);
  }
  return retval;
}

/*-----------------------------------------------------------------
  firkLsStackedATimes

  Generates the stacked stage-system matrix-vector product
    z_i = M v_i - h sum_j a_ij J v_j.
  J v_j is computed with SUNMatMatvec on the saved Jacobian when its
  matrix provides it, and with the jtimes routine (the user's or the
  difference quotient) otherwise (design Listing 8.4).
  -----------------------------------------------------------------*/
int firkLsStackedATimes(void* firkode_mem, N_Vector v_stk, N_Vector z_stk)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  const struct FIRKodeTableMem* T;
  N_Vector v_i, z_i, mv;
  int i, j, s, retval;
  sunbooleantype use_matvec;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  T = firk_mem->tbl;
  s = firk_mem->s;

  /* J v from the saved Jacobian if its matrix provides the product, then
     from jtimes (the user's or the difference quotient) */
  use_matvec = (firkls_mem->savedJ != NULL) &&
               (firkls_mem->savedJ->ops->matvec != NULL);

  /* Jacobian-vector products for each stage */
  for (j = 0; j < s; j++)
  {
    v_i = firkStage(firk_mem, v_stk, j);
    if (use_matvec)
    {
      retval = SUNMatMatvec(firkls_mem->savedJ, v_i, firkls_mem->Jv[j]);
      if (retval != 0) { return retval; }
    }
    else
    {
      retval = firkls_mem->jtimes(v_i, firkls_mem->Jv[j], firk_mem->tn,
                                  firkls_mem->ycur, firkls_mem->fcur,
                                  firkls_mem->jt_data, firkls_mem->ytemp);
      firkls_mem->njtimes++;
      if (retval != 0) { return retval; }
    }
  }

  /* z_i = M v_i - h sum_j a_ij J v_j */
  for (i = 0; i < s; i++)
  {
    v_i = firkStage(firk_mem, v_stk, i);
    z_i = firkStage(firk_mem, z_stk, i);

    if (firk_mem->mass_set)
    {
      retval = firk_mem->mmult(firk_mem, v_i, z_i);
      if (retval != 0) { return retval; }
      mv = z_i;
    }
    else { mv = v_i; }

    firk_mem->cvals[0] = ONE;
    firk_mem->Xvecs[0] = mv;
    for (j = 0; j < s; j++)
    {
      firk_mem->cvals[j + 1] = -firk_mem->h * T->A[i * s + j];
      firk_mem->Xvecs[j + 1] = firkls_mem->Jv[j];
    }
    retval = N_VLinearCombination(s + 1, firk_mem->cvals, firk_mem->Xvecs, z_i);
    if (retval != 0) { return -1; }
  }

  return 0;
}

/*-----------------------------------------------------------------
  firkLsStackedPSolve

  Block-diagonal preconditioner for the stacked stage system (design
  Listing 8.4). A direct block solver solves each stage block exactly
  with M - gamma J; a matrix-iterative one applies one inner solve to
  the tolerance eps_L ||r_i|| (R6). A matrix-free block solver
  contributes only the user's preconditioner, applied on the sides of
  the user's pretype: the user's Krylov iteration never runs here
  (D2, H37); without a preconditioner a block is the identity.
  Recoverable failures of a block are accepted and counted in
  nblockfails; an unrecoverable one returns -1, which SPFGMR reports
  as SUNLS_PSOLVE_FAIL_UNREC.
  -----------------------------------------------------------------*/
int firkLsStackedPSolve(void* firkode_mem, N_Vector r_stk, N_Vector z_stk,
                        SUNDIALS_MAYBE_UNUSED sunrealtype tol,
                        SUNDIALS_MAYBE_UNUSED int lr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  N_Vector r_i, z_i;
  int i, pretype, retval;
  long int nbf0;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  pretype = (firkls_mem->psolve != NULL) ? firkLsUserPretype(firkls_mem->LS)
                                         : SUN_PREC_NONE;
  nbf0    = firkls_mem->nblockfails;

  for (i = 0; i < firk_mem->s; i++)
  {
    r_i = firkStage(firk_mem, r_stk, i);
    z_i = firkStage(firk_mem, z_stk, i);

    if (firkls_mem->matrixbased)
    {
      N_VScale(ONE, r_i, z_i);
      retval = firkLsSolveBlock(firk_mem, z_i, z_i, FIRKLS_BLOCK_PRECOND, 0);
      if (retval != 0) { return -1; }
    }
    else if (pretype != SUN_PREC_NONE)
    {
      retval = firkLsUserPSolveBlock(firk_mem, firkls_mem, r_i, z_i, pretype);
      if (retval < 0) { return -1; }
      if (retval > 0) { firkls_mem->nblockfails++; }
      firkls_mem->nblocksolves++;
    }
    else { N_VScale(ONE, r_i, z_i); }
  }

  SUNLogExtraDebug(FIRK_LOGGER, "block-solve",
                   "blocks = %i, accepted-failures = %li", firk_mem->s,
                   firkls_mem->nblockfails - nbf0);
  (void)nbf0; /* used only by the log record */
  return 0;
}

/*-----------------------------------------------------------------
  firkLsUserPretype

  The sides on which the user's Krylov solver applies the
  preconditioner. The SUNDIALS Krylov solvers are recognized by their
  solve routine, so that a wrapper of one (which may report the
  wrapped solver's ID) is never mistaken for it; PCG applies its
  preconditioner as SUN_PREC_LEFT whatever its pretype. Any other
  solver is assumed to precondition on the left (R19).
  -----------------------------------------------------------------*/
static int firkLsUserPretype(SUNLinearSolver LS)
{
  if (LS->ops->solve == SUNLinSolSolve_SPGMR)
  {
    return ((SUNLinearSolverContent_SPGMR)LS->content)->pretype;
  }
  if (LS->ops->solve == SUNLinSolSolve_SPFGMR)
  {
    return ((SUNLinearSolverContent_SPFGMR)LS->content)->pretype;
  }
  if (LS->ops->solve == SUNLinSolSolve_SPBCGS)
  {
    return ((SUNLinearSolverContent_SPBCGS)LS->content)->pretype;
  }
  if (LS->ops->solve == SUNLinSolSolve_SPTFQMR)
  {
    return ((SUNLinearSolverContent_SPTFQMR)LS->content)->pretype;
  }
  if (LS->ops->solve == SUNLinSolSolve_PCG)
  {
    return (((SUNLinearSolverContent_PCG)LS->content)->pretype == SUN_PREC_NONE)
             ? SUN_PREC_NONE
             : SUN_PREC_LEFT;
  }
  return SUN_PREC_LEFT;
}

/*-----------------------------------------------------------------
  firkLsUserPSolveBlock

  Applies the user's preconditioner to one stage block, z = P^{-1} r,
  with the tolerance eps_L ||r||_WRMS nrmfac; SUN_PREC_BOTH applies
  the left and then the right factor. Returns 0, +1 after a
  recoverable failure of a call, or -1.
  -----------------------------------------------------------------*/
static int firkLsUserPSolveBlock(FIRKodeMem firk_mem, FIRKLsMem firkls_mem,
                                 N_Vector r, N_Vector z, int pretype)
{
  sunrealtype delta;
  int retval, failed = 0;

  delta = firkls_mem->eplifac * N_VWrmsNorm(r, firk_mem->ewt) *
          firkls_mem->nrmfac;

  if (pretype == SUN_PREC_BOTH)
  {
    retval = firkls_mem->psolve(firk_mem->tn, firkls_mem->ycur, firkls_mem->fcur,
                                r, firkls_mem->ytemp, firk_mem->gamma, delta,
                                SUN_PREC_LEFT, firkls_mem->P_data);
    firkls_mem->nps++;
    if (retval < 0) { return -1; }
    failed  = (retval > 0);
    r       = firkls_mem->ytemp;
    pretype = SUN_PREC_RIGHT;
  }

  retval = firkls_mem->psolve(firk_mem->tn, firkls_mem->ycur, firkls_mem->fcur,
                              r, z, firk_mem->gamma, delta, pretype,
                              firkls_mem->P_data);
  firkls_mem->nps++;
  if (retval < 0) { return -1; }
  return (failed || retval > 0) ? 1 : 0;
}

/*-----------------------------------------------------------------
  firkLsDQJac

  Wrapper for the dense and band difference quotient Jacobian
  approximation routines.
  ---------------------------------------------------------------*/
int firkLsDQJac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac,
                void* firkode_mem, N_Vector tmp1, N_Vector tmp2,
                SUNDIALS_MAYBE_UNUSED N_Vector tmp3)
{
  FIRKodeMem firk_mem;
  int retval;

  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRKLS_MEM_NULL, __LINE__, __func__, __FILE__,
                     MSG_LS_FIRKMEM_NULL);
    return FIRKLS_MEM_NULL;
  }
  firk_mem = (FIRKodeMem)firkode_mem;

  if (Jac == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_LMEM_NULL, __LINE__, __func__, __FILE__,
                     "SUNMatrix is NULL");
    return FIRKLS_LMEM_NULL;
  }

  if (y->ops->nvcloneempty == NULL || y->ops->nvwrmsnorm == NULL ||
      y->ops->nvlinearsum == NULL || y->ops->nvdestroy == NULL ||
      y->ops->nvscale == NULL || y->ops->nvgetarraypointer == NULL ||
      y->ops->nvsetarraypointer == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSG_LS_BAD_NVECTOR);
    return FIRKLS_ILL_INPUT;
  }

  if (SUNMatGetID(Jac) == SUNMATRIX_DENSE)
  {
    retval = firkLsDenseDQJac(t, y, fy, Jac, firk_mem, tmp1);
  }
  else if (SUNMatGetID(Jac) == SUNMATRIX_BAND)
  {
    retval = firkLsBandDQJac(t, y, fy, Jac, firk_mem, tmp1, tmp2);
  }
  else
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "unrecognized matrix type for firkLsDQJac");
    retval = FIRKLS_ILL_INPUT;
  }
  return retval;
}

/*-----------------------------------------------------------------
  firkLsDenseDQJac

  Dense difference quotient approximation to the Jacobian of f. The
  columns perturb y in place; each component is restored after its
  evaluation, whatever its outcome (H2).
  -----------------------------------------------------------------*/
int firkLsDenseDQJac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac,
                     FIRKodeMem firk_mem, N_Vector tmp1)
{
  sunrealtype fnorm, minInc, inc, inc_inv, yjsaved, srur;
  sunrealtype *y_data, *ewt_data;
  N_Vector ftemp, jthCol;
  sunindextype j, N;
  FIRKLsMem firkls_mem;
  int retval = 0;

  firkls_mem = (FIRKLsMem)firk_mem->lmem;

  N = SUNDenseMatrix_Columns(Jac);

  ftemp = tmp1;

  jthCol = N_VCloneEmpty(tmp1);
  if (jthCol == NULL) { return FIRKLS_MEM_FAIL; }

  ewt_data = N_VGetArrayPointer(firk_mem->ewt);
  y_data   = N_VGetArrayPointer(y);

  srur   = SUNRsqrt(firk_mem->uround);
  fnorm  = N_VWrmsNorm(fy, firk_mem->ewt);
  minInc = (fnorm != ZERO) ? (MIN_INC_MULT * SUNRabs(firk_mem->h) *
                              firk_mem->uround * (sunrealtype)N * fnorm)
                           : ONE;

  for (j = 0; j < N; j++)
  {
    N_VSetArrayPointer(SUNDenseMatrix_Column(Jac, j), jthCol);

    yjsaved = y_data[j];
    inc     = SUNMAX(srur * SUNRabs(yjsaved), minInc / ewt_data[j]);
    y_data[j] += inc;

    retval = firk_mem->f(t, y, ftemp, firk_mem->user_data);
    firkls_mem->nfeDQ++;

    /* restore the component, also after a failure: y is the solution (H2) */
    y_data[j] = yjsaved;
    if (retval != 0) { break; }

    inc_inv = ONE / inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fy, jthCol);
  }

  N_VSetArrayPointer(NULL, jthCol);
  N_VDestroy(jthCol);

  return retval;
}

/*-----------------------------------------------------------------
  firkLsBandDQJac

  Banded difference quotient approximation to the Jacobian of f. The
  groups perturb a copy of y, which is never modified (H2).
  -----------------------------------------------------------------*/
int firkLsBandDQJac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac,
                    FIRKodeMem firk_mem, N_Vector tmp1, N_Vector tmp2)
{
  N_Vector ftemp, ytemp;
  sunrealtype fnorm, minInc, inc, inc_inv, srur;
  sunrealtype *col_j, *ewt_data, *fy_data, *ftemp_data;
  sunrealtype *y_data, *ytemp_data;
  sunindextype group, i, j, width, ngroups, i1, i2;
  sunindextype N, mupper, mlower;
  FIRKLsMem firkls_mem;
  int retval = 0;

  firkls_mem = (FIRKLsMem)firk_mem->lmem;

  N      = SUNBandMatrix_Columns(Jac);
  mupper = SUNBandMatrix_UpperBandwidth(Jac);
  mlower = SUNBandMatrix_LowerBandwidth(Jac);

  ftemp = tmp1;
  ytemp = tmp2;

  ewt_data   = N_VGetArrayPointer(firk_mem->ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);

  N_VScale(ONE, y, ytemp);

  srur   = SUNRsqrt(firk_mem->uround);
  fnorm  = N_VWrmsNorm(fy, firk_mem->ewt);
  minInc = (fnorm != ZERO) ? (MIN_INC_MULT * SUNRabs(firk_mem->h) *
                              firk_mem->uround * (sunrealtype)N * fnorm)
                           : ONE;

  width   = mlower + mupper + 1;
  ngroups = SUNMIN(width, N);

  for (group = 1; group <= ngroups; group++)
  {
    for (j = group - 1; j < N; j += width)
    {
      inc = SUNMAX(srur * SUNRabs(y_data[j]), minInc / ewt_data[j]);
      ytemp_data[j] += inc;
    }

    retval = firk_mem->f(t, ytemp, ftemp, firk_mem->user_data);
    firkls_mem->nfeDQ++;
    if (retval != 0) { break; }

    for (j = group - 1; j < N; j += width)
    {
      ytemp_data[j] = y_data[j];
      col_j         = SUNBandMatrix_Column(Jac, j);
      inc           = SUNMAX(srur * SUNRabs(y_data[j]), minInc / ewt_data[j]);
      inc_inv       = ONE / inc;
      i1            = SUNMAX(0, j - mupper);
      i2            = SUNMIN(j + mlower, N - 1);
      for (i = i1; i <= i2; i++)
      {
        SM_COLUMN_ELEMENT_B(col_j, i, j) = inc_inv * (ftemp_data[i] - fy_data[i]);
      }
    }
  }

  return retval;
}

/*-----------------------------------------------------------------
  firkLsDQJtimes

  Difference quotient approximation to J v:
    Jv = [f(y + v sig) - f(y)] / sig,  sig = 1 / ||v||_WRMS,
  and J v = 0 for v = 0.
  -----------------------------------------------------------------*/
int firkLsDQJtimes(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y,
                   N_Vector fy, void* firkode_mem, N_Vector work)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  sunrealtype sig, siginv;
  int iter, retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  /* J 0 = 0 without a division (H7) */
  sig = N_VWrmsNorm(v, firk_mem->ewt);
  if (sig == ZERO)
  {
    N_VConst(ZERO, Jv);
    return 0;
  }
  sig = ONE / sig;

  for (iter = 0; iter < MAX_DQITERS; iter++)
  {
    N_VLinearSum(sig, v, ONE, y, work);

    retval = firkls_mem->jtimesRhs(t, work, Jv, firk_mem->user_data);
    firkls_mem->nfeDQ++;
    if (retval == 0) { break; }
    if (retval < 0) { return -1; }

    sig *= PT25;
  }

  if (retval > 0) { return +1; }

  siginv = ONE / sig;
  N_VLinearSum(siginv, Jv, -siginv, fy, Jv);

  return 0;
}

/*-----------------------------------------------------------------
  firkLsLinSys

  Sets up the block linear system A = M - gamma J.
  -----------------------------------------------------------------*/
static int firkLsLinSys(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix A,
                        SUNMatrix M, sunbooleantype jok, sunbooleantype* jcur,
                        sunrealtype gamma, void* firkode_mem, N_Vector vtemp1,
                        N_Vector vtemp2, N_Vector vtemp3)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (jok && firkls_mem->jac_saved)
  {
    /* use saved copy of J */
    *jcur = SUNFALSE;

    retval = SUNMatCopy(firkls_mem->savedJ, A);
    if (retval)
    {
      firkProcessError(firk_mem, FIRKLS_SUNMAT_FAIL, __LINE__, __func__,
                       __FILE__, MSG_LS_SUNMAT_FAILED);
      firkls_mem->last_flag = FIRKLS_SUNMAT_FAIL;
      return firkls_mem->last_flag;
    }
  }
  else
  {
    /* call jac() routine to update J */
    *jcur = SUNTRUE;

    if (SUNLinSolGetType(firkls_mem->LS) == SUNLINEARSOLVER_DIRECT)
    {
      retval = SUNMatZero(A);
      if (retval)
      {
        firkProcessError(firk_mem, FIRKLS_SUNMAT_FAIL, __LINE__, __func__,
                         __FILE__, MSG_LS_SUNMAT_FAILED);
        firkls_mem->last_flag = FIRKLS_SUNMAT_FAIL;
        return firkls_mem->last_flag;
      }
    }

    retval = firkls_mem->jac(t, y, fy, A, firkls_mem->J_data, vtemp1, vtemp2,
                             vtemp3);
    if (retval < 0)
    {
      firkProcessError(firk_mem, FIRKLS_JACFUNC_UNRECVR, __LINE__, __func__,
                       __FILE__, MSG_LS_JACFUNC_FAILED);
      firkls_mem->last_flag = FIRKLS_JACFUNC_UNRECVR;
      return -1;
    }
    if (retval > 0)
    {
      firkls_mem->last_flag = FIRKLS_JACFUNC_RECVR;
      return 1;
    }

    retval = SUNMatCopy(A, firkls_mem->savedJ);
    if (retval)
    {
      firkProcessError(firk_mem, FIRKLS_SUNMAT_FAIL, __LINE__, __func__,
                       __FILE__, MSG_LS_SUNMAT_FAILED);
      firkls_mem->last_flag = FIRKLS_SUNMAT_FAIL;
      return firkls_mem->last_flag;
    }
    firkls_mem->jac_saved = SUNTRUE;

    /* the products of the stacked operator use the saved Jacobian: set them
       up once per Jacobian evaluation */
    if (firkls_mem->savedJ->ops->matvecsetup)
    {
      retval = SUNMatMatvecSetup(firkls_mem->savedJ);
      if (retval)
      {
        firkProcessError(firk_mem, FIRKLS_SUNMAT_FAIL, __LINE__, __func__,
                         __FILE__, MSG_LS_SUNMAT_FAILED);
        firkls_mem->last_flag = FIRKLS_SUNMAT_FAIL;
        return firkls_mem->last_flag;
      }
    }
  }

  /* A = M - gamma J */
  if (M == NULL) { retval = SUNMatScaleAddI(-gamma, A); }
  else { retval = SUNMatScaleAdd(-gamma, A, M); }
  if (retval)
  {
    firkProcessError(firk_mem, FIRKLS_SUNMAT_FAIL, __LINE__, __func__, __FILE__,
                     MSG_LS_SUNMAT_FAILED);
    firkls_mem->last_flag = FIRKLS_SUNMAT_FAIL;
    return firkls_mem->last_flag;
  }

  return FIRKLS_SUCCESS;
}

/*-----------------------------------------------------------------
  firkLsSetupStacked / firkLsFreeStacked

  Create/free the internal flexible GMRES solver and work vectors
  for the stacked stage system (only used for s >= 2).
  -----------------------------------------------------------------*/
static int firkLsSetupStacked(FIRKodeMem firk_mem, FIRKLsMem firkls_mem)
{
  int s, maxl, retval;

  s = firk_mem->s;

  if (firkls_mem->SPFGMR != NULL && firkls_mem->s_alloc == s &&
      !firkls_mem->stk_stale)
  {
    return FIRKLS_SUCCESS;
  }
  firkLsFreeStacked(firkls_mem);
  firkls_mem->stk_stale = SUNFALSE;

  if (s < 2)
  {
    firkls_mem->s_alloc = s;
    return FIRKLS_SUCCESS;
  }

  maxl             = firkls_mem->maxl_set ? firkls_mem->maxl : SUNMAX(5, 3 * s);
  firkls_mem->maxl = maxl;

  firkls_mem->SPFGMR = SUNLinSol_SPFGMR(firk_mem->g_stk, SUN_PREC_RIGHT, maxl,
                                        firk_mem->sunctx);
  if (firkls_mem->SPFGMR == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSG_LS_MEM_FAIL);
    return FIRKLS_MEM_FAIL;
  }

  retval = SUNLinSol_SPFGMRSetMaxRestarts(firkls_mem->SPFGMR, firkls_mem->maxrs);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__, __FILE__,
                     "Error in calling SUNLinSol_SPFGMRSetMaxRestarts");
    return FIRKLS_SUNLS_FAIL;
  }

  retval = SUNLinSol_SPFGMRSetGSType(firkls_mem->SPFGMR, firkls_mem->gstype);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__, __FILE__,
                     "Error in calling SUNLinSol_SPFGMRSetGSType");
    return FIRKLS_SUNLS_FAIL;
  }

  retval = SUNLinSolSetATimes(firkls_mem->SPFGMR, firk_mem, firkLsStackedATimes);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__, __FILE__,
                     "Error in calling SUNLinSolSetATimes");
    return FIRKLS_SUNLS_FAIL;
  }

  retval = SUNLinSolSetPreconditioner(firkls_mem->SPFGMR, firk_mem, NULL,
                                      firkLsStackedPSolve);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__, __FILE__,
                     "Error in calling SUNLinSolSetPreconditioner");
    return FIRKLS_SUNLS_FAIL;
  }

  retval = SUNLinSolInitialize(firkls_mem->SPFGMR);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__, __FILE__,
                     "Error in calling SUNLinSolInitialize");
    return FIRKLS_SUNLS_FAIL;
  }

  firkls_mem->Jv = firkCloneVectorArray(s, firkls_mem->ytemp);
  if (firkls_mem->Jv == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSG_LS_MEM_FAIL);
    return FIRKLS_MEM_FAIL;
  }

  firkls_mem->s_alloc = s;

  return FIRKLS_SUCCESS;
}

static void firkLsFreeStacked(FIRKLsMem firkls_mem)
{
  if (firkls_mem->SPFGMR != NULL)
  {
    SUNLinSolFree(firkls_mem->SPFGMR);
    firkls_mem->SPFGMR = NULL;
  }
  if (firkls_mem->Jv != NULL)
  {
    N_VDestroyVectorArray(firkls_mem->Jv, firkls_mem->s_alloc);
    firkls_mem->Jv = NULL;
  }
  firkls_mem->s_alloc = 0;
}

/*-----------------------------------------------------------------
  firkLsInitialize

  Performs the remaining initializations of the linear solver
  interface; called from the initial setup.
  -----------------------------------------------------------------*/
int firkLsInitialize(FIRKodeMem firk_mem)
{
  FIRKLsMem firkls_mem;
  FIRKLsMassMem firkls_massmem;
  int retval;

  if (firk_mem->lmem == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_LMEM_NULL, __LINE__, __func__, __FILE__,
                     MSG_LS_LMEM_NULL);
    return FIRKLS_LMEM_NULL;
  }
  firkls_mem = (FIRKLsMem)firk_mem->lmem;

  /* the Jacobian and mass linear solvers must both be matrix-based or both
     be matrix-free, since the block system M - gamma J is assembled from the
     two SUNMatrix objects, and the matrices must share a matrix type */
  if (firk_mem->mass_set && firk_mem->mass_mem != NULL)
  {
    firkls_massmem = (FIRKLsMassMem)firk_mem->mass_mem;
    if ((firkls_mem->A == NULL) != (firkls_massmem->M == NULL))
    {
      firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                       "Cannot combine NULL and non-NULL system and mass "
                       "matrices");
      firkls_mem->last_flag = FIRKLS_ILL_INPUT;
      return FIRKLS_ILL_INPUT;
    }
  }
  if (firk_mem->mass_set && firk_mem->mass_mem != NULL && firkls_mem->A != NULL)
  {
    firkls_massmem = (FIRKLsMassMem)firk_mem->mass_mem;
    if (firkls_massmem->M != NULL &&
        (firkls_massmem->M->ops->getid == NULL ||
         firkls_mem->A->ops->getid == NULL ||
         SUNMatGetID(firkls_massmem->M) != SUNMatGetID(firkls_mem->A)))
    {
      firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                       MSG_LS_MATRIX_TYPES);
      firkls_mem->last_flag = FIRKLS_ILL_INPUT;
      return FIRKLS_ILL_INPUT;
    }
  }

  if (firkls_mem->A != NULL)
  {
    /* matrix-based case */
    if (firkls_mem->user_linsys) { firkls_mem->A_data = firk_mem->user_data; }
    else
    {
      firkls_mem->linsys = firkLsLinSys;
      firkls_mem->A_data = firk_mem;

      if (!firkls_mem->user_jac)
      {
        retval = 0;
        if (firkls_mem->A->ops->getid)
        {
          if ((SUNMatGetID(firkls_mem->A) == SUNMATRIX_DENSE) ||
              (SUNMatGetID(firkls_mem->A) == SUNMATRIX_BAND))
          {
            firkls_mem->jac    = firkLsDQJac;
            firkls_mem->J_data = firk_mem;
          }
          else { retval++; }
        }
        else { retval++; }
        if (retval)
        {
          firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__,
                           __FILE__,
                           "No Jacobian constructor available for SUNMatrix "
                           "type");
          firkls_mem->last_flag = FIRKLS_ILL_INPUT;
          return FIRKLS_ILL_INPUT;
        }
      }
      else { firkls_mem->J_data = firk_mem->user_data; }

      /* the stacked operator takes J v from the saved Jacobian when its
         matrix provides SUNMatMatvec, else from a user jtimes, else from
         difference quotients, with this warning (R7) */
      if (firk_mem->s > 1 && firkls_mem->jtimesDQ &&
          firkls_mem->A->ops->matvec == NULL)
      {
        firkProcessError(firk_mem, FIRK_WARNING, __LINE__, __func__, __FILE__,
                         "The SUNMatrix does not support SUNMatMatvec and no "
                         "Jacobian-vector product routine was supplied: the "
                         "stage system uses difference-quotient products");
      }

      if (firkls_mem->savedJ == NULL)
      {
        firkls_mem->savedJ = SUNMatClone(firkls_mem->A);
        if (firkls_mem->savedJ == NULL)
        {
          firkProcessError(firk_mem, FIRKLS_MEM_FAIL, __LINE__, __func__,
                           __FILE__, MSG_LS_MEM_FAIL);
          firkls_mem->last_flag = FIRKLS_MEM_FAIL;
          return FIRKLS_MEM_FAIL;
        }
        firkls_mem->jac_saved = SUNFALSE;
      }
    }
  }
  else
  {
    /* matrix-free case */
    firkls_mem->user_jac = SUNFALSE;
    firkls_mem->jac      = NULL;
    firkls_mem->J_data   = NULL;

    firkls_mem->user_linsys = SUNFALSE;
    firkls_mem->linsys      = NULL;
    firkls_mem->A_data      = NULL;
  }

  /* the counters are zeroed when the solver is attached and by ReInit, not
     here: Reset keeps them (H3) */

  /* Jacobian-vector product related fields */
  if (firkls_mem->jtimesDQ)
  {
    firkls_mem->jtsetup = NULL;
    firkls_mem->jtimes  = firkLsDQJtimes;
    firkls_mem->jt_data = firk_mem;
  }
  else { firkls_mem->jt_data = firk_mem->user_data; }

  /* the linearization point defaults */
  firkls_mem->ycur = firk_mem->yn;
  firkls_mem->fcur = firk_mem->fn;

  /* if A is NULL and psetup is not present, lsetup is not needed */
  firk_mem->lsetup = firkLsSetup;
  if ((firkls_mem->A == NULL) && (firkls_mem->pset == NULL))
  {
    firk_mem->lsetup = NULL;
  }

  /* matrix-embedded linear solvers do not use lsetup */
  if (SUNLinSolGetType(firkls_mem->LS) == SUNLINEARSOLVER_MATRIX_EMBEDDED)
  {
    firk_mem->lsetup = NULL;
  }

  /* stacked solver for the coupled stage system */
  retval = firkLsSetupStacked(firk_mem, firkls_mem);
  if (retval != FIRKLS_SUCCESS)
  {
    firkls_mem->last_flag = retval;
    return retval;
  }

  firkls_mem->last_flag = SUNLinSolInitialize(firkls_mem->LS);
  return firkls_mem->last_flag;
}

int firkLsReInitialize(FIRKodeMem firk_mem)
{
  FIRKLsMem firkls_mem;

  if (firk_mem->lmem == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_LMEM_NULL, __LINE__, __func__, __FILE__,
                     MSG_LS_LMEM_NULL);
    return FIRKLS_LMEM_NULL;
  }
  firkls_mem = (FIRKLsMem)firk_mem->lmem;

  firkLsInitializeCounters(firkls_mem);

  return FIRKLS_SUCCESS;
}

/*-----------------------------------------------------------------
  firkLsJacSchedule

  The FIRKLS options that the PERIODIC reuse policy reads (design
  decision D9): whether msbj steps passed since the last Jacobian
  evaluation (SetJacEvalFrequency), and the change in gamma below which
  a convergence failure with a stale Jacobian is blamed on it
  (SetDeltaGammaMaxBadJac).
  -----------------------------------------------------------------*/
void firkLsJacSchedule(FIRKodeMem firk_mem, sunbooleantype* due,
                       sunrealtype* dgmax_jbad)
{
  FIRKLsMem firkls_mem;

  *due        = SUNTRUE;
  *dgmax_jbad = FIRKLS_DGMAX;
  if (firk_mem->lmem == NULL) { return; }
  firkls_mem  = (FIRKLsMem)firk_mem->lmem;
  *due        = (firk_mem->nst >= firkls_mem->nstlj + firkls_mem->msbj);
  *dgmax_jbad = firkls_mem->dgmax_jbad;
}

/*-----------------------------------------------------------------
  firkLsSetup

  Forms A = M - gamma J at (tn, yn, fn), with a new Jacobian if jbad,
  calls the block solver's setup routine and records the setup in
  gammap, nsetups and nstlp, whatever its outcome (design Listing
  8.4).
  -----------------------------------------------------------------*/
int firkLsSetup(FIRKodeMem firk_mem, sunbooleantype jbad, sunbooleantype* jcur)
{
  FIRKLsMem firkls_mem;
  int retval;

  if (firk_mem->lmem == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_LMEM_NULL, __LINE__, __func__, __FILE__,
                     MSG_LS_LMEM_NULL);
    return FIRKLS_LMEM_NULL;
  }
  firkls_mem = (FIRKLsMem)firk_mem->lmem;

  SUNDIALS_MARK_FUNCTION_BEGIN(FIRK_PROFILER);
  retval = firkLsSetupMatrix(firk_mem, firkls_mem, jbad, jcur);

  firk_mem->gammap = firk_mem->gamma;
  firk_mem->nsetups++;
  firk_mem->nstlp = firk_mem->nst;

  SUNLogDebug(FIRK_LOGGER, "lsetup",
              "jbad = %i, jcur = %i, gamma = " SUN_FORMAT_G ", retval = %i",
              (int)jbad, (int)*jcur, firk_mem->gamma, retval);
  SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
  return retval;
}

static int firkLsSetupMatrix(FIRKodeMem firk_mem, FIRKLsMem firkls_mem,
                             sunbooleantype jbad, sunbooleantype* jcur)
{
  FIRKLsMassMem firkls_massmem;
  SUNMatrix M;
  int retval;

  if (SUNLinSolGetType(firkls_mem->LS) == SUNLINEARSOLVER_MATRIX_EMBEDDED)
  {
    firkls_mem->last_flag = FIRKLS_SUCCESS;
    return firkls_mem->last_flag;
  }

  /* linearization point */
  firkls_mem->ycur = firk_mem->yn;
  firkls_mem->fcur = firk_mem->fn;

  firkls_mem->jbad = jbad;

  if (firkls_mem->A != NULL)
  {
    M = NULL;
    if (firk_mem->mass_set && firk_mem->mass_mem != NULL)
    {
      firkls_massmem = (FIRKLsMassMem)firk_mem->mass_mem;
      M              = firkls_massmem->M;
    }

    retval = firkls_mem->linsys(firk_mem->tn, firk_mem->yn, firk_mem->fn,
                                firkls_mem->A, M, !(firkls_mem->jbad), jcur,
                                firk_mem->gamma, firkls_mem->A_data,
                                firk_mem->tempv1, firk_mem->tempv2,
                                firk_mem->tempv3);

    if (*jcur)
    {
      firkls_mem->nje++;
      firkls_mem->nstlj = firk_mem->nst;
      firkls_mem->tnlj  = firk_mem->tn;
    }

    if (retval != FIRKLS_SUCCESS)
    {
      if (firkls_mem->user_linsys)
      {
        if (retval < 0)
        {
          firkProcessError(firk_mem, FIRKLS_JACFUNC_UNRECVR, __LINE__, __func__,
                           __FILE__, MSG_LS_JACFUNC_FAILED);
          firkls_mem->last_flag = FIRKLS_JACFUNC_UNRECVR;
          return -1;
        }
        firkls_mem->last_flag = FIRKLS_JACFUNC_RECVR;
        return 1;
      }
      return retval;
    }
  }
  else
  {
    /* matrix-free: jcur = jbad unless the user's psetup reports otherwise */
    firk_mem->jcur = firkls_mem->jbad;
  }

  /* call the block solver's setup routine, which calls the user's psetup */
  firkls_mem->last_flag = SUNLinSolSetup(firkls_mem->LS, firkls_mem->A);

  /* whatever psetup reports is the currency of the Jacobian data (H29) */
  if (firkls_mem->A == NULL)
  {
    *jcur = firk_mem->jcur;
    if (*jcur)
    {
      firkls_mem->nstlj = firk_mem->nst;
      firkls_mem->tnlj  = firk_mem->tn;
    }
  }

  return firkls_mem->last_flag;
}

/*-----------------------------------------------------------------
  firkLsMapSolveFlag

  Translates a SUNLinSolSolve return value into the 0/+1/-1
  convention of the integrator (design Listing 8.4). A reduced but
  unconverged residual is accepted only where accept_reduced holds:
  at the first Newton iteration (H18).
  -----------------------------------------------------------------*/
static int firkLsMapSolveFlag(FIRKodeMem firk_mem, FIRKLsMem firkls_mem,
                              int retval, sunbooleantype accept_reduced)
{
  firkls_mem->last_flag = retval;

  switch (retval)
  {
  case SUN_SUCCESS: return 0;
  case SUNLS_RES_REDUCED: return accept_reduced ? 0 : 1;
  case SUNLS_RECOV_FAILURE:
  case SUNLS_CONV_FAIL:
  case SUNLS_ATIMES_FAIL_REC:
  case SUNLS_PSET_FAIL_REC:
  case SUNLS_PSOLVE_FAIL_REC:
  case SUNLS_PACKAGE_FAIL_REC:
  case SUNLS_QRFACT_FAIL:
  case SUNLS_LUFACT_FAIL: return 1;
  case SUN_ERR_ARG_CORRUPT:
  case SUN_ERR_ARG_INCOMPATIBLE:
  case SUN_ERR_MEM_FAIL:
  case SUNLS_ATIMES_NULL:
  case SUNLS_PSOLVE_NULL:
  case SUNLS_GS_FAIL:
  case SUNLS_QRSOL_FAIL: return -1;
  case SUN_ERR_EXT_FAIL:
    /* an unrecoverable failure of an external package (R8) */
    firkProcessError(firk_mem, SUN_ERR_EXT_FAIL, __LINE__, __func__, __FILE__,
                     "Failure in SUNLinSol external package");
    return -1;
  case SUNLS_ATIMES_FAIL_UNREC:
    firkProcessError(firk_mem, SUNLS_ATIMES_FAIL_UNREC, __LINE__, __func__,
                     __FILE__, MSG_LS_JTIMES_FAILED);
    return -1;
  case SUNLS_PSET_FAIL_UNREC:
    firkProcessError(firk_mem, SUNLS_PSET_FAIL_UNREC, __LINE__, __func__,
                     __FILE__, MSG_LS_PSET_FAILED);
    return -1;
  case SUNLS_PSOLVE_FAIL_UNREC:
    firkProcessError(firk_mem, SUNLS_PSOLVE_FAIL_UNREC, __LINE__, __func__,
                     __FILE__, MSG_LS_PSOLVE_FAILED);
    return -1;
  default:
    firkProcessError(firk_mem, retval, __LINE__, __func__, __FILE__,
                     "Unrecognized error return value from SUNLinSolSolve");
    return -1;
  }
}

/* For a preconditioner application in the stacked Krylov iteration, an
   iterative block solver that only reduced the residual (or ran out of
   iterations) still provides a useful approximate solution, which the
   flexible GMRES iteration can use, so such outcomes are not reported as
   failures.

   The tolerances of the Newton system and of the error estimate bound the
   WRMS norm of the solution x, which has the units of y. Without a mass
   matrix b has the same units, and b and the Krylov residual are tested
   against these tolerances directly, as in CVODE. With a mass matrix b has
   the units of M y, so for entries of M of order m these tests would accept
   a correction of size deltar/m. For an iterative block solver one solve
   with M then converts b to the units of y: x0 = M^{-1} b is the right-hand
   side of the same solve for the equivalent problem y' = M^{-1} f without a
   mass matrix, so its norm takes the place of ||b|| in the small right-hand
   side test of the Newton system and in the tolerance of the error
   estimate, and x0 takes the place of b as the correction returned on the
   first Newton iteration. The ratio sysfac = ||x0|| / ||b|| converts the
   Krylov residual: a matrix-free solver is given the block system
   multiplied by sysfac (the operator in firkLsATimes, the right-hand side
   here and the preconditioner in firkLsPSolve, which leaves the iterates
   unchanged), whose residual has the units of y and whose left-preconditioned
   residual is unchanged, so the tolerance holds for either preconditioning
   side. For a matrix-based iterative solver the tolerance is divided by
   sysfac instead, assuming that its test applies to the residual itself. */

/*-----------------------------------------------------------------
  firkLsJtSetup

  Calls the user's jtsetup, if attached, at (tn, yn, fn): once per
  Newton iteration, and once per estimate solve with a matrix-free
  block solver, whose operator uses jtimes (H6).
  -----------------------------------------------------------------*/
static int firkLsJtSetup(FIRKodeMem firk_mem, FIRKLsMem firkls_mem)
{
  if (firkls_mem->jtsetup == NULL) { return 0; }
  firkls_mem->last_flag = firkls_mem->jtsetup(firk_mem->tn, firk_mem->yn,
                                              firk_mem->fn, firkls_mem->jt_data);
  firkls_mem->njtsetup++;
  if (firkls_mem->last_flag != 0)
  {
    firkProcessError(firk_mem, firkls_mem->last_flag, __LINE__, __func__,
                     __FILE__, MSG_LS_JTSETUP_FAILED);
    return (firkls_mem->last_flag < 0) ? -1 : 1;
  }
  return 0;
}

/*-----------------------------------------------------------------
  firkLsSolveBlock

  Solves the block system (M - gamma J) x = b with the user's linear
  solver for the given purpose, setting tolerances and scaling
  vectors for iterative solvers; k is the Newton iteration. x may be
  b; b may be overwritten.
  -----------------------------------------------------------------*/
static int firkLsSolveBlockBody(FIRKodeMem firk_mem, FIRKLsMem firkls_mem,
                                N_Vector b, N_Vector x, firkLsBlockSolve purpose,
                                int k, sunbooleantype* begun, const char** status)
{
  N_Vector weight;
  sunrealtype bnorm  = ZERO;
  sunrealtype xnorm  = ZERO;
  sunrealtype sysfac = ONE;
  sunrealtype deltar, delta, w_mean;
  sunbooleantype convert = SUNFALSE;
  int nli_inc, retval, flag;

  weight = firk_mem->ewt;

  /* iterative block solver: tolerance and, for the Newton system from its
     second iteration, the small right-hand side shortcut */
  if (firkls_mem->iterative)
  {
    bnorm = N_VWrmsNorm(b, weight);

    /* with a mass matrix, b in the units of y (see above) */
    xnorm   = bnorm;
    convert = firk_mem->mass_set && (purpose != FIRKLS_BLOCK_PRECOND) &&
              (bnorm > ZERO);
    if (convert)
    {
      N_VScale(ONE, b, firkls_mem->x);
      retval = firkLsMassSolve(firk_mem, firkls_mem->x, bnorm, SUNFALSE);
      if (retval != 0)
      {
        *status = "failed mass solve";
        return (retval < 0) ? -1 : 1;
      }
      xnorm = N_VWrmsNorm(firkls_mem->x, weight);
    }

    if (purpose == FIRKLS_BLOCK_NEWTON)
    {
      /* Newton system: tolerance relative to the nonlinear solver tolerance,
         as in CVODE; from the second iteration a right-hand side below it (in
         the units of y) gives the correction 0 (CVLS's shortcut, design
         Listing 8.4) */
      deltar = firkls_mem->eplifac * firk_mem->nlscoef;

      SUNLogInfo(FIRK_LOGGER, "begin-linear-solve",
                 "purpose = newton, iterative = 1, b-norm = " SUN_FORMAT_G
                 ", x-norm = " SUN_FORMAT_G ", b-tol = " SUN_FORMAT_G
                 ", res-tol = " SUN_FORMAT_G,
                 bnorm, xnorm, deltar, deltar * firkls_mem->nrmfac);
      *begun = SUNTRUE;

      if (k >= 2 && xnorm <= deltar)
      {
        N_VConst(ZERO, x);
        firkls_mem->last_flag = FIRKLS_SUCCESS;
        firkls_mem->nblocksolves++;
        *status = "success small rhs";
        return firkls_mem->last_flag;
      }
    }
    else
    {
      /* error estimate: the solution is the estimate itself, so the solve is
         not skipped for a small right-hand side, and the tolerance is the
         fixed accuracy FIRKLS_EPSEST of the estimate, tightened to one
         relative to b when b (in the units of y) is smaller, so that a
         small estimate is resolved instead of returned as zero; the Newton
         tolerance, which follows the relative tolerance, plays no part.
         Preconditioner application: the right-hand side is a Krylov basis
         vector of arbitrary scale, so the tolerance is relative to it. */
      deltar = (purpose == FIRKLS_BLOCK_ESTIMATE)
                 ? firkls_mem->eplifac * SUNMIN(FIRKLS_EPSEST, xnorm)
                 : firkls_mem->eplifac * bnorm;

      /* no shortcut here (H17): the estimate is the solution itself, and a
         zero right-hand side is solved exactly with the zero guess */

      /* a preconditioner block logs one record, from firkLsSolveBlock */
      if (purpose == FIRKLS_BLOCK_ESTIMATE)
      {
        SUNLogInfo(FIRK_LOGGER, "begin-linear-solve",
                   "purpose = estimate, iterative = 1, b-norm = " SUN_FORMAT_G
                   ", x-norm = " SUN_FORMAT_G ", b-tol = " SUN_FORMAT_G
                   ", res-tol = " SUN_FORMAT_G,
                   bnorm, xnorm, deltar, deltar * firkls_mem->nrmfac);
        *begun = SUNTRUE;
      }
    }
    delta = deltar * firkls_mem->nrmfac;

    /* with a mass matrix, the Krylov residual in the units of y (see above) */
    if (convert)
    {
      sysfac = xnorm / bnorm;
      if (SUNLinSolGetType(firkls_mem->LS) != SUNLINEARSOLVER_ITERATIVE)
      {
        delta /= sysfac;
        sysfac = ONE;
      }
    }
  }
  else
  {
    delta = ZERO;
    if (purpose != FIRKLS_BLOCK_PRECOND)
    {
      SUNLogInfo(FIRK_LOGGER, "begin-linear-solve", "purpose = %s, iterative = 0",
                 (purpose == FIRKLS_BLOCK_NEWTON) ? "newton" : "estimate");
      *begun = SUNTRUE;
    }
  }

  /* linearization point for ATimes/PSolve */
  firkls_mem->ycur = firk_mem->yn;
  firkls_mem->fcur = firk_mem->fn;

  if (firkls_mem->LS->ops->setscalingvectors)
  {
    retval = SUNLinSolSetScalingVectors(firkls_mem->LS, weight, weight);
    if (retval != SUN_SUCCESS)
    {
      firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__,
                       __FILE__, "Error in calling SUNLinSolSetScalingVectors");
      firkls_mem->last_flag = FIRKLS_SUNLS_FAIL;
      return firkls_mem->last_flag;
    }
  }
  else if (firkls_mem->iterative)
  {
    N_VConst(ONE, firkls_mem->x);
    w_mean = N_VWrmsNorm(weight, firkls_mem->x);
    delta /= w_mean;
  }

  /* the user's jtsetup: once per Newton iteration of the one-stage method
     that reaches the solver, and for an estimate solve with a matrix-free
     block solver (H6) */
  if (purpose == FIRKLS_BLOCK_NEWTON ||
      (purpose == FIRKLS_BLOCK_ESTIMATE && !firkls_mem->matrixbased))
  {
    retval = firkLsJtSetup(firk_mem, firkls_mem);
    if (retval != 0) { return retval; }
  }

  N_VConst(ZERO, firkls_mem->x);
  retval = SUNLinSolSetZeroGuess(firkls_mem->LS, SUNTRUE);
  if (retval != SUN_SUCCESS) { return -1; }

  /* the operator and the user's preconditioner use gamma for the Newton
     system and gammap, the gamma of the last setup, for the estimate and a
     preconditioner block, whose matrix-based operator is that of the setup
     (design Table 8.6); without any setup, gamma */
  firkls_mem->gamma_op = (purpose == FIRKLS_BLOCK_NEWTON ||
                          firk_mem->gammap == ZERO)
                           ? firk_mem->gamma
                           : firk_mem->gammap;

  /* a matrix-free solver solves the system multiplied by sysfac */
  if (sysfac != ONE) { N_VScale(sysfac, b, b); }
  firkls_mem->sysfac = sysfac;
  retval = SUNLinSolSolve(firkls_mem->LS, firkls_mem->A, firkls_mem->x, b, delta);
  firkls_mem->sysfac = ONE;
  N_VScale(ONE, firkls_mem->x, x);
  firkls_mem->nblocksolves++;

  /* the s = 1 Newton correction from the matrix of the last setup, scaled
     for the change of gamma since then, as in CVODE; never the estimate or
     a preconditioner block */
  if (purpose == FIRKLS_BLOCK_NEWTON && firkls_mem->scalesol &&
      firk_mem->gamma != firk_mem->gammap)
  {
    N_VScale(TWO / (ONE + firk_mem->gamma / firk_mem->gammap), x, x);
  }

  nli_inc = 0;
  if (firkls_mem->iterative && firkls_mem->LS->ops->numiters)
  {
    nli_inc = SUNLinSolNumIters(firkls_mem->LS);
  }
  firkls_mem->nli += nli_inc;

  /* a preconditioner block accepts a recoverable failure and counts it */
  flag = firkLsMapSolveFlag(firk_mem, firkls_mem, retval,
                            purpose == FIRKLS_BLOCK_NEWTON && k == 1);
  if (purpose == FIRKLS_BLOCK_PRECOND && flag > 0)
  {
    firkls_mem->nblockfails++;
    *status = "accepted preconditioner failure";
    return 0;
  }
  if (retval != SUN_SUCCESS) { firkls_mem->ncfl++; }

  return flag;
}

/*-----------------------------------------------------------------
  firkLsSolveBlock

  The profiling region and log records of a block solve around
  firkLsSolveBlockBody: a Newton or estimate solve is one
  begin/end-linear-solve region, closed on every exit path once
  opened (H31), an estimate within an element of the list region
  "estimate"; a preconditioner block has no record of its own.
  -----------------------------------------------------------------*/
int firkLsSolveBlock(FIRKodeMem firk_mem, N_Vector b, N_Vector x,
                     firkLsBlockSolve purpose, int k)
{
  FIRKLsMem firkls_mem;
  sunbooleantype begun = SUNFALSE;
  const char* status   = NULL;
  long int nli0, nps0;
  int retval;

  if (firk_mem->lmem == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_LMEM_NULL, __LINE__, __func__, __FILE__,
                     MSG_LS_LMEM_NULL);
    return FIRKLS_LMEM_NULL;
  }
  firkls_mem = (FIRKLsMem)firk_mem->lmem;

  SUNDIALS_MARK_FUNCTION_BEGIN(FIRK_PROFILER);
  /* an attempt may filter the estimate twice: a list region */
  SUNLogInfoIf(purpose == FIRKLS_BLOCK_ESTIMATE, FIRK_LOGGER,
               "begin-estimate-list", "");
  nli0   = firkls_mem->nli;
  nps0   = firkls_mem->nps;
  retval = firkLsSolveBlockBody(firk_mem, firkls_mem, b, x, purpose, k, &begun,
                                &status);
  if (status == NULL) { status = (retval == 0) ? "success" : "failed"; }

  SUNLogInfoIf(begun, FIRK_LOGGER, "end-linear-solve",
               "status = %s, retval = %i, iters = %li, p-solves = %li", status,
               retval, firkls_mem->nli - nli0, firkls_mem->nps - nps0);
  SUNLogInfoIf(purpose == FIRKLS_BLOCK_ESTIMATE, FIRK_LOGGER,
               "end-estimate-list", "");
  SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
  (void)begun; /* used only by the log records */
  (void)nli0;
  (void)nps0;
  return retval;
}

/*-----------------------------------------------------------------
  firkLsSolveStacked

  Solves the stacked Newton system (I_s (x) M - h A (x) J) x = b in
  place with the internal flexible GMRES iteration, preconditioned
  block-wise by the user's linear solver. For s = 1 the block solve
  is exact and used directly (see firkLsSolveBlock for its
  treatment of a mass matrix). k is the Newton iteration; x may be
  b_stk, and b_stk may be overwritten.

  The tolerance deltar bounds the WRMS norm of the correction x, which
  has the units of y. Without a mass matrix the residual b of the stage
  system has the same units, and b and the Krylov residual are tested
  against deltar directly, as in CVODE. With a mass matrix b has the
  units of M y (b_i = -(M Z_i - h sum_j a_ij F_j)), so for entries of M
  of order m both tests would accept a correction of size deltar/m.
  As for the block solves (see firkLsSolveBlock), one solve with M
  per stage converts b to the units of y: x0 = (I_s (x) M)^{-1} b is the
  right-hand side of the stacked Newton system of the equivalent
  problem y' = M^{-1} f without a mass matrix, and the ratio
  ||b|| / ||x0|| converts the Krylov tolerance to the units of b.
  Scaling the tolerance is exact here since the flexible GMRES
  iteration is preconditioned on the right only. There is no
  small-right-hand-side shortcut: a right-hand side already below the
  tolerance ends the Krylov iteration at its zero initial guess, so
  the correction is zero (design Listing 8.4).
  -----------------------------------------------------------------*/
static int firkLsSolveStackedBody(FIRKodeMem firk_mem, FIRKLsMem firkls_mem,
                                  N_Vector b_stk, N_Vector x, int k,
                                  sunbooleantype* begun, const char** status)
{
  N_Vector x_i;
  sunrealtype bnorm, xnorm, deltar, delta;
  sunbooleantype convert;
  int i, nli_inc, retval;

  /* single stage: the block solve is the exact Newton solve */
  if (firk_mem->s == 1)
  {
    return firkLsSolveBlock(firk_mem, firkStage(firk_mem, b_stk, 0),
                            firkStage(firk_mem, x, 0), FIRKLS_BLOCK_NEWTON, k);
  }

  /* a stacked solver made stale by a setter is created again here, where
     no solve is in progress (H4) */
  if (firkls_mem->stk_stale || firkls_mem->SPFGMR == NULL)
  {
    if (firkLsSetupStacked(firk_mem, firkls_mem) != FIRKLS_SUCCESS)
    {
      return -1;
    }
  }

  /* the user's jtsetup, exactly once per Newton iteration (H6) */
  retval = firkLsJtSetup(firk_mem, firkls_mem);
  if (retval != 0) { return retval; }

  /* the tolerance eps_stk eps_nls sqrt(s) nrmfac, with the nrmfac of the
     block system (design Table 8.6), which SetLSNormFactor may set */
  deltar = firkls_mem->eps_stk * firk_mem->nlscoef;
  delta  = deltar * SUNRsqrt((sunrealtype)firk_mem->s) * firkls_mem->nrmfac;
  bnorm  = N_VWrmsNorm(b_stk, firk_mem->ewt_stk);

  /* linearization point for ATimes/PSolve */
  firkls_mem->ycur = firk_mem->yn;
  firkls_mem->fcur = firk_mem->fn;

  /* with a mass matrix, b in the units of y (see above); each stage block
     is solved to a tolerance relative to the norm of all of b */
  xnorm   = bnorm;
  convert = firk_mem->mass_set && (bnorm > ZERO);
  if (convert)
  {
    for (i = 0; i < firk_mem->s; i++)
    {
      x_i = firkStage(firk_mem, x, i);
      N_VScale(ONE, firkStage(firk_mem, b_stk, i), x_i);
      retval = firkLsMassSolve(firk_mem, x_i, bnorm, SUNFALSE);
      if (retval != 0)
      {
        *status = "failed mass solve";
        return (retval < 0) ? -1 : 1;
      }
    }
    xnorm = N_VWrmsNorm(x, firk_mem->ewt_stk);
    if (xnorm > ZERO) { delta *= bnorm / xnorm; }
  }

  SUNLogInfo(FIRK_LOGGER, "begin-linear-solve",
             "purpose = newton, stacked = 1, iterative = 1, b-norm "
             "= " SUN_FORMAT_G ", x-norm = " SUN_FORMAT_G
             ", b-tol = " SUN_FORMAT_G ", res-tol = " SUN_FORMAT_G,
             bnorm, xnorm, deltar, delta);
  *begun = SUNTRUE;

  retval = SUNLinSolSetScalingVectors(firkls_mem->SPFGMR, firk_mem->ewt_stk,
                                      firk_mem->ewt_stk);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__, __FILE__,
                     "Error in calling SUNLinSolSetScalingVectors");
    firkls_mem->last_flag = FIRKLS_SUNLS_FAIL;
    return firkls_mem->last_flag;
  }

  /* the solution in the caller's vector, from the zero guess (the
     conversion above only used it for the norm) */
  N_VConst(ZERO, x);
  retval = SUNLinSolSetZeroGuess(firkls_mem->SPFGMR, SUNTRUE);
  if (retval != SUN_SUCCESS) { return -1; }

  retval = SUNLinSolSolve(firkls_mem->SPFGMR, NULL, x, b_stk, delta);

  nli_inc = SUNLinSolNumIters(firkls_mem->SPFGMR);
  firkls_mem->nsli += nli_inc;
  if (retval != SUN_SUCCESS) { firkls_mem->nscf++; }

  return firkLsMapSolveFlag(firk_mem, firkls_mem, retval, k == 1);
}

/*-----------------------------------------------------------------
  firkLsSolveStacked

  The profiling region and the begin/end-linear-solve region of the
  Newton solve around firkLsSolveStackedBody; the region is closed
  on every exit path once opened (H31). One stage logs as a block
  solve.
  -----------------------------------------------------------------*/
int firkLsSolveStacked(FIRKodeMem firk_mem, N_Vector b_stk, N_Vector x, int k)
{
  FIRKLsMem firkls_mem;
  sunbooleantype begun = SUNFALSE;
  const char* status   = NULL;
  long int nsli0, nps0, nbs0;
  int retval;

  if (firk_mem->lmem == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_LMEM_NULL, __LINE__, __func__, __FILE__,
                     MSG_LS_LMEM_NULL);
    return FIRKLS_LMEM_NULL;
  }
  firkls_mem = (FIRKLsMem)firk_mem->lmem;

  SUNDIALS_MARK_FUNCTION_BEGIN(FIRK_PROFILER);
  nsli0  = firkls_mem->nsli;
  nps0   = firkls_mem->nps;
  nbs0   = firkls_mem->nblocksolves;
  retval = firkLsSolveStackedBody(firk_mem, firkls_mem, b_stk, x, k, &begun,
                                  &status);
  if (status == NULL) { status = (retval == 0) ? "success" : "failed"; }

  SUNLogInfoIf(begun, FIRK_LOGGER, "end-linear-solve",
               "status = %s, retval = %i, iters = %li, p-solves = %li, "
               "block-solves = %li, resnorm = " SUN_FORMAT_G,
               status, retval, firkls_mem->nsli - nsli0, firkls_mem->nps - nps0,
               firkls_mem->nblocksolves - nbs0,
               (firkls_mem->SPFGMR != NULL) ? SUNLinSolResNorm(firkls_mem->SPFGMR)
                                            : ZERO);
  SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
  (void)begun; /* used only by the log records */
  (void)nsli0;
  (void)nps0;
  (void)nbs0;
  return retval;
}

/*-----------------------------------------------------------------
  firkLsFree

  Frees memory associated with the FIRKLS interface.
  -----------------------------------------------------------------*/
int firkLsFree(FIRKodeMem firk_mem)
{
  FIRKLsMem firkls_mem;

  if (firk_mem == NULL) { return FIRKLS_SUCCESS; }
  if (firk_mem->lmem == NULL) { return FIRKLS_SUCCESS; }
  firkls_mem = (FIRKLsMem)firk_mem->lmem;

  firkLsFreeStacked(firkls_mem);

  /* detach the interface routines from the user's solver, which outlives
     this memory (as the mass interface does) */
  if (firkls_mem->LS && firkls_mem->LS->ops)
  {
    if (firkls_mem->LS->ops->setatimes)
    {
      SUNLinSolSetATimes(firkls_mem->LS, NULL, NULL);
    }
    if (firkls_mem->LS->ops->setpreconditioner)
    {
      SUNLinSolSetPreconditioner(firkls_mem->LS, NULL, NULL, NULL);
    }
  }

  if (firkls_mem->ytemp)
  {
    N_VDestroy(firkls_mem->ytemp);
    firkls_mem->ytemp = NULL;
  }
  if (firkls_mem->x)
  {
    N_VDestroy(firkls_mem->x);
    firkls_mem->x = NULL;
  }
  if (firkls_mem->savedJ)
  {
    SUNMatDestroy(firkls_mem->savedJ);
    firkls_mem->savedJ = NULL;
  }

  firkls_mem->ycur = NULL;
  firkls_mem->fcur = NULL;
  firkls_mem->A    = NULL;

  free(firk_mem->lmem);
  firk_mem->lmem = NULL;

  firk_mem->linit      = NULL;
  firk_mem->lreinit    = NULL;
  firk_mem->lsetup     = NULL;
  firk_mem->lsolve_stk = NULL;
  firk_mem->lsolve_blk = NULL;
  firk_mem->lfree      = NULL;

  return FIRKLS_SUCCESS;
}

/*-----------------------------------------------------------------
  firkLsInitializeCounters
  -----------------------------------------------------------------*/
int firkLsInitializeCounters(FIRKLsMem firkls_mem)
{
  firkls_mem->nje          = 0;
  firkls_mem->nfeDQ        = 0;
  firkls_mem->nstlj        = 0;
  firkls_mem->npe          = 0;
  firkls_mem->nli          = 0;
  firkls_mem->nps          = 0;
  firkls_mem->ncfl         = 0;
  firkls_mem->njtsetup     = 0;
  firkls_mem->njtimes      = 0;
  firkls_mem->nblocksolves = 0;
  firkls_mem->nblockfails  = 0;
  firkls_mem->nsli         = 0;
  firkls_mem->nscf         = 0;
  firkls_mem->tnlj         = ZERO;
  return 0;
}

/*---------------------------------------------------------------
  firkLs_AccessLMem

  Unpacks the firk_mem and lmem structures from a void* pointer.
  ---------------------------------------------------------------*/
int firkLs_AccessLMem(void* firkode_mem, const char* fname,
                      FIRKodeMem* firk_mem, FIRKLsMem* firkls_mem)
{
  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRKLS_MEM_NULL, __LINE__, fname, __FILE__,
                     MSG_LS_FIRKMEM_NULL);
    return FIRKLS_MEM_NULL;
  }
  *firk_mem = (FIRKodeMem)firkode_mem;
  if ((*firk_mem)->lmem == NULL)
  {
    firkProcessError(*firk_mem, FIRKLS_LMEM_NULL, __LINE__, fname, __FILE__,
                     MSG_LS_LMEM_NULL);
    return FIRKLS_LMEM_NULL;
  }
  *firkls_mem = (FIRKLsMem)(*firk_mem)->lmem;
  return FIRKLS_SUCCESS;
}

/*---------------------------------------------------------------
  firkLs_AccessLMemExpert

  As firkLs_AccessLMem, for the setters of EXPERT options: after
  the memory check, BEGINNER mode rejects the call with the
  integrator code FIRK_ILL_INPUT (design decision D8), before the
  linear solver memory is checked.
  ---------------------------------------------------------------*/
int firkLs_AccessLMemExpert(void* firkode_mem, const char* fname,
                            FIRKodeMem* firk_mem, FIRKLsMem* firkls_mem)
{
  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRKLS_MEM_NULL, __LINE__, fname, __FILE__,
                     MSG_LS_FIRKMEM_NULL);
    return FIRKLS_MEM_NULL;
  }
  if (!firkExpertAllowed((FIRKodeMem)firkode_mem, fname))
  {
    return FIRK_ILL_INPUT;
  }
  return firkLs_AccessLMem(firkode_mem, fname, firk_mem, firkls_mem);
}

/*---------------------------------------------------------------
  firkLsSetExpertDefaults

  Sets the options of the linear solver interface that only EXPERT
  mode may change to their defaults (design decision D8): the
  Jacobian evaluation frequency, the bad-Jacobian threshold, and the
  options of the stacked FGMRES solver (dimension, restarts,
  Gram-Schmidt type). Used by FIRKodeSetLinearSolver and by
  FIRKodeSetMode.
  ---------------------------------------------------------------*/
int firkLsSetExpertDefaults(FIRKodeMem firk_mem)
{
  FIRKLsMem firkls_mem = (FIRKLsMem)firk_mem->lmem;

  firkls_mem->msbj       = FIRKLS_MSBJ;
  firkls_mem->dgmax_jbad = FIRKLS_DGMAX;
  firkls_mem->eps_stk    = FIRKLS_EPLIN;
  firkls_mem->maxrs      = FIRKLS_MAXRS;
  firkls_mem->gstype     = SUN_MODIFIED_GS;
  if (firkls_mem->maxl_set)
  {
    /* the stacked solver is created again with the default dimension at
       the next stacked solve: it may be live now (H4) */
    firkls_mem->maxl_set  = SUNFALSE;
    firkls_mem->maxl      = 0;
    firkls_mem->stk_stale = SUNTRUE;
  }
  if (firkls_mem->SPFGMR != NULL &&
      (SUNLinSol_SPFGMRSetMaxRestarts(firkls_mem->SPFGMR, firkls_mem->maxrs) !=
         SUN_SUCCESS ||
       SUNLinSol_SPFGMRSetGSType(firkls_mem->SPFGMR, firkls_mem->gstype) !=
         SUN_SUCCESS))
  {
    firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__, __FILE__,
                     "Error in calling the SPFGMR option setters");
    return FIRKLS_SUNLS_FAIL;
  }
  return FIRKLS_SUCCESS;
}

/*===============================================================
  FIRKLS mass matrix interface (constant mass matrices only)
  ===============================================================*/

/*---------------------------------------------------------------
  FIRKodeSetMassLinearSolver

  Attaches a SUNLinearSolver (and optionally a SUNMatrix) for the
  mass matrix M in M y' = f(t, y). Only time-independent mass
  matrices are supported: M is constructed and factored once, at
  the first call to FIRKodeEvolve after FIRKodeInit/ReInit.
  ---------------------------------------------------------------*/
int FIRKodeSetMassLinearSolver(void* firkode_mem, SUNLinearSolver LS,
                               SUNMatrix M, sunbooleantype time_dep)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval, LSType;
  sunbooleantype iterative;   /* is the solver iterative?    */
  sunbooleantype matrixbased; /* is a matrix structure used? */

  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRKLS_MEM_NULL, __LINE__, __func__, __FILE__,
                     MSG_LS_FIRKMEM_NULL);
    return FIRKLS_MEM_NULL;
  }
  if (LS == NULL)
  {
    firkProcessError(NULL, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "LS must be non-NULL");
    return FIRKLS_ILL_INPUT;
  }
  firk_mem = (FIRKodeMem)firkode_mem;

  if (firk_mem->MallocDone == SUNFALSE)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MALLOC);
    return FIRKLS_ILL_INPUT;
  }

  if (time_dep)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "FIRKODE only supports time-independent mass matrices");
    return FIRKLS_ILL_INPUT;
  }

  /* test if solver is compatible with the LS interface */
  if ((LS->ops->gettype == NULL) || (LS->ops->solve == NULL))
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "LS object is missing a required operation");
    return FIRKLS_ILL_INPUT;
  }

  LSType = SUNLinSolGetType(LS);

  iterative   = (LSType != SUNLINEARSOLVER_DIRECT);
  matrixbased = ((LSType != SUNLINEARSOLVER_ITERATIVE) &&
                 (LSType != SUNLINEARSOLVER_MATRIX_EMBEDDED));

  /* test if vector is compatible with the LS interface */
  if ((firk_mem->tempv1->ops->nvconst == NULL) ||
      (firk_mem->tempv1->ops->nvwrmsnorm == NULL) ||
      (firk_mem->tempv1->ops->nvgetlength == NULL))
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSG_LS_BAD_NVECTOR);
    return FIRKLS_ILL_INPUT;
  }

  /* ensure that M is NULL when LS is matrix-embedded */
  if ((LSType == SUNLINEARSOLVER_MATRIX_EMBEDDED) && (M != NULL))
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Incompatible inputs: matrix-embedded LS requires NULL "
                     "matrix");
    return FIRKLS_ILL_INPUT;
  }

  /* check for compatible LS type, matrix and "atimes" support */
  if (iterative)
  {
    if (!matrixbased && (LSType != SUNLINEARSOLVER_MATRIX_EMBEDDED) &&
        (LS->ops->setatimes == NULL))
    {
      firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                       "Incompatible inputs: iterative LS must support ATimes "
                       "routine");
      return FIRKLS_ILL_INPUT;
    }
    if (matrixbased && (M == NULL))
    {
      firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                       "Incompatible inputs: matrix-iterative LS requires "
                       "non-NULL matrix");
      return FIRKLS_ILL_INPUT;
    }
  }
  else if (M == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Incompatible inputs: direct LS requires non-NULL matrix");
    return FIRKLS_ILL_INPUT;
  }

  /* free any existing mass matrix interface memory */
  if (firk_mem->mfree != NULL) { firk_mem->mfree(firk_mem); }

  /* allocate the interface memory */
  firkls_mem = (FIRKLsMassMem)malloc(sizeof(struct FIRKLsMassMemRec));
  if (firkls_mem == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSG_LS_MEM_FAIL);
    return FIRKLS_MEM_FAIL;
  }
  memset(firkls_mem, 0, sizeof(struct FIRKLsMassMemRec));

  firkls_mem->LS          = LS;
  firkls_mem->iterative   = iterative;
  firkls_mem->matrixbased = matrixbased;

  firkls_mem->mass    = NULL;
  firkls_mem->M_data  = firk_mem->user_data;
  firkls_mem->mtsetup = NULL;
  firkls_mem->mtimes  = NULL;
  firkls_mem->mt_data = NULL;

  firkls_mem->pset   = NULL;
  firkls_mem->psolve = NULL;
  firkls_mem->P_data = firk_mem->user_data;

  firkls_mem->eplifac   = FIRKLS_EPLIN;
  firkls_mem->nrmfac    = ZERO;
  firkls_mem->last_flag = FIRKLS_SUCCESS;

  firkls_mem->nmsetups = 0;
  firkls_mem->nmtsetup = 0;
  firkls_mem->nmmults  = 0;
  firkls_mem->nmsolves = 0;
  firkls_mem->nmpe     = 0;
  firkls_mem->nmli     = 0;
  firkls_mem->nmps     = 0;
  firkls_mem->nmcfl    = 0;

  /* if LS supports ATimes, attach the FIRKLs routine */
  if (LS->ops->setatimes)
  {
    retval = SUNLinSolSetATimes(LS, firk_mem, firkLsMTimes);
    if (retval != SUN_SUCCESS)
    {
      firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__,
                       __FILE__, "Error in calling SUNLinSolSetATimes");
      free(firkls_mem);
      return FIRKLS_SUNLS_FAIL;
    }
  }

  /* if LS supports preconditioning, initialize pset/psolve to NULL */
  if (LS->ops->setpreconditioner)
  {
    retval = SUNLinSolSetPreconditioner(LS, firk_mem, NULL, NULL);
    if (retval != SUN_SUCCESS)
    {
      firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__,
                       __FILE__, "Error in calling SUNLinSolSetPreconditioner");
      free(firkls_mem);
      return FIRKLS_SUNLS_FAIL;
    }
  }

  /* store M and, for direct solvers, create M_lu to hold its factorization
     so that M itself stays available for matrix-vector products */
  firkls_mem->M    = M;
  firkls_mem->M_lu = M;
  if (M != NULL && !iterative)
  {
    firkls_mem->M_lu = SUNMatClone(M);
    if (firkls_mem->M_lu == NULL)
    {
      firkProcessError(firk_mem, FIRKLS_MEM_FAIL, __LINE__, __func__, __FILE__,
                       MSG_LS_MEM_FAIL);
      free(firkls_mem);
      return FIRKLS_MEM_FAIL;
    }
  }

  /* work vector for the mass solve */
  firkls_mem->x = N_VClone(firk_mem->tempv1);
  if (firkls_mem->x == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSG_LS_MEM_FAIL);
    if (M != NULL && !iterative) { SUNMatDestroy(firkls_mem->M_lu); }
    free(firkls_mem);
    return FIRKLS_MEM_FAIL;
  }

  /* default norm conversion factor for iterative solvers */
  if (iterative)
  {
    firkls_mem->nrmfac = SUNRsqrt((sunrealtype)N_VGetLength(firkls_mem->x));
  }

  /* attach the interface to the integrator */
  firk_mem->minit          = firkLsMassInitialize;
  firk_mem->mreinit        = firkLsMassReInitialize;
  firk_mem->msetup         = firkLsMassSetup;
  firk_mem->mmult          = firkLsMassMult;
  firk_mem->msolve         = firkLsMassSolveDeriv;
  firk_mem->mfree          = firkLsMassFree;
  firk_mem->mass_mem       = firkls_mem;
  firk_mem->mass_set       = SUNTRUE;
  firk_mem->mass_iterative = iterative;

  return FIRKLS_SUCCESS;
}

/*---------------------------------------------------------------
  Mass matrix optional inputs
  ---------------------------------------------------------------*/

int FIRKodeSetMassFn(void* firkode_mem, FIRKLsMassFn mass)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (mass == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Mass-matrix routine must be non-NULL");
    return FIRKLS_ILL_INPUT;
  }
  if (firkls_mem->M == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__,
                     __FILE__, "Mass-matrix routine cannot be supplied for NULL SUNMatrix");
    return FIRKLS_ILL_INPUT;
  }

  firkls_mem->mass        = mass;
  firkls_mem->M_data      = firk_mem->user_data;
  firkls_mem->msetup_done = SUNFALSE;

  return FIRKLS_SUCCESS;
}

int FIRKodeSetMassTimes(void* firkode_mem, FIRKLsMassTimesSetupFn msetup,
                        FIRKLsMassTimesVecFn mtimes, void* mtimes_data)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (mtimes == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "non-NULL mtimes function must be supplied");
    return FIRKLS_ILL_INPUT;
  }

  firkls_mem->mtsetup     = msetup;
  firkls_mem->mtimes      = mtimes;
  firkls_mem->mt_data     = mtimes_data;
  firkls_mem->msetup_done = SUNFALSE;

  return FIRKLS_SUCCESS;
}

int FIRKodeSetMassPreconditioner(void* firkode_mem, FIRKLsMassPrecSetupFn psetup,
                                 FIRKLsMassPrecSolveFn psolve)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  SUNPSetupFn firkls_mpsetup;
  SUNPSolveFn firkls_mpsolve;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (firkls_mem->LS->ops->setpreconditioner == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "SUNLinearSolver object does not support user-supplied "
                     "preconditioning");
    return FIRKLS_ILL_INPUT;
  }

  firkls_mem->pset        = psetup;
  firkls_mem->psolve      = psolve;
  firkls_mem->P_data      = firk_mem->user_data;
  firkls_mem->msetup_done = SUNFALSE;

  firkls_mpsetup = (psetup == NULL) ? NULL : firkLsMPSetup;
  firkls_mpsolve = (psolve == NULL) ? NULL : firkLsMPSolve;
  retval = SUNLinSolSetPreconditioner(firkls_mem->LS, firk_mem, firkls_mpsetup,
                                      firkls_mpsolve);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__, __FILE__,
                     "Error in calling SUNLinSolSetPreconditioner");
    return FIRKLS_SUNLS_FAIL;
  }

  return FIRKLS_SUCCESS;
}

int FIRKodeSetMassEpsLin(void* firkode_mem, sunrealtype eplifac)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  firkls_mem->eplifac = (eplifac <= ZERO) ? FIRKLS_EPLIN : eplifac;

  return FIRKLS_SUCCESS;
}

int FIRKodeSetMassLSNormFactor(void* firkode_mem, sunrealtype nrmfac)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (nrmfac > ZERO) { firkls_mem->nrmfac = nrmfac; }
  else if (nrmfac < ZERO)
  {
    if (firkls_mem->x->ops->nvdotprod == NULL)
    {
      firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                       "N_VDotProd unimplemented (required for "
                       "FIRKodeSetMassLSNormFactor)");
      return FIRKLS_ILL_INPUT;
    }
    N_VConst(ONE, firkls_mem->x);
    firkls_mem->nrmfac = SUNRsqrt(N_VDotProd(firkls_mem->x, firkls_mem->x));
  }
  else
  {
    firkls_mem->nrmfac = SUNRsqrt((sunrealtype)N_VGetLength(firkls_mem->x));
  }

  return FIRKLS_SUCCESS;
}

/*---------------------------------------------------------------
  Mass matrix optional outputs
  ---------------------------------------------------------------*/

int FIRKodeGetCurrentMassMatrix(void* firkode_mem, SUNMatrix* M_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *M_ptr = firkls_mem->M;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumMassSetups(void* firkode_mem, long int* nmsetups_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nmsetups_ptr = firkls_mem->nmsetups;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumMassMultSetups(void* firkode_mem, long int* nmvsetups_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nmvsetups_ptr = firkls_mem->nmtsetup;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumMassMult(void* firkode_mem, long int* nmvevals_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nmvevals_ptr = firkls_mem->nmmults;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumMassSolves(void* firkode_mem, long int* nmsolves_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nmsolves_ptr = firkls_mem->nmsolves;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumMassPrecEvals(void* firkode_mem, long int* nmpevals_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nmpevals_ptr = firkls_mem->nmpe;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumMassPrecSolves(void* firkode_mem, long int* nmpsolves_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nmpsolves_ptr = firkls_mem->nmps;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumMassIters(void* firkode_mem, long int* nmiters_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nmiters_ptr = firkls_mem->nmli;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumMassConvFails(void* firkode_mem, long int* nmcfails_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nmcfails_ptr = firkls_mem->nmcfl;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetLastMassFlag(void* firkode_mem, long int* flag_ptr)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *flag_ptr = firkls_mem->last_flag;
  return FIRKLS_SUCCESS;
}

/*---------------------------------------------------------------
  firkLsMTimes

  ATimes routine for the mass linear solver: z = M v, using the
  user's mtimes routine when supplied and SUNMatMatvec otherwise.
  ---------------------------------------------------------------*/
int firkLsMTimes(void* firkode_mem, N_Vector v, N_Vector z)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (firkls_mem->mtimes)
  {
    retval = firkls_mem->mtimes(v, z, firk_mem->tcur, firkls_mem->mt_data);
    firkls_mem->nmmults++;
    if (retval != 0)
    {
      firkProcessError(firk_mem, retval, __LINE__, __func__, __FILE__,
                       "Error in user mass matrix-vector product routine");
    }
    return retval;
  }
  else if (firkls_mem->M && firkls_mem->M->ops->matvec)
  {
    retval = SUNMatMatvec(firkls_mem->M, v, z);
    firkls_mem->nmmults++;
    if (retval != 0)
    {
      firkProcessError(firk_mem, retval, __LINE__, __func__, __FILE__,
                       "Error in SUNMatrix mass matrix-vector product routine");
    }
    return retval;
  }

  firkProcessError(firk_mem, -1, __LINE__, __func__, __FILE__,
                   "Missing mass matrix-vector product routine");
  return -1;
}

/*---------------------------------------------------------------
  firkLsMPSetup / firkLsMPSolve

  Interface the user's mass matrix preconditioner routines with
  the mass linear solver. Since the mass matrix is constant, the
  setup routine is called at most once.
  ---------------------------------------------------------------*/
int firkLsMPSetup(void* firkode_mem)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (firkls_mem->nmpe > 0) { return 0; }

  retval = firkls_mem->pset(firk_mem->tcur, firkls_mem->P_data);
  firkls_mem->nmpe++;
  return retval;
}

int firkLsMPSolve(void* firkode_mem, N_Vector r, N_Vector z, sunrealtype tol,
                  int lr)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  retval = firkls_mem->psolve(firk_mem->tcur, r, z, tol, lr, firkls_mem->P_data);
  firkls_mem->nmps++;
  return retval;
}

/*---------------------------------------------------------------
  firkLsMassInitialize

  Checks the mass matrix configuration and initializes the mass
  linear solver. Called from firkInitialSetup.
  ---------------------------------------------------------------*/
int firkLsMassInitialize(FIRKodeMem firk_mem)
{
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firk_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  /* the counters are zeroed when the solver is attached and by ReInit, not
     here: Reset keeps them (H3) */

  /* matrix-based mass system */
  if (firkls_mem->M != NULL)
  {
    if (firkls_mem->mass == NULL)
    {
      firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                       "Missing user-provided mass-matrix routine");
      firkls_mem->last_flag = FIRKLS_ILL_INPUT;
      return firkls_mem->last_flag;
    }
    if ((firkls_mem->mtimes == NULL) && (firkls_mem->M->ops->matvec == NULL))
    {
      firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                       "No available mass matrix-vector product routine");
      firkls_mem->last_flag = FIRKLS_ILL_INPUT;
      return firkls_mem->last_flag;
    }
  }

  /* matrix-free mass system */
  if ((firkls_mem->M == NULL) && (firkls_mem->mtimes == NULL) &&
      (SUNLinSolGetType(firkls_mem->LS) != SUNLINEARSOLVER_MATRIX_EMBEDDED))
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__,
                     __FILE__, "Missing user-provided mass matrix-vector product routine");
    firkls_mem->last_flag = FIRKLS_ILL_INPUT;
    return firkls_mem->last_flag;
  }

  /* the setup routine has nothing to do for a matrix-free system without
     preconditioner or mtimes setup routines, or for a matrix-embedded LS */
  if (((firkls_mem->M == NULL) && (firkls_mem->pset == NULL) &&
       (firkls_mem->mtsetup == NULL)) ||
      (SUNLinSolGetType(firkls_mem->LS) == SUNLINEARSOLVER_MATRIX_EMBEDDED))
  {
    firk_mem->msetup = NULL;
  }
  else { firk_mem->msetup = firkLsMassSetup; }

  /* the solver of the constant mass matrix is initialized and set up once:
     an initialization after the setup could discard it (H28) */
  firkls_mem->last_flag = FIRKLS_SUCCESS;
  if (!firkls_mem->msetup_done)
  {
    firkls_mem->last_flag = SUNLinSolInitialize(firkls_mem->LS);
  }
  return firkls_mem->last_flag;
}

/*---------------------------------------------------------------
  firkLsMassReInitialize

  Zeroes the mass counters (ReInit) and has the mass matrix set up
  anew for the new problem; Reset keeps both (H3, H28).
  ---------------------------------------------------------------*/
int firkLsMassReInitialize(FIRKodeMem firk_mem)
{
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firk_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  firkls_mem->nmsetups = 0;
  firkls_mem->nmtsetup = 0;
  firkls_mem->nmmults  = 0;
  firkls_mem->nmsolves = 0;
  firkls_mem->nmpe     = 0;
  firkls_mem->nmli     = 0;
  firkls_mem->nmps     = 0;
  firkls_mem->nmcfl    = 0;

  firkls_mem->msetup_done = SUNFALSE;
  return FIRKLS_SUCCESS;
}

/*---------------------------------------------------------------
  firkLsMassSetup

  Constructs the (constant) mass matrix and sets up the mass
  linear solver. Called from firkInitialSetup; returns at once when
  the mass matrix is already set up (H28).
  ---------------------------------------------------------------*/
int firkLsMassSetup(FIRKodeMem firk_mem, sunrealtype t, N_Vector vtemp1,
                    N_Vector vtemp2, N_Vector vtemp3)
{
  FIRKLsMassMem firkls_mem;
  sunbooleantype call_lssetup;
  int retval;

  retval = firkLs_AccessMassMem(firk_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  /* the mass matrix is constant: set up once after the mass solver or a
     mass routine is attached and after ReInit, not again after Reset (H28) */
  if (firkls_mem->msetup_done)
  {
    firkls_mem->last_flag = FIRKLS_SUCCESS;
    return firkls_mem->last_flag;
  }

  /* user-supplied mass-times setup routine */
  if (firkls_mem->mtsetup)
  {
    firkls_mem->last_flag = firkls_mem->mtsetup(t, firkls_mem->mt_data);
    firkls_mem->nmtsetup++;
    if (firkls_mem->last_flag != 0)
    {
      firkProcessError(firk_mem, firkls_mem->last_flag, __LINE__, __func__,
                       __FILE__, MSG_LS_MTSETUP_FAILED);
      return firkls_mem->last_flag;
    }
  }

  if (firkls_mem->M == NULL)
  {
    /* matrix-free: only call the LS setup if a preconditioner setup exists */
    call_lssetup = (firkls_mem->pset != NULL);
  }
  else
  {
    /* clear the mass matrix for direct solvers */
    if (!firkls_mem->iterative)
    {
      retval = SUNMatZero(firkls_mem->M);
      if (retval)
      {
        firkProcessError(firk_mem, FIRKLS_SUNMAT_FAIL, __LINE__, __func__,
                         __FILE__, MSG_LS_SUNMAT_FAILED);
        firkls_mem->last_flag = FIRKLS_SUNMAT_FAIL;
        return firkls_mem->last_flag;
      }
    }

    /* call the user routine to fill the mass matrix */
    retval = firkls_mem->mass(t, firkls_mem->M, firkls_mem->M_data, vtemp1,
                              vtemp2, vtemp3);
    if (retval < 0)
    {
      firkProcessError(firk_mem, FIRKLS_MASSFUNC_UNRECVR, __LINE__, __func__,
                       __FILE__, MSG_LS_MASSFUNC_FAILED);
      firkls_mem->last_flag = FIRKLS_MASSFUNC_UNRECVR;
      return -1;
    }
    if (retval > 0)
    {
      firkls_mem->last_flag = FIRKLS_MASSFUNC_RECVR;
      return 1;
    }

    /* copy M into M_lu for factorization by direct solvers */
    if (!firkls_mem->iterative)
    {
      retval = SUNMatCopy(firkls_mem->M, firkls_mem->M_lu);
      if (retval)
      {
        firkProcessError(firk_mem, FIRKLS_SUNMAT_FAIL, __LINE__, __func__,
                         __FILE__, MSG_LS_SUNMAT_FAILED);
        firkls_mem->last_flag = FIRKLS_SUNMAT_FAIL;
        return firkls_mem->last_flag;
      }
    }

    /* matrix-vector product setup when the SUNMatrix provides one and the
       user did not supply mtimes */
    if ((firkls_mem->mtimes == NULL) && (firkls_mem->M->ops->matvecsetup))
    {
      retval = SUNMatMatvecSetup(firkls_mem->M);
      firkls_mem->nmtsetup++;
      if (retval)
      {
        firkProcessError(firk_mem, FIRKLS_SUNMAT_FAIL, __LINE__, __func__,
                         __FILE__, MSG_LS_SUNMAT_FAILED);
        firkls_mem->last_flag = FIRKLS_SUNMAT_FAIL;
        return firkls_mem->last_flag;
      }
    }

    call_lssetup = SUNTRUE;
  }

  if (call_lssetup)
  {
    firkls_mem->last_flag = SUNLinSolSetup(firkls_mem->LS, firkls_mem->M_lu);
    firkls_mem->nmsetups++;
  }
  firkls_mem->msetup_done = (firkls_mem->last_flag == 0);

  return firkls_mem->last_flag;
}

/*---------------------------------------------------------------
  firkLsMassMult

  Computes Mv = M v for the integrator (residuals, error estimate
  and stacked matrix-vector products).
  ---------------------------------------------------------------*/
int firkLsMassMult(FIRKodeMem firk_mem, N_Vector v, N_Vector Mv)
{
  int retval;
  SUNDIALS_MARK_FUNCTION_BEGIN(FIRK_PROFILER);
  retval = firkLsMTimes(firk_mem, v, Mv);
  SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
  return retval;
}

/*---------------------------------------------------------------
  firkLsMassSolveBody (firkLsMassSolve)

  Solves M x = b in place: sets the tolerance and scaling vectors,
  calls the mass linear solver and accumulates statistics. The
  return value is 0 on success, 1 for a recoverable failure and
  -1 for an unrecoverable failure (design Table 8.6). A residual that
  was reduced without reaching the tolerance is accepted when
  accept_reduced holds (a derivative of y, which only feeds
  heuristics), and is a recoverable failure otherwise (the units
  conversion of an enclosing solve, which then fails recoverably).
  ---------------------------------------------------------------*/
static int firkLsMassSolveBody(FIRKodeMem firk_mem, FIRKLsMassMem firkls_mem,
                               N_Vector b, sunrealtype delta,
                               sunbooleantype accept_reduced)
{
  sunrealtype ewt_mean;
  int nli_inc, retval;

  if (firkls_mem->LS->ops->setscalingvectors)
  {
    retval = SUNLinSolSetScalingVectors(firkls_mem->LS, firk_mem->ewt,
                                        firk_mem->ewt);
    if (retval != SUN_SUCCESS)
    {
      firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__,
                       __FILE__, "Error in call to SUNLinSolSetScalingVectors");
      firkls_mem->last_flag = FIRKLS_SUNLS_FAIL;
      return -1;
    }
  }
  else if (firkls_mem->iterative)
  {
    /* the solver cannot scale: approximate the WRMS tolerance in the 2-norm
       by assuming all weights equal their RMS value */
    N_VConst(ONE, firkls_mem->x);
    ewt_mean = N_VWrmsNorm(firk_mem->ewt, firkls_mem->x);
    delta /= ewt_mean;
  }

  N_VConst(ZERO, firkls_mem->x);
  retval = SUNLinSolSetZeroGuess(firkls_mem->LS, SUNTRUE);
  if (retval != SUN_SUCCESS) { return -1; }

  retval = SUNLinSolSolve(firkls_mem->LS, firkls_mem->M_lu, firkls_mem->x, b,
                          delta);
  N_VScale(ONE, firkls_mem->x, b);
  firkls_mem->nmsolves++;

  nli_inc = 0;
  if (firkls_mem->iterative && firkls_mem->LS->ops->numiters)
  {
    nli_inc = SUNLinSolNumIters(firkls_mem->LS);
  }
  firkls_mem->nmli += nli_inc;
  if (retval != SUN_SUCCESS) { firkls_mem->nmcfl++; }

  firkls_mem->last_flag = retval;

  switch (retval)
  {
  case SUN_SUCCESS: return 0;
  case SUNLS_RES_REDUCED: return accept_reduced ? 0 : 1;
  case SUNLS_RECOV_FAILURE:
  case SUNLS_CONV_FAIL:
  case SUNLS_ATIMES_FAIL_REC:
  case SUNLS_PSET_FAIL_REC:
  case SUNLS_PSOLVE_FAIL_REC:
  case SUNLS_PACKAGE_FAIL_REC:
  case SUNLS_QRFACT_FAIL:
  case SUNLS_LUFACT_FAIL: return 1;
  case SUN_ERR_ARG_CORRUPT:
  case SUN_ERR_ARG_INCOMPATIBLE:
  case SUN_ERR_MEM_FAIL:
  case SUNLS_ATIMES_NULL:
  case SUNLS_PSOLVE_NULL:
  case SUNLS_GS_FAIL:
  case SUNLS_QRSOL_FAIL: return -1;
  case SUN_ERR_EXT_FAIL:
    firkProcessError(firk_mem, SUN_ERR_EXT_FAIL, __LINE__, __func__, __FILE__,
                     "Failure in SUNLinSol external package");
    return -1;
  case SUNLS_ATIMES_FAIL_UNREC:
    firkProcessError(firk_mem, SUNLS_ATIMES_FAIL_UNREC, __LINE__, __func__,
                     __FILE__, MSG_LS_MTIMES_FAILED);
    return -1;
  case SUNLS_PSET_FAIL_UNREC:
    firkProcessError(firk_mem, SUNLS_PSET_FAIL_UNREC, __LINE__, __func__,
                     __FILE__, MSG_LS_PSET_FAILED);
    return -1;
  case SUNLS_PSOLVE_FAIL_UNREC:
    firkProcessError(firk_mem, SUNLS_PSOLVE_FAIL_UNREC, __LINE__, __func__,
                     __FILE__, MSG_LS_PSOLVE_FAILED);
    return -1;
  default:
    firkProcessError(firk_mem, retval, __LINE__, __func__, __FILE__,
                     "Unrecognized error return value from SUNLinSolSolve");
    return -1;
  }
}

/*-----------------------------------------------------------------
  firkLsMassSolve (see the comment above firkLsMassSolveBody)

  The profiling region and the log records of a mass solve around
  firkLsMassSolveBody (H31): a list region, since one region of the
  integrator may contain several mass solves (s in a conversion).
  -----------------------------------------------------------------*/
int firkLsMassSolve(FIRKodeMem firk_mem, N_Vector b, sunrealtype tol,
                    sunbooleantype accept_reduced)
{
  FIRKLsMassMem firkls_mem;
  sunrealtype delta;
  long int nmli0;
  int retval;

  retval = firkLs_AccessMassMem(firk_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  SUNDIALS_MARK_FUNCTION_BEGIN(FIRK_PROFILER);

  /* residual tolerance for iterative solvers (in the 2-norm) */
  delta = firkls_mem->iterative ? firkls_mem->eplifac * tol * firkls_mem->nrmfac
                                : ZERO;
  SUNLogInfo(FIRK_LOGGER, "begin-mass-linear-solve-list",
             "iterative = %i, b-tol = " SUN_FORMAT_G
             ", res-tol = " SUN_FORMAT_G,
             (int)firkls_mem->iterative, tol, delta);
  nmli0  = firkls_mem->nmli;
  retval = firkLsMassSolveBody(firk_mem, firkls_mem, b, delta, accept_reduced);
  SUNLogInfo(FIRK_LOGGER, "end-mass-linear-solve-list",
             "status = %s, retval = %i, last-flag = %li, iters = %li",
             (retval == 0) ? "success" : "failed", retval,
             (long int)firkls_mem->last_flag, firkls_mem->nmli - nmli0);
  SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
  (void)nmli0; /* used only by the log records */
  return retval;
}

/*---------------------------------------------------------------
  firkLsMassSolveDeriv

  The mass solve of the integrator (msolve): a derivative of y, for
  the initial step and rootfinding, which accepts a reduced residual.
  ---------------------------------------------------------------*/
int firkLsMassSolveDeriv(FIRKodeMem firk_mem, N_Vector b, sunrealtype tol)
{
  return firkLsMassSolve(firk_mem, b, tol, SUNTRUE);
}

/*---------------------------------------------------------------
  firkLsMassFree

  Frees the mass matrix interface memory and detaches it from the
  integrator.
  ---------------------------------------------------------------*/
int firkLsMassFree(FIRKodeMem firk_mem)
{
  FIRKLsMassMem firkls_mem;

  if (firk_mem == NULL) { return FIRKLS_SUCCESS; }
  if (firk_mem->mass_mem == NULL) { return FIRKLS_SUCCESS; }
  firkls_mem = (FIRKLsMassMem)firk_mem->mass_mem;

  /* detach the interface routines from the LS object */
  if (firkls_mem->LS && firkls_mem->LS->ops)
  {
    if (firkls_mem->LS->ops->setatimes)
    {
      SUNLinSolSetATimes(firkls_mem->LS, NULL, NULL);
    }
    if (firkls_mem->LS->ops->setpreconditioner)
    {
      SUNLinSolSetPreconditioner(firkls_mem->LS, NULL, NULL, NULL);
    }
  }

  if (firkls_mem->x)
  {
    N_VDestroy(firkls_mem->x);
    firkls_mem->x = NULL;
  }
  if (!firkls_mem->iterative && firkls_mem->M_lu)
  {
    SUNMatDestroy(firkls_mem->M_lu);
  }
  firkls_mem->M_lu = NULL;
  firkls_mem->M    = NULL;

  free(firk_mem->mass_mem);
  firk_mem->mass_mem = NULL;

  firk_mem->minit          = NULL;
  firk_mem->mreinit        = NULL;
  firk_mem->msetup         = NULL;
  firk_mem->mmult          = NULL;
  firk_mem->msolve         = NULL;
  firk_mem->mfree          = NULL;
  firk_mem->mass_set       = SUNFALSE;
  firk_mem->mass_iterative = SUNFALSE;

  return FIRKLS_SUCCESS;
}

/*---------------------------------------------------------------
  firkLs_AccessMassMem

  Unpacks the firk_mem and mass_mem structures from a void*
  pointer.
  ---------------------------------------------------------------*/
int firkLs_AccessMassMem(void* firkode_mem, const char* fname,
                         FIRKodeMem* firk_mem, FIRKLsMassMem* firkls_mem)
{
  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRKLS_MEM_NULL, __LINE__, fname, __FILE__,
                     MSG_LS_FIRKMEM_NULL);
    return FIRKLS_MEM_NULL;
  }
  *firk_mem = (FIRKodeMem)firkode_mem;
  if ((*firk_mem)->mass_mem == NULL)
  {
    firkProcessError(*firk_mem, FIRKLS_MASSMEM_NULL, __LINE__, fname, __FILE__,
                     MSG_LS_MASSMEM_NULL);
    return FIRKLS_MASSMEM_NULL;
  }
  *firkls_mem = (FIRKLsMassMem)(*firk_mem)->mass_mem;
  return FIRKLS_SUCCESS;
}
