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

#include <nvector/nvector_manyvector.h>
#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_spfgmr.h>
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
                              int retval, int curiter);

/* Purposes of a solve with the block matrix M - gamma J; they differ in the
   tolerance given to an iterative block solver and in how its failures are
   treated */
typedef enum
{
  FIRKLS_BLOCK_NEWTON,   /* Newton system of the single-stage method */
  FIRKLS_BLOCK_ESTIMATE, /* filtered error estimate */
  FIRKLS_BLOCK_PRECOND   /* preconditioner application in the stacked solve */
} firkLsBlockSolve;

static int firkLsSolveBlockImpl(FIRKodeMem firk_mem, N_Vector b,
                                N_Vector weight, firkLsBlockSolve kind);

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

  /* set defaults for Jacobian-related fields */
  if (A != NULL)
  {
    firkls_mem->jacDQ  = SUNTRUE;
    firkls_mem->jac    = firkLsDQJac;
    firkls_mem->J_data = firk_mem;
  }
  else
  {
    firkls_mem->jacDQ  = SUNFALSE;
    firkls_mem->jac    = NULL;
    firkls_mem->J_data = NULL;
  }

  firkls_mem->jtimesDQ = SUNTRUE;
  firkls_mem->jtsetup  = NULL;
  firkls_mem->jtimes   = firkLsDQJtimes;
  firkls_mem->jt_f     = firk_mem->f;
  firkls_mem->jt_data  = firk_mem;

  firkls_mem->user_linsys = SUNFALSE;
  firkls_mem->linsys      = firkLsLinSys;
  firkls_mem->A_data      = firk_mem;

  /* set defaults for preconditioner-related fields */
  firkls_mem->pset   = NULL;
  firkls_mem->psolve = NULL;
  firkls_mem->pfree  = NULL;
  firkls_mem->P_data = firk_mem->user_data;

  /* initialize counters */
  firkLsInitializeCounters(firkls_mem);

  /* set default values for the rest of the LS parameters */
  firkls_mem->msbj         = FIRKLS_MSBJ;
  firkls_mem->jbad         = SUNTRUE;
  firkls_mem->dgmax_jbad   = FIRKLS_DGMAX;
  firkls_mem->eplifac      = FIRKLS_EPLIN;
  firkls_mem->eplifac_stk  = FIRKLS_EPLIN;
  firkls_mem->maxl_stk     = 0; /* set from s in firkLsInitialize */
  firkls_mem->maxl_stk_set = SUNFALSE;
  firkls_mem->maxrs_stk    = FIRKLS_MAXRS;
  firkls_mem->last_flag    = FIRKLS_SUCCESS;

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
  firkls_mem->LS_stk  = NULL;
  firkls_mem->x_stk   = NULL;
  firkls_mem->Jv      = NULL;
  firkls_mem->s_alloc = 0;

  /* attach linear solver memory to integrator memory */
  firk_mem->lmem = firkls_mem;

  return FIRKLS_SUCCESS;
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
    firkls_mem->jacDQ  = SUNFALSE;
    firkls_mem->jac    = jac;
    firkls_mem->J_data = firk_mem->user_data;
  }
  else
  {
    firkls_mem->jacDQ  = SUNTRUE;
    firkls_mem->jac    = firkLsDQJac;
    firkls_mem->J_data = firk_mem;
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

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
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

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
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

  firkls_mem->pset   = pset;
  firkls_mem->psolve = psolve;

  if (firkls_mem->LS->ops->setpreconditioner == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "SUNLinearSolver object does not support user-supplied "
                     "preconditioning");
    return FIRKLS_ILL_INPUT;
  }

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
    firkls_mem->jtimesDQ = SUNTRUE;
    firkls_mem->jtsetup  = NULL;
    firkls_mem->jtimes   = firkLsDQJtimes;
    firkls_mem->jt_f     = firk_mem->f;
    firkls_mem->jt_data  = firk_mem;
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

  firkls_mem->jt_f = (jtimesRhsFn != NULL) ? jtimesRhsFn : firk_mem->f;

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
  }
  else
  {
    firkls_mem->user_linsys = SUNFALSE;
    firkls_mem->linsys      = firkLsLinSys;
    firkls_mem->A_data      = firk_mem;
  }

  return FIRKLS_SUCCESS;
}

int FIRKodeSetStageSolverMaxl(void* firkode_mem, int maxl)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (maxl <= 0)
  {
    firkls_mem->maxl_stk_set = SUNFALSE;
    firkls_mem->maxl_stk     = 0;
  }
  else
  {
    firkls_mem->maxl_stk_set = SUNTRUE;
    firkls_mem->maxl_stk     = maxl;
  }

  /* force re-creation of the stacked solver */
  firkLsFreeStacked(firkls_mem);

  return FIRKLS_SUCCESS;
}

int FIRKodeSetStageSolverMaxRestarts(void* firkode_mem, int maxrs)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  firkls_mem->maxrs_stk = (maxrs < 0) ? FIRKLS_MAXRS : maxrs;
  if (firkls_mem->LS_stk != NULL)
  {
    if (SUNLinSol_SPFGMRSetMaxRestarts(firkls_mem->LS_stk,
                                       firkls_mem->maxrs_stk) != SUN_SUCCESS)
    {
      firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__, __FILE__,
                       "Error in calling SUNLinSol_SPFGMRSetMaxRestarts");
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

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  if (eplifac < ZERO)
  {
    firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "eplifac < 0 illegal.");
    return FIRKLS_ILL_INPUT;
  }

  firkls_mem->eplifac_stk = (eplifac == ZERO) ? FIRKLS_EPLIN : eplifac;

  return FIRKLS_SUCCESS;
}

/*===============================================================
  Optional Get routines
  ===============================================================*/

int FIRKodeGetJac(void* firkode_mem, SUNMatrix* J)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *J = firkls_mem->savedJ;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetJacTime(void* firkode_mem, sunrealtype* t_J)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *t_J = firkls_mem->tnlj;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetJacNumSteps(void* firkode_mem, long int* nst_J)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nst_J = firkls_mem->nstlj;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumJacEvals(void* firkode_mem, long int* njevals)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *njevals = firkls_mem->nje;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumLinRhsEvals(void* firkode_mem, long int* nfevalsLS)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nfevalsLS = firkls_mem->nfeDQ;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumPrecEvals(void* firkode_mem, long int* npevals)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *npevals = firkls_mem->npe;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumPrecSolves(void* firkode_mem, long int* npsolves)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *npsolves = firkls_mem->nps;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumLinIters(void* firkode_mem, long int* nliters)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nliters = firkls_mem->nli;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumLinConvFails(void* firkode_mem, long int* nlcfails)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nlcfails = firkls_mem->ncfl;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumJTSetupEvals(void* firkode_mem, long int* njtsetups)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *njtsetups = firkls_mem->njtsetup;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumJtimesEvals(void* firkode_mem, long int* njvevals)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *njvevals = firkls_mem->njtimes;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumBlockSolves(void* firkode_mem, long int* nbsolves)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nbsolves = firkls_mem->nbs;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumStageLinIters(void* firkode_mem, long int* nsliters)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nsliters = firkls_mem->nsli;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumStageLinConvFails(void* firkode_mem, long int* nslcfails)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nslcfails = firkls_mem->nscf;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetLinSolveStats(void* firkode_mem, long int* njevals,
                            long int* nfevalsLS, long int* nliters,
                            long int* nlcfails, long int* npevals,
                            long int* npsolves, long int* njtsetups,
                            long int* njtimes)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  *njevals   = firkls_mem->nje;
  *nfevalsLS = firkls_mem->nfeDQ;
  *nliters   = firkls_mem->nli;
  *nlcfails  = firkls_mem->ncfl;
  *npevals   = firkls_mem->npe;
  *npsolves  = firkls_mem->nps;
  *njtsetups = firkls_mem->njtsetup;
  *njtimes   = firkls_mem->njtimes;

  return FIRKLS_SUCCESS;
}

int FIRKodeGetLastLinFlag(void* firkode_mem, long int* flag)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  int retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *flag = firkls_mem->last_flag;
  return FIRKLS_SUCCESS;
}

char* FIRKodeGetLinReturnFlagName(long int flag)
{
  char* name = (char*)malloc(30 * sizeof(char));
  if (name == NULL) { return NULL; }

  switch (flag)
  {
  case FIRKLS_SUCCESS: sprintf(name, "FIRKLS_SUCCESS"); break;
  case FIRKLS_MEM_NULL: sprintf(name, "FIRKLS_MEM_NULL"); break;
  case FIRKLS_LMEM_NULL: sprintf(name, "FIRKLS_LMEM_NULL"); break;
  case FIRKLS_ILL_INPUT: sprintf(name, "FIRKLS_ILL_INPUT"); break;
  case FIRKLS_MEM_FAIL: sprintf(name, "FIRKLS_MEM_FAIL"); break;
  case FIRKLS_PMEM_NULL: sprintf(name, "FIRKLS_PMEM_NULL"); break;
  case FIRKLS_MASSMEM_NULL: sprintf(name, "FIRKLS_MASSMEM_NULL"); break;
  case FIRKLS_JACFUNC_UNRECVR: sprintf(name, "FIRKLS_JACFUNC_UNRECVR"); break;
  case FIRKLS_JACFUNC_RECVR: sprintf(name, "FIRKLS_JACFUNC_RECVR"); break;
  case FIRKLS_MASSFUNC_UNRECVR: sprintf(name, "FIRKLS_MASSFUNC_UNRECVR"); break;
  case FIRKLS_MASSFUNC_RECVR: sprintf(name, "FIRKLS_MASSFUNC_RECVR"); break;
  case FIRKLS_SUNMAT_FAIL: sprintf(name, "FIRKLS_SUNMAT_FAIL"); break;
  case FIRKLS_SUNLS_FAIL: sprintf(name, "FIRKLS_SUNLS_FAIL"); break;
  default: sprintf(name, "NONE");
  }

  return name;
}

/*=================================================================
  FIRKLS private functions
  =================================================================*/

/*-----------------------------------------------------------------
  firkLsATimes

  Generates the block matrix-vector product z = (M - gamma J) v for
  the user's iterative block solver. J v is obtained from the
  jtimes routine.
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
    N_VLinearSum(ONE, firkls_mem->ytemp, -firk_mem->gamma, z, z);
  }
  else { N_VLinearSum(ONE, v, -firk_mem->gamma, z, z); }

  return 0;
}

/*---------------------------------------------------------------
  firkLsPSetup

  Interfaces the user's psetup routine with the block solver.
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
  return retval;
}

/*-----------------------------------------------------------------
  firkLsPSolve

  Interfaces the user's psolve routine with the block solver.
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
                              r, z, firk_mem->gamma, tol, lr, firkls_mem->P_data);
  firkls_mem->nps++;
  return retval;
}

/*-----------------------------------------------------------------
  firkLsStackedATimes

  Generates the stacked stage-system matrix-vector product
    z_i = M v_i - h sum_j a_ij J v_j.
  J v_j is computed with SUNMatMatvec on the saved Jacobian when
  available, and with the jtimes routine otherwise.
  -----------------------------------------------------------------*/
int firkLsStackedATimes(void* firkode_mem, N_Vector v_stk, N_Vector z_stk)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  FIRKodeTable T;
  N_Vector v_i, z_i, mv;
  int i, j, s, retval;
  sunbooleantype use_matvec;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  T = firk_mem->T;
  s = firk_mem->s;

  use_matvec = (firkls_mem->savedJ != NULL) && firkls_mem->jtimesDQ &&
               (firkls_mem->savedJ->ops->matvec != NULL);

  /* Jacobian-vector products for each stage */
  for (j = 0; j < s; j++)
  {
    v_i = N_VGetSubvector_ManyVector(v_stk, (sunindextype)j);
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
    v_i = N_VGetSubvector_ManyVector(v_stk, (sunindextype)i);
    z_i = N_VGetSubvector_ManyVector(z_stk, (sunindextype)i);

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
      firk_mem->cvals[j + 1] = -firk_mem->h * T->A[i][j];
      firk_mem->Xvecs[j + 1] = firkls_mem->Jv[j];
    }
    retval = N_VLinearCombination(s + 1, firk_mem->cvals, firk_mem->Xvecs, z_i);
    if (retval != 0) { return -1; }
  }

  return 0;
}

/*-----------------------------------------------------------------
  firkLsStackedPSolve

  Block-diagonal preconditioner for the stacked stage system. With a
  matrix-based block solver each stage block is solved (inexactly)
  with the user's linear solver on (M - gamma J). With a matrix-free
  iterative block solver the user's preconditioner is applied to each
  stage block instead (the identity if none was supplied): nesting a
  Krylov iteration inside the preconditioner would make its
  convergence test act on arbitrarily scaled Krylov vectors, and the
  outer flexible GMRES iteration performs the Krylov work anyway.
  -----------------------------------------------------------------*/
int firkLsStackedPSolve(void* firkode_mem, N_Vector r_stk, N_Vector z_stk,
                        sunrealtype tol, SUNDIALS_MAYBE_UNUSED int lr)
{
  FIRKodeMem firk_mem;
  FIRKLsMem firkls_mem;
  N_Vector r_i, z_i;
  int i, retval;

  retval = firkLs_AccessLMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  for (i = 0; i < firk_mem->s; i++)
  {
    r_i = N_VGetSubvector_ManyVector(r_stk, (sunindextype)i);
    z_i = N_VGetSubvector_ManyVector(z_stk, (sunindextype)i);

    if (firkls_mem->matrixbased)
    {
      N_VScale(ONE, r_i, z_i);
      retval = firkLsSolveBlockImpl(firk_mem, z_i, firk_mem->ewt,
                                    FIRKLS_BLOCK_PRECOND);
      if (retval != 0) { return retval; }
    }
    else if (firkls_mem->psolve != NULL)
    {
      retval = firkls_mem->psolve(firk_mem->tn, firkls_mem->ycur,
                                  firkls_mem->fcur, r_i, z_i, firk_mem->gamma,
                                  tol, SUN_PREC_LEFT, firkls_mem->P_data);
      firkls_mem->nps++;
      if (retval != 0) { return retval; }
    }
    else { N_VScale(ONE, r_i, z_i); }
  }

  return 0;
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

  Dense difference quotient approximation to the Jacobian of f.
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
    if (retval != 0) { break; }

    y_data[j] = yjsaved;

    inc_inv = ONE / inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fy, jthCol);
  }

  N_VSetArrayPointer(NULL, jthCol);
  N_VDestroy(jthCol);

  return retval;
}

/*-----------------------------------------------------------------
  firkLsBandDQJac

  Banded difference quotient approximation to the Jacobian of f.
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
    Jv = [f(y + v sig) - f(y)] / sig,  sig = 1 / ||v||_WRMS.
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

  sig = ONE / N_VWrmsNorm(v, firk_mem->ewt);

  for (iter = 0; iter < MAX_DQITERS; iter++)
  {
    N_VLinearSum(sig, v, ONE, y, work);

    retval = firkls_mem->jt_f(t, work, Jv, firk_mem->user_data);
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

  if (jok)
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

  if (firkls_mem->LS_stk != NULL && firkls_mem->s_alloc == s)
  {
    return FIRKLS_SUCCESS;
  }
  firkLsFreeStacked(firkls_mem);

  if (s < 2)
  {
    firkls_mem->s_alloc = s;
    return FIRKLS_SUCCESS;
  }

  maxl = firkls_mem->maxl_stk_set ? firkls_mem->maxl_stk : SUNMAX(5, 3 * s);
  firkls_mem->maxl_stk = maxl;

  firkls_mem->LS_stk = SUNLinSol_SPFGMR(firk_mem->zpred_stk, SUN_PREC_RIGHT,
                                        maxl, firk_mem->sunctx);
  if (firkls_mem->LS_stk == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSG_LS_MEM_FAIL);
    return FIRKLS_MEM_FAIL;
  }

  retval = SUNLinSol_SPFGMRSetMaxRestarts(firkls_mem->LS_stk,
                                          firkls_mem->maxrs_stk);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__, __FILE__,
                     "Error in calling SUNLinSol_SPFGMRSetMaxRestarts");
    return FIRKLS_SUNLS_FAIL;
  }

  retval = SUNLinSolSetATimes(firkls_mem->LS_stk, firk_mem, firkLsStackedATimes);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__, __FILE__,
                     "Error in calling SUNLinSolSetATimes");
    return FIRKLS_SUNLS_FAIL;
  }

  retval = SUNLinSolSetPreconditioner(firkls_mem->LS_stk, firk_mem, NULL,
                                      firkLsStackedPSolve);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__, __FILE__,
                     "Error in calling SUNLinSolSetPreconditioner");
    return FIRKLS_SUNLS_FAIL;
  }

  retval = SUNLinSolInitialize(firkls_mem->LS_stk);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__, __FILE__,
                     "Error in calling SUNLinSolInitialize");
    return FIRKLS_SUNLS_FAIL;
  }

  firkls_mem->x_stk = N_VClone(firk_mem->zpred_stk);
  firkls_mem->Jv    = N_VCloneVectorArray(s, firkls_mem->ytemp);
  if (firkls_mem->x_stk == NULL || firkls_mem->Jv == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSG_LS_MEM_FAIL);
    return FIRKLS_MEM_FAIL;
  }

  firkls_mem->nrmfac_stk =
    SUNRsqrt((sunrealtype)N_VGetLength(firk_mem->zpred_stk));
  firkls_mem->s_alloc = s;

  return FIRKLS_SUCCESS;
}

static void firkLsFreeStacked(FIRKLsMem firkls_mem)
{
  if (firkls_mem->LS_stk != NULL)
  {
    SUNLinSolFree(firkls_mem->LS_stk);
    firkls_mem->LS_stk = NULL;
  }
  if (firkls_mem->x_stk != NULL)
  {
    N_VDestroy(firkls_mem->x_stk);
    firkls_mem->x_stk = NULL;
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

      if (firkls_mem->jacDQ)
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

      /* the stacked operator needs J v: via SUNMatMatvec or a jtimes fn */
      if (firk_mem->s > 1 && firkls_mem->jtimesDQ &&
          firkls_mem->A->ops->matvec == NULL)
      {
        firkProcessError(firk_mem, FIRKLS_ILL_INPUT, __LINE__, __func__, __FILE__,
                         "The SUNMatrix must support SUNMatMatvec, or a "
                         "Jacobian-vector product routine must be supplied");
        firkls_mem->last_flag = FIRKLS_ILL_INPUT;
        return FIRKLS_ILL_INPUT;
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
      }
    }
  }
  else
  {
    /* matrix-free case */
    firkls_mem->jacDQ  = SUNFALSE;
    firkls_mem->jac    = NULL;
    firkls_mem->J_data = NULL;

    firkls_mem->user_linsys = SUNFALSE;
    firkls_mem->linsys      = NULL;
    firkls_mem->A_data      = NULL;
  }

  /* reset counters */
  firkLsInitializeCounters(firkls_mem);

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
  firkLsSetup

  Conditionally updates the Jacobian, forms A = M - gamma J and
  calls the block solver's setup routine.
  -----------------------------------------------------------------*/
int firkLsSetup(FIRKodeMem firk_mem, int convfail, N_Vector ypred,
                N_Vector fpred, sunbooleantype* jcurPtr, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3)
{
  FIRKLsMem firkls_mem;
  FIRKLsMassMem firkls_massmem;
  SUNMatrix M;
  sunrealtype dgamma;
  int retval;

  if (firk_mem->lmem == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_LMEM_NULL, __LINE__, __func__, __FILE__,
                     MSG_LS_LMEM_NULL);
    return FIRKLS_LMEM_NULL;
  }
  firkls_mem = (FIRKLsMem)firk_mem->lmem;

  if (SUNLinSolGetType(firkls_mem->LS) == SUNLINEARSOLVER_MATRIX_EMBEDDED)
  {
    firkls_mem->last_flag = FIRKLS_SUCCESS;
    return firkls_mem->last_flag;
  }

  /* linearization point */
  firkls_mem->ycur = ypred;
  firkls_mem->fcur = fpred;

  /* use nst, gamma/gammap, and convfail to set the J/P evaluation flag */
  dgamma           = (firk_mem->gammap != ZERO)
                       ? SUNRabs((firk_mem->gamma / firk_mem->gammap) - ONE)
                       : ONE;
  firkls_mem->jbad = (firk_mem->firststage) ||
                     (firk_mem->nst >= firkls_mem->nstlj + firkls_mem->msbj) ||
                     ((convfail == FIRK_FAIL_BAD_J) &&
                      (dgamma < firkls_mem->dgmax_jbad)) ||
                     (convfail == FIRK_FAIL_OTHER);

  if (firkls_mem->A != NULL)
  {
    M = NULL;
    if (firk_mem->mass_set && firk_mem->mass_mem != NULL)
    {
      firkls_massmem = (FIRKLsMassMem)firk_mem->mass_mem;
      M              = firkls_massmem->M;
    }

    retval = firkls_mem->linsys(firk_mem->tn, ypred, fpred, firkls_mem->A, M,
                                !(firkls_mem->jbad), jcurPtr, firk_mem->gamma,
                                firkls_mem->A_data, vtemp1, vtemp2, vtemp3);

    if (*jcurPtr)
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
    /* matrix-free case, set jcur to jbad */
    *jcurPtr = firkls_mem->jbad;
  }

  /* call the block solver's setup routine */
  firkls_mem->last_flag = SUNLinSolSetup(firkls_mem->LS, firkls_mem->A);

  if (firkls_mem->A == NULL)
  {
    if (*jcurPtr)
    {
      firkls_mem->npe++;
      firkls_mem->nstlj = firk_mem->nst;
      firkls_mem->tnlj  = firk_mem->tn;
    }
    if (firkls_mem->jbad) { *jcurPtr = SUNTRUE; }
  }

  return firkls_mem->last_flag;
}

/*-----------------------------------------------------------------
  firkLsMapSolveFlag

  Translates a SUNLinSolSolve return value into the 0/+1/-1
  convention of the integrator.
  -----------------------------------------------------------------*/
static int firkLsMapSolveFlag(FIRKodeMem firk_mem, FIRKLsMem firkls_mem,
                              int retval, int curiter)
{
  firkls_mem->last_flag = retval;

  switch (retval)
  {
  case SUN_SUCCESS: return 0;
  case SUNLS_RES_REDUCED:
    /* allow reduction but not solution on the first Newton iteration */
    return (curiter == 0) ? 0 : 1;
  case SUNLS_CONV_FAIL:
  case SUNLS_ATIMES_FAIL_REC:
  case SUNLS_PSOLVE_FAIL_REC:
  case SUNLS_PACKAGE_FAIL_REC:
  case SUNLS_QRFACT_FAIL:
  case SUNLS_LUFACT_FAIL: return 1;
  case SUN_ERR_ARG_CORRUPT:
  case SUN_ERR_ARG_INCOMPATIBLE:
  case SUN_ERR_MEM_FAIL:
  case SUNLS_GS_FAIL:
  case SUNLS_QRSOL_FAIL: return -1;
  case SUN_ERR_EXT_FAIL:
    firkProcessError(firk_mem, SUN_ERR_EXT_FAIL, __LINE__, __func__, __FILE__,
                     "Failure in SUNLinSol external package");
    return -1;
  case SUNLS_ATIMES_FAIL_UNREC:
    firkProcessError(firk_mem, SUNLS_ATIMES_FAIL_UNREC, __LINE__, __func__,
                     __FILE__, MSG_LS_JTIMES_FAILED);
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
  firkLsSolveBlock

  Solves the block system (M - gamma J) x = b of the filtered error
  estimate in place with the user's linear solver, setting
  tolerances and scaling vectors for iterative solvers.
  -----------------------------------------------------------------*/
int firkLsSolveBlock(FIRKodeMem firk_mem, N_Vector b, N_Vector weight)
{
  return firkLsSolveBlockImpl(firk_mem, b, weight, FIRKLS_BLOCK_ESTIMATE);
}

/* For a preconditioner application in the stacked Krylov iteration, an
   iterative block solver that only reduced the residual (or ran out of
   iterations) still provides a useful approximate solution, which the
   flexible GMRES iteration can use, so such outcomes are not reported as
   failures. */
static int firkLsSolveBlockImpl(FIRKodeMem firk_mem, N_Vector b,
                                N_Vector weight, firkLsBlockSolve kind)
{
  FIRKLsMem firkls_mem;
  sunrealtype bnorm = ZERO;
  sunrealtype deltar, delta, w_mean;
  int curiter, nli_inc, retval;

  if (firk_mem->lmem == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_LMEM_NULL, __LINE__, __func__, __FILE__,
                     MSG_LS_LMEM_NULL);
    return FIRKLS_LMEM_NULL;
  }
  firkls_mem = (FIRKLsMem)firk_mem->lmem;

  curiter = 0;
  if (firk_mem->NLS != NULL)
  {
    (void)SUNNonlinSolGetCurIter(firk_mem->NLS, &curiter);
  }

  /* iterative block solver: tolerance and small right-hand side shortcut */
  if (firkls_mem->iterative)
  {
    bnorm = N_VWrmsNorm(b, weight);

    if (kind == FIRKLS_BLOCK_NEWTON)
    {
      /* Newton system: tolerance relative to the nonlinear solver tolerance,
         as in CVODE; a smaller correction is negligible after the first
         iteration */
      deltar = firkls_mem->eplifac * firk_mem->nlscoef;

      SUNLogInfo(FIRK_LOGGER, "begin-block-linear-solve",
                 "iterative = 1, b-norm = " SUN_FORMAT_G
                 ", b-tol = " SUN_FORMAT_G ", res-tol = " SUN_FORMAT_G,
                 bnorm, deltar, deltar * firkls_mem->nrmfac);

      if (bnorm <= deltar)
      {
        if (curiter > 0) { N_VConst(ZERO, b); }
        firkls_mem->last_flag = FIRKLS_SUCCESS;
        firkls_mem->nbs++;
        SUNLogInfo(FIRK_LOGGER, "end-block-linear-solve",
                   "status = success small rhs");
        return firkls_mem->last_flag;
      }
    }
    else
    {
      /* error estimate: the solution is the estimate itself, so the solve is
         not skipped for a small right-hand side, and the tolerance of the
         Newton system is tightened to one relative to b when b is smaller
         than the nonlinear solver tolerance, so that a small estimate is
         resolved instead of returned as zero. Preconditioner application:
         the right-hand side is a Krylov basis vector of arbitrary scale, so
         the tolerance is relative to it. */
      deltar = (kind == FIRKLS_BLOCK_ESTIMATE)
                 ? firkls_mem->eplifac * SUNMIN(firk_mem->nlscoef, bnorm)
                 : firkls_mem->eplifac * bnorm;

      SUNLogInfo(FIRK_LOGGER, "begin-block-linear-solve",
                 "iterative = 1, %s, b-norm = " SUN_FORMAT_G
                 ", res-tol = " SUN_FORMAT_G,
                 (kind == FIRKLS_BLOCK_ESTIMATE) ? "error estimate"
                                                 : "preconditioner application",
                 bnorm, deltar * firkls_mem->nrmfac);

      if (bnorm == ZERO)
      {
        firkls_mem->last_flag = FIRKLS_SUCCESS;
        firkls_mem->nbs++;
        return firkls_mem->last_flag;
      }
    }
    delta = deltar * firkls_mem->nrmfac;
  }
  else
  {
    delta = ZERO;
    SUNLogInfo(FIRK_LOGGER, "begin-block-linear-solve", "iterative = 0");
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

  N_VConst(ZERO, firkls_mem->x);
  retval = SUNLinSolSetZeroGuess(firkls_mem->LS, SUNTRUE);
  if (retval != SUN_SUCCESS) { return -1; }

  retval = SUNLinSolSolve(firkls_mem->LS, firkls_mem->A, firkls_mem->x, b, delta);
  N_VScale(ONE, firkls_mem->x, b);
  firkls_mem->nbs++;

  nli_inc = 0;
  if (firkls_mem->iterative && firkls_mem->LS->ops->numiters)
  {
    nli_inc = SUNLinSolNumIters(firkls_mem->LS);
  }
  firkls_mem->nli += nli_inc;

  if (kind == FIRKLS_BLOCK_PRECOND &&
      (retval == SUNLS_RES_REDUCED || retval == SUNLS_CONV_FAIL))
  {
    firkls_mem->last_flag = retval;
    SUNLogInfo(FIRK_LOGGER, "end-block-linear-solve",
               "status = inexact preconditioner application, retval = %i, "
               "iters = %i",
               retval, nli_inc);
    return 0;
  }
  if (retval != SUN_SUCCESS) { firkls_mem->ncfl++; }

  SUNLogInfoIf(retval == SUN_SUCCESS, FIRK_LOGGER, "end-block-linear-solve",
               "status = success, iters = %i", nli_inc);
  SUNLogInfoIf(retval != SUN_SUCCESS, FIRK_LOGGER, "end-block-linear-solve",
               "status = failed, retval = %i, iters = %i", retval, nli_inc);

  return firkLsMapSolveFlag(firk_mem, firkls_mem, retval, curiter);
}

/*-----------------------------------------------------------------
  firkLsSolveStacked

  Solves the stacked Newton system (I_s (x) M - h A (x) J) x = b in
  place with the internal flexible GMRES iteration, preconditioned
  block-wise by the user's linear solver. For s = 1 the block solve
  is exact and used directly.
  -----------------------------------------------------------------*/
int firkLsSolveStacked(FIRKodeMem firk_mem, N_Vector b_stk)
{
  FIRKLsMem firkls_mem;
  sunrealtype bnorm, deltar, delta;
  int curiter, nli_inc, retval;

  if (firk_mem->lmem == NULL)
  {
    firkProcessError(firk_mem, FIRKLS_LMEM_NULL, __LINE__, __func__, __FILE__,
                     MSG_LS_LMEM_NULL);
    return FIRKLS_LMEM_NULL;
  }
  firkls_mem = (FIRKLsMem)firk_mem->lmem;

  /* single stage: the block solve is the exact Newton solve */
  if (firk_mem->s == 1)
  {
    return firkLsSolveBlockImpl(firk_mem, N_VGetSubvector_ManyVector(b_stk, 0),
                                firk_mem->ewt, FIRKLS_BLOCK_NEWTON);
  }

  curiter = 0;
  (void)SUNNonlinSolGetCurIter(firk_mem->NLS, &curiter);

  deltar = firkls_mem->eplifac_stk * firk_mem->nlscoef;
  bnorm  = N_VWrmsNorm(b_stk, firk_mem->ewt_stk);

  SUNLogInfo(FIRK_LOGGER, "begin-stage-linear-solve",
             "b-norm = " SUN_FORMAT_G ", b-tol = " SUN_FORMAT_G
             ", res-tol = " SUN_FORMAT_G,
             bnorm, deltar, deltar * firkls_mem->nrmfac_stk);

  if (bnorm <= deltar)
  {
    if (curiter > 0) { N_VConst(ZERO, b_stk); }
    firkls_mem->last_flag = FIRKLS_SUCCESS;
    SUNLogInfo(FIRK_LOGGER, "end-stage-linear-solve",
               "status = success small rhs");
    return firkls_mem->last_flag;
  }
  delta = deltar * firkls_mem->nrmfac_stk;

  /* linearization point for ATimes/PSolve */
  firkls_mem->ycur = firk_mem->yn;
  firkls_mem->fcur = firk_mem->fn;

  retval = SUNLinSolSetScalingVectors(firkls_mem->LS_stk, firk_mem->ewt_stk,
                                      firk_mem->ewt_stk);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRKLS_SUNLS_FAIL, __LINE__, __func__, __FILE__,
                     "Error in calling SUNLinSolSetScalingVectors");
    firkls_mem->last_flag = FIRKLS_SUNLS_FAIL;
    return firkls_mem->last_flag;
  }

  N_VConst(ZERO, firkls_mem->x_stk);
  retval = SUNLinSolSetZeroGuess(firkls_mem->LS_stk, SUNTRUE);
  if (retval != SUN_SUCCESS) { return -1; }

  /* user-provided jtsetup routine, called once per Newton solve */
  if (firkls_mem->jtsetup)
  {
    firkls_mem->last_flag = firkls_mem->jtsetup(firk_mem->tn, firk_mem->yn,
                                                firk_mem->fn,
                                                firkls_mem->jt_data);
    firkls_mem->njtsetup++;
    if (firkls_mem->last_flag != 0)
    {
      firkProcessError(firk_mem, firkls_mem->last_flag, __LINE__, __func__,
                       __FILE__, MSG_LS_JTSETUP_FAILED);
      return firkls_mem->last_flag;
    }
  }

  retval = SUNLinSolSolve(firkls_mem->LS_stk, NULL, firkls_mem->x_stk, b_stk,
                          delta);
  N_VScale(ONE, firkls_mem->x_stk, b_stk);

  nli_inc = SUNLinSolNumIters(firkls_mem->LS_stk);
  firkls_mem->nsli += nli_inc;
  if (retval != SUN_SUCCESS) { firkls_mem->nscf++; }

  SUNLogInfoIf(retval == SUN_SUCCESS, FIRK_LOGGER, "end-stage-linear-solve",
               "status = success, iters = %i, res-norm = " SUN_FORMAT_G,
               nli_inc, SUNLinSolResNorm(firkls_mem->LS_stk));
  SUNLogInfoIf(retval != SUN_SUCCESS, FIRK_LOGGER,
               "end-stage-linear-solve", "status = failed, retval = %i, iters = %i, res-norm = " SUN_FORMAT_G,
               retval, nli_inc, SUNLinSolResNorm(firkls_mem->LS_stk));

  return firkLsMapSolveFlag(firk_mem, firkls_mem, retval, curiter);
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

  if (firkls_mem->pfree) { firkls_mem->pfree(firk_mem); }

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
  firkls_mem->nje      = 0;
  firkls_mem->nfeDQ    = 0;
  firkls_mem->nstlj    = 0;
  firkls_mem->npe      = 0;
  firkls_mem->nli      = 0;
  firkls_mem->nps      = 0;
  firkls_mem->ncfl     = 0;
  firkls_mem->njtsetup = 0;
  firkls_mem->njtimes  = 0;
  firkls_mem->nbs      = 0;
  firkls_mem->nsli     = 0;
  firkls_mem->nscf     = 0;
  firkls_mem->tnlj     = ZERO;
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
  firkls_mem->pfree  = NULL;
  firkls_mem->P_data = firk_mem->user_data;

  firkls_mem->eplifac   = FIRKLS_EPLIN;
  firkls_mem->nrmfac    = ZERO;
  firkls_mem->last_flag = FIRKLS_SUCCESS;

  firkls_mem->nmsetups  = 0;
  firkls_mem->nmvsetups = 0;
  firkls_mem->nmvevals  = 0;
  firkls_mem->nmsolves  = 0;
  firkls_mem->nmpe      = 0;
  firkls_mem->nmli      = 0;
  firkls_mem->nmps      = 0;
  firkls_mem->nmcfails  = 0;

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
  firk_mem->minit    = firkLsMassInitialize;
  firk_mem->msetup   = firkLsMassSetup;
  firk_mem->mmult    = firkLsMassMult;
  firk_mem->msolve   = firkLsMassSolve;
  firk_mem->mfree    = firkLsMassFree;
  firk_mem->mass_mem = firkls_mem;
  firk_mem->mass_set = SUNTRUE;

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

  firkls_mem->mass   = mass;
  firkls_mem->M_data = firk_mem->user_data;

  return FIRKLS_SUCCESS;
}

int FIRKodeSetMassTimes(void* firkode_mem, FIRKLsMassTimesSetupFn mtsetup,
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

  firkls_mem->mtsetup = mtsetup;
  firkls_mem->mtimes  = mtimes;
  firkls_mem->mt_data = mtimes_data;

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

  firkls_mem->pset   = psetup;
  firkls_mem->psolve = psolve;
  firkls_mem->P_data = firk_mem->user_data;

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

int FIRKodeGetCurrentMassMatrix(void* firkode_mem, SUNMatrix* M)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *M = firkls_mem->M;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumMassSetups(void* firkode_mem, long int* nmsetups)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nmsetups = firkls_mem->nmsetups;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumMassMultSetups(void* firkode_mem, long int* nmvsetups)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nmvsetups = firkls_mem->nmvsetups;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumMassMult(void* firkode_mem, long int* nmvevals)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nmvevals = firkls_mem->nmvevals;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumMassSolves(void* firkode_mem, long int* nmsolves)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nmsolves = firkls_mem->nmsolves;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumMassPrecEvals(void* firkode_mem, long int* nmpevals)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nmpevals = firkls_mem->nmpe;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumMassPrecSolves(void* firkode_mem, long int* nmpsolves)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nmpsolves = firkls_mem->nmps;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumMassIters(void* firkode_mem, long int* nmiters)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nmiters = firkls_mem->nmli;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetNumMassConvFails(void* firkode_mem, long int* nmcfails)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *nmcfails = firkls_mem->nmcfails;
  return FIRKLS_SUCCESS;
}

int FIRKodeGetLastMassFlag(void* firkode_mem, long int* flag)
{
  FIRKodeMem firk_mem;
  FIRKLsMassMem firkls_mem;
  int retval;

  retval = firkLs_AccessMassMem(firkode_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }
  *flag = firkls_mem->last_flag;
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
    firkls_mem->nmvevals++;
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
    firkls_mem->nmvevals++;
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

  firkls_mem->nmsetups  = 0;
  firkls_mem->nmvsetups = 0;
  firkls_mem->nmvevals  = 0;
  firkls_mem->nmsolves  = 0;
  firkls_mem->nmpe      = 0;
  firkls_mem->nmli      = 0;
  firkls_mem->nmps      = 0;
  firkls_mem->nmcfails  = 0;

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

  firkls_mem->last_flag = SUNLinSolInitialize(firkls_mem->LS);
  return firkls_mem->last_flag;
}

/*---------------------------------------------------------------
  firkLsMassSetup

  Constructs the (constant) mass matrix and sets up the mass
  linear solver. Called once from firkInitialSetup; subsequent
  calls return immediately.
  ---------------------------------------------------------------*/
int firkLsMassSetup(FIRKodeMem firk_mem, sunrealtype t, N_Vector vtemp1,
                    N_Vector vtemp2, N_Vector vtemp3)
{
  FIRKLsMassMem firkls_mem;
  sunbooleantype call_lssetup;
  int retval;

  retval = firkLs_AccessMassMem(firk_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  /* the mass matrix is constant: set up at most once */
  if (firkls_mem->nmsetups > 0 || firkls_mem->nmvsetups > 0)
  {
    firkls_mem->last_flag = FIRKLS_SUCCESS;
    return firkls_mem->last_flag;
  }

  /* user-supplied mass-times setup routine */
  if (firkls_mem->mtsetup)
  {
    firkls_mem->last_flag = firkls_mem->mtsetup(t, firkls_mem->mt_data);
    firkls_mem->nmvsetups++;
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
      firkls_mem->nmvsetups++;
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

  return firkls_mem->last_flag;
}

/*---------------------------------------------------------------
  firkLsMassMult

  Computes Mv = M v for the integrator (residuals, error estimate
  and stacked matrix-vector products).
  ---------------------------------------------------------------*/
int firkLsMassMult(FIRKodeMem firk_mem, N_Vector v, N_Vector Mv)
{
  return firkLsMTimes(firk_mem, v, Mv);
}

/*---------------------------------------------------------------
  firkLsMassSolve

  Solves M x = b in place: sets the tolerance and scaling vectors,
  calls the mass linear solver and accumulates statistics. The
  return value is 0 on success, 1 for a recoverable failure and
  -1 for an unrecoverable failure.
  ---------------------------------------------------------------*/
int firkLsMassSolve(FIRKodeMem firk_mem, N_Vector b, sunrealtype tol)
{
  FIRKLsMassMem firkls_mem;
  sunrealtype delta, ewt_mean;
  int nli_inc, retval;

  retval = firkLs_AccessMassMem(firk_mem, __func__, &firk_mem, &firkls_mem);
  if (retval != FIRKLS_SUCCESS) { return retval; }

  /* residual tolerance for iterative solvers (in the 2-norm) */
  delta = firkls_mem->iterative ? firkls_mem->eplifac * tol * firkls_mem->nrmfac
                                : ZERO;

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
  if (retval != SUN_SUCCESS) { firkls_mem->nmcfails++; }

  firkls_mem->last_flag = retval;

  switch (retval)
  {
  case SUN_SUCCESS: return 0;
  case SUNLS_RES_REDUCED:
  case SUNLS_CONV_FAIL:
  case SUNLS_ATIMES_FAIL_REC:
  case SUNLS_PSOLVE_FAIL_REC:
  case SUNLS_PACKAGE_FAIL_REC:
  case SUNLS_QRFACT_FAIL:
  case SUNLS_LUFACT_FAIL: return 1;
  case SUN_ERR_ARG_CORRUPT:
  case SUN_ERR_ARG_INCOMPATIBLE:
  case SUN_ERR_MEM_FAIL:
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
  firkls_mem->ycur = NULL;

  if (firkls_mem->pfree) { firkls_mem->pfree(firk_mem); }

  free(firk_mem->mass_mem);
  firk_mem->mass_mem = NULL;

  firk_mem->minit    = NULL;
  firk_mem->msetup   = NULL;
  firk_mem->mmult    = NULL;
  firk_mem->msolve   = NULL;
  firk_mem->mfree    = NULL;
  firk_mem->mass_set = SUNFALSE;

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
