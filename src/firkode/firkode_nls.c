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
 * Interface between FIRKODE and the SUNNonlinearSolver Newton
 * iteration on the stacked stage system.
 *
 * With the stage increments z_i = Y_i - y_n stacked in a ManyVector
 * Z, the nonlinear system is
 *   G(Z)_i = M z_i - h sum_j a_ij f(t_n + c_j h, y_n + z_j) = 0.
 * The Newton solver iterates on the correction Zcor to the predicted
 * increments Zpred. The convergence test follows RADAU5: with the
 * contraction estimate theta_k = ||dZ_k|| / ||dZ_{k-1}||, the
 * iteration is accepted when theta/(1-theta) ||dZ_k|| <= tol, and
 * abandoned early when it diverges or is predicted not to converge
 * within the allowed number of iterations.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_manyvector.h>

#include "firkode_impl.h"
#include "sundials_macros.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* private functions */
static int firkNlsResidual(N_Vector ycor, N_Vector res, void* firkode_mem);
static int firkNlsLSetup(sunbooleantype jbad, sunbooleantype* jcur,
                         void* firkode_mem);
static int firkNlsLSolve(N_Vector delta, void* firkode_mem);
static int firkNlsConvTest(SUNNonlinearSolver NLS, N_Vector ycor, N_Vector del,
                           sunrealtype tol, N_Vector ewt, void* firkode_mem);
static SUNErrCode firkNlsNorm(N_Vector delta, N_Vector ewt, sunrealtype* delnrm,
                              void* firkode_mem);
static SUNErrCode firkNlsGetUpdateNorm(sunrealtype* delnrm, void* firkode_mem);
static SUNErrCode firkNlsGetConvRate(sunrealtype* crate, void* firkode_mem);

/* -----------------------------------------------------------------------------
 * Package-internal functions
 * ---------------------------------------------------------------------------*/

/*---------------------------------------------------------------
  firkNlsAttach

  Attaches a (root-finding type) nonlinear solver operating on the
  stacked stage vector to the FIRKODE memory.
  ---------------------------------------------------------------*/
int firkNlsAttach(FIRKodeMem firk_mem, SUNNonlinearSolver NLS)
{
  int retval;

  if (NLS == NULL)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "NLS must be non-NULL");
    return FIRK_ILL_INPUT;
  }

  if (NLS->ops->gettype == NULL || NLS->ops->solve == NULL ||
      NLS->ops->setsysfn == NULL)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "NLS does not support required operations");
    return FIRK_ILL_INPUT;
  }

  if (SUNNonlinSolGetType(NLS) != SUNNONLINEARSOLVER_ROOTFIND)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "FIRKODE requires a root-finding type nonlinear solver");
    return FIRK_ILL_INPUT;
  }

  /* free any existing nonlinear solver */
  if ((firk_mem->NLS != NULL) && (firk_mem->ownNLS))
  {
    (void)SUNNonlinSolFree(firk_mem->NLS);
  }

  firk_mem->NLS    = NLS;
  firk_mem->ownNLS = SUNFALSE;

  retval = SUNNonlinSolSetSysFn(NLS, firkNlsResidual);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Setting nonlinear system function failed");
    return FIRK_ILL_INPUT;
  }

  retval = SUNNonlinSolSetConvTestFn(NLS, firkNlsConvTest, firk_mem);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Setting convergence test function failed");
    return FIRK_ILL_INPUT;
  }

  retval = SUNNonlinSolSetNormFn(NLS, firkNlsNorm, firk_mem);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Setting convergence-test norm function failed");
    return FIRK_ILL_INPUT;
  }

  retval = SUNNonlinSolSetGetUpdateNormFn(NLS, firkNlsGetUpdateNorm, firk_mem);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Setting update-norm getter failed");
    return FIRK_ILL_INPUT;
  }

  retval = SUNNonlinSolSetGetConvRateFn(NLS, firkNlsGetConvRate, firk_mem);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Setting convergence-rate getter failed");
    return FIRK_ILL_INPUT;
  }

  retval = SUNNonlinSolSetMaxIters(NLS, firk_mem->maxcor);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Setting maximum number of nonlinear iterations failed");
    return FIRK_ILL_INPUT;
  }

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkNlsInit

  Sets the linear solver wrapper functions and initializes the
  nonlinear solver; called from the initial setup.
  ---------------------------------------------------------------*/
int firkNlsInit(FIRKodeMem firk_mem)
{
  int retval;

  if (firk_mem->lsetup)
  {
    retval = SUNNonlinSolSetLSetupFn(firk_mem->NLS, firkNlsLSetup);
  }
  else { retval = SUNNonlinSolSetLSetupFn(firk_mem->NLS, NULL); }
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Setting the linear solver setup function failed");
    return FIRK_NLS_INIT_FAIL;
  }

  retval = SUNNonlinSolSetLSolveFn(firk_mem->NLS, firkNlsLSolve);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Setting linear solver solve function failed");
    return FIRK_NLS_INIT_FAIL;
  }

  retval = SUNNonlinSolSetMaxIters(firk_mem->NLS, firk_mem->maxcor);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Setting maximum number of nonlinear iterations failed");
    return FIRK_NLS_INIT_FAIL;
  }

  retval = SUNNonlinSolInitialize(firk_mem->NLS);
  if (retval != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGFIRK_NLS_INIT_FAIL);
    return FIRK_NLS_INIT_FAIL;
  }

  return FIRK_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Private functions
 * ---------------------------------------------------------------------------*/

/*---------------------------------------------------------------
  firkNlsLSetup

  Wrapper calling the linear solver setup routine at the
  linearization point (tn, yn, f(tn,yn)).
  ---------------------------------------------------------------*/
static int firkNlsLSetup(sunbooleantype jbad, sunbooleantype* jcur,
                         void* firkode_mem)
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

  if (jbad) { firk_mem->convfail = FIRK_FAIL_BAD_J; }

  retval = firk_mem->lsetup(firk_mem, firk_mem->convfail, firk_mem->yn,
                            firk_mem->fn, &(firk_mem->jcur), firk_mem->tempv1,
                            firk_mem->tempv2, firk_mem->tempv3);
  firk_mem->nsetups++;

  *jcur = firk_mem->jcur;

  firk_mem->gamrat = ONE;
  firk_mem->gammap = firk_mem->gamma;
  firk_mem->delnrm = ZERO;
  firk_mem->nstlp  = firk_mem->nst;

  if (retval < 0) { return FIRK_LSETUP_FAIL; }
  if (retval > 0) { return SUN_NLS_CONV_RECVR; }

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkNlsLSolve

  Wrapper solving the stacked Newton system in place.
  ---------------------------------------------------------------*/
static int firkNlsLSolve(N_Vector delta, void* firkode_mem)
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

  retval = firk_mem->lsolve_stk(firk_mem, delta);

  if (retval < 0) { return FIRK_LSOLVE_FAIL; }
  if (retval > 0) { return SUN_NLS_CONV_RECVR; }

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkNlsConvTest

  RADAU5-style convergence test on the stacked Newton update.
  ---------------------------------------------------------------*/
static int firkNlsConvTest(SUNNonlinearSolver NLS,
                           SUNDIALS_MAYBE_UNUSED N_Vector ycor, N_Vector delta,
                           sunrealtype tol, N_Vector ewt, void* firkode_mem)
{
  FIRKodeMem firk_mem;
  int m, k, retval;
  sunrealtype eta_k, dyth, qnewt;

  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRK_MEM_NULL, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MEM);
    return FIRK_MEM_NULL;
  }
  firk_mem = (FIRKodeMem)firkode_mem;

  if (firkNlsNorm(delta, ewt, &firk_mem->delnrm, firkode_mem) != SUN_SUCCESS)
  {
    firkProcessError(firk_mem, FIRK_NLS_FAIL, __LINE__, __func__, __FILE__,
                     MSGFIRK_NLS_FAIL);
    return FIRK_NLS_FAIL;
  }

  /* m = number of completed iterations before this one; k = m + 1 */
  retval = SUNNonlinSolGetCurIter(NLS, &m);
  if (retval != SUN_SUCCESS) { return FIRK_MEM_NULL; }
  k = m + 1;

  if (m == 0)
  {
    /* reuse the contraction estimate of the previous solve */
    eta_k = SUNRpowerR(SUNMAX(firk_mem->eta_nls, firk_mem->uround),
                       RATE_REUSE_EXP);
  }
  else
  {
    firk_mem->theta = firk_mem->delnrm / SUNMAX(firk_mem->delp, firk_mem->uround);

    /* divergence */
    if (firk_mem->theta >= firk_mem->theta_max)
    {
      firk_mem->nls_reason  = FIRK_NLS_REASON_DIVERGE;
      firk_mem->eta_nlsfail = firk_mem->hadapt_mem->etacf;
      SUNLogDebug(FIRK_LOGGER, "nonlinear-diverge", "theta = " SUN_FORMAT_G,
                  firk_mem->theta);
      return SUN_NLS_CONV_RECVR;
    }

    eta_k = firk_mem->theta / (ONE - firk_mem->theta);

    /* predicted failure to converge within maxcor iterations */
    if (k < firk_mem->maxcor)
    {
      dyth = eta_k * firk_mem->delnrm *
             SUNRpowerI(firk_mem->theta, firk_mem->maxcor - 1 - k) / tol;
      if (dyth >= ONE)
      {
        qnewt = SUNMAX(QNEWT_MIN, SUNMIN(QNEWT_MAX, dyth));
        firk_mem->eta_nlsfail =
          PREDICT_FAIL_FAC *
          SUNRpowerR(qnewt, -ONE / (SUN_RCONST(4.0) +
                                    (sunrealtype)(firk_mem->maxcor - 1 - k)));
        firk_mem->nls_reason = FIRK_NLS_REASON_PREDICT;
        SUNLogDebug(FIRK_LOGGER, "nonlinear-predicted-failure",
                    "theta = " SUN_FORMAT_G ", eta = " SUN_FORMAT_G,
                    firk_mem->theta, firk_mem->eta_nlsfail);
        return SUN_NLS_CONV_RECVR;
      }
    }
  }

  firk_mem->eta_nls = eta_k;
  firk_mem->delp    = firk_mem->delnrm;

  SUNLogDebug(FIRK_LOGGER, "nonlinear-convtest",
              "iter = %i, delnrm = " SUN_FORMAT_G ", theta = " SUN_FORMAT_G
              ", test = " SUN_FORMAT_G,
              k, firk_mem->delnrm, firk_mem->theta, eta_k * firk_mem->delnrm);

  if (eta_k * firk_mem->delnrm <= tol) { return SUN_SUCCESS; }

  return SUN_NLS_CONTINUE;
}

static SUNErrCode firkNlsNorm(N_Vector delta, N_Vector ewt, sunrealtype* delnrm,
                              SUNDIALS_MAYBE_UNUSED void* firkode_mem)
{
  *delnrm = N_VWrmsNorm(delta, ewt);
  return SUN_SUCCESS;
}

static SUNErrCode firkNlsGetUpdateNorm(sunrealtype* delnrm, void* firkode_mem)
{
  FIRKodeMem firk_mem;

  if (firkode_mem == NULL) { return SUN_ERR_ARG_CORRUPT; }
  firk_mem = (FIRKodeMem)firkode_mem;

  *delnrm = firk_mem->delnrm;
  return SUN_SUCCESS;
}

static SUNErrCode firkNlsGetConvRate(sunrealtype* crate, void* firkode_mem)
{
  FIRKodeMem firk_mem;

  if (firkode_mem == NULL) { return SUN_ERR_ARG_CORRUPT; }
  firk_mem = (FIRKodeMem)firkode_mem;

  *crate = firk_mem->theta;
  return SUN_SUCCESS;
}

/*---------------------------------------------------------------
  firkNlsResidual

  Evaluates the stacked stage residual
    res_i = M (Zpred_i + ycor_i) - h sum_j a_ij f(tn + c_j h, yn + Z_j).
  ---------------------------------------------------------------*/
static int firkNlsResidual(N_Vector ycor, N_Vector res, void* firkode_mem)
{
  FIRKodeMem firk_mem;
  FIRKodeTable T;
  N_Vector res_i, mz;
  int i, j, s, retval;

  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRK_MEM_NULL, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MEM);
    return FIRK_MEM_NULL;
  }
  firk_mem = (FIRKodeMem)firkode_mem;
  T        = firk_mem->T;
  s        = firk_mem->s;

  /* update the stage increments and evaluate the stage right-hand sides */
  for (i = 0; i < s; i++)
  {
    N_VLinearSum(ONE, firk_mem->Zpred[i], ONE,
                 N_VGetSubvector_ManyVector(ycor, (sunindextype)i),
                 firk_mem->Z[i]);
    N_VLinearSum(ONE, firk_mem->yn, ONE, firk_mem->Z[i], firk_mem->tempv1);
    retval = firk_mem->f(firk_mem->tn + T->c[i] * firk_mem->h, firk_mem->tempv1,
                         firk_mem->F[i], firk_mem->user_data);
    firk_mem->nfe++;
    if (retval < 0) { return FIRK_RHSFUNC_FAIL; }
    if (retval > 0) { return RHSFUNC_RECVR; }
  }

  /* assemble the residual for each stage */
  for (i = 0; i < s; i++)
  {
    res_i = N_VGetSubvector_ManyVector(res, (sunindextype)i);

    if (firk_mem->mass_set)
    {
      retval = firk_mem->mmult(firk_mem, firk_mem->Z[i], res_i);
      if (retval != 0) { return FIRK_MASSMULT_FAIL; }
      mz = res_i;
    }
    else { mz = firk_mem->Z[i]; }

    firk_mem->cvals[0] = ONE;
    firk_mem->Xvecs[0] = mz;
    for (j = 0; j < s; j++)
    {
      firk_mem->cvals[j + 1] = -firk_mem->h * T->A[i][j];
      firk_mem->Xvecs[j + 1] = firk_mem->F[j];
    }
    retval = N_VLinearCombination(s + 1, firk_mem->cvals, firk_mem->Xvecs, res_i);
    if (retval != 0) { return FIRK_VECTOROP_ERR; }
  }

  return FIRK_SUCCESS;
}
