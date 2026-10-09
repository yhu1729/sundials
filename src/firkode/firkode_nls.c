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
 * The simplified Newton iteration on the stage system (design
 * decision D3, Section 4.4, Listing 8.3).
 *
 * With the stage increments Z_i = Y_i - y_n, the stage system is
 *   G_i(Z) = M Z_i - h sum_j a_ij f(t_n + c_j h, y_n + Z_j) = 0,
 * and each iteration k = 1, ..., m solves
 *   (I_s (x) M - h A (x) J) delta = -G(Z),  Z <- Z + delta,
 * with J and the block preconditioner set up by the step driver.
 *
 * The convergence test is RADAU5's, in the WRMS norm over all sN
 * components: with the ratio of successive update norms, the rate
 * theta (from the third iteration the geometric mean of the last two
 * ratios) and eta = theta / (1 - theta), the iteration
 *   - diverges if theta >= theta_max,
 *   - is predicted to fail if eta theta^(m-1-k) ||delta|| > eps,
 *   - converges if eta ||delta|| <= eps after the update,
 * where the first two tests precede the update. The driver never
 * changes h: it returns its status and the step factor nls_hfac for
 * firkHandleNFlag.
 * -----------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "firkode_impl.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/*---------------------------------------------------------------
  firkNewton

  Solves the stage system from the predictor in Z. Returns
  FIRK_NLS_SUCCESS, a positive FIRK_NLS_* status of a recoverable
  failure (with nls_hfac set), or a negative FIRKODE return code.
  Sets nni_last (the iterations of this solve, each one linear solve,
  all counted in nni), theta (the last rate) and eta_nls (the last
  convergence factor, the start of the next solve).
  ---------------------------------------------------------------*/
int firkNewton(FIRKodeMem firk_mem)
{
  const struct FIRKodeTableMem* T;
  N_Vector b_i, mz;
  sunrealtype eta, theta, dyno, dynold, ratio, ratio_prev, dyth, qnewt, hfac;
  int i, j, k, m, s, status, retval;

  SUNDIALS_MARK_FUNCTION_BEGIN(FIRK_PROFILER);

  T = firk_mem->tbl;
  s = firk_mem->s;
  m = firk_mem->maxcor;

  /* the first iteration reuses the convergence factor of the last solve */
  eta        = SUNRpowerR(SUNMAX(firk_mem->eta_nls, firk_mem->uround),
                          firk_mem->rate_exp);
  theta      = firk_mem->theta_jac;
  dyno       = ZERO;
  dynold     = ZERO;
  ratio_prev = ZERO;
  status     = FIRK_NLS_MAXITER;
  hfac       = HFAC_NEWTON_FAIL;

  SUNLogInfo(FIRK_LOGGER, "begin-nonlinear-solve", "tol = " SUN_FORMAT_G,
             firk_mem->nlscoef);

  for (k = 1; k <= m; k++)
  {
    SUNLogInfo(FIRK_LOGGER, "begin-iterations-list", "");

    /* stage right-hand sides F_i = f(tn + c_i h, yn + Z_i) */
    retval = FIRK_SUCCESS;
    for (i = 0; i < s; i++)
    {
      N_VLinearSum(ONE, firk_mem->yn, ONE, firk_mem->Z[i], firk_mem->tempv1);
      retval = firkRhs(firk_mem, firk_mem->tn + T->c[i] * firk_mem->h,
                       firk_mem->tempv1, firk_mem->F[i]);
      if (retval != FIRK_SUCCESS) { break; }
    }
    if (retval < 0)
    {
      status = retval;
      break;
    }
    if (retval > 0)
    {
      status = FIRK_NLS_RHS_RECVR;
      hfac   = firk_mem->hadapt_mem->etacf;
      break;
    }

    /* right-hand side of the Newton system, b_i = -G_i(Z) */
    for (i = 0; i < s; i++)
    {
      b_i = firk_mem->G[i];
      if (firk_mem->mass_set)
      {
        if (firk_mem->mmult(firk_mem, firk_mem->Z[i], b_i) != 0)
        {
          retval = FIRK_MASSMULT_FAIL;
          break;
        }
        mz = b_i;
      }
      else { mz = firk_mem->Z[i]; }

      firk_mem->cvals[0] = -ONE;
      firk_mem->Xvecs[0] = mz;
      for (j = 0; j < s; j++)
      {
        firk_mem->cvals[j + 1] = firk_mem->h * T->A[i * s + j];
        firk_mem->Xvecs[j + 1] = firk_mem->F[j];
      }
      if (N_VLinearCombination(s + 1, firk_mem->cvals, firk_mem->Xvecs, b_i) != 0)
      {
        retval = FIRK_VECTOROP_ERR;
        break;
      }
    }
    if (retval != FIRK_SUCCESS)
    {
      status = retval;
      break;
    }

    /* the linear solve; every solve is one iteration in nni */
    if (s == 1)
    {
      retval = firk_mem->lsolve_blk(firk_mem, firk_mem->G[0],
                                    firk_mem->delta[0], FIRKLS_BLOCK_NEWTON, k);
    }
    else
    {
      retval = firk_mem->lsolve_stk(firk_mem, firk_mem->g_stk,
                                    firk_mem->delta_stk, k);
    }
    firk_mem->nni++;
    if (retval > 0)
    {
      status = FIRK_NLS_LSOLVE_RECVR;
      hfac   = firk_mem->hadapt_mem->etacf;
      break;
    }
    if (retval < 0)
    {
      status = FIRK_LSOLVE_FAIL;
      break;
    }

    /* a non-finite update norm is a divergence */
    dyno = N_VWrmsNorm(firk_mem->delta_stk, firk_mem->ewt_stk);
    if (!isfinite(dyno))
    {
      status = FIRK_NLS_DIVERGED;
      hfac   = HFAC_NEWTON_FAIL;
      break;
    }

    /* rate, divergence and predicted failure, before the update */
    if (k >= 2 && k < m)
    {
      ratio      = dyno / dynold;
      theta      = (k == 2) ? ratio : SUNRsqrt(ratio * ratio_prev);
      ratio_prev = ratio;
      if (theta >= firk_mem->theta_max)
      {
        status = FIRK_NLS_DIVERGED;
        hfac   = HFAC_NEWTON_FAIL;
        break;
      }
      eta  = theta / (ONE - theta);
      dyth = eta * dyno * SUNRpowerI(theta, m - 1 - k) / firk_mem->nlscoef;
      if (dyth >= ONE)
      {
        qnewt = SUNMAX(QNEWT_MIN, SUNMIN(QNEWT_MAX, dyth));
        hfac  = PREDICT_FAIL_FAC *
               SUNRpowerR(qnewt,
                          -ONE / (SUN_RCONST(4.0) + (sunrealtype)(m - 1 - k)));
        status = FIRK_NLS_PREDICTED_FAIL;
        break;
      }
    }
    dynold = SUNMAX(dyno, firk_mem->uround);

    /* the update */
    for (i = 0; i < s; i++)
    {
      N_VLinearSum(ONE, firk_mem->Z[i], ONE, firk_mem->delta[i], firk_mem->Z[i]);
    }

    SUNLogDebug(FIRK_LOGGER, "nonlinear-iterate",
                "cur-iter = %i, total-iters = %li, update-norm = " SUN_FORMAT_G
                ", theta = " SUN_FORMAT_G ", eta = " SUN_FORMAT_G,
                k, firk_mem->nni, dyno, theta, eta);

    /* the stop test, after the update */
    if (eta * dyno <= firk_mem->nlscoef)
    {
      status = FIRK_NLS_SUCCESS;
      break;
    }

    SUNLogInfoIf(k < m, FIRK_LOGGER, "end-iterations-list", "status = continue");
  }

  firk_mem->eta_nls  = eta;
  firk_mem->theta    = theta;
  firk_mem->nni_last = (k > m) ? m : k;
  firk_mem->nls_hfac = hfac;

  SUNLogInfoIf(status == FIRK_NLS_SUCCESS, FIRK_LOGGER, "end-iterations-list",
               "status = success");
  SUNLogInfoIf(status != FIRK_NLS_SUCCESS, FIRK_LOGGER, "end-iterations-list",
               "status = failure, retval = %i", status);
  SUNLogDebug(FIRK_LOGGER, "newton-exit",
              "nls-status = %i, nls-hfac = " SUN_FORMAT_G, status, hfac);
  SUNLogInfoIf(status == FIRK_NLS_SUCCESS, FIRK_LOGGER, "end-nonlinear-solve",
               "status = success, iters = %li", firk_mem->nni_last);
  SUNLogInfoIf(status != FIRK_NLS_SUCCESS, FIRK_LOGGER,
               "end-nonlinear-solve", "status = failed, retval = %i, iters = %li, hfac = " SUN_FORMAT_G,
               status, firk_mem->nni_last, hfac);

  if (status != FIRK_NLS_SUCCESS) { firk_mem->nnf++; }
  else { firk_mem->jcur = SUNFALSE; }

  SUNDIALS_MARK_FUNCTION_END(FIRK_PROFILER);
  return status;
}
