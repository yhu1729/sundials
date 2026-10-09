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
 * This is the implementation file for the collocation polynomial of
 * FIRKODE: the stage predictor and dense output.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <firkode/firkode.h>
#include <sundials/sundials_types.h>

#include "firkode_impl.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/*---------------------------------------------------------------
  firkPredict

  Predicts the stage increments, either as zero or by extrapolating
  the collocation polynomial of the previous step:
    Z_i = u_prev(tn + c_i h) - yn
            = sum_j (L_j(theta_i) - d_j) Zprev_j,  theta_i = 1 + c_i h/hold.
  ---------------------------------------------------------------*/
int firkPredict(FIRKodeMem firk_mem)
{
  int i, j, m, s, retval;
  sunrealtype theta, tpow;
  const struct FIRKodeTableMem* T;

  s = firk_mem->s;
  T = firk_mem->tbl;

  if (firk_mem->predictor == FIRK_PREDICT_TRIVIAL || !firk_mem->have_prev ||
      firk_mem->hold == ZERO)
  {
    for (i = 0; i < s; i++) { N_VConst(ZERO, firk_mem->Z[i]); }
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
      for (j = 0; j < s; j++) { firk_mem->wtmp[j] += tpow * T->P[m * s + j]; }
    }
    retval = N_VLinearCombination(s, firk_mem->wtmp, firk_mem->Zprev,
                                  firk_mem->Z[i]);
    if (retval != 0) { return FIRK_VECTOROP_ERR; }
    SUNLogExtraDebugVec(FIRK_LOGGER, "predictor", firk_mem->Z[i], "Z_%i(:) =", i);
  }

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkDenseEval

  Evaluates the k-th derivative of the collocation polynomial of
  the last completed step at time t (checked by FIRKodeGetDky):
    u(tcur + (theta - 1) hold) = yn + sum_j (L_j(theta) - d_j) Zprev_j,
  with L_j the Lagrange basis on the nodes {0, c_1, ..., c_s}.
  ---------------------------------------------------------------*/
int firkDenseEval(FIRKodeMem firk_mem, sunrealtype t, int k, N_Vector dky)
{
  const struct FIRKodeTableMem* T;
  sunrealtype theta, tpow, fac;
  int i, j, m, s, ier;

  T     = firk_mem->tbl;
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
      for (j = 0; j < s; j++)
      {
        firk_mem->cvals[j + 1] += tpow * T->P[m * s + j];
      }
    }
    ier = N_VLinearCombination(s + 1, firk_mem->cvals, firk_mem->Xvecs, dky);
    if (ier != 0) { return FIRK_VECTOROP_ERR; }
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
    for (j = 0; j < s; j++)
    {
      firk_mem->cvals[j] += fac * tpow * T->P[m * s + j];
    }
  }
  ier = N_VLinearCombination(s, firk_mem->cvals, firk_mem->Xvecs, dky);
  if (ier != 0) { return FIRK_VECTOROP_ERR; }
  N_VScale(SUNRpowerI(firk_mem->hold, -k), dky, dky);

  return FIRK_SUCCESS;
}
