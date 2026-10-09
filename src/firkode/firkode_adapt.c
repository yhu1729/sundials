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
 * This is the implementation file for the FIRKODE time step
 * adaptivity utilities: the SUNAdaptController glue and the step
 * size heuristics applied to the controller's proposal.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include <sundials/sundials_types.h>

#include "firkode_impl.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)

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
  firkHminEff

  The smallest step size of an adaptive attempt: the user's hmin, or
  RADAU5's 10 u |tn| if larger (decision R18), below which a step
  could not resolve tn. A failure at this size ends the integration
  through the failure returns of the error test or the stage solve.
  ---------------------------------------------------------------*/
sunrealtype firkHminEff(FIRKodeMem firk_mem)
{
  return SUNMAX(firk_mem->hmin,
                HMIN_REL_ROUNDOFF * firk_mem->uround * SUNRabs(firk_mem->tn));
}

/*---------------------------------------------------------------
  firkAdapt

  Step size selection for the next attempt (design Listing 8.6),
  called once per attempt, before the controller's history update:
  sets firk_mem->eta and firk_mem->keep_setup. With
  kappa = safety min(1, (1 + 2m) / (k_nls + 2m)) (the Newton count
  safety factor of RADAU5) and the I-controller proposal
  h_I = kappa |h| dsm^(-1/(p+1)), the new step is
    BEGINNER, accepted: min(safety |h_C|, h_I)  (RADAU5's Gustafsson
                        step with its I-controller minimum),
    EXPERT, accepted:   kappa |h_C|,
    rejected:           h_I in BEGINNER (RADAU5 uses the predictive
                        controller on accepted steps only), and
                        min(kappa |h_C|, h_I) in EXPERT, so that a
                        rejected step always shrinks,
  where h_C is the controller's proposal; then the growth is capped by
  etamax (etamx1 on the first step, growth after an accepted one, 1
  after a rejection: RADAU5's HNEW = MIN(HNEW, H) after REJECT), the
  reduction by etamin, an accepted step whose new size is within the
  dead band [lbound, ubound] and whose Newton rate was at most
  theta_jac keeps h and its setup, and hmin and hmax apply last.
  ---------------------------------------------------------------*/
int firkAdapt(FIRKodeMem firk_mem, sunrealtype dsm, sunbooleantype accepted)
{
  int retval;
  sunrealtype hcur, habs, kappa, h_I, h_C, hnew;
  FIRKodeHAdaptMem hadapt_mem = firk_mem->hadapt_mem;

  hcur                 = firk_mem->h;
  habs                 = SUNRabs(hcur);
  firk_mem->keep_setup = SUNFALSE;

  if (firk_mem->fixedstep || hadapt_mem->hcontroller == NULL)
  {
    firk_mem->eta = ONE;
    return FIRK_SUCCESS;
  }

  /* guard pow(0, negative) (H13) */
  dsm = SUNMAX(dsm, DSM_FLOOR);

  /* the I-controller proposal with the Newton count safety factor */
  kappa = hadapt_mem->safety;
  if (hadapt_mem->newt_safety)
  {
    kappa *= SUNMIN(ONE, (TWO * (sunrealtype)firk_mem->maxcor + ONE) /
                           (TWO * (sunrealtype)firk_mem->maxcor +
                            (sunrealtype)firk_mem->nni_last));
  }
  h_I = kappa * habs * SUNRpowerR(dsm, -ONE / ((sunrealtype)hadapt_mem->p + ONE));

  if (firk_mem->mode == FIRK_BEGINNER && !accepted) { hnew = h_I; }
  else
  {
    retval = SUNAdaptController_EstimateStep(hadapt_mem->hcontroller, hcur,
                                             hadapt_mem->p, dsm, &h_C);
    if (retval != SUN_SUCCESS)
    {
      firkProcessError(firk_mem, FIRK_CONTROLLER_ERR, __LINE__, __func__,
                       __FILE__, "SUNAdaptController_EstimateStep failure.");
      return FIRK_CONTROLLER_ERR;
    }
    h_C = SUNRabs(h_C);
    if (firk_mem->mode == FIRK_BEGINNER)
    {
      hnew = SUNMIN(hadapt_mem->safety * h_C, h_I);
    }
    else if (accepted) { hnew = kappa * h_C; }
    else { hnew = SUNMIN(kappa * h_C, h_I); }
  }

  SUNLogDebug(FIRK_LOGGER, "new-step-before-bounds",
              "h_I = " SUN_FORMAT_G ", h_acc_controller = " SUN_FORMAT_G, h_I,
              hnew);

  /* bounds on the step size growth and reduction */
  hnew = SUNMIN(hnew, hadapt_mem->etamax * habs);
  hnew = SUNMAX(hnew, hadapt_mem->etamin * habs);

  SUNLogDebug(FIRK_LOGGER, "new-step-after-max-min-bounds",
              "h_acc = " SUN_FORMAT_G, hnew);

  /* dead band: keep h, and with it the setup, since the reuse policies see
     gamma == gammap, when the Newton iteration contracted well; keep_setup
     records the decision for the tests */
  if (accepted && firk_mem->theta <= firk_mem->theta_jac &&
      hnew > habs * hadapt_mem->lbound * ONEMSM &&
      hnew < habs * hadapt_mem->ubound * ONEPSM)
  {
    hnew                 = habs;
    firk_mem->keep_setup = SUNTRUE;
  }

  /* hmin (with RADAU5's 10 u |tn|) and hmax apply last */
  firk_mem->eta = hnew / habs;
  firk_mem->eta = SUNMAX(firk_mem->eta, firkHminEff(firk_mem) / habs);
  firk_mem->eta /= SUNMAX(ONE, habs * firk_mem->hmax_inv * firk_mem->eta);

  SUNLogDebug(FIRK_LOGGER, "new-step-eta", "eta = " SUN_FORMAT_G, firk_mem->eta);

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkReusePolicy

  Decides at the start of a step attempt whether to set up the
  linear solver and whether to evaluate the Jacobian (design
  decision D9); FIRKLS executes the decision. Both policies evaluate
  J and set up on the first step, and when the Newton rate of the
  last solve exceeded theta_jbad (H1).

  FIRK_POLICY_RADAU5, the default, is RADAU5's: after an accepted
  step J is evaluated unless the Newton rate was at most theta_jac,
  after a failure only if it is not current (RADAU5's CALJAC), and
  the factorization is redone whenever J or gamma changed, so that a
  step kept in the dead band, or a fixed step, reuses it (H30).

  FIRK_POLICY_PERIODIC is CVODE's schedule: a setup on the first
  step, after any failure, msbp steps after the last setup and when
  gamma changed by more than dgmax_lsetup; a Jacobian on the first
  step and msbj steps after the last one. After a convergence failure
  the Jacobian is evaluated only if it was not current at the failed
  attempt and gamma was then within dgmax_jbad of its setup (CVLS's
  attribution of the failure to the Jacobian), so that J is never
  evaluated twice at the same (tn, yn) (H8). A Newton rate above
  theta_jbad on the last solve forces both (H1).
  ---------------------------------------------------------------*/
void firkReusePolicy(FIRKodeMem firk_mem, int nflag, sunbooleantype* jbad,
                     sunbooleantype* need_setup)
{
  sunbooleantype first, slow, due;
  sunrealtype dgmax_jbad;

  *jbad       = SUNFALSE;
  *need_setup = SUNFALSE;
  if (firk_mem->lsetup == NULL) { return; }

  /* the first step, and a slow Newton iteration on the last accepted step
     (H1), ask for a Jacobian at the first attempt of a step; after a failure
     a Jacobian is evaluated only if the current point has none (H8) */
  first = (nflag == FIRST_CALL);
  slow  = first && (firk_mem->theta_jbad < ONE) &&
         (firk_mem->theta > firk_mem->theta_jbad);

  if (firk_mem->reuse_policy == FIRK_POLICY_RADAU5)
  {
    if (first)
    {
      *jbad = firk_mem->firststage || slow ||
              (firk_mem->theta > firk_mem->theta_jac);
    }
    else { *jbad = !firk_mem->jac_current; }
    *need_setup = *jbad || (firk_mem->gammap == ZERO) ||
                  (firk_mem->gamma != firk_mem->gammap);
    return;
  }

  firkLsJacSchedule(firk_mem, &due, &dgmax_jbad);

  *need_setup = !first || (firk_mem->firststage) || (firk_mem->gammap == ZERO) ||
                (firk_mem->nst >= firk_mem->nstlp + firk_mem->msbp) ||
                (SUNRabs(firk_mem->gamrat - ONE) > firk_mem->dgmax_lsetup) ||
                slow;
  if (!*need_setup) { return; }

  *jbad = (firk_mem->firststage && !firk_mem->jac_current) || due || slow ||
          ((nflag == PREV_CONV_FAIL) && firk_mem->jac_stale_fail &&
           (firk_mem->dgamma_fail < dgmax_jbad));
}

/*---------------------------------------------------------------
  firkReplaceAdaptController

  Replaces the current SUNAdaptController; a NULL input installs
  the default, the implicit Gustafsson controller with k1 = k2 = 1
  (RADAU5's, design Table 5.9).
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
    /* RADAU5's predictive controller: both error exponents 1 / (p + 1) */
    retval = SUNAdaptController_SetParams_ImpGus(C, ONE, ONE);
    if (retval != SUN_SUCCESS)
    {
      (void)SUNAdaptController_Destroy(C);
      firkProcessError(firk_mem, FIRK_CONTROLLER_ERR, __LINE__, __func__,
                       __FILE__, "SUNAdaptController_SetParams_ImpGus failure");
      return FIRK_CONTROLLER_ERR;
    }
    firk_mem->hadapt_mem->owncontroller = SUNTRUE;
  }
  else { firk_mem->hadapt_mem->owncontroller = take_ownership; }

  firk_mem->hadapt_mem->hcontroller = C;

  return FIRK_SUCCESS;
}
