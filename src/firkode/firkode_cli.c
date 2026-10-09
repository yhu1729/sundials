/* -----------------------------------------------------------------------------
 * Programmer(s): Yifan Hu @ UMBC
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * Command-line control over the optional inputs of FIRKODE (design Section
 * 9.4): FIRKodeSetOptions reads "firkode.<key> <value>" pairs, one key per
 * scalar setter named after it in snake case, through the helpers of
 * sundials_cli.h as CVODE and ARKODE do. A key of an EXPERT option is
 * rejected in BEGINNER mode like its setter; "firkode.mode" is applied
 * before every other key of the same call.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"
#include "firkode_impl.h"
#include "sundials_cli.h"

/*---------------------------------------------------------------
  Setters of the options whose values are names
  ---------------------------------------------------------------*/

static int firkSetModeByName(void* firkode_mem, const char* name)
{
  if (strcmp(name, "beginner") == 0 || strcmp(name, "BEGINNER") == 0)
  {
    return FIRKodeSetMode(firkode_mem, FIRK_BEGINNER);
  }
  if (strcmp(name, "expert") == 0 || strcmp(name, "EXPERT") == 0)
  {
    return FIRKodeSetMode(firkode_mem, FIRK_EXPERT);
  }
  return FIRK_ILL_INPUT;
}

static int firkSetReusePolicyByName(void* firkode_mem, const char* name)
{
  if (strcmp(name, "radau5") == 0 || strcmp(name, "RADAU5") == 0)
  {
    return FIRKodeSetReusePolicy(firkode_mem, FIRK_POLICY_RADAU5);
  }
  if (strcmp(name, "periodic") == 0 || strcmp(name, "PERIODIC") == 0)
  {
    return FIRKodeSetReusePolicy(firkode_mem, FIRK_POLICY_PERIODIC);
  }
  return FIRK_ILL_INPUT;
}

static int firkSetStageSolverGSTypeByName(void* firkode_mem, const char* name)
{
  if (strcmp(name, "modified") == 0)
  {
    return FIRKodeSetStageSolverGSType(firkode_mem, SUN_MODIFIED_GS);
  }
  if (strcmp(name, "classical") == 0)
  {
    return FIRKodeSetStageSolverGSType(firkode_mem, SUN_CLASSICAL_GS);
  }
  return FIRK_ILL_INPUT;
}

/*---------------------------------------------------------------
  The keys, one per scalar setter
  ---------------------------------------------------------------*/

static const struct sunKeyIntPair int_pairs[] =
  {{"max_conv_fails", FIRKodeSetMaxConvFails},
   {"max_err_test_fails", FIRKodeSetMaxErrTestFails},
   {"max_hnil_warns", FIRKodeSetMaxHnilWarns},
   {"max_nonlin_iters", FIRKodeSetMaxNonlinIters},
   {"num_stages", FIRKodeSetNumStages},
   {"order", FIRKodeSetOrder},
   {"predictor_method", FIRKodeSetPredictorMethod},
   {"small_num_efails", FIRKodeSetSmallNumEFails},
   {"error_refilter", FIRKodeSetErrorRefilter},
   {"interpolate_stop_time", FIRKodeSetInterpolateStopTime},
   {"newton_count_safety", FIRKodeSetNewtonCountSafety},
   {"linear_solution_scaling", FIRKodeSetLinearSolutionScaling},
   {"stage_solver_maxl", FIRKodeSetStageSolverMaxl},
   {"stage_solver_max_restarts", FIRKodeSetStageSolverMaxRestarts}};
static const int num_int_keys = sizeof(int_pairs) / sizeof(*int_pairs);

static const struct sunKeyLongPair long_pairs[] =
  {{"max_num_steps", FIRKodeSetMaxNumSteps},
   {"jac_eval_frequency", FIRKodeSetJacEvalFrequency},
   {"lsetup_frequency", FIRKodeSetLSetupFrequency}};
static const int num_long_keys = sizeof(long_pairs) / sizeof(*long_pairs);

static const struct sunKeyRealPair real_pairs[] =
  {{"init_step", FIRKodeSetInitStep},
   {"min_step", FIRKodeSetMinStep},
   {"max_step", FIRKodeSetMaxStep},
   {"fixed_step", FIRKodeSetFixedStep},
   {"stop_time", FIRKodeSetStopTime},
   {"error_bias", FIRKodeSetErrorBias},
   {"safety_factor", FIRKodeSetSafetyFactor},
   {"max_growth", FIRKodeSetMaxGrowth},
   {"max_first_growth", FIRKodeSetMaxFirstGrowth},
   {"min_reduction", FIRKodeSetMinReduction},
   {"max_efail_growth", FIRKodeSetMaxEFailGrowth},
   {"max_cfail_growth", FIRKodeSetMaxCFailGrowth},
   {"nonlin_conv_coef", FIRKodeSetNonlinConvCoef},
   {"nonlin_divergence_rate", FIRKodeSetNonlinDivergenceRate},
   {"nonlin_rate_exponent", FIRKodeSetNonlinRateExponent},
   {"jac_reuse_conv_rate", FIRKodeSetJacReuseConvRate},
   {"jac_bad_conv_rate", FIRKodeSetJacBadConvRate},
   {"delta_gamma_max_lsetup", FIRKodeSetDeltaGammaMaxLSetup},
   {"delta_gamma_max_bad_jac", FIRKodeSetDeltaGammaMaxBadJac},
   {"eps_lin", FIRKodeSetEpsLin},
   {"ls_norm_factor", FIRKodeSetLSNormFactor},
   {"stage_eps_lin", FIRKodeSetStageEpsLin},
   {"mass_eps_lin", FIRKodeSetMassEpsLin},
   {"mass_ls_norm_factor", FIRKodeSetMassLSNormFactor}};
static const int num_real_keys = sizeof(real_pairs) / sizeof(*real_pairs);

static const struct sunKeyTwoRealPair tworeal_pairs[] =
  {{"scalar_tolerances", FIRKodeSStolerances},
   {"fixed_step_bounds", FIRKodeSetFixedStepBounds}};
static const int num_tworeal_keys = sizeof(tworeal_pairs) /
                                    sizeof(*tworeal_pairs);

static const struct sunKeyCharPair char_pairs[] =
  {{"reuse_policy", firkSetReusePolicyByName},
   {"adapt_controller_by_name", FIRKodeSetAdaptControllerByName},
   {"stage_solver_gs_type", firkSetStageSolverGSTypeByName}};
static const int num_char_keys = sizeof(char_pairs) / sizeof(*char_pairs);

static const struct sunKeyActionPair action_pairs[] =
  {{"clear_stop_time", FIRKodeClearStopTime},
   {"no_inactive_root_warn", FIRKodeSetNoInactiveRootWarn}};
static const int num_action_keys = sizeof(action_pairs) / sizeof(*action_pairs);

/* the number of values a key takes: 0 to 2, or -1 for an unknown key */
static int firkOptionValues(const char* key)
{
  int j;
  for (j = 0; j < num_int_keys; j++)
  {
    if (strcmp(key, int_pairs[j].key) == 0) { return 1; }
  }
  for (j = 0; j < num_long_keys; j++)
  {
    if (strcmp(key, long_pairs[j].key) == 0) { return 1; }
  }
  for (j = 0; j < num_real_keys; j++)
  {
    if (strcmp(key, real_pairs[j].key) == 0) { return 1; }
  }
  for (j = 0; j < num_char_keys; j++)
  {
    if (strcmp(key, char_pairs[j].key) == 0) { return 1; }
  }
  for (j = 0; j < num_tworeal_keys; j++)
  {
    if (strcmp(key, tworeal_pairs[j].key) == 0) { return 2; }
  }
  for (j = 0; j < num_action_keys; j++)
  {
    if (strcmp(key, action_pairs[j].key) == 0) { return 0; }
  }
  if (strcmp(key, "mode") == 0) { return 1; }
  return -1;
}

/*---------------------------------------------------------------
  FIRKodeSetOptions

  Sets FIRKODE options from strings. File-based options are not
  supported.
  ---------------------------------------------------------------*/
int FIRKodeSetOptions(void* firkode_mem, const char* firkid,
                      const char* file_name, int argc, char* argv[])
{
  FIRKodeMem firk_mem;
  const char* id;
  const char* key;
  char* prefix;
  size_t offset;
  int idx, j, nvals, retval;
  sunbooleantype arg_used;

  if (firkode_mem == NULL)
  {
    firkProcessError(NULL, FIRK_MEM_NULL, __LINE__, __func__, __FILE__,
                     MSGFIRK_NO_MEM);
    return FIRK_MEM_NULL;
  }
  firk_mem = (FIRKodeMem)firkode_mem;

  if (file_name != NULL && strlen(file_name) > 0)
  {
    firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "file-based options are not supported.");
    return FIRK_ILL_INPUT;
  }
  if (argc <= 0 || argv == NULL) { return FIRK_SUCCESS; }

  /* the prefix "<firkid>." */
  id     = (firkid != NULL && strlen(firkid) > 0) ? firkid : "firkode";
  offset = strlen(id) + 1;
  prefix = (char*)malloc(offset + 1);
  if (prefix == NULL)
  {
    firkProcessError(firk_mem, FIRK_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSGFIRK_MEM_FAIL);
    return FIRK_MEM_FAIL;
  }
  strcpy(prefix, id);
  strcat(prefix, ".");

  /* the mode first, wherever it is, so that it unlocks the EXPERT keys
     beside it */
  for (idx = 1; idx < argc; idx++)
  {
    if (strncmp(argv[idx], prefix, offset) != 0 ||
        strcmp(argv[idx] + offset, "mode") != 0)
    {
      continue;
    }
    if (idx + 1 >= argc) { retval = FIRK_ILL_INPUT; }
    else { retval = firkSetModeByName(firkode_mem, argv[idx + 1]); }
    if (retval != FIRK_SUCCESS)
    {
      firkProcessError(firk_mem, retval, __LINE__, __func__, __FILE__,
                       "error setting key: mode");
      free(prefix);
      return retval;
    }
    idx++;
  }

  for (idx = 1; idx < argc; idx++)
  {
    /* only "<firkid>." arguments */
    if (strncmp(argv[idx], prefix, offset) != 0) { continue; }
    key   = argv[idx] + offset;
    nvals = firkOptionValues(key);
    if (nvals < 0)
    {
      firkProcessError(firk_mem, FIRK_WARNING, __LINE__, __func__, __FILE__,
                       "key %s was not handled", argv[idx]);
      continue;
    }
    if (idx + nvals >= argc)
    {
      firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                       "key %s needs %i value(s)", argv[idx], nvals);
      free(prefix);
      return FIRK_ILL_INPUT;
    }

    /* the mode is set */
    if (strcmp(key, "mode") == 0)
    {
      idx++;
      continue;
    }

    arg_used = SUNFALSE;
    j        = 0;
    retval   = sunCheckAndSetIntArgs(firkode_mem, &idx, argv, offset, int_pairs,
                                     num_int_keys, &arg_used, &j);
    if (retval == FIRK_SUCCESS && !arg_used)
    {
      retval = sunCheckAndSetLongArgs(firkode_mem, &idx, argv, offset,
                                      long_pairs, num_long_keys, &arg_used, &j);
    }
    if (retval == FIRK_SUCCESS && !arg_used)
    {
      retval = sunCheckAndSetRealArgs(firkode_mem, &idx, argv, offset,
                                      real_pairs, num_real_keys, &arg_used, &j);
    }
    if (retval == FIRK_SUCCESS && !arg_used)
    {
      retval = sunCheckAndSetTwoRealArgs(firkode_mem, &idx, argv, offset,
                                         tworeal_pairs, num_tworeal_keys,
                                         &arg_used, &j);
    }
    if (retval == FIRK_SUCCESS && !arg_used)
    {
      retval = sunCheckAndSetCharArgs(firkode_mem, &idx, argv, offset,
                                      char_pairs, num_char_keys, &arg_used, &j);
    }
    if (retval == FIRK_SUCCESS && !arg_used)
    {
      retval = sunCheckAndSetActionArgs(firkode_mem, &idx, argv, offset,
                                        action_pairs, num_action_keys,
                                        &arg_used, &j);
    }
    if (retval != FIRK_SUCCESS)
    {
      firkProcessError(firk_mem, retval, __LINE__, __func__, __FILE__,
                       "error setting key: %s", key);
      free(prefix);
      return retval;
    }
  }

  free(prefix);
  return FIRK_SUCCESS;
}
