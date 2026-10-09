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
 * Unit test for the FIRKODE modes (design decision D8 and Section 9.1): the
 * default is BEGINNER, FIRKodeSetMode accepts FIRK_BEGINNER and FIRK_EXPERT
 * until the first step and again after FIRKodeReInit, rejects other values and
 * a change after a step with FIRK_ILL_INPUT, FIRKodeGetMode reports the mode,
 * and FIRKodeWriteParameters prints the mode first. In BEGINNER mode every
 * setter of an EXPERT option returns FIRK_ILL_INPUT with the message
 * "<function> requires FIRKodeSetMode(mem, FIRK_EXPERT).", before and after a
 * linear solver is attached; in EXPERT mode it succeeds; and returning to
 * BEGINNER mode restores every default. Pass criterion: every return value,
 * reported mode, message and parameter listing matches.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_context.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)t;
  (void)user_data;
  N_VScale(-SUN_RCONST(1.0), y, ydot);
  return 0;
}

static int check(int ok, const char* what)
{
  if (!ok)
  {
    printf("  FAIL: %s\n", what);
    return 1;
  }
  return 0;
}

/* the last error reported through the context's error handlers */
static char last_func[128];
static char last_msg[256];

static void capture(int line, const char* func, const char* file, const char* msg,
                    SUNErrCode err_code, void* err_user_data, SUNContext sunctx)
{
  (void)line;
  (void)file;
  (void)err_code;
  (void)err_user_data;
  (void)sunctx;
  snprintf(last_func, sizeof(last_func), "%s", func);
  snprintf(last_msg, sizeof(last_msg), "%s", msg);
}

/* the setters of EXPERT options (design Table 9.1), with valid arguments */
#define NUM_EXPERT 32
static const char* expert_names[NUM_EXPERT] =
  {"FIRKodeSetMaxErrTestFails",
   "FIRKodeSetMaxConvFails",
   "FIRKodeSetPredictorMethod",
   "FIRKodeSetErrorRefilter",
   "FIRKodeSetAdaptController",
   "FIRKodeSetSafetyFactor",
   "FIRKodeSetErrorBias",
   "FIRKodeSetMaxGrowth",
   "FIRKodeSetMinReduction",
   "FIRKodeSetFixedStepBounds",
   "FIRKodeSetMaxFirstGrowth",
   "FIRKodeSetMaxEFailGrowth",
   "FIRKodeSetSmallNumEFails",
   "FIRKodeSetMaxCFailGrowth",
   "FIRKodeSetNewtonCountSafety",
   "FIRKodeSetMaxNonlinIters",
   "FIRKodeSetNonlinConvCoef",
   "FIRKodeSetNonlinDivergenceRate",
   "FIRKodeSetJacBadConvRate",
   "FIRKodeSetLSetupFrequency",
   "FIRKodeSetDeltaGammaMaxLSetup",
   "FIRKodeSetNumStages",
   "FIRKodeSetOrder",
   "FIRKodeSetJacEvalFrequency",
   "FIRKodeSetDeltaGammaMaxBadJac",
   "FIRKodeSetStageSolverMaxl",
   "FIRKodeSetStageSolverMaxRestarts",
   "FIRKodeSetStageEpsLin",
   "FIRKodeSetAdaptControllerByName",
   "FIRKodeSetNonlinRateExponent",
   "FIRKodeSetJacReuseConvRate",
   "FIRKodeSetReusePolicy"};

static int call_expert(void* mem, int i)
{
  switch (i)
  {
  case 0: return FIRKodeSetMaxErrTestFails(mem, 5);
  case 1: return FIRKodeSetMaxConvFails(mem, 5);
  case 2: return FIRKodeSetPredictorMethod(mem, FIRK_PREDICT_TRIVIAL);
  case 3: return FIRKodeSetErrorRefilter(mem, SUNFALSE);
  case 4: return FIRKodeSetAdaptController(mem, NULL);
  case 5: return FIRKodeSetSafetyFactor(mem, SUN_RCONST(0.8));
  case 6: return FIRKodeSetErrorBias(mem, SUN_RCONST(1.2));
  case 7: return FIRKodeSetMaxGrowth(mem, SUN_RCONST(5.0));
  case 8: return FIRKodeSetMinReduction(mem, SUN_RCONST(0.1));
  case 9:
    return FIRKodeSetFixedStepBounds(mem, SUN_RCONST(1.0), SUN_RCONST(1.5));
  case 10: return FIRKodeSetMaxFirstGrowth(mem, SUN_RCONST(100.0));
  case 11: return FIRKodeSetMaxEFailGrowth(mem, SUN_RCONST(0.2));
  case 12: return FIRKodeSetSmallNumEFails(mem, 3);
  case 13: return FIRKodeSetMaxCFailGrowth(mem, SUN_RCONST(0.2));
  case 14: return FIRKodeSetNewtonCountSafety(mem, SUNFALSE);
  case 15: return FIRKodeSetMaxNonlinIters(mem, 5);
  case 16: return FIRKodeSetNonlinConvCoef(mem, SUN_RCONST(0.05));
  case 17: return FIRKodeSetNonlinDivergenceRate(mem, SUN_RCONST(0.9));
  case 18: return FIRKodeSetJacBadConvRate(mem, SUN_RCONST(0.5));
  case 19: return FIRKodeSetLSetupFrequency(mem, 10);
  case 20: return FIRKodeSetDeltaGammaMaxLSetup(mem, SUN_RCONST(0.2));
  case 21: return FIRKodeSetNumStages(mem, 3);
  case 22: return FIRKodeSetOrder(mem, 5);
  case 23: return FIRKodeSetJacEvalFrequency(mem, 30);
  case 24: return FIRKodeSetDeltaGammaMaxBadJac(mem, SUN_RCONST(0.1));
  case 25: return FIRKodeSetStageSolverMaxl(mem, 12);
  case 26: return FIRKodeSetStageSolverMaxRestarts(mem, 2);
  case 27: return FIRKodeSetStageEpsLin(mem, SUN_RCONST(0.01));
  case 28: return FIRKodeSetAdaptControllerByName(mem, "PI");
  case 29: return FIRKodeSetNonlinRateExponent(mem, SUN_RCONST(0.7));
  case 30: return FIRKodeSetJacReuseConvRate(mem, SUN_RCONST(0.01));
  case 31: return FIRKodeSetReusePolicy(mem, FIRK_POLICY_PERIODIC);
  default: return -1;
  }
}

/* in BEGINNER mode every EXPERT setter is rejected with the mode message */
static int check_rejections(void* mem, const char* when)
{
  int i, fails = 0;
  char expected[256];

  for (i = 0; i < NUM_EXPERT; i++)
  {
    last_func[0] = last_msg[0] = '\0';
    snprintf(expected, sizeof(expected),
             "%s requires FIRKodeSetMode(mem, FIRK_EXPERT).", expert_names[i]);
    if (call_expert(mem, i) != FIRK_ILL_INPUT ||
        strcmp(last_func, expert_names[i]) != 0 ||
        strcmp(last_msg, expected) != 0)
    {
      printf("  FAIL: %s in BEGINNER mode %s: function '%s', message '%s'\n",
             expert_names[i], when, last_func, last_msg);
      fails++;
    }
  }
  return fails;
}

/* the complete output of FIRKodeWriteParameters */
static int write_parameters(void* mem, char* text, size_t n)
{
  FILE* fp = tmpfile();
  size_t len;

  if (fp == NULL || FIRKodeWriteParameters(mem, fp) != FIRK_SUCCESS)
  {
    if (fp) { fclose(fp); }
    return -1;
  }
  rewind(fp);
  len       = fread(text, 1, n - 1, fp);
  text[len] = '\0';
  fclose(fp);
  return 0;
}

/* a memory block with a linear solver attached, ready to integrate */
static void* new_mem(SUNContext sunctx, N_Vector y, SUNLinearSolver LS,
                     SUNMatrix A)
{
  void* mem = FIRKodeCreate(sunctx);
  if (mem == NULL || FIRKodeInit(mem, rhs, SUN_RCONST(0.0), y) ||
      FIRKodeSStolerances(mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-10)) ||
      FIRKodeSetLinearSolver(mem, LS, A))
  {
    return NULL;
  }
  return mem;
}

/* the first parameter line written by FIRKodeWriteParameters */
static int mode_line(void* mem, char* line, int n)
{
  FILE* fp = tmpfile();
  int retval;

  if (fp == NULL) { return -1; }
  retval = FIRKodeWriteParameters(mem, fp);
  rewind(fp);
  if (retval != FIRK_SUCCESS || fgets(line, n, fp) == NULL ||
      fgets(line, n, fp) == NULL)
  {
    retval = -1;
  }
  fclose(fp);
  return retval;
}

int main(void)
{
  SUNContext sunctx  = NULL;
  void* mem          = NULL;
  N_Vector y         = NULL;
  SUNMatrix A        = NULL;
  SUNLinearSolver LS = NULL;
  sunrealtype t;
  char line[128];
  int mode = -1, fails = 0;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  /* NULL memory */
  fails += check(FIRKodeSetMode(NULL, FIRK_EXPERT) == FIRK_MEM_NULL,
                 "SetMode(NULL) returns FIRK_MEM_NULL");
  fails += check(FIRKodeGetMode(NULL, &mode) == FIRK_MEM_NULL,
                 "GetMode(NULL) returns FIRK_MEM_NULL");

  mem = FIRKodeCreate(sunctx);
  if (mem == NULL) { return 1; }

  /* default and changes before Init */
  fails += check(FIRKodeGetMode(mem, &mode) == FIRK_SUCCESS &&
                   mode == FIRK_BEGINNER,
                 "the default mode is BEGINNER");
  fails += check(FIRKodeSetMode(mem, FIRK_EXPERT) == FIRK_SUCCESS,
                 "SetMode(EXPERT) before Init");
  fails += check(FIRKodeGetMode(mem, &mode) == FIRK_SUCCESS && mode == FIRK_EXPERT,
                 "GetMode reports EXPERT");
  fails += check(FIRKodeSetMode(mem, 2) == FIRK_ILL_INPUT,
                 "SetMode rejects an unknown mode");
  fails += check(FIRKodeSetMode(mem, -1) == FIRK_ILL_INPUT,
                 "SetMode rejects a negative mode");
  fails += check(FIRKodeGetMode(mem, &mode) == FIRK_SUCCESS && mode == FIRK_EXPERT,
                 "a rejected SetMode keeps the mode");
  fails += check(mode_line(mem, line, (int)sizeof(line)) == FIRK_SUCCESS &&
                   strcmp(line, "  Mode = EXPERT\n") == 0,
                 "WriteParameters prints the mode first");

  /* changes after Init and before the first step */
  y = N_VNew_Serial(1, sunctx);
  N_VConst(SUN_RCONST(1.0), y);
  A  = SUNDenseMatrix(1, 1, sunctx);
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (FIRKodeInit(mem, rhs, SUN_RCONST(0.0), y) ||
      FIRKodeSStolerances(mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-10)) ||
      FIRKodeSetLinearSolver(mem, LS, A))
  {
    printf("  FAIL: setup failed\n");
    return 1;
  }
  fails += check(FIRKodeSetMode(mem, FIRK_BEGINNER) == FIRK_SUCCESS,
                 "SetMode(BEGINNER) after Init");
  fails += check(mode_line(mem, line, (int)sizeof(line)) == FIRK_SUCCESS &&
                   strcmp(line, "  Mode = BEGINNER\n") == 0,
                 "WriteParameters prints BEGINNER");

  /* the mode is fixed once a step has been attempted */
  if (FIRKodeEvolve(mem, SUN_RCONST(1.0), y, &t, FIRK_ONE_STEP) < 0)
  {
    printf("  FAIL: FIRKodeEvolve failed\n");
    return 1;
  }
  fails += check(FIRKodeSetMode(mem, FIRK_EXPERT) == FIRK_ILL_INPUT,
                 "SetMode after a step returns FIRK_ILL_INPUT");
  fails += check(FIRKodeGetMode(mem, &mode) == FIRK_SUCCESS &&
                   mode == FIRK_BEGINNER,
                 "a rejected SetMode after a step keeps the mode");

  /* Reset keeps the step history, ReInit starts afresh */
  fails += check(FIRKodeReset(mem, t, y) == FIRK_SUCCESS &&
                   FIRKodeSetMode(mem, FIRK_EXPERT) == FIRK_ILL_INPUT,
                 "SetMode after Reset returns FIRK_ILL_INPUT");
  fails += check(FIRKodeReInit(mem, SUN_RCONST(0.0), y) == FIRK_SUCCESS &&
                   FIRKodeSetMode(mem, FIRK_EXPERT) == FIRK_SUCCESS,
                 "SetMode after ReInit");

  FIRKodeFree(&mem);

  /* EXPERT options in BEGINNER mode, before and after a linear solver is
     attached, and in EXPERT mode */
  SUNContext_PushErrHandler(sunctx, capture, NULL);
  mem = FIRKodeCreate(sunctx);
  if (mem == NULL || FIRKodeInit(mem, rhs, SUN_RCONST(0.0), y) ||
      FIRKodeSStolerances(mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-10)))
  {
    return 1;
  }
  fails += check_rejections(mem, "before SetLinearSolver");
  if (FIRKodeSetLinearSolver(mem, LS, A)) { return 1; }
  fails += check_rejections(mem, "after SetLinearSolver");
  if (FIRKodeSetMode(mem, FIRK_EXPERT)) { return 1; }
  {
    int i;
    for (i = 0; i < NUM_EXPERT; i++)
    {
      if (call_expert(mem, i) != FIRK_SUCCESS)
      {
        printf("  FAIL: %s in EXPERT mode: %s\n", expert_names[i], last_msg);
        fails++;
      }
    }
  }

  /* returning to BEGINNER mode restores every default */
  {
    static char changed[8192], restored[8192], fresh[8192];
    void* mem2 = new_mem(sunctx, y, LS, A);
    if (mem2 == NULL || write_parameters(mem, changed, sizeof(changed)) ||
        FIRKodeSetMode(mem, FIRK_BEGINNER) ||
        write_parameters(mem, restored, sizeof(restored)) ||
        write_parameters(mem2, fresh, sizeof(fresh)))
    {
      printf("  FAIL: parameter listings could not be written\n");
      return 1;
    }
    fails += check(strcmp(changed, fresh) != 0,
                   "the EXPERT options change the parameter listing");
    fails += check(strcmp(restored, fresh) == 0,
                   "SetMode(BEGINNER) restores the defaults");
    if (strcmp(restored, fresh) != 0)
    {
      printf("--- restored:\n%s--- fresh:\n%s", restored, fresh);
    }
    FIRKodeFree(&mem2);
  }
  FIRKodeFree(&mem);

  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("All tests passed\n");
  return 0;
}
