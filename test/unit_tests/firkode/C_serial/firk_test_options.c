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
 * Unit test of FIRKodeSetOptions (design Section 9.4): integer, long, real,
 * two-real, named and action keys reach their setters (white box); an EXPERT
 * key is rejected in BEGINNER mode, and firkode.mode is applied before every
 * other key of the call wherever it stands; file names, unknown names, a
 * missing value and a FIRKLS key before the linear solver are errors; an
 * unknown key only warns; arguments without the prefix, or with another
 * one, are ignored.
 *
 * Pass criterion: every return value and every option set matches.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "firkode/firkode.h"
#include "firkode/firkode_impl.h"
#include "firkode/firkode_ls.h"
#include "firkode/firkode_ls_impl.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_spgmr.h"

#define ZERO SUN_RCONST(0.0)

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
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

#define NARGS(a) ((int)(sizeof(a) / sizeof((a)[0])))

/* a modifiable copy of an argument list, as main() receives one */
static char** argv_copy(const char* const* src, int n)
{
  int i;
  char** argv = (char**)malloc((size_t)(n + 1) * sizeof(char*));
  if (argv == NULL) { return NULL; }
  for (i = 0; i < n; i++)
  {
    argv[i] = (char*)malloc(strlen(src[i]) + 1);
    if (argv[i] == NULL) { return NULL; }
    strcpy(argv[i], src[i]);
  }
  argv[n] = NULL;
  return argv;
}

static void argv_free(char** argv, int n)
{
  int i;
  for (i = 0; i < n; i++) { free(argv[i]); }
  free(argv);
}

/* FIRKodeSetOptions on a copy of src */
static int set_options(void* mem, const char* id, const char* const* src, int n)
{
  int retval;
  char** argv = argv_copy(src, n);
  if (argv == NULL) { return -999; }
  retval = FIRKodeSetOptions(mem, id, NULL, n, argv);
  argv_free(argv, n);
  return retval;
}

int main(void)
{
  SUNContext sunctx = NULL;
  N_Vector y;
  SUNLinearSolver LS;
  void* mem;
  FIRKodeMem fm;
  int fails = 0, mode = -1;

  const char* plain_src[]        = {"prog",
                                    "firkode.max_num_steps",
                                    "1234",
                                    "firkode.scalar_tolerances",
                                    "1e-5",
                                    "1e-9",
                                    "firkode.max_hnil_warns",
                                    "4",
                                    "firkode.stop_time",
                                    "2.5",
                                    "firkode.unknown_key",
                                    "3",
                                    "other.max_num_steps",
                                    "7",
                                    "firkode",
                                    "max_num_steps"};
  const char* expert_alone_src[] = {"prog", "firkode.reuse_policy", "periodic"};
  const char* expert_mode_after_src[] =
    {"prog", "firkode.reuse_policy", "periodic", "firkode.jac_reuse_conv_rate",
     "0.01", "firkode.mode",         "expert"};
  const char* bad_name_src[] = {"prog", "firkode.reuse_policy", "bogus"};
  const char* bad_gs_src[]  = {"prog", "firkode.stage_solver_gs_type", "bogus"};
  const char* missing_src[] = {"prog", "firkode.max_num_steps"};
  const char* action_src[]  = {"prog", "firkode.clear_stop_time"};
  const char* ls_keys_src[] = {"prog",      "firkode.eps_lin",
                               "0.1",       "firkode.stage_solver_gs_type",
                               "classical", "firkode.stage_solver_maxl",
                               "4"};
  const char* custom_src[]  = {"prog", "myode.max_nonlin_iters", "6",
                               "firkode.max_nonlin_iters", "5"};
  const char* controller_src[] = {"prog", "firkode.adapt_controller_by_name",
                                  "Soderlind"};

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }
  y = N_VNew_Serial(3, sunctx);
  N_VConst(SUN_RCONST(1.0), y);

  fails +=
    check(set_options(NULL, NULL, plain_src, NARGS(plain_src)) == FIRK_MEM_NULL,
          "a NULL memory is FIRK_MEM_NULL");

  mem = FIRKodeCreate(sunctx);
  if (mem == NULL || FIRKodeInit(mem, f, ZERO, y)) { return 1; }
  fm = (FIRKodeMem)mem;

  fails += check(FIRKodeSetOptions(mem, NULL, "options.txt", 0, NULL) ==
                   FIRK_ILL_INPUT,
                 "a file name is FIRK_ILL_INPUT");
  fails += check(FIRKodeSetOptions(mem, NULL, NULL, 0, NULL) == FIRK_SUCCESS,
                 "no arguments succeed");

  /* BEGINNER keys of every kind; unknown and foreign arguments */
  fails += check(set_options(mem, "firkode", plain_src, NARGS(plain_src)) ==
                   FIRK_SUCCESS,
                 "BEGINNER keys succeed, an unknown key only warns");
  fails += check(fm->mxstep == 1234, "max_num_steps reaches its setter");
  fails += check(fm->reltol == SUN_RCONST(1.0e-5) &&
                   fm->Sabstol == SUN_RCONST(1.0e-9),
                 "scalar_tolerances reaches its setter");
  fails += check(fm->mxhnil == 4, "max_hnil_warns reaches its setter");
  fails += check(fm->tstopset && fm->tstop == SUN_RCONST(2.5),
                 "stop_time reaches its setter");

  /* an action key */
  fails += check(set_options(mem, NULL, action_src, NARGS(action_src)) ==
                     FIRK_SUCCESS &&
                   !fm->tstopset,
                 "clear_stop_time reaches its setter");

  /* EXPERT keys and the mode */
  fails += check(set_options(mem, NULL, expert_alone_src,
                             NARGS(expert_alone_src)) == FIRK_ILL_INPUT,
                 "an EXPERT key is rejected in BEGINNER mode");
  fails += check(set_options(mem, NULL, expert_mode_after_src,
                             NARGS(expert_mode_after_src)) == FIRK_SUCCESS,
                 "firkode.mode after the EXPERT keys unlocks them");
  FIRKodeGetMode(mem, &mode);
  fails += check(mode == FIRK_EXPERT, "mode reaches its setter");
  fails += check(fm->reuse_policy == FIRK_POLICY_PERIODIC,
                 "reuse_policy reaches its setter by name");
  fails += check(fm->theta_jac == SUN_RCONST(0.01),
                 "jac_reuse_conv_rate reaches its setter");
  fails += check(set_options(mem, NULL, controller_src, NARGS(controller_src)) ==
                   FIRK_SUCCESS,
                 "adapt_controller_by_name reaches its setter");

  /* errors */
  fails += check(set_options(mem, NULL, bad_name_src, NARGS(bad_name_src)) ==
                   FIRK_ILL_INPUT,
                 "an unknown name is FIRK_ILL_INPUT");
  fails += check(set_options(mem, NULL, missing_src, NARGS(missing_src)) ==
                   FIRK_ILL_INPUT,
                 "a missing value is FIRK_ILL_INPUT");
  fails += check(set_options(mem, NULL, ls_keys_src, NARGS(ls_keys_src)) ==
                   FIRKLS_LMEM_NULL,
                 "a FIRKLS key before the linear solver is FIRKLS_LMEM_NULL");

  /* FIRKLS keys with the linear solver */
  LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx);
  if (FIRKodeSetLinearSolver(mem, LS, NULL)) { return 1; }
  fails += check(set_options(mem, NULL, ls_keys_src, NARGS(ls_keys_src)) ==
                   FIRK_SUCCESS,
                 "FIRKLS keys with the linear solver succeed");
  fails += check(((FIRKLsMem)fm->lmem)->eplifac == SUN_RCONST(0.1) &&
                   ((FIRKLsMem)fm->lmem)->gstype == SUN_CLASSICAL_GS &&
                   ((FIRKLsMem)fm->lmem)->maxl == 4,
                 "eps_lin, stage_solver_gs_type and stage_solver_maxl reach "
                 "their setters");
  fails += check(set_options(mem, NULL, bad_gs_src, NARGS(bad_gs_src)) ==
                   FIRK_ILL_INPUT,
                 "an unknown Gram-Schmidt name is FIRK_ILL_INPUT");

  /* another prefix */
  fails += check(set_options(mem, "myode", custom_src, NARGS(custom_src)) ==
                     FIRK_SUCCESS &&
                   fm->maxcor == 6,
                 "a custom prefix selects its keys only");

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
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
