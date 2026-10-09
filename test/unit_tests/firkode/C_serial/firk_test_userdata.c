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
 * Unit test of FIRKodeSetUserData (design hazard H5) on
 * M y' = M diag(lambda) y with M = diag(1, 2, 3, 4), lambda = -1, -10, -100,
 * -1000, y(0) = 1, user error weights and the root of y_1 - 1/2, in three
 * configurations of the linear solvers:
 *
 * - matrix-free: SPGMR with jtsetup, jtimes, psetup and psolve; PCG mass
 *   solver with mtimes (own data) and a mass preconditioner;
 * - dense with a Jacobian function and a mass-matrix function;
 * - dense with a linear-system function and a mass-matrix function.
 *
 * SetUserData is called after every routine is attached, again between
 * two Evolve calls, and a third time before FIRKodeReInit, after which the
 * setup-time callbacks (the Jacobian, linear-system, mass and
 * preconditioner setups) run again; every callback that receives user_data
 * checks that it is the pointer last set, and mtimes that it receives its
 * own data.
 *
 * Pass criterion: no callback sees another pointer, and each callback of a
 * configuration runs after the first SetUserData, the solve-time ones after
 * the second, and all of them after the third.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunlinsol/sunlinsol_pcg.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ZERO SUN_RCONST(0.0)
#define HALF SUN_RCONST(0.5)
#define ONE  SUN_RCONST(1.0)
#define NEQ  4
#define RTOL SUN_RCONST(1.0e-6)

static const sunrealtype lambda[NEQ] = {SUN_RCONST(-1.0), SUN_RCONST(-10.0),
                                        SUN_RCONST(-100.0), SUN_RCONST(-1000.0)};

static sunrealtype mdiag(int i) { return (sunrealtype)(i + 1); }

/* the callbacks, their record, and the pointer they must receive */
enum
{
  CB_F,
  CB_EWT,
  CB_ROOT,
  CB_JAC,
  CB_LINSYS,
  CB_JTSETUP,
  CB_JTIMES,
  CB_PSETUP,
  CB_PSOLVE,
  CB_MASS,
  CB_MPSETUP,
  CB_MPSOLVE,
  CB_MTIMES,
  NCB
};

static const char* cb_names[NCB] = {"f",      "efun",        "root",
                                    "jac",    "linsys",      "jtsetup",
                                    "jtimes", "psetup",      "psolve",
                                    "mass",   "mass psetup", "mass psolve",
                                    "mtimes"};
static void* expected = NULL; /* the pointer last set               */
static int mt_data    = 0;    /* the data given to SetMassTimes     */
static int phase      = 0;    /* the SetUserData calls before, 0 to 2 */
static long int seen[3][NCB];
static long int bad[NCB];

static void record(int cb, void* user_data)
{
  seen[phase][cb]++;
  if (user_data != ((cb == CB_MTIMES) ? (void*)&mt_data : expected))
  {
    bad[cb]++;
  }
}

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  int i;
  (void)t;
  record(CB_F, user_data);
  for (i = 0; i < NEQ; i++)
  {
    NV_Ith_S(ydot, i) = mdiag(i) * lambda[i] * NV_Ith_S(y, i);
  }
  return 0;
}

static int ewt(N_Vector y, N_Vector w, void* user_data)
{
  int i;
  record(CB_EWT, user_data);
  for (i = 0; i < NEQ; i++)
  {
    NV_Ith_S(w, i) = ONE / (RTOL * SUNRabs(NV_Ith_S(y, i)) + RTOL);
  }
  return 0;
}

static int g(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)
{
  (void)t;
  record(CB_ROOT, user_data);
  gout[0] = NV_Ith_S(y, 0) - HALF;
  return 0;
}

static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int i;
  (void)t;
  (void)y;
  (void)fy;
  (void)tmp1;
  (void)tmp2;
  (void)tmp3;
  record(CB_JAC, user_data);
  for (i = 0; i < NEQ; i++) { SM_ELEMENT_D(J, i, i) = mdiag(i) * lambda[i]; }
  return 0;
}

static int linsys(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix A,
                  SUNMatrix M, sunbooleantype jok, sunbooleantype* jcur,
                  sunrealtype gamma, void* user_data, N_Vector tmp1,
                  N_Vector tmp2, N_Vector tmp3)
{
  int i;
  (void)t;
  (void)y;
  (void)fy;
  (void)M;
  (void)tmp1;
  (void)tmp2;
  (void)tmp3;
  record(CB_LINSYS, user_data);
  SUNMatZero(A);
  for (i = 0; i < NEQ; i++)
  {
    SM_ELEMENT_D(A, i, i) = mdiag(i) * (ONE - gamma * lambda[i]);
  }
  *jcur = !jok;
  return 0;
}

static int jtsetup(sunrealtype t, N_Vector y, N_Vector fy, void* user_data)
{
  (void)t;
  (void)y;
  (void)fy;
  record(CB_JTSETUP, user_data);
  return 0;
}

static int jtimes(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y,
                  N_Vector fy, void* user_data, N_Vector tmp)
{
  int i;
  (void)t;
  (void)y;
  (void)fy;
  (void)tmp;
  record(CB_JTIMES, user_data);
  for (i = 0; i < NEQ; i++)
  {
    NV_Ith_S(Jv, i) = mdiag(i) * lambda[i] * NV_Ith_S(v, i);
  }
  return 0;
}

static int psetup(sunrealtype t, N_Vector y, N_Vector fy, sunbooleantype jok,
                  sunbooleantype* jcur, sunrealtype gamma, void* user_data)
{
  (void)t;
  (void)y;
  (void)fy;
  (void)gamma;
  record(CB_PSETUP, user_data);
  *jcur = !jok;
  return 0;
}

/* the exact solve with M - gamma J */
static int psolve(sunrealtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z,
                  sunrealtype gamma, sunrealtype delta, int lr, void* user_data)
{
  int i;
  (void)t;
  (void)y;
  (void)fy;
  (void)delta;
  (void)lr;
  record(CB_PSOLVE, user_data);
  for (i = 0; i < NEQ; i++)
  {
    NV_Ith_S(z, i) = NV_Ith_S(r, i) / (mdiag(i) * (ONE - gamma * lambda[i]));
  }
  return 0;
}

static int mass(sunrealtype t, SUNMatrix M, void* user_data, N_Vector tmp1,
                N_Vector tmp2, N_Vector tmp3)
{
  int i;
  (void)t;
  (void)tmp1;
  (void)tmp2;
  (void)tmp3;
  record(CB_MASS, user_data);
  SUNMatZero(M);
  for (i = 0; i < NEQ; i++) { SM_ELEMENT_D(M, i, i) = mdiag(i); }
  return 0;
}

static int mtimes(N_Vector v, N_Vector Mv, sunrealtype t, void* mtimes_data)
{
  int i;
  (void)t;
  record(CB_MTIMES, mtimes_data);
  for (i = 0; i < NEQ; i++) { NV_Ith_S(Mv, i) = mdiag(i) * NV_Ith_S(v, i); }
  return 0;
}

static int mpsetup(sunrealtype t, void* user_data)
{
  (void)t;
  record(CB_MPSETUP, user_data);
  return 0;
}

static int mpsolve(sunrealtype t, N_Vector r, N_Vector z, sunrealtype delta,
                   int lr, void* user_data)
{
  int i;
  (void)t;
  (void)delta;
  (void)lr;
  record(CB_MPSOLVE, user_data);
  for (i = 0; i < NEQ; i++) { NV_Ith_S(z, i) = NV_Ith_S(r, i) / mdiag(i); }
  return 0;
}

/* integrate to tout through the root */
static int evolve(void* mem, sunrealtype tout, N_Vector y)
{
  int flag;
  sunrealtype tret;
  do {
    flag = FIRKodeEvolve(mem, tout, y, &tret, FIRK_NORMAL);
  }
  while (flag == FIRK_ROOT_RETURN);
  return flag;
}

/* config 0: matrix-free; 1: dense with jac; 2: dense with linsys */
static int run(SUNContext sunctx, int config)
{
  static const char* config_names[3] = {"matrix-free", "dense, jac",
                                        "dense, linsys"};
  /* the callbacks of each configuration; the ones that run at every step */
  static const int used[3][NCB]    = {{1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1},
                                      {1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                                      {1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0}};
  static const int every_step[NCB] = {1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1};
  int fails = 0, cb, ud1 = 1, ud2 = 2, ud3 = 3;
  void* mem = NULL;
  N_Vector y;
  SUNMatrix A = NULL, M = NULL;
  SUNLinearSolver LS, MLS;

  for (cb = 0; cb < NCB; cb++)
  {
    seen[0][cb] = seen[1][cb] = seen[2][cb] = 0;
    bad[cb]                                 = 0;
  }
  expected = NULL;
  phase    = 0;

  y = N_VNew_Serial(NEQ, sunctx);
  N_VConst(ONE, y);
  if (config == 0)
  {
    LS  = SUNLinSol_SPGMR(y, SUN_PREC_LEFT, 0, sunctx);
    MLS = SUNLinSol_PCG(y, SUN_PREC_LEFT, 0, sunctx);
  }
  else
  {
    A   = SUNDenseMatrix(NEQ, NEQ, sunctx);
    M   = SUNDenseMatrix(NEQ, NEQ, sunctx);
    LS  = SUNLinSol_Dense(y, A, sunctx);
    MLS = SUNLinSol_Dense(y, M, sunctx);
  }
  mem = FIRKodeCreate(sunctx);
  if (FIRKodeInit(mem, f, ZERO, y) || FIRKodeWFtolerances(mem, ewt) ||
      FIRKodeSetLinearSolver(mem, LS, A) ||
      FIRKodeSetMassLinearSolver(mem, MLS, M, SUNFALSE) ||
      FIRKodeRootInit(mem, 1, g))
  {
    return 1;
  }
  if (config == 0)
  {
    if (FIRKodeSetJacTimes(mem, jtsetup, jtimes) ||
        FIRKodeSetPreconditioner(mem, psetup, psolve) ||
        FIRKodeSetMassTimes(mem, NULL, mtimes, &mt_data) ||
        FIRKodeSetMassPreconditioner(mem, mpsetup, mpsolve))
    {
      return 1;
    }
  }
  else
  {
    if (((config == 1) ? FIRKodeSetJacFn(mem, jac)
                       : FIRKodeSetLinSysFn(mem, linsys)) ||
        FIRKodeSetMassFn(mem, mass))
    {
      return 1;
    }
  }

  /* after every routine is attached */
  expected = &ud1;
  if (FIRKodeSetUserData(mem, &ud1)) { return 1; }
  if (evolve(mem, HALF, y) != FIRK_SUCCESS)
  {
    printf("  FAIL (%s): the first Evolve failed\n", config_names[config]);
    fails++;
  }

  /* between two Evolve calls */
  expected = &ud2;
  phase    = 1;
  if (FIRKodeSetUserData(mem, &ud2)) { return 1; }
  if (evolve(mem, ONE, y) != FIRK_SUCCESS)
  {
    printf("  FAIL (%s): the second Evolve failed\n", config_names[config]);
    fails++;
  }

  /* before ReInit, after which every setup runs again */
  expected = &ud3;
  phase    = 2;
  N_VConst(ONE, y);
  if (FIRKodeSetUserData(mem, &ud3) || FIRKodeReInit(mem, ZERO, y))
  {
    return 1;
  }
  if (evolve(mem, HALF, y) != FIRK_SUCCESS)
  {
    printf("  FAIL (%s): the Evolve after ReInit failed\n", config_names[config]);
    fails++;
  }

  for (cb = 0; cb < NCB; cb++)
  {
    if (bad[cb] > 0)
    {
      printf("  FAIL (%s): H5: %s received another pointer in %ld of %ld "
             "calls\n",
             config_names[config], cb_names[cb], bad[cb],
             seen[0][cb] + seen[1][cb] + seen[2][cb]);
      fails++;
    }
    if (used[config][cb] && seen[0][cb] == 0)
    {
      printf("  FAIL (%s): %s was not called after the first SetUserData\n",
             config_names[config], cb_names[cb]);
      fails++;
    }
    if (used[config][cb] && every_step[cb] && seen[1][cb] == 0)
    {
      printf("  FAIL (%s): %s was not called after the second SetUserData\n",
             config_names[config], cb_names[cb]);
      fails++;
    }
    if (used[config][cb] && seen[2][cb] == 0)
    {
      printf("  FAIL (%s): %s was not called after the third SetUserData "
             "and ReInit\n",
             config_names[config], cb_names[cb]);
      fails++;
    }
  }

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNLinSolFree(MLS);
  if (A) { SUNMatDestroy(A); }
  if (M) { SUNMatDestroy(M); }
  N_VDestroy(y);
  return fails;
}

int main(void)
{
  SUNContext sunctx = NULL;
  int fails         = 0, config;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }
  for (config = 0; config < 3; config++) { fails += run(sunctx, config); }
  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("All tests passed\n");
  return 0;
}
