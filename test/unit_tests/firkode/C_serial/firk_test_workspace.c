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
 * Unit test of the vector budget of design Table 8.7 (hazard H27): the
 * vectors of length N that FIRKODE holds, counted as the live clones of y0
 * (a clone of a stacked vector clones each of its stages), after the first
 * step:
 *
 *   11 + s (2 maxl + 10) for s >= 2, 16 for s = 1,
 *   + 2 with a mass matrix, + 1 with vector absolute tolerances,
 *
 * for SPGMR and dense block solvers, the default and a user maxl, a mass
 * matrix and the one-stage method; the user's solvers clone another
 * template and are not counted. FIRKodeFree releases every vector.
 *
 * Pass criterion: every count equals its budget.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>

#include "firkode/firkode.h"
#include "firkode/firkode_impl.h" /* firkSetNumStagesAny */
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define NEQ  4

static const sunrealtype lambda[NEQ] = {SUN_RCONST(-1.0), SUN_RCONST(-10.0),
                                        SUN_RCONST(-100.0), SUN_RCONST(-1000.0)};

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  int i;
  (void)t;
  (void)user_data;
  for (i = 0; i < NEQ; i++) { NV_Ith_S(ydot, i) = lambda[i] * NV_Ith_S(y, i); }
  return 0;
}

static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int i;
  (void)t;
  (void)y;
  (void)fy;
  (void)user_data;
  (void)tmp1;
  (void)tmp2;
  (void)tmp3;
  for (i = 0; i < NEQ; i++) { SM_ELEMENT_D(J, i, i) = lambda[i]; }
  return 0;
}

static int mass(sunrealtype t, SUNMatrix M, void* user_data, N_Vector tmp1,
                N_Vector tmp2, N_Vector tmp3)
{
  int i;
  (void)t;
  (void)user_data;
  (void)tmp1;
  (void)tmp2;
  (void)tmp3;
  SUNMatZero(M);
  for (i = 0; i < NEQ; i++) { SM_ELEMENT_D(M, i, i) = ONE; }
  return 0;
}

/* the live clones of y0 */
static long int nlive = 0;

static void destroy_counting(N_Vector v)
{
  nlive--;
  N_VDestroy_Serial(v);
}

static N_Vector clone_counting(N_Vector w)
{
  N_Vector v = N_VClone_Serial(w); /* copies the operations of w */
  if (v == NULL) { return NULL; }
  v->ops->nvdestroy = destroy_counting;
  nlive++;
  return v;
}

struct config
{
  const char* name;
  int s;
  int dense;
  int maxl; /* 0: the default max(5, 3s) */
  int mass;
  int vtol;
};

static int run(SUNContext sunctx, const struct config* c)
{
  int fails = 0, maxl;
  long int budget;
  sunrealtype tret;
  N_Vector y, tmpl, abstol = NULL;
  SUNMatrix A = NULL, M   = NULL;
  SUNLinearSolver LS, MLS = NULL;
  void* mem;

  /* y0 counts its clones; the user's objects clone another template */
  y               = N_VNew_Serial(NEQ, sunctx);
  tmpl            = N_VNew_Serial(NEQ, sunctx);
  y->ops->nvclone = clone_counting;
  nlive           = 0;
  N_VConst(ONE, y);
  if (c->dense)
  {
    A  = SUNDenseMatrix(NEQ, NEQ, sunctx);
    LS = SUNLinSol_Dense(tmpl, A, sunctx);
  }
  else { LS = SUNLinSol_SPGMR(tmpl, SUN_PREC_NONE, 0, sunctx); }
  if (c->mass)
  {
    M   = SUNDenseMatrix(NEQ, NEQ, sunctx);
    MLS = SUNLinSol_Dense(tmpl, M, sunctx);
  }

  mem = FIRKodeCreate(sunctx);
  if (FIRKodeInit(mem, f, ZERO, y)) { return 1; }
  if (c->vtol)
  {
    abstol = N_VNew_Serial(NEQ, sunctx);
    N_VConst(SUN_RCONST(1.0e-6), abstol);
    if (FIRKodeSVtolerances(mem, SUN_RCONST(1.0e-6), abstol)) { return 1; }
  }
  else if (FIRKodeSStolerances(mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-6)))
  {
    return 1;
  }
  if (FIRKodeSetLinearSolver(mem, LS, A) ||
      (c->dense && FIRKodeSetJacFn(mem, jac)) ||
      (c->mass && (FIRKodeSetMassLinearSolver(mem, MLS, M, SUNFALSE) ||
                   FIRKodeSetMassFn(mem, mass))) ||
      (c->s != 3 && firkSetNumStagesAny(mem, c->s)) ||
      (c->maxl > 0 && (FIRKodeSetMode(mem, FIRK_EXPERT) ||
                       FIRKodeSetStageSolverMaxl(mem, c->maxl))))
  {
    return 1;
  }
  if (FIRKodeEvolve(mem, ONE, y, &tret, FIRK_ONE_STEP) != FIRK_SUCCESS)
  {
    printf("  FAIL (%s): the first step failed\n", c->name);
    return 1;
  }

  maxl   = (c->maxl > 0) ? c->maxl : SUNMAX(5, 3 * c->s);
  budget = (c->s == 1) ? 16 : 11 + (long int)c->s * (2 * maxl + 10);
  if (c->mass) { budget += 2; }
  if (c->vtol) { budget += 1; }
  if (nlive != budget)
  {
    printf("  FAIL (%s): H27: %ld vectors, the budget is %ld\n", c->name, nlive,
           budget);
    fails++;
  }

  FIRKodeFree(&mem);
  if (nlive != 0)
  {
    printf("  FAIL (%s): %ld vectors left after FIRKodeFree\n", c->name, nlive);
    fails++;
  }

  SUNLinSolFree(LS);
  if (MLS) { SUNLinSolFree(MLS); }
  if (A) { SUNMatDestroy(A); }
  if (M) { SUNMatDestroy(M); }
  if (abstol) { N_VDestroy(abstol); }
  N_VDestroy(tmpl);
  N_VDestroy(y);
  return fails;
}

int main(void)
{
  SUNContext sunctx                    = NULL;
  int fails                            = 0, k;
  static const struct config configs[] = {{"SPGMR", 3, 0, 0, 0, 0},
                                          {"SPGMR, maxl = 4", 3, 0, 4, 0, 0},
                                          {"dense", 3, 1, 0, 0, 0},
                                          {"dense, mass matrix", 3, 1, 0, 1, 0},
                                          {"SPGMR, vector tolerances", 3, 0, 0,
                                           0, 1},
                                          {"dense, one stage", 1, 1, 0, 0, 0},
                                          {"dense, one stage, mass matrix", 1,
                                           1, 0, 1, 0}};

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }
  for (k = 0; k < (int)(sizeof(configs) / sizeof(configs[0])); k++)
  {
    fails += run(sunctx, &configs[k]);
  }
  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("All tests passed\n");
  return 0;
}
