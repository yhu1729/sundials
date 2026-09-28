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
 * Dahlquist test for FIRKODE: one fixed step of the s-stage Radau IIA method
 * applied to the linear system y' = L y must reproduce y1 = R(h L) y0 with the
 * stability function R(z) = 1 + z b^T (I - z A)^{-1} 1, evaluated here by
 * solving the s*d x s*d stage system directly. Real eigenvalues (d = 1)
 * including the stiff limit (L-stability), and complex conjugate pairs
 * (d = 2 rotation systems) are tested. For a linear problem the Newton
 * iteration must converge in at most two iterations.
 *
 * Usage: firk_test_dahlquist [stages]
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"
#include "firkode/firkode_tables.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* integrator tolerances: fixed steps are used, so these only control the
   Newton stopping criterion and must be representable in the precision */
#if defined(SUNDIALS_SINGLE_PRECISION)
#define ITOL SUN_RCONST(1.0e-6)
#else
#define ITOL SUN_RCONST(1.0e-12)
#endif

/* problem data: y' = L y with L a d x d matrix (row-major, d <= 2) */
typedef struct
{
  int d;
  sunrealtype L[4];
} ProblemData;

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  ProblemData* pd = (ProblemData*)user_data;
  sunrealtype* yd = N_VGetArrayPointer(y);
  sunrealtype* fd = N_VGetArrayPointer(ydot);
  int i, j, d = pd->d;
  (void)t;
  for (i = 0; i < d; i++)
  {
    fd[i] = ZERO;
    for (j = 0; j < d; j++) { fd[i] += pd->L[i * d + j] * yd[j]; }
  }
  return 0;
}

static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  ProblemData* pd = (ProblemData*)user_data;
  int i, j, d = pd->d;
  (void)t;
  (void)y;
  (void)fy;
  (void)tmp1;
  (void)tmp2;
  (void)tmp3;
  for (i = 0; i < d; i++)
  {
    for (j = 0; j < d; j++) { SM_ELEMENT_D(J, i, j) = pd->L[i * d + j]; }
  }
  return 0;
}

/* Solve the dense n x n system a x = b in place (partial pivoting) */
static int dense_solve(sunrealtype* a, sunrealtype* b, int n)
{
  int i, j, k, p;
  sunrealtype amax, tmp;
  for (k = 0; k < n; k++)
  {
    p    = k;
    amax = SUNRabs(a[k * n + k]);
    for (i = k + 1; i < n; i++)
    {
      if (SUNRabs(a[i * n + k]) > amax)
      {
        amax = SUNRabs(a[i * n + k]);
        p    = i;
      }
    }
    if (amax == ZERO) { return -1; }
    if (p != k)
    {
      for (j = 0; j < n; j++)
      {
        tmp          = a[k * n + j];
        a[k * n + j] = a[p * n + j];
        a[p * n + j] = tmp;
      }
      tmp  = b[k];
      b[k] = b[p];
      b[p] = tmp;
    }
    for (i = k + 1; i < n; i++)
    {
      tmp = a[i * n + k] / a[k * n + k];
      for (j = k; j < n; j++) { a[i * n + j] -= tmp * a[k * n + j]; }
      b[i] -= tmp * b[k];
    }
  }
  for (i = n - 1; i >= 0; i--)
  {
    for (j = i + 1; j < n; j++) { b[i] -= a[i * n + j] * b[j]; }
    b[i] /= a[i * n + i];
  }
  return 0;
}

/* Reference solution of one Radau IIA step for y' = L y: solve the stage
   system (I - h A (x) L) Y = 1 (x) y0, then y1 = y0 + h sum_i b_i L Y_i */
static int reference_step(FIRKodeTable T, ProblemData* pd, sunrealtype h,
                          const sunrealtype* y0, sunrealtype* y1)
{
  int s = T->s, d = pd->d, n = s * d, i, j, k, l, retval;
  sunrealtype* K = (sunrealtype*)calloc((size_t)(n * n), sizeof(sunrealtype));
  sunrealtype* Y = (sunrealtype*)calloc((size_t)n, sizeof(sunrealtype));

  for (i = 0; i < s; i++)
  {
    for (k = 0; k < d; k++)
    {
      K[(i * d + k) * n + (i * d + k)] = ONE;
      Y[i * d + k]                     = y0[k];
      for (j = 0; j < s; j++)
      {
        for (l = 0; l < d; l++)
        {
          K[(i * d + k) * n + (j * d + l)] -= h * T->A[i][j] * pd->L[k * d + l];
        }
      }
    }
  }
  retval = dense_solve(K, Y, n);
  for (k = 0; k < d; k++)
  {
    y1[k] = y0[k];
    for (i = 0; i < s; i++)
    {
      for (l = 0; l < d; l++)
      {
        y1[k] += h * T->b[i] * pd->L[k * d + l] * Y[i * d + l];
      }
    }
  }
  free(K);
  free(Y);
  return retval;
}

static int run_case(SUNContext sunctx, int s, ProblemData* pd,
                    const sunrealtype* y0, sunrealtype h, sunrealtype tol,
                    const char* label)
{
  void* mem;
  N_Vector y;
  SUNMatrix A;
  SUNLinearSolver LS;
  FIRKodeTable T;
  sunrealtype t, ref[2], err, nrm;
  long int nni, nnf, nst;
  int k, retval, fails = 0;

  y = N_VNew_Serial(pd->d, sunctx);
  for (k = 0; k < pd->d; k++) { NV_Ith_S(y, k) = y0[k]; }

  mem    = FIRKodeCreate(sunctx);
  retval = FIRKodeInit(mem, f, ZERO, y);
  retval += FIRKodeSetUserData(mem, pd);
  retval += FIRKodeSetNumStages(mem, s);
  retval += FIRKodeSStolerances(mem, ITOL, ITOL);
  retval += FIRKodeSetFixedStep(mem, h);
  retval += FIRKodeSetNonlinConvCoef(mem, SUN_RCONST(1.0e-2));
  A  = SUNDenseMatrix(pd->d, pd->d, sunctx);
  LS = SUNLinSol_Dense(y, A, sunctx);
  retval += FIRKodeSetLinearSolver(mem, LS, A);
  retval += FIRKodeSetJacFn(mem, Jac);
  retval += FIRKodeSetStageEpsLin(mem, SUN_RCONST(1.0e-3));
  if (retval != 0)
  {
    printf("  FAIL (%s): setup returned %d\n", label, retval);
    return 1;
  }

  retval = FIRKodeEvolve(mem, h, y, &t, FIRK_ONE_STEP);
  if (retval != FIRK_SUCCESS)
  {
    printf("  FAIL (%s): FIRKodeEvolve returned %d\n", label, retval);
    fails++;
  }

  FIRKodeGetCurrentTable(mem, &T);
  reference_step(T, pd, h, y0, ref);

  err = ZERO;
  nrm = ZERO;
  for (k = 0; k < pd->d; k++)
  {
    err = SUNMAX(err, SUNRabs(NV_Ith_S(y, k) - ref[k]));
    nrm = SUNMAX(nrm, SUNRabs(ref[k]));
  }
  if (err > tol * (ONE + nrm))
  {
    printf("  FAIL (%s): |y1 - R(hL) y0| = " SUN_FORMAT_E " > " SUN_FORMAT_E "\n",
           label, err, tol * (ONE + nrm));
    fails++;
  }
  if (SUNRabs(t - h) > SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF * h)
  {
    printf("  FAIL (%s): returned t = " SUN_FORMAT_E "\n", label, t);
    fails++;
  }

  FIRKodeGetNumSteps(mem, &nst);
  FIRKodeGetNonlinSolvStats(mem, &nni, &nnf);
  if (nst != 1 || nni > 2 || nnf != 0)
  {
    printf("  FAIL (%s): nst = %li, nni = %li, nnf = %li\n", label, nst, nni,
           nnf);
    fails++;
  }

  printf("  %-28s |y1 - ref| = " SUN_FORMAT_E "  |y1| = " SUN_FORMAT_E
         "  nni = %li\n",
         label, err, SUNRabs(NV_Ith_S(y, 0)), nni);

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  return fails;
}

int main(int argc, char* argv[])
{
  SUNContext sunctx;
  ProblemData pd;
  sunrealtype y0[2], h, tol, y1;
  int s, i, fails = 0;
  char label[64];
  const sunrealtype lambdas[5]  = {SUN_RCONST(-1.0), SUN_RCONST(-10.0),
                                   SUN_RCONST(-1.0e3), SUN_RCONST(-1.0e6),
                                   SUN_RCONST(-1.0e8)};
  const sunrealtype pairs[3][2] = {{SUN_RCONST(-1.0), SUN_RCONST(10.0)},
                                   {SUN_RCONST(-100.0), SUN_RCONST(1000.0)},
                                   {SUN_RCONST(0.0), SUN_RCONST(1.0)}};

  s = (argc > 1) ? atoi(argv[1]) : 3;
  h = SUN_RCONST(0.1);
#if defined(SUNDIALS_SINGLE_PRECISION)
  tol = SUN_RCONST(1.0e-4);
#else
  tol = SUN_RCONST(1.0e5) * SUN_UNIT_ROUNDOFF;
#endif

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  printf("Dahlquist test, Radau IIA with s = %d, h = " SUN_FORMAT_G "\n", s, h);

  /* real eigenvalues */
  pd.d  = 1;
  y0[0] = ONE;
  for (i = 0; i < 5; i++)
  {
#if defined(SUNDIALS_SINGLE_PRECISION)
    if (i >= 3) { break; }
#endif
    pd.L[0] = lambdas[i];
    snprintf(label, sizeof(label), "lambda = " SUN_FORMAT_G, lambdas[i]);
    fails += run_case(sunctx, s, &pd, y0, h, tol, label);
  }

  /* L-stability: |R(z)| -> 0 as z -> -infinity */
  {
    void* mem;
    N_Vector y         = N_VNew_Serial(1, sunctx);
    SUNMatrix A        = SUNDenseMatrix(1, 1, sunctx);
    SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
    sunrealtype t;
    pd.d           = 1;
    pd.L[0]        = SUN_RCONST(-1.0e8);
    NV_Ith_S(y, 0) = ONE;
    mem            = FIRKodeCreate(sunctx);
    FIRKodeInit(mem, f, ZERO, y);
    FIRKodeSetUserData(mem, &pd);
    FIRKodeSetNumStages(mem, s);
    FIRKodeSStolerances(mem, ITOL, ITOL);
    FIRKodeSetFixedStep(mem, h);
    FIRKodeSetLinearSolver(mem, LS, A);
    FIRKodeSetJacFn(mem, Jac);
    FIRKodeEvolve(mem, h, y, &t, FIRK_ONE_STEP);
    y1 = SUNRabs(NV_Ith_S(y, 0));
    if (y1 > SUN_RCONST(1.0e-5))
    {
      printf("  FAIL (L-stability): |R(-1e7)| = " SUN_FORMAT_E "\n", y1);
      fails++;
    }
    else { printf("  L-stability: |R(-1e7)| = " SUN_FORMAT_E "\n", y1); }
    FIRKodeFree(&mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    N_VDestroy(y);
  }

  /* complex conjugate eigenvalue pairs alpha +/- i beta */
  pd.d  = 2;
  y0[0] = ONE;
  y0[1] = SUN_RCONST(0.5);
  for (i = 0; i < 3; i++)
  {
    pd.L[0] = pairs[i][0];
    pd.L[1] = -pairs[i][1];
    pd.L[2] = pairs[i][1];
    pd.L[3] = pairs[i][0];
    snprintf(label, sizeof(label), "alpha,beta = " SUN_FORMAT_G "," SUN_FORMAT_G,
             pairs[i][0], pairs[i][1]);
    fails += run_case(sunctx, s, &pd, y0, h, tol, label);
  }

  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("\nAll tests passed\n");
  return 0;
}
