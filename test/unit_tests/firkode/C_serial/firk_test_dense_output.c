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
 * Dense output test for FIRKODE: the collocation polynomial of an s-stage
 * Radau IIA method reproduces polynomial solutions of degree <= s exactly,
 * so for y' = p'(t) with deg p <= s the values and derivatives returned by
 * FIRKodeGetDky inside a step must match p to roundoff. A polynomial of degree
 * s+1 must not be reproduced exactly (negative control), and requests outside
 * the last step must be rejected.
 *
 * Usage: firk_test_dense_output [stages]
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* integrator tolerances: fixed steps are used, so these only control the
   Newton stopping criterion and must be representable in the precision */
#if defined(SUNDIALS_SINGLE_PRECISION)
#define ITOL    SUN_RCONST(1.0e-5)
#define NLSCOEF SUN_RCONST(1.0e-1)
#else
#define ITOL    SUN_RCONST(1.0e-10)
#define NLSCOEF SUN_RCONST(1.0e-3)
#endif

/* y = (p(t), q(t)) with p of degree deg and q of degree deg-1, coefficients
   p_k = (k+1)/(k+2), q_k = 1/(k+1) */
typedef struct
{
  int deg;
} ProblemData;

static sunrealtype poly(int deg, sunrealtype t, int k, sunrealtype (*coef)(int))
{
  /* k-th derivative of sum_m coef(m) t^m */
  int m, i;
  sunrealtype sum = ZERO, fac;
  for (m = k; m <= deg; m++)
  {
    fac = ONE;
    for (i = m; i > m - k; i--) { fac *= (sunrealtype)i; }
    sum += coef(m) * fac * SUNRpowerI(t, m - k);
  }
  return sum;
}

static sunrealtype pcoef(int m)
{
  return (sunrealtype)(m + 1) / (sunrealtype)(m + 2);
}

static sunrealtype qcoef(int m) { return ONE / (sunrealtype)(m + 1); }

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  ProblemData* pd = (ProblemData*)user_data;
  (void)y;
  NV_Ith_S(ydot, 0) = poly(pd->deg, t, 1, pcoef);
  NV_Ith_S(ydot, 1) = poly(pd->deg - 1, t, 1, qcoef);
  return 0;
}

static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  (void)t;
  (void)y;
  (void)fy;
  (void)user_data;
  (void)tmp1;
  (void)tmp2;
  (void)tmp3;
  SUNMatZero(J);
  return 0;
}

static int run(SUNContext sunctx, int s, int deg, sunrealtype* maxerr0,
               sunrealtype* maxerr1, int* range_ok)
{
  void* mem;
  N_Vector y, dky;
  SUNMatrix A;
  SUNLinearSolver LS;
  ProblemData pd;
  sunrealtype t0 = SUN_RCONST(0.3), h = SUN_RCONST(0.7), t, tq, e;
  const sunrealtype theta[5] = {SUN_RCONST(0.137), SUN_RCONST(0.42),
                                SUN_RCONST(0.61), SUN_RCONST(0.83),
                                SUN_RCONST(0.97)};
  int i, retval;

  pd.deg = deg;

  y              = N_VNew_Serial(2, sunctx);
  dky            = N_VClone(y);
  NV_Ith_S(y, 0) = poly(deg, t0, 0, pcoef);
  NV_Ith_S(y, 1) = poly(deg - 1, t0, 0, qcoef);

  mem    = FIRKodeCreate(sunctx);
  retval = FIRKodeInit(mem, f, t0, y);
  retval += FIRKodeSetUserData(mem, &pd);
  retval += FIRKodeSetNumStages(mem, s);
  retval += FIRKodeSStolerances(mem, ITOL, ITOL);
  retval += FIRKodeSetFixedStep(mem, h);
  retval += FIRKodeSetNonlinConvCoef(mem, NLSCOEF);
  A  = SUNDenseMatrix(2, 2, sunctx);
  LS = SUNLinSol_Dense(y, A, sunctx);
  retval += FIRKodeSetLinearSolver(mem, LS, A);
  retval += FIRKodeSetJacFn(mem, Jac);
  if (retval != 0) { return 1; }

  /* requests before the first step must fail (except the current state) */
  *range_ok = 1;
  if (FIRKodeGetDky(mem, t0 - h, 0, dky) != FIRK_BAD_T) { *range_ok = 0; }
  if (FIRKodeGetDky(mem, t0, 0, dky) != FIRK_SUCCESS) { *range_ok = 0; }

  retval = FIRKodeEvolve(mem, t0 + h, y, &t, FIRK_ONE_STEP);
  if (retval != FIRK_SUCCESS) { return 1; }

  *maxerr0 = ZERO;
  *maxerr1 = ZERO;
  for (i = 0; i < 5; i++)
  {
    tq = t0 + theta[i] * h;

    retval = FIRKodeGetDky(mem, tq, 0, dky);
    if (retval != FIRK_SUCCESS) { return 1; }
    e = SUNRabs(NV_Ith_S(dky, 0) - poly(deg, tq, 0, pcoef)) /
        (ONE + SUNRabs(poly(deg, tq, 0, pcoef)));
    *maxerr0 = SUNMAX(*maxerr0, e);
    e        = SUNRabs(NV_Ith_S(dky, 1) - poly(deg - 1, tq, 0, qcoef)) /
        (ONE + SUNRabs(poly(deg - 1, tq, 0, qcoef)));
    *maxerr0 = SUNMAX(*maxerr0, e);

    retval = FIRKodeGetDky(mem, tq, 1, dky);
    if (retval != FIRK_SUCCESS) { return 1; }
    e = SUNRabs(NV_Ith_S(dky, 0) - poly(deg, tq, 1, pcoef)) /
        (ONE + SUNRabs(poly(deg, tq, 1, pcoef)));
    *maxerr1 = SUNMAX(*maxerr1, e);
    e        = SUNRabs(NV_Ith_S(dky, 1) - poly(deg - 1, tq, 1, qcoef)) /
        (ONE + SUNRabs(poly(deg - 1, tq, 1, qcoef)));
    *maxerr1 = SUNMAX(*maxerr1, e);
  }

  /* out-of-range requests and illegal derivative orders are rejected */
  if (FIRKodeGetDky(mem, t0 - h, 0, dky) != FIRK_BAD_T) { *range_ok = 0; }
  if (FIRKodeGetDky(mem, t0 + SUN_RCONST(2.0) * h, 0, dky) != FIRK_BAD_T)
  {
    *range_ok = 0;
  }
  if (FIRKodeGetDky(mem, t, s + 1, dky) != FIRK_BAD_K) { *range_ok = 0; }
  if (FIRKodeGetDky(mem, t, s, dky) != FIRK_SUCCESS) { *range_ok = 0; }

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  N_VDestroy(dky);
  return 0;
}

int main(int argc, char* argv[])
{
  SUNContext sunctx;
  int s, range_ok, fails = 0;
  sunrealtype err0, err1, tol0, tol1;

  s = (argc > 1) ? atoi(argv[1]) : 3;
#if defined(SUNDIALS_SINGLE_PRECISION)
  tol0 = SUN_RCONST(1.0e-4);
  tol1 = SUN_RCONST(1.0e-3);
#else
  tol0 = SUN_RCONST(1.0e4) * SUN_UNIT_ROUNDOFF;
  tol1 = SUN_RCONST(1.0e5) * SUN_UNIT_ROUNDOFF;
#endif

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  printf("Dense output test, Radau IIA with s = %d stages\n", s);

  /* polynomial of degree s: reproduced exactly */
  if (run(sunctx, s, s, &err0, &err1, &range_ok))
  {
    printf("  FAIL: integration failed\n");
    return 1;
  }
  printf("  degree %d: value error = " SUN_FORMAT_E
         ", derivative error = " SUN_FORMAT_E "\n",
         s, err0, err1);
  if (err0 > tol0)
  {
    printf("  FAIL: dense output value error exceeds " SUN_FORMAT_E "\n", tol0);
    fails++;
  }
  if (err1 > tol1)
  {
    printf("  FAIL: dense output derivative error exceeds " SUN_FORMAT_E "\n",
           tol1);
    fails++;
  }
  if (!range_ok)
  {
    printf("  FAIL: out-of-range or illegal requests were not rejected\n");
    fails++;
  }

  /* polynomial of degree s+1: not reproduced exactly */
  if (run(sunctx, s, s + 1, &err0, &err1, &range_ok))
  {
    printf("  FAIL: integration failed\n");
    return 1;
  }
  printf("  degree %d: value error = " SUN_FORMAT_E
         ", derivative error = " SUN_FORMAT_E "\n",
         s + 1, err0, err1);
  if (err0 < SUN_RCONST(1.0e-8))
  {
    printf("  FAIL: degree %d polynomial unexpectedly reproduced\n", s + 1);
    fails++;
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
