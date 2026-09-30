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
 * Unit test for the FIRKODE Radau IIA tables: order conditions (B, C),
 * collocation structure, stiff accuracy, L-stability (R(inf) = 0), the
 * derived quantities (A^{-1}, d, e, P, gamma0), agreement between the
 * closed-form and numerically computed tables, and FIRKodeTable_CheckOrder.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "firkode/firkode_tables.h"
#include "firkode/firkode_tables_impl.h"
#include "sundials/sundials_math.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* Tolerance scaling: base * s^2 * unit roundoff */
static sunrealtype tolerance(int s)
{
  return SUN_RCONST(1000.0) * (sunrealtype)(s * s) * SUN_UNIT_ROUNDOFF;
}

static int check(sunbooleantype ok, const char* what, int s, sunrealtype err,
                 sunrealtype tol)
{
  if (!ok)
  {
    printf("  FAIL: %s (s = %d): error = " SUN_FORMAT_E ", tol = " SUN_FORMAT_E
           "\n",
           what, s, err, tol);
    return 1;
  }
  return 0;
}

/* max_i |x_i - y_i| */
static sunrealtype max_diff(const sunrealtype* x, const sunrealtype* y, int n)
{
  int i;
  sunrealtype err = ZERO;
  for (i = 0; i < n; i++) { err = SUNMAX(err, SUNRabs(x[i] - y[i])); }
  return err;
}

/* Smallest |pivot| / largest |pivot| in an LU factorization of B, used to
   test (near) singularity of small matrices */
static sunrealtype lu_pivot_ratio(sunrealtype** B, int s)
{
  int i, j, k, p;
  sunrealtype pmin = SUN_BIG_REAL, pmax = ZERO, amax, tmp;
  sunrealtype** a = (sunrealtype**)malloc((size_t)s * sizeof(sunrealtype*));
  for (i = 0; i < s; i++)
  {
    a[i] = (sunrealtype*)malloc((size_t)s * sizeof(sunrealtype));
    for (j = 0; j < s; j++) { a[i][j] = B[i][j]; }
  }
  for (k = 0; k < s; k++)
  {
    p    = k;
    amax = SUNRabs(a[k][k]);
    for (i = k + 1; i < s; i++)
    {
      if (SUNRabs(a[i][k]) > amax)
      {
        amax = SUNRabs(a[i][k]);
        p    = i;
      }
    }
    if (p != k)
    {
      for (j = 0; j < s; j++)
      {
        tmp     = a[k][j];
        a[k][j] = a[p][j];
        a[p][j] = tmp;
      }
    }
    pmin = SUNMIN(pmin, amax);
    pmax = SUNMAX(pmax, amax);
    if (amax == ZERO) { break; }
    for (i = k + 1; i < s; i++)
    {
      a[i][k] /= a[k][k];
      for (j = k + 1; j < s; j++) { a[i][j] -= a[i][k] * a[k][j]; }
    }
  }
  for (i = 0; i < s; i++) { free(a[i]); }
  free(a);
  return (pmax > ZERO) ? pmin / pmax : ZERO;
}

static int test_table(int s)
{
  int i, j, k, fails = 0, q, p, retval;
  sunrealtype tol, err, sum, cpow, scale, gamma0_ref, ratio;
  sunrealtype* x;
  sunrealtype** B;
  FIRKodeTable T, Tc, Td;

  /* Reference shifts for even s (no real eigenvalue of A^{-1}) and for s = 3
     from scripts/firkode_radau_tables.py; for odd s the real eigenvalue of
     A^{-1} is verified directly below */
  const sunrealtype gamma0_even[FIRK_MAX_STAGES + 1] =
    {ZERO,
     ONE,
     SUN_RCONST(0.5),
     SUN_RCONST(0.2748888295956773677478286),
     SUN_RCONST(0.208890675278273738551678),
     ZERO,
     SUN_RCONST(0.1334999853021741882617498),
     ZERO,
     SUN_RCONST(0.09833377347648450252097062),
     ZERO};
  /* RADAU5 error estimate constants DD1..DD3 = e/gamma0 for s = 3 */
  const sunrealtype dd3[3] = {SUN_RCONST(-10.048809399827415562),
                              SUN_RCONST(1.3821427331607488958),
                              SUN_RCONST(-0.33333333333333333333)};

  tol = tolerance(s);
  printf("Testing Radau IIA table with s = %d stages (tol = " SUN_FORMAT_E
         ")\n",
         s, tol);

  T = FIRKodeTable_RadauIIA(s);
  if (T == NULL)
  {
    printf("  FAIL: FIRKodeTable_RadauIIA(%d) returned NULL\n", s);
    return 1;
  }

  /* metadata */
  fails += check(T->s == s, "stage count", s, ZERO, ZERO);
  fails += check(T->q == 2 * s - 1, "method order", s, ZERO, ZERO);
  fails += check(T->p == s, "estimate order", s, ZERO, ZERO);
  fails += check(T->stiffly_accurate, "stiffly accurate flag", s, ZERO, ZERO);

  /* nodes: increasing in (0, 1], last node exactly 1 */
  for (i = 0; i < s; i++)
  {
    fails += check(T->c[i] > ZERO && T->c[i] <= ONE, "node range", s, T->c[i],
                   ONE);
    if (i > 0)
    {
      fails += check(T->c[i] > T->c[i - 1], "node ordering", s, T->c[i], ZERO);
    }
  }
  fails += check(T->c[s - 1] == ONE, "last node equals one", s,
                 SUNRabs(T->c[s - 1] - ONE), ZERO);

  /* row sums A 1 = c and b^T 1 = 1 */
  err = ZERO;
  for (i = 0; i < s; i++)
  {
    sum = ZERO;
    for (j = 0; j < s; j++) { sum += T->A[i][j]; }
    err = SUNMAX(err, SUNRabs(sum - T->c[i]));
  }
  fails += check(err <= tol, "row sums of A", s, err, tol);
  sum = ZERO;
  for (i = 0; i < s; i++) { sum += T->b[i]; }
  fails += check(SUNRabs(sum - ONE) <= tol, "sum of b", s, SUNRabs(sum - ONE),
                 tol);

  /* stiff accuracy: b equals the last row of A */
  err = max_diff(T->b, T->A[s - 1], s);
  fails += check(err <= tol, "stiff accuracy", s, err, tol);

  /* B(2s-1): sum_i b_i c_i^{k-1} = 1/k */
  err = ZERO;
  for (k = 1; k <= 2 * s - 1; k++)
  {
    sum = ZERO;
    for (i = 0; i < s; i++) { sum += T->b[i] * SUNRpowerI(T->c[i], k - 1); }
    err = SUNMAX(err, SUNRabs(sum - ONE / (sunrealtype)k));
  }
  fails += check(err <= tol, "quadrature conditions B(2s-1)", s, err, tol);

  /* C(s): sum_j a_ij c_j^{k-1} = c_i^k / k */
  err = ZERO;
  for (k = 1; k <= s; k++)
  {
    for (i = 0; i < s; i++)
    {
      sum = ZERO;
      for (j = 0; j < s; j++)
      {
        sum += T->A[i][j] * SUNRpowerI(T->c[j], k - 1);
      }
      cpow = SUNRpowerI(T->c[i], k) / (sunrealtype)k;
      err  = SUNMAX(err, SUNRabs(sum - cpow));
    }
  }
  fails += check(err <= tol, "collocation conditions C(s)", s, err, tol);

  /* A^{-1} A = I */
  err = ZERO;
  for (i = 0; i < s; i++)
  {
    for (j = 0; j < s; j++)
    {
      sum = ZERO;
      for (k = 0; k < s; k++) { sum += T->Ainv[i][k] * T->A[k][j]; }
      err = SUNMAX(err, SUNRabs(sum - ((i == j) ? ONE : ZERO)));
    }
  }
  fails += check(err <= tol, "A^{-1} A = I", s, err, tol);

  /* d = b^T A^{-1} equals the last unit vector for stiffly accurate tables */
  err = ZERO;
  for (j = 0; j < s; j++)
  {
    err = SUNMAX(err, SUNRabs(T->d[j] - ((j == s - 1) ? ONE : ZERO)));
  }
  fails += check(err <= tol, "d = b^T A^{-1} = e_s", s, err, tol);

  /* L-stability: R(inf) = 1 - b^T A^{-1} 1 = 0 */
  x = (sunrealtype*)calloc((size_t)s, sizeof(sunrealtype));
  for (i = 0; i < s; i++)
  {
    for (j = 0; j < s; j++) { x[i] += T->Ainv[i][j]; }
  }
  sum = ONE;
  for (i = 0; i < s; i++) { sum -= T->b[i] * x[i]; }
  fails += check(SUNRabs(sum) <= tol, "R(inf) = 0", s, SUNRabs(sum), tol);
  free(x);

  /* error estimate weights: sum_i (e A)_i c_i^{k-1} = -gamma0 delta_{k1},
     i.e. the embedded weights bhat = b - gamma0 l(0) satisfy the quadrature
     conditions with an explicit first stage */
  err = ZERO;
  for (k = 1; k <= s; k++)
  {
    sum = ZERO;
    for (i = 0; i < s; i++)
    {
      cpow = ZERO;
      for (j = 0; j < s; j++) { cpow += T->e[j] * T->A[j][i]; }
      sum += cpow * SUNRpowerI(T->c[i], k - 1);
    }
    if (k == 1) { sum += T->gamma0; }
    err = SUNMAX(err, SUNRabs(sum));
  }
  fails += check(err <= tol, "error estimate weights e", s, err, tol);
  cpow = ((s % 2) == 1) ? -ONE : ONE;
  err  = SUNRabs(T->e[s - 1] - cpow * T->gamma0 / (sunrealtype)s);
  fails += check(err <= tol, "e_s = (-1)^s gamma0/s", s, err, tol);

  /* dense output coefficients: L_i(c_j) = delta_ij, L_i(1) = delta_is; the
     monomial coefficients grow with s so the tolerance is scaled by |P| */
  err   = ZERO;
  scale = ONE;
  for (i = 0; i < s; i++)
  {
    for (j = 0; j < s; j++) { scale = SUNMAX(scale, SUNRabs(T->P[i][j])); }
  }
  for (i = 0; i < s; i++)
  {
    for (j = 0; j < s; j++)
    {
      sum  = ZERO;
      cpow = T->c[j];
      for (k = 0; k < s; k++)
      {
        sum += T->P[k][i] * cpow;
        cpow *= T->c[j];
      }
      err = SUNMAX(err, SUNRabs(sum - ((i == j) ? ONE : ZERO)));
    }
  }
  fails += check(err <= tol * scale, "dense output coefficients P", s, err,
                 tol * scale);

  /* gamma0 */
  if (s % 2 == 1 && s >= 3)
  {
    /* A^{-1} - (1/gamma0) I must be singular */
    B = (sunrealtype**)malloc((size_t)s * sizeof(sunrealtype*));
    for (i = 0; i < s; i++)
    {
      B[i] = (sunrealtype*)malloc((size_t)s * sizeof(sunrealtype));
      for (j = 0; j < s; j++)
      {
        B[i][j] = T->Ainv[i][j] - ((i == j) ? ONE / T->gamma0 : ZERO);
      }
    }
    ratio = lu_pivot_ratio(B, s);
    fails += check(ratio <= SUN_RCONST(1.0e6) * tol,
                   "1/gamma0 is a real eigenvalue of A^{-1}", s, ratio,
                   SUN_RCONST(1.0e6) * tol);
    for (i = 0; i < s; i++) { free(B[i]); }
    free(B);
  }
  gamma0_ref = gamma0_even[s];
  if (gamma0_ref > ZERO)
  {
    err = SUNRabs(T->gamma0 - gamma0_ref);
    fails += check(err <= tol, "gamma0 reference value", s, err, tol);
  }

  /* RADAU5 error estimate constants for s = 3 */
  if (s == 3)
  {
    err = ZERO;
    for (i = 0; i < 3; i++)
    {
      err = SUNMAX(err, SUNRabs(T->e[i] / T->gamma0 - dd3[i]));
    }
    fails += check(err <= SUN_RCONST(10.0) * tol, "RADAU5 DD constants", s, err,
                   SUN_RCONST(10.0) * tol);
  }

  /* closed-form and numerically computed tables agree */
  Tc = firkTable_RadauIIA_Computed(s);
  if (Tc == NULL)
  {
    printf("  FAIL: firkTable_RadauIIA_Computed(%d) returned NULL\n", s);
    fails++;
  }
  else
  {
    err = max_diff(T->c, Tc->c, s);
    for (i = 0; i < s; i++)
    {
      err = SUNMAX(err, max_diff(T->A[i], Tc->A[i], s));
    }
    fails += check(err <= SUN_RCONST(10.0) * SUN_UNIT_ROUNDOFF,
                   "closed-form vs computed table", s, err,
                   SUN_RCONST(10.0) * SUN_UNIT_ROUNDOFF);
    FIRKodeTable_Free(Tc);
  }

  /* recomputing the derived quantities from the stored c, A, b reproduces
     them to working precision (relative to the size of A^{-1}) */
  Td = FIRKodeTable_Copy(T);
  if (Td == NULL || firkTable_Derive(Td) != 0)
  {
    printf("  FAIL: firkTable_Derive(%d) failed\n", s);
    fails++;
  }
  else
  {
    scale = ZERO;
    for (i = 0; i < s; i++)
    {
      for (j = 0; j < s; j++) { scale = SUNMAX(scale, SUNRabs(T->Ainv[i][j])); }
    }
    err = max_diff(T->d, Td->d, s);
    err = SUNMAX(err, max_diff(T->e, Td->e, s));
    for (i = 0; i < s; i++)
    {
      err = SUNMAX(err, max_diff(T->Ainv[i], Td->Ainv[i], s));
      err = SUNMAX(err, max_diff(T->P[i], Td->P[i], s));
    }
    fails += check(err <= tol * scale, "firkTable_Derive reproducibility", s,
                   err, tol * scale);
  }
  FIRKodeTable_Free(Td);

  /* order check (in single precision the order conditions can only be
     resolved for s <= 3) */
  retval = FIRKodeTable_CheckOrder(T, &q, &p, NULL);
#if defined(SUNDIALS_SINGLE_PRECISION)
  if (s <= 3)
#endif
  {
    fails += check(retval == 0 && q == 2 * s - 1 && p == s,
                   "FIRKodeTable_CheckOrder", s, (sunrealtype)q, (sunrealtype)p);
  }

  FIRKodeTable_Free(T);

  if (fails == 0) { printf("  passed\n"); }
  return fails;
}

int main(void)
{
  int s, smax, fails = 0;
  FIRKodeTable T;

#if defined(SUNDIALS_SINGLE_PRECISION)
  smax = 4;
#else
  smax = FIRK_MAX_STAGES;
#endif

  for (s = 1; s <= smax; s++) { fails += test_table(s); }

  /* invalid stage counts are rejected */
  T = FIRKodeTable_RadauIIA(0);
  if (T != NULL)
  {
    printf("  FAIL: FIRKodeTable_RadauIIA(0) did not return NULL\n");
    FIRKodeTable_Free(T);
    fails++;
  }
  T = FIRKodeTable_RadauIIA(FIRK_MAX_STAGES + 1);
  if (T != NULL)
  {
    printf("  FAIL: FIRKodeTable_RadauIIA(%d) did not return NULL\n",
           FIRK_MAX_STAGES + 1);
    FIRKodeTable_Free(T);
    fails++;
  }

  /* the return value of FIRKodeTable_CheckOrder distinguishes a stored order
     that is too high (failure) from one that is too low (warning) */
  T = FIRKodeTable_RadauIIA(2);
  if (T != NULL)
  {
    int q, p, retval;
    T->q   = 5;
    retval = FIRKodeTable_CheckOrder(T, &q, &p, NULL);
    if (retval != -1 || q != 3)
    {
      printf("  FAIL: CheckOrder with stored order 5 returned %d, q = %d\n",
             retval, q);
      fails++;
    }
    T->q   = 1;
    retval = FIRKodeTable_CheckOrder(T, &q, &p, NULL);
    if (retval != 1 || q != 3)
    {
      printf("  FAIL: CheckOrder with stored order 1 returned %d, q = %d\n",
             retval, q);
      fails++;
    }
    FIRKodeTable_Free(T);
  }

  /* copying a table and writing it does not crash */
  T = FIRKodeTable_RadauIIA(3);
  if (T != NULL)
  {
    FIRKodeTable Tcopy = FIRKodeTable_Copy(T);
    if (Tcopy == NULL)
    {
      printf("  FAIL: FIRKodeTable_Copy returned NULL\n");
      fails++;
    }
    else
    {
      FIRKodeTable_Write(Tcopy, stdout);
      FIRKodeTable_Free(Tcopy);
    }
    FIRKodeTable_Free(T);
  }

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("\nAll tests passed\n");
  return 0;
}
