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
 * This is the implementation file for FIRKODE fully implicit
 * Runge-Kutta (Radau IIA collocation) tables.
 *
 * All internal computations are carried out in long double so that
 * the resulting tables are accurate to working precision for every
 * supported sunrealtype (float, double, long double).
 * -----------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <firkode/firkode_tables.h>
#include <sundials/sundials_math.h>

#include "firkode_tables_impl.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

typedef long double firk_ld;

#define LD_ZERO 0.0L
#define LD_ONE  1.0L
#define LD_TWO  2.0L

/* Number of bisection iterations used to locate the Radau nodes */
#define FIRK_NODE_BISECT_ITERS 120

/* Grid points per stage used to bracket the Radau nodes */
#define FIRK_NODE_GRID_PER_STAGE 40

/* Highest order tested by FIRKodeTable_CheckOrder */
#define FIRK_MAX_CHECK_ORDER (2 * FIRK_MAX_STAGES + 1)

/*---------------------------------------------------------------
  Real shift gamma0 for each stage count (derived with
  scripts/firkode_radau_tables.py). For odd s, gamma0 = 1/U1 with
  U1 the unique real eigenvalue of A^{-1}. For even s, A^{-1} has
  no real eigenvalue and gamma0 = 1/Re(lambda) with lambda the
  smallest-modulus eigenvalue of A^{-1}, which is the real shift
  minimizing the stiff-limit deviation of the block-preconditioned
  Newton matrix from the identity.
  ---------------------------------------------------------------*/
static const firk_ld firk_gamma0[FIRK_MAX_STAGES + 1] =
  {LD_ZERO,
   LD_ONE,
   0.5L,
   0.2748888295956773677478286L,
   0.208890675278273738551678L,
   0.1590658444274691204779152L,
   0.1334999853021741882617498L,
   0.1118964653000350759359053L,
   0.09833377347648450252097062L,
   0.08630100242868910058599241L};

/*===============================================================
  Long double dense linear algebra helpers (s <= FIRK_MAX_STAGES)
  ===============================================================*/

static firk_ld** ld_matrix_alloc(int s)
{
  int i;
  firk_ld** a = (firk_ld**)calloc((size_t)s, sizeof(firk_ld*));
  if (a == NULL) { return NULL; }
  for (i = 0; i < s; i++)
  {
    a[i] = (firk_ld*)calloc((size_t)s, sizeof(firk_ld));
    if (a[i] == NULL)
    {
      while (--i >= 0) { free(a[i]); }
      free(a);
      return NULL;
    }
  }
  return a;
}

static void ld_matrix_free(firk_ld** a, int s)
{
  int i;
  if (a == NULL) { return; }
  for (i = 0; i < s; i++) { free(a[i]); }
  free(a);
}

/* In-place LU factorization with partial pivoting; returns 0 on
   success and -1 if a zero pivot is encountered */
static int ld_lu_factor(firk_ld** a, int s, int* piv)
{
  int i, j, k, p;
  firk_ld amax, tmp;

  for (k = 0; k < s; k++)
  {
    p    = k;
    amax = fabsl(a[k][k]);
    for (i = k + 1; i < s; i++)
    {
      if (fabsl(a[i][k]) > amax)
      {
        amax = fabsl(a[i][k]);
        p    = i;
      }
    }
    if (amax == LD_ZERO) { return -1; }
    piv[k] = p;
    if (p != k)
    {
      for (j = 0; j < s; j++)
      {
        tmp     = a[k][j];
        a[k][j] = a[p][j];
        a[p][j] = tmp;
      }
    }
    for (i = k + 1; i < s; i++)
    {
      a[i][k] /= a[k][k];
      for (j = k + 1; j < s; j++) { a[i][j] -= a[i][k] * a[k][j]; }
    }
  }
  return 0;
}

/* Solve LU x = b in place using the factors from ld_lu_factor */
static void ld_lu_solve(firk_ld** a, int s, const int* piv, firk_ld* b)
{
  int i, j;
  firk_ld tmp;

  for (i = 0; i < s; i++)
  {
    if (piv[i] != i)
    {
      tmp       = b[i];
      b[i]      = b[piv[i]];
      b[piv[i]] = tmp;
    }
    for (j = 0; j < i; j++) { b[i] -= a[i][j] * b[j]; }
  }
  for (i = s - 1; i >= 0; i--)
  {
    for (j = i + 1; j < s; j++) { b[i] -= a[i][j] * b[j]; }
    b[i] /= a[i][i];
  }
}

/*===============================================================
  Radau IIA nodes and collocation coefficients
  ===============================================================*/

/* Evaluate P_s(x) - P_{s-1}(x) with the Legendre three-term
   recurrence; the interior Radau IIA nodes are its roots in (-1,1) */
static firk_ld firk_radau_poly(int s, firk_ld x)
{
  int k;
  firk_ld pkm1 = LD_ONE, pk = x, pkp1;

  if (s == 1) { return x - LD_ONE; }
  for (k = 1; k < s; k++)
  {
    pkp1 = ((LD_TWO * (firk_ld)k + LD_ONE) * x * pk - (firk_ld)k * pkm1) /
           ((firk_ld)k + LD_ONE);
    pkm1 = pk;
    pk   = pkp1;
  }
  return pk - pkm1;
}

/* Compute the Radau IIA nodes c_1 < ... < c_s = 1 in (0,1]. The
   s-1 interior roots are bracketed on a uniform grid in [-1,1) and
   refined by bisection. Returns 0 on success, -1 if the expected
   number of roots is not found. */
static int firk_radau_nodes(int s, firk_ld* c)
{
  int j, k, n, nroots;
  firk_ld a, b, m, fa, fb, fm;

  if (s == 1)
  {
    c[0] = LD_ONE;
    return 0;
  }

  n      = FIRK_NODE_GRID_PER_STAGE * s;
  nroots = 0;
  a      = -LD_ONE;
  fa     = firk_radau_poly(s, a);
  for (j = 1; j < n; j++)
  {
    b  = -LD_ONE + LD_TWO * (firk_ld)j / (firk_ld)n;
    fb = firk_radau_poly(s, b);
    if (fa == LD_ZERO)
    {
      if (nroots >= s - 1) { return -1; }
      c[nroots++] = (a + LD_ONE) / LD_TWO;
    }
    else if (fa * fb < LD_ZERO)
    {
      if (nroots >= s - 1) { return -1; }
      for (k = 0; k < FIRK_NODE_BISECT_ITERS; k++)
      {
        m  = (a + b) / LD_TWO;
        fm = firk_radau_poly(s, m);
        if (fm == LD_ZERO) { break; }
        if (fa * fm < LD_ZERO)
        {
          b  = m;
          fb = fm;
        }
        else
        {
          a  = m;
          fa = fm;
        }
      }
      m           = (a + b) / LD_TWO;
      c[nroots++] = (m + LD_ONE) / LD_TWO;
      b           = -LD_ONE + LD_TWO * (firk_ld)j / (firk_ld)n;
      fb          = firk_radau_poly(s, b);
    }
    a  = b;
    fa = fb;
  }
  if (nroots != s - 1) { return -1; }
  c[s - 1] = LD_ONE;
  return 0;
}

/* Monomial coefficients r[0..s-1] of the Lagrange basis polynomial
   l_j(tau) = prod_{k != j} (tau - c_k)/(c_j - c_k) on the nodes c */
static void firk_lagrange_coeffs(int s, const firk_ld* c, int j, firk_ld* r)
{
  int k, m, deg;
  firk_ld scale;

  for (m = 0; m < s; m++) { r[m] = LD_ZERO; }
  r[0] = LD_ONE;
  deg  = 0;
  for (k = 0; k < s; k++)
  {
    if (k == j) { continue; }
    scale = LD_ONE / (c[j] - c[k]);
    /* multiply by (tau - c_k) */
    for (m = deg; m >= 0; m--)
    {
      r[m + 1] += r[m];
      r[m] = -c[k] * r[m];
    }
    deg++;
    for (m = 0; m <= deg; m++) { r[m] *= scale; }
  }
}

/* Collocation coefficients a_ij = int_0^{c_i} l_j(tau) dtau */
static void firk_collocation_coeffs(int s, const firk_ld* c, firk_ld* r,
                                    firk_ld** A)
{
  int i, j, m;
  firk_ld sum, cpow;

  for (j = 0; j < s; j++)
  {
    firk_lagrange_coeffs(s, c, j, r);
    for (i = 0; i < s; i++)
    {
      sum  = LD_ZERO;
      cpow = c[i];
      for (m = 0; m < s; m++)
      {
        sum += r[m] * cpow / ((firk_ld)m + LD_ONE);
        cpow *= c[i];
      }
      A[i][j] = sum;
    }
  }
}

/*---------------------------------------------------------------
  Derived quantities from (c, A, b, gamma0), all in long double:
    Ainv = A^{-1}
    d    = b^T A^{-1}
    e    = -gamma0 * l(0)^T A^{-1}, l_i(0) = prod_{k!=i} (-c_k)/(c_i-c_k)
    P    = monomial coefficients of L_i(theta) = theta l_i(theta)/c_i
  Returns 0 on success, -1 on failure (singular A or allocation).
  ---------------------------------------------------------------*/
static int firk_derive_ld(int s, const firk_ld* c, firk_ld** A,
                          const firk_ld* b, firk_ld gamma0, firk_ld** Ainv,
                          firk_ld* d, firk_ld* e, firk_ld** P)
{
  int i, j, k, retval;
  int* piv      = NULL;
  firk_ld* work = NULL;
  firk_ld* l0   = NULL;
  firk_ld** LU  = NULL;
  retval        = -1;

  piv  = (int*)calloc((size_t)s, sizeof(int));
  work = (firk_ld*)calloc((size_t)s, sizeof(firk_ld));
  l0   = (firk_ld*)calloc((size_t)s, sizeof(firk_ld));
  LU   = ld_matrix_alloc(s);
  if (piv == NULL || work == NULL || l0 == NULL || LU == NULL) { goto cleanup; }

  /* Ainv via LU factorization of A */
  for (i = 0; i < s; i++)
  {
    for (j = 0; j < s; j++) { LU[i][j] = A[i][j]; }
  }
  if (ld_lu_factor(LU, s, piv) != 0) { goto cleanup; }
  for (j = 0; j < s; j++)
  {
    for (i = 0; i < s; i++) { work[i] = (i == j) ? LD_ONE : LD_ZERO; }
    ld_lu_solve(LU, s, piv, work);
    for (i = 0; i < s; i++) { Ainv[i][j] = work[i]; }
  }

  /* d = b^T A^{-1} */
  for (j = 0; j < s; j++)
  {
    d[j] = LD_ZERO;
    for (i = 0; i < s; i++) { d[j] += b[i] * Ainv[i][j]; }
  }

  /* Lagrange basis values at 0 and dense output coefficients */
  for (i = 0; i < s; i++)
  {
    firk_lagrange_coeffs(s, c, i, work);
    l0[i] = work[0];
    for (k = 0; k < s; k++) { P[k][i] = work[k] / c[i]; }
  }

  /* error estimate weights */
  for (j = 0; j < s; j++)
  {
    e[j] = LD_ZERO;
    for (i = 0; i < s; i++) { e[j] += l0[i] * Ainv[i][j]; }
    e[j] *= -gamma0;
  }

  retval = 0;

cleanup:
  free(piv);
  free(work);
  free(l0);
  ld_matrix_free(LU, s);
  return retval;
}

/* Store long double data into the sunrealtype table fields */
static void firk_store_table(FIRKodeTable T, const firk_ld* c, firk_ld** A,
                             const firk_ld* b, firk_ld gamma0, firk_ld** Ainv,
                             const firk_ld* d, const firk_ld* e, firk_ld** P)
{
  int i, j, s = T->s;

  for (i = 0; i < s; i++)
  {
    T->c[i] = (sunrealtype)c[i];
    T->b[i] = (sunrealtype)b[i];
    T->d[i] = (sunrealtype)d[i];
    T->e[i] = (sunrealtype)e[i];
    for (j = 0; j < s; j++)
    {
      T->A[i][j]    = (sunrealtype)A[i][j];
      T->Ainv[i][j] = (sunrealtype)Ainv[i][j];
      T->P[i][j]    = (sunrealtype)P[i][j];
    }
  }
  T->gamma0 = (sunrealtype)gamma0;
}

/* Fill the closed-form Radau IIA coefficients for s = 1, 2, 3 */
static void firk_radau_closed_form(int s, firk_ld* c, firk_ld** A)
{
  firk_ld r6;

  switch (s)
  {
  case 1:
    c[0]    = LD_ONE;
    A[0][0] = LD_ONE;
    break;
  case 2:
    c[0]    = LD_ONE / 3.0L;
    c[1]    = LD_ONE;
    A[0][0] = 5.0L / 12.0L;
    A[0][1] = -LD_ONE / 12.0L;
    A[1][0] = 3.0L / 4.0L;
    A[1][1] = LD_ONE / 4.0L;
    break;
  default: /* s == 3 */
    r6      = sqrtl(6.0L);
    c[0]    = (4.0L - r6) / 10.0L;
    c[1]    = (4.0L + r6) / 10.0L;
    c[2]    = LD_ONE;
    A[0][0] = (88.0L - 7.0L * r6) / 360.0L;
    A[0][1] = (296.0L - 169.0L * r6) / 1800.0L;
    A[0][2] = (-2.0L + 3.0L * r6) / 225.0L;
    A[1][0] = (296.0L + 169.0L * r6) / 1800.0L;
    A[1][1] = (88.0L + 7.0L * r6) / 360.0L;
    A[1][2] = (-2.0L - 3.0L * r6) / 225.0L;
    A[2][0] = (16.0L - r6) / 36.0L;
    A[2][1] = (16.0L + r6) / 36.0L;
    A[2][2] = LD_ONE / 9.0L;
    break;
  }
}

/* Build a Radau IIA table, using the closed-form coefficients for
   s <= 3 when closed_form is true and numerically computed nodes
   and collocation coefficients otherwise */
static FIRKodeTable firk_radau_build(int s, sunbooleantype closed_form)
{
  int i, retval;
  FIRKodeTable T = NULL;
  firk_ld *c = NULL, *b = NULL, *d = NULL, *e = NULL, *work = NULL;
  firk_ld **A = NULL, **Ainv = NULL, **P = NULL;

  if (s < 1 || s > FIRK_MAX_STAGES) { return NULL; }

  T = FIRKodeTable_Alloc(s);
  if (T == NULL) { return NULL; }

  c      = (firk_ld*)calloc((size_t)s, sizeof(firk_ld));
  b      = (firk_ld*)calloc((size_t)s, sizeof(firk_ld));
  d      = (firk_ld*)calloc((size_t)s, sizeof(firk_ld));
  e      = (firk_ld*)calloc((size_t)s, sizeof(firk_ld));
  work   = (firk_ld*)calloc((size_t)s, sizeof(firk_ld));
  A      = ld_matrix_alloc(s);
  Ainv   = ld_matrix_alloc(s);
  P      = ld_matrix_alloc(s);
  retval = -1;
  if (c == NULL || b == NULL || d == NULL || e == NULL || work == NULL ||
      A == NULL || Ainv == NULL || P == NULL)
  {
    goto cleanup;
  }

  if (closed_form && s <= 3) { firk_radau_closed_form(s, c, A); }
  else
  {
    if (firk_radau_nodes(s, c) != 0) { goto cleanup; }
    firk_collocation_coeffs(s, c, work, A);
  }
  for (i = 0; i < s; i++) { b[i] = A[s - 1][i]; }

  if (firk_derive_ld(s, c, A, b, firk_gamma0[s], Ainv, d, e, P) != 0)
  {
    goto cleanup;
  }

  firk_store_table(T, c, A, b, firk_gamma0[s], Ainv, d, e, P);
  T->q                = 2 * s - 1;
  T->p                = s;
  T->stiffly_accurate = SUNTRUE;
  retval              = 0;

cleanup:
  free(c);
  free(b);
  free(d);
  free(e);
  free(work);
  ld_matrix_free(A, s);
  ld_matrix_free(Ainv, s);
  ld_matrix_free(P, s);
  if (retval != 0)
  {
    FIRKodeTable_Free(T);
    return NULL;
  }
  return T;
}

/*===============================================================
  Exported routines
  ===============================================================*/

/*---------------------------------------------------------------
  Routine to allocate an empty table structure with s stages
  ---------------------------------------------------------------*/
FIRKodeTable FIRKodeTable_Alloc(int s)
{
  int i;
  FIRKodeTable T;

  if (s < 1) { return NULL; }

  T = (FIRKodeTable)malloc(sizeof(struct FIRKodeTableMem));
  if (T == NULL) { return NULL; }

  T->s                = s;
  T->q                = 0;
  T->p                = 0;
  T->c                = NULL;
  T->A                = NULL;
  T->b                = NULL;
  T->Ainv             = NULL;
  T->d                = NULL;
  T->e                = NULL;
  T->gamma0           = ZERO;
  T->P                = NULL;
  T->stiffly_accurate = SUNFALSE;

  T->c    = (sunrealtype*)calloc((size_t)s, sizeof(sunrealtype));
  T->b    = (sunrealtype*)calloc((size_t)s, sizeof(sunrealtype));
  T->d    = (sunrealtype*)calloc((size_t)s, sizeof(sunrealtype));
  T->e    = (sunrealtype*)calloc((size_t)s, sizeof(sunrealtype));
  T->A    = (sunrealtype**)calloc((size_t)s, sizeof(sunrealtype*));
  T->Ainv = (sunrealtype**)calloc((size_t)s, sizeof(sunrealtype*));
  T->P    = (sunrealtype**)calloc((size_t)s, sizeof(sunrealtype*));
  if (T->c == NULL || T->b == NULL || T->d == NULL || T->e == NULL ||
      T->A == NULL || T->Ainv == NULL || T->P == NULL)
  {
    FIRKodeTable_Free(T);
    return NULL;
  }

  for (i = 0; i < s; i++)
  {
    T->A[i]    = (sunrealtype*)calloc((size_t)s, sizeof(sunrealtype));
    T->Ainv[i] = (sunrealtype*)calloc((size_t)s, sizeof(sunrealtype));
    T->P[i]    = (sunrealtype*)calloc((size_t)s, sizeof(sunrealtype));
    if (T->A[i] == NULL || T->Ainv[i] == NULL || T->P[i] == NULL)
    {
      FIRKodeTable_Free(T);
      return NULL;
    }
  }

  return T;
}

/*---------------------------------------------------------------
  Routine to construct the Radau IIA table with s stages
  ---------------------------------------------------------------*/
FIRKodeTable FIRKodeTable_RadauIIA(int s)
{
  return firk_radau_build(s, SUNTRUE);
}

/*---------------------------------------------------------------
  Private routine to construct the Radau IIA table numerically
  ---------------------------------------------------------------*/
FIRKodeTable firkTable_RadauIIA_Computed(int s)
{
  return firk_radau_build(s, SUNFALSE);
}

/*---------------------------------------------------------------
  Private routine to (re)compute the derived quantities of a table
  from its c, A, b, and gamma0 entries
  ---------------------------------------------------------------*/
int firkTable_Derive(FIRKodeTable T)
{
  int i, j, s, retval;
  firk_ld *c = NULL, *b = NULL, *d = NULL, *e = NULL;
  firk_ld **A = NULL, **Ainv = NULL, **P = NULL;

  if (T == NULL || T->s < 1) { return -1; }
  s = T->s;

  c      = (firk_ld*)calloc((size_t)s, sizeof(firk_ld));
  b      = (firk_ld*)calloc((size_t)s, sizeof(firk_ld));
  d      = (firk_ld*)calloc((size_t)s, sizeof(firk_ld));
  e      = (firk_ld*)calloc((size_t)s, sizeof(firk_ld));
  A      = ld_matrix_alloc(s);
  Ainv   = ld_matrix_alloc(s);
  P      = ld_matrix_alloc(s);
  retval = -1;
  if (c == NULL || b == NULL || d == NULL || e == NULL || A == NULL ||
      Ainv == NULL || P == NULL)
  {
    goto cleanup;
  }

  for (i = 0; i < s; i++)
  {
    c[i] = (firk_ld)T->c[i];
    b[i] = (firk_ld)T->b[i];
    for (j = 0; j < s; j++) { A[i][j] = (firk_ld)T->A[i][j]; }
  }

  if (firk_derive_ld(s, c, A, b, (firk_ld)T->gamma0, Ainv, d, e, P) != 0)
  {
    goto cleanup;
  }
  firk_store_table(T, c, A, b, (firk_ld)T->gamma0, Ainv, d, e, P);
  retval = 0;

cleanup:
  free(c);
  free(b);
  free(d);
  free(e);
  ld_matrix_free(A, s);
  ld_matrix_free(Ainv, s);
  ld_matrix_free(P, s);
  return retval;
}

/*---------------------------------------------------------------
  Routine to copy a table structure
  ---------------------------------------------------------------*/
FIRKodeTable FIRKodeTable_Copy(FIRKodeTable T)
{
  int i, j, s;
  FIRKodeTable Tcopy;

  if (T == NULL) { return NULL; }
  s     = T->s;
  Tcopy = FIRKodeTable_Alloc(s);
  if (Tcopy == NULL) { return NULL; }

  Tcopy->q                = T->q;
  Tcopy->p                = T->p;
  Tcopy->gamma0           = T->gamma0;
  Tcopy->stiffly_accurate = T->stiffly_accurate;
  for (i = 0; i < s; i++)
  {
    Tcopy->c[i] = T->c[i];
    Tcopy->b[i] = T->b[i];
    Tcopy->d[i] = T->d[i];
    Tcopy->e[i] = T->e[i];
    for (j = 0; j < s; j++)
    {
      Tcopy->A[i][j]    = T->A[i][j];
      Tcopy->Ainv[i][j] = T->Ainv[i][j];
      Tcopy->P[i][j]    = T->P[i][j];
    }
  }
  return Tcopy;
}

/*---------------------------------------------------------------
  Routine to free a table structure
  ---------------------------------------------------------------*/
void FIRKodeTable_Free(FIRKodeTable T)
{
  int i;

  if (T == NULL) { return; }
  if (T->A != NULL)
  {
    for (i = 0; i < T->s; i++) { free(T->A[i]); }
    free(T->A);
  }
  if (T->Ainv != NULL)
  {
    for (i = 0; i < T->s; i++) { free(T->Ainv[i]); }
    free(T->Ainv);
  }
  if (T->P != NULL)
  {
    for (i = 0; i < T->s; i++) { free(T->P[i]); }
    free(T->P);
  }
  free(T->c);
  free(T->b);
  free(T->d);
  free(T->e);
  free(T);
}

/*---------------------------------------------------------------
  Routine to print a table structure
  ---------------------------------------------------------------*/
void FIRKodeTable_Write(FIRKodeTable T, FILE* outfile)
{
  int i, j;

  if (T == NULL || outfile == NULL) { return; }

  fprintf(outfile, "  stages = %i,  q = %i,  p = %i\n", T->s, T->q, T->p);
  fprintf(outfile, "  A = \n");
  for (i = 0; i < T->s; i++)
  {
    fprintf(outfile, "      ");
    for (j = 0; j < T->s; j++)
    {
      fprintf(outfile, SUN_FORMAT_E "  ", T->A[i][j]);
    }
    fprintf(outfile, "\n");
  }
  fprintf(outfile, "  c = ");
  for (i = 0; i < T->s; i++) { fprintf(outfile, SUN_FORMAT_E "  ", T->c[i]); }
  fprintf(outfile, "\n");
  fprintf(outfile, "  b = ");
  for (i = 0; i < T->s; i++) { fprintf(outfile, SUN_FORMAT_E "  ", T->b[i]); }
  fprintf(outfile, "\n");
  fprintf(outfile, "  e = ");
  for (i = 0; i < T->s; i++) { fprintf(outfile, SUN_FORMAT_E "  ", T->e[i]); }
  fprintf(outfile, "\n");
  fprintf(outfile, "  gamma0 = " SUN_FORMAT_E "\n", T->gamma0);
}

/*---------------------------------------------------------------
  Routine to determine the order of a table from the simplifying
  assumptions B(q), C(eta), and D(zeta): the method has order at
  least min(qB, 2 eta + 2, eta + zeta + 1). Returns 0 if the
  computed order matches T->q, 1 if it is higher (a warning), and
  -1 if it is lower or the table is invalid. The error estimate
  order p is reported as stored in the table.
  ---------------------------------------------------------------*/
int FIRKodeTable_CheckOrder(FIRKodeTable T, int* q, int* p, FILE* outfile)
{
  int i, j, k, s, qB, eta, zeta, qcomp;
  sunrealtype tol, sum, cpow, resid;
  sunbooleantype ok;

  if (T == NULL || q == NULL || p == NULL) { return -1; }
  if (T->s < 1 || T->c == NULL || T->A == NULL || T->b == NULL) { return -1; }
  s = T->s;

  /* The smallest residual of a violated condition among the supported
     Radau IIA tables is about 1e-10 (B(18) for s = 9), so this tolerance
     separates satisfied from violated conditions in double precision for
     all s <= FIRK_MAX_STAGES, and in single precision for s <= 3. */
  tol = SUN_RCONST(1.0e3) * SUN_UNIT_ROUNDOFF * (sunrealtype)s;

  /* B(k): sum_i b_i c_i^{k-1} = 1/k */
  qB = 0;
  for (k = 1; k <= FIRK_MAX_CHECK_ORDER; k++)
  {
    sum = ZERO;
    for (i = 0; i < s; i++) { sum += T->b[i] * SUNRpowerI(T->c[i], k - 1); }
    resid = SUNRabs(sum - ONE / (sunrealtype)k);
    if (resid > tol) { break; }
    qB = k;
  }

  /* C(k): sum_j a_ij c_j^{k-1} = c_i^k / k for all i */
  eta = 0;
  for (k = 1; k <= FIRK_MAX_CHECK_ORDER; k++)
  {
    ok = SUNTRUE;
    for (i = 0; i < s; i++)
    {
      sum = ZERO;
      for (j = 0; j < s; j++)
      {
        sum += T->A[i][j] * SUNRpowerI(T->c[j], k - 1);
      }
      cpow  = SUNRpowerI(T->c[i], k) / (sunrealtype)k;
      resid = SUNRabs(sum - cpow);
      if (resid > tol) { ok = SUNFALSE; }
    }
    if (!ok) { break; }
    eta = k;
  }

  /* D(k): sum_i b_i c_i^{k-1} a_ij = b_j (1 - c_j^k) / k for all j */
  zeta = 0;
  for (k = 1; k <= FIRK_MAX_CHECK_ORDER; k++)
  {
    ok = SUNTRUE;
    for (j = 0; j < s; j++)
    {
      sum = ZERO;
      for (i = 0; i < s; i++)
      {
        sum += T->b[i] * SUNRpowerI(T->c[i], k - 1) * T->A[i][j];
      }
      cpow  = T->b[j] * (ONE - SUNRpowerI(T->c[j], k)) / (sunrealtype)k;
      resid = SUNRabs(sum - cpow);
      if (resid > tol) { ok = SUNFALSE; }
    }
    if (!ok) { break; }
    zeta = k;
  }

  qcomp = qB;
  if (2 * eta + 2 < qcomp) { qcomp = 2 * eta + 2; }
  if (eta + zeta + 1 < qcomp) { qcomp = eta + zeta + 1; }

  *q = qcomp;
  *p = T->p;

  if (outfile != NULL)
  {
    fprintf(outfile, "  FIRKodeTable_CheckOrder: B(%i), C(%i), D(%i)", qB, eta,
            zeta);
    fprintf(outfile, " -> q = %i (table q = %i), p = %i\n", qcomp, T->q, T->p);
  }

  /* an order below the stored one is a failure, above it a warning */
  if (qcomp < T->q) { return -1; }
  return (qcomp > T->q) ? 1 : 0;
}
