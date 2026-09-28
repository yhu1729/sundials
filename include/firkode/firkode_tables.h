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
 * This is the header file for FIRKODE fully implicit Runge-Kutta
 * (Radau IIA collocation) table structures.
 * -----------------------------------------------------------------*/

#ifndef _FIRKODE_TABLES_H
#define _FIRKODE_TABLES_H

#include <stdio.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Maximum number of stages for which tables can be constructed */
#define FIRK_MAX_STAGES 9

/*---------------------------------------------------------------
  Types : struct FIRKodeTableMem, FIRKodeTable

  A fully implicit Runge-Kutta table with s stages together with
  the derived quantities used by the FIRKODE integrator:

    A, b, c  : Butcher table coefficients (A is a full s x s array)
    Ainv     : A^{-1}
    d        : d = b^T A^{-1}  (d = e_s for stiffly accurate tables)
    e        : error estimate weights; with the stage increments
               z_i = Y_i - y_n, the embedded solution satisfies
               yhat - y = gamma0 h f(t_n,y_n) + sum_i e_i M z_i
    gamma0   : real shift used for the block preconditioner and
               the filtered error estimate (M - h gamma0 J)^{-1}
    P        : dense output coefficients, P[k][i] is the
               coefficient of theta^(k+1) in the Lagrange basis
               polynomial L_i(theta) on the nodes {0, c_1, ..., c_s},
               so that u(t_n + theta h) = y_n + sum_i L_i(theta) z_i
  ---------------------------------------------------------------*/

struct FIRKodeTableMem
{
  int s;                           /* number of stages                */
  int q;                           /* method order of accuracy        */
  int p;                           /* error estimate order            */
  sunrealtype* c;                  /* stage nodes                     */
  sunrealtype** A;                 /* Butcher table coefficients      */
  sunrealtype* b;                  /* solution weights                */
  sunrealtype** Ainv;              /* inverse of A                    */
  sunrealtype* d;                  /* b^T A^{-1}                      */
  sunrealtype* e;                  /* error estimate weights          */
  sunrealtype gamma0;              /* real shift for preconditioning  */
  sunrealtype** P;                 /* dense output coefficients       */
  sunbooleantype stiffly_accurate; /* b == last row of A           */
};

typedef _SUNDIALS_STRUCT_ FIRKodeTableMem* FIRKodeTable;

/* Utility routines to construct/copy/free/output/check tables */

SUNDIALS_EXPORT FIRKodeTable FIRKodeTable_RadauIIA(int s);
SUNDIALS_EXPORT FIRKodeTable FIRKodeTable_Alloc(int s);
SUNDIALS_EXPORT FIRKodeTable FIRKodeTable_Copy(FIRKodeTable T);
SUNDIALS_EXPORT void FIRKodeTable_Free(FIRKodeTable T);
SUNDIALS_EXPORT void FIRKodeTable_Write(FIRKodeTable T, FILE* outfile);
SUNDIALS_EXPORT int FIRKodeTable_CheckOrder(FIRKodeTable T, int* q, int* p,
                                            FILE* outfile);

#ifdef __cplusplus
}
#endif

#endif
