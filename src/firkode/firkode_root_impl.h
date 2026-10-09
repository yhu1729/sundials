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
 * Implementation header file for FIRKODE's rootfinding module
 * (adapted from ARKODE's).
 * -----------------------------------------------------------------*/

#ifndef _FIRKODE_ROOT_IMPL_H
#define _FIRKODE_ROOT_IMPL_H

#include <stdio.h>

#include <firkode/firkode.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*---------------------------------------------------------------
  Types : struct FIRKodeRootMemRec, FIRKodeRootMem
  ---------------------------------------------------------------*/
typedef struct FIRKodeRootMemRec
{
  FIRKRootFn gfun;         /* function g for roots sought                  */
  int nrtfn;               /* number of components of g                    */
  int* iroots;             /* array for root information                   */
  int* rootdir;            /* array specifying direction of zero-crossing  */
  sunrealtype tlo;         /* nearest endpoint of interval in root search  */
  sunrealtype thi;         /* farthest endpoint of interval in root search */
  sunrealtype trout;       /* t value returned by rootfinding routine      */
  sunrealtype* glo;        /* saved array of g values at t = tlo           */
  sunrealtype* ghi;        /* saved array of g values at t = thi           */
  sunrealtype* grout;      /* array of g values at t = trout               */
  sunrealtype ttol;        /* tolerance on root location                   */
  int irfnd;               /* flag showing whether last step had a root    */
  long int nge;            /* counter for g evaluations                    */
  sunbooleantype* gactive; /* array with active/inactive event functions   */
  int mxgnull;             /* num. warning messages about possible g==0    */
  void* root_data;         /* pointer to user_data                         */

}* FIRKodeRootMem;

/*---------------------------------------------------------------
  Rootfinding routines
  ---------------------------------------------------------------*/

int firkRootFree(void* firkode_mem);
int firkPrintRootMem(void* firkode_mem, FILE* outfile);
int firkRootCheck1(void* firkode_mem);
int firkRootCheck2(void* firkode_mem);
int firkRootCheck3(void* firkode_mem, sunrealtype tout, int itask);
int firkRootfind(void* firkode_mem);

#ifdef __cplusplus
}
#endif

#endif
