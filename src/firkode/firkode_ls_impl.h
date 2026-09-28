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
 * Implementation header file for FIRKODE's linear solver interface.
 * -----------------------------------------------------------------*/

#ifndef _FIRKLS_IMPL_H
#define _FIRKLS_IMPL_H

#include <firkode/firkode_ls.h>

#include "firkode_impl.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*-----------------------------------------------------------------
  FIRKLS solver constants

  FIRKLS_MSBJ   maximum number of steps between Jacobian and/or
                preconditioner evaluations
  FIRKLS_DGMAX  maximum change in gamma between Jacobian and/or
                preconditioner evaluations
  FIRKLS_EPLIN  default value for factor by which the tolerance on
                the nonlinear iteration is multiplied to get a
                tolerance on the linear iteration
  FIRKLS_MAXRS  default number of restarts of the stage-system
                flexible GMRES iteration
  -----------------------------------------------------------------*/
#define FIRKLS_MSBJ  51
#define FIRKLS_DGMAX SUN_RCONST(0.2)
#define FIRKLS_EPLIN SUN_RCONST(0.05)
#define FIRKLS_MAXRS 1

/*-----------------------------------------------------------------
  Types : FIRKLsMemRec, FIRKLsMem
  -----------------------------------------------------------------*/
typedef struct FIRKLsMemRec
{
  /* Linear solver type information */
  sunbooleantype iterative;   /* is the block solver iterative?      */
  sunbooleantype matrixbased; /* is a matrix structure used?         */

  /* Jacobian construction & storage */
  sunbooleantype jacDQ;   /* SUNTRUE if using internal DQ Jac approx. */
  FIRKLsJacFn jac;        /* Jacobian routine to be called            */
  void* J_data;           /* user data is passed to jac               */
  sunbooleantype jbad;    /* heuristic suggestion for pset            */
  sunrealtype dgmax_jbad; /* if convfail = FAIL_BAD_J and the gamma
                             ratio |gamma/gammap-1| < dgmax_jbad then
                             J is bad                                 */

  /* Iterative block solver tolerance */
  sunrealtype eplifac; /* nonlinear -> linear tol scaling factor   */
  sunrealtype nrmfac;  /* integrator -> LS norm conversion factor  */

  /* Block linear solver, matrix and vector objects/pointers */
  SUNLinearSolver LS; /* user's N x N linear solver                */
  SUNMatrix A;        /* A = M - gamma * df/dy                     */
  SUNMatrix savedJ;   /* savedJ = old Jacobian                     */
  N_Vector ytemp;     /* temp vector passed to jtimes and psolve   */
  N_Vector x;         /* temp vector used by the block solve       */
  N_Vector ycur;      /* linearization point y                     */
  N_Vector fcur;      /* fcur = f(tn, ycur)                        */

  /* Stacked stage-system solver (flexible GMRES on ManyVectors) */
  SUNLinearSolver LS_stk;      /* internal SPFGMR solver, NULL if s == 1 */
  N_Vector x_stk;              /* stacked solution vector                */
  N_Vector* Jv;                /* J v_j products, one per stage          */
  int s_alloc;                 /* number of stages the above hold        */
  int maxl_stk;                /* Krylov subspace dimension              */
  int maxrs_stk;               /* maximum number of restarts             */
  sunbooleantype maxl_stk_set; /* user set maxl_stk                 */
  sunrealtype eplifac_stk;     /* nonlinear -> stacked linear tol factor */
  sunrealtype nrmfac_stk;      /* stacked norm conversion factor         */

  /* Statistics and associated parameters */
  long int msbj;     /* max num steps between jac/pset calls         */
  long int nje;      /* nje = no. of calls to jac                    */
  long int nfeDQ;    /* no. of calls to f due to DQ Jacobian or J*v  */
  long int nstlj;    /* nstlj = nst at last jac/pset call            */
  long int npe;      /* npe = total number of pset calls             */
  long int nli;      /* nli = total number of block linear iters     */
  long int nps;      /* nps = total number of psolve calls           */
  long int ncfl;     /* ncfl = total number of block conv failures   */
  long int njtsetup; /* njtsetup = total number of calls to jtsetup  */
  long int njtimes;  /* njtimes = total number of calls to jtimes    */
  long int nbs;      /* nbs = total number of block solves           */
  long int nsli;     /* nsli = total number of stacked linear iters  */
  long int nscf;     /* nscf = total number of stacked conv failures */
  sunrealtype tnlj;  /* tnlj = t_n at last jac/pset call             */

  /* Preconditioner computation
   * (a) user-provided:
   *     - P_data == user_data
   *     - pfree == NULL (the user deallocates memory for user_data)
   * (b) internal preconditioner module
   *     - P_data == firkode_mem
   *     - pfree == set by the prec. module and called in FIRKodeFree */
  FIRKLsPrecSetupFn pset;
  FIRKLsPrecSolveFn psolve;
  int (*pfree)(FIRKodeMem firk_mem);
  void* P_data;

  /* Jacobian times vector computation
   * (a) jtimes function provided by the user:
   *     - jt_data == user_data
   *     - jtimesDQ == SUNFALSE
   * (b) internal jtimes
   *     - jt_data == firkode_mem
   *     - jtimesDQ == SUNTRUE */
  sunbooleantype jtimesDQ;
  FIRKLsJacTimesSetupFn jtsetup;
  FIRKLsJacTimesVecFn jtimes;
  FIRKRhsFn jt_f;
  void* jt_data;

  /* Linear system setup function
   * (a) user-provided linsys function:
   *     - user_linsys = SUNTRUE
   *     - A_data      = user_data
   * (b) internal linsys function:
   *     - user_linsys = SUNFALSE
   *     - A_data      = firkode_mem */
  sunbooleantype user_linsys;
  FIRKLsLinSysFn linsys;
  void* A_data;

  int last_flag; /* last error flag returned by any function */

}* FIRKLsMem;

/*-----------------------------------------------------------------
  Types : FIRKLsMassMemRec, FIRKLsMassMem
  -----------------------------------------------------------------*/
typedef struct FIRKLsMassMemRec
{
  /* Linear solver type information */
  sunbooleantype iterative;   /* is the solver iterative?    */
  sunbooleantype matrixbased; /* is a matrix structure used? */

  /* Mass matrix construction & storage */
  FIRKLsMassFn mass; /* user-provided mass matrix routine        */
  SUNMatrix M;       /* mass matrix structure                    */
  SUNMatrix M_lu;    /* factored mass matrix structure           */
  void* M_data;      /* user data pointer passed to mass         */

  /* Iterative solver tolerance */
  sunrealtype eplifac; /* nonlinear -> linear tol scaling factor */
  sunrealtype nrmfac;  /* integrator -> LS norm conversion factor */

  /* Statistics and associated parameters */
  long int nmsetups;  /* total number of mass matrix setup calls  */
  long int nmvsetups; /* total number of mass-vector setup calls  */
  long int nmvevals;  /* total number of mass-vector products     */
  long int nmsolves;  /* total number of mass solves              */
  long int nmpe;      /* total number of mass preconditioner setups */
  long int nmli;      /* total number of mass linear iterations   */
  long int nmps;      /* total number of mass psolve calls        */
  long int nmcfails;  /* total number of mass solver conv fails   */

  /* Linear solver, matrix and vector objects/pointers */
  SUNLinearSolver LS; /* generic linear solver object            */
  N_Vector x;         /* temp vector used by the mass solve      */
  N_Vector ycur;      /* current solution vector                 */

  /* Preconditioner computation */
  FIRKLsMassPrecSetupFn pset;
  FIRKLsMassPrecSolveFn psolve;
  int (*pfree)(FIRKodeMem firk_mem);
  void* P_data;

  /* Mass matrix times vector computation */
  FIRKLsMassTimesSetupFn mtsetup;
  FIRKLsMassTimesVecFn mtimes;
  void* mt_data;

  int last_flag; /* last error flag returned by any function */

}* FIRKLsMassMem;

/*-----------------------------------------------------------------
  Prototypes of internal functions
  -----------------------------------------------------------------*/

/* Interface routines called by the user's block SUNLinearSolver */
int firkLsATimes(void* firkode_mem, N_Vector v, N_Vector z);
int firkLsPSetup(void* firkode_mem);
int firkLsPSolve(void* firkode_mem, N_Vector r, N_Vector z, sunrealtype tol,
                 int lr);

/* Interface routines called by the internal stacked solver */
int firkLsStackedATimes(void* firkode_mem, N_Vector v_stk, N_Vector z_stk);
int firkLsStackedPSolve(void* firkode_mem, N_Vector r_stk, N_Vector z_stk,
                        sunrealtype tol, int lr);

/* Difference quotient approximation for Jacobian and J*v */
int firkLsDQJac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac,
                void* firkode_mem, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int firkLsDenseDQJac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac,
                     FIRKodeMem firk_mem, N_Vector tmp1);
int firkLsBandDQJac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac,
                    FIRKodeMem firk_mem, N_Vector tmp1, N_Vector tmp2);
int firkLsDQJtimes(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y,
                   N_Vector fy, void* firkode_mem, N_Vector work);

/* Generic linear solver interface functions */
int firkLsInitialize(FIRKodeMem firk_mem);
int firkLsReInitialize(FIRKodeMem firk_mem);
int firkLsSetup(FIRKodeMem firk_mem, int convfail, N_Vector ypred,
                N_Vector fpred, sunbooleantype* jcurPtr, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3);
int firkLsSolveBlock(FIRKodeMem firk_mem, N_Vector b, N_Vector weight);
int firkLsSolveStacked(FIRKodeMem firk_mem, N_Vector b_stk);
int firkLsFree(FIRKodeMem firk_mem);

/* Auxiliary functions */
int firkLsInitializeCounters(FIRKLsMem firkls_mem);
int firkLs_AccessLMem(void* firkode_mem, const char* fname,
                      FIRKodeMem* firk_mem, FIRKLsMem* firkls_mem);

/* Mass matrix interface functions */
int firkLsMassInitialize(FIRKodeMem firk_mem);
int firkLsMassSetup(FIRKodeMem firk_mem, sunrealtype t, N_Vector vtemp1,
                    N_Vector vtemp2, N_Vector vtemp3);
int firkLsMassMult(FIRKodeMem firk_mem, N_Vector v, N_Vector Mv);
int firkLsMassSolve(FIRKodeMem firk_mem, N_Vector b, sunrealtype tol);
int firkLsMassFree(FIRKodeMem firk_mem);
int firkLsMTimes(void* firkode_mem, N_Vector v, N_Vector z);
int firkLsMPSetup(void* firkode_mem);
int firkLsMPSolve(void* firkode_mem, N_Vector r, N_Vector z, sunrealtype tol,
                  int lr);
int firkLs_AccessMassMem(void* firkode_mem, const char* fname,
                         FIRKodeMem* firk_mem, FIRKLsMassMem* firkls_mem);

/*-----------------------------------------------------------------
  Error Messages
  -----------------------------------------------------------------*/

#define MSG_LS_FIRKMEM_NULL "Integrator memory is NULL."
#define MSG_LS_MEM_FAIL     "A memory request failed."
#define MSG_LS_BAD_NVECTOR  "A required vector operation is not implemented."
#define MSG_LS_LMEM_NULL    "Linear solver memory is NULL."
#define MSG_LS_MASSMEM_NULL "Mass matrix solver memory is NULL."
#define MSG_LS_BAD_SIZES \
  "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSG_LS_JACFUNC_FAILED \
  "The Jacobian routine failed in an unrecoverable manner."
#define MSG_LS_MASSFUNC_FAILED \
  "The mass matrix routine failed in an unrecoverable manner."
#define MSG_LS_SUNMAT_FAILED \
  "A SUNMatrix routine failed in an unrecoverable manner."
#define MSG_LS_JTSETUP_FAILED \
  "The Jacobian x vector setup routine failed in an unrecoverable manner."
#define MSG_LS_JTIMES_FAILED \
  "The Jacobian x vector routine failed in an unrecoverable manner."
#define MSG_LS_MTSETUP_FAILED \
  "The mass matrix x vector setup routine failed in an unrecoverable manner."
#define MSG_LS_MTIMES_FAILED \
  "The mass matrix x vector routine failed in an unrecoverable manner."
#define MSG_LS_PSOLVE_FAILED \
  "The preconditioner solve routine failed in an unrecoverable manner."
#define MSG_LS_MATRIX_TYPES \
  "The Jacobian and mass matrix objects must have the same SUNMatrix type."

#ifdef __cplusplus
}
#endif

#endif
