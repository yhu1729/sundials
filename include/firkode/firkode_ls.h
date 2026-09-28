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
 * This is the header file for FIRKODE's linear solver interface
 * (FIRKLS). The user attaches an N x N SUNLinearSolver for the
 * shifted systems (M - gamma J); FIRKODE applies it block-wise as
 * the preconditioner of an internal flexible GMRES iteration on the
 * coupled stage system and reuses it for the error estimate.
 * -----------------------------------------------------------------*/

#ifndef _FIRKLS_H
#define _FIRKLS_H

#include <sundials/sundials_direct.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*=================================================================
  FIRKLS Constants
  =================================================================*/

#define FIRKLS_SUCCESS          0
#define FIRKLS_MEM_NULL         -1
#define FIRKLS_LMEM_NULL        -2
#define FIRKLS_ILL_INPUT        -3
#define FIRKLS_MEM_FAIL         -4
#define FIRKLS_PMEM_NULL        -5
#define FIRKLS_MASSMEM_NULL     -6
#define FIRKLS_JACFUNC_UNRECVR  -7
#define FIRKLS_JACFUNC_RECVR    -8
#define FIRKLS_MASSFUNC_UNRECVR -9
#define FIRKLS_MASSFUNC_RECVR   -10
#define FIRKLS_SUNMAT_FAIL      -11
#define FIRKLS_SUNLS_FAIL       -12

/*=================================================================
  FIRKLS user-supplied function prototypes

  All Jacobian-related functions are called at the single
  linearization point (t, y, fy) = (tn, yn, f(tn,yn)) of the
  current step with gamma = h * gamma0.
  =================================================================*/

typedef int (*FIRKLsJacFn)(sunrealtype t, N_Vector y, N_Vector fy,
                           SUNMatrix Jac, void* user_data, N_Vector tmp1,
                           N_Vector tmp2, N_Vector tmp3);

typedef int (*FIRKLsPrecSetupFn)(sunrealtype t, N_Vector y, N_Vector fy,
                                 sunbooleantype jok, sunbooleantype* jcurPtr,
                                 sunrealtype gamma, void* user_data);

typedef int (*FIRKLsPrecSolveFn)(sunrealtype t, N_Vector y, N_Vector fy,
                                 N_Vector r, N_Vector z, sunrealtype gamma,
                                 sunrealtype delta, int lr, void* user_data);

typedef int (*FIRKLsJacTimesSetupFn)(sunrealtype t, N_Vector y, N_Vector fy,
                                     void* user_data);

typedef int (*FIRKLsJacTimesVecFn)(N_Vector v, N_Vector Jv, sunrealtype t,
                                   N_Vector y, N_Vector fy, void* user_data,
                                   N_Vector tmp);

typedef int (*FIRKLsLinSysFn)(sunrealtype t, N_Vector y, N_Vector fy,
                              SUNMatrix A, SUNMatrix M, sunbooleantype jok,
                              sunbooleantype* jcur, sunrealtype gamma,
                              void* user_data, N_Vector tmp1, N_Vector tmp2,
                              N_Vector tmp3);

typedef int (*FIRKLsMassFn)(sunrealtype t, SUNMatrix M, void* user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

typedef int (*FIRKLsMassTimesSetupFn)(sunrealtype t, void* mtimes_data);

typedef int (*FIRKLsMassTimesVecFn)(N_Vector v, N_Vector Mv, sunrealtype t,
                                    void* mtimes_data);

typedef int (*FIRKLsMassPrecSetupFn)(sunrealtype t, void* user_data);

typedef int (*FIRKLsMassPrecSolveFn)(sunrealtype t, N_Vector r, N_Vector z,
                                     sunrealtype delta, int lr, void* user_data);

/*=================================================================
  FIRKLS Exported functions
  =================================================================*/

SUNDIALS_EXPORT int FIRKodeSetLinearSolver(void* firkode_mem,
                                           SUNLinearSolver LS, SUNMatrix A);

SUNDIALS_EXPORT int FIRKodeSetMassLinearSolver(void* firkode_mem,
                                               SUNLinearSolver LS, SUNMatrix M,
                                               sunbooleantype time_dep);

/*-----------------------------------------------------------------
  Optional inputs to the FIRKLS linear solver interface
  -----------------------------------------------------------------*/

SUNDIALS_EXPORT int FIRKodeSetJacFn(void* firkode_mem, FIRKLsJacFn jac);
SUNDIALS_EXPORT int FIRKodeSetJacEvalFrequency(void* firkode_mem, long int msbj);
SUNDIALS_EXPORT int FIRKodeSetDeltaGammaMaxBadJac(void* firkode_mem,
                                                  sunrealtype dgmax_jbad);
SUNDIALS_EXPORT int FIRKodeSetEpsLin(void* firkode_mem, sunrealtype eplifac);
SUNDIALS_EXPORT int FIRKodeSetLSNormFactor(void* firkode_mem, sunrealtype nrmfac);
SUNDIALS_EXPORT int FIRKodeSetPreconditioner(void* firkode_mem,
                                             FIRKLsPrecSetupFn pset,
                                             FIRKLsPrecSolveFn psolve);
SUNDIALS_EXPORT int FIRKodeSetJacTimes(void* firkode_mem,
                                       FIRKLsJacTimesSetupFn jtsetup,
                                       FIRKLsJacTimesVecFn jtimes);
SUNDIALS_EXPORT int FIRKodeSetLinSysFn(void* firkode_mem, FIRKLsLinSysFn linsys);

/* Options for the internal flexible GMRES iteration on the coupled
   stage system (only used when the method has more than one stage) */

SUNDIALS_EXPORT int FIRKodeSetStageSolverMaxl(void* firkode_mem, int maxl);
SUNDIALS_EXPORT int FIRKodeSetStageSolverMaxRestarts(void* firkode_mem,
                                                     int maxrs);
SUNDIALS_EXPORT int FIRKodeSetStageEpsLin(void* firkode_mem, sunrealtype eplifac);

/* Mass matrix optional inputs */

SUNDIALS_EXPORT int FIRKodeSetMassFn(void* firkode_mem, FIRKLsMassFn mass);
SUNDIALS_EXPORT int FIRKodeSetMassTimes(void* firkode_mem,
                                        FIRKLsMassTimesSetupFn msetup,
                                        FIRKLsMassTimesVecFn mtimes,
                                        void* mtimes_data);
SUNDIALS_EXPORT int FIRKodeSetMassPreconditioner(void* firkode_mem,
                                                 FIRKLsMassPrecSetupFn psetup,
                                                 FIRKLsMassPrecSolveFn psolve);
SUNDIALS_EXPORT int FIRKodeSetMassEpsLin(void* firkode_mem, sunrealtype eplifac);
SUNDIALS_EXPORT int FIRKodeSetMassLSNormFactor(void* firkode_mem,
                                               sunrealtype nrmfac);

/*-----------------------------------------------------------------
  Optional outputs from the FIRKLS linear solver interface
  -----------------------------------------------------------------*/

SUNDIALS_EXPORT int FIRKodeGetJac(void* firkode_mem,
                                  SUNMatrix* J); // nb::rv_policy::reference
SUNDIALS_EXPORT int FIRKodeGetJacTime(void* firkode_mem, sunrealtype* t_J);
SUNDIALS_EXPORT int FIRKodeGetJacNumSteps(void* firkode_mem, long int* nst_J);
SUNDIALS_EXPORT int FIRKodeGetNumJacEvals(void* firkode_mem, long int* njevals);
SUNDIALS_EXPORT int FIRKodeGetNumPrecEvals(void* firkode_mem, long int* npevals);
SUNDIALS_EXPORT int FIRKodeGetNumPrecSolves(void* firkode_mem,
                                            long int* npsolves);
SUNDIALS_EXPORT int FIRKodeGetNumLinIters(void* firkode_mem, long int* nliters);
SUNDIALS_EXPORT int FIRKodeGetNumLinConvFails(void* firkode_mem,
                                              long int* nlcfails);
SUNDIALS_EXPORT int FIRKodeGetNumJTSetupEvals(void* firkode_mem,
                                              long int* njtsetups);
SUNDIALS_EXPORT int FIRKodeGetNumJtimesEvals(void* firkode_mem,
                                             long int* njvevals);
SUNDIALS_EXPORT int FIRKodeGetNumLinRhsEvals(void* firkode_mem,
                                             long int* nfevalsLS);
SUNDIALS_EXPORT int FIRKodeGetNumBlockSolves(void* firkode_mem,
                                             long int* nbsolves);
SUNDIALS_EXPORT int FIRKodeGetNumStageLinIters(void* firkode_mem,
                                               long int* nsliters);
SUNDIALS_EXPORT int FIRKodeGetNumStageLinConvFails(void* firkode_mem,
                                                   long int* nslcfails);
SUNDIALS_EXPORT int FIRKodeGetLinSolveStats(
  void* firkode_mem, long int* njevals, long int* nfevalsLS, long int* nliters,
  long int* nlcfails, long int* npevals, long int* npsolves,
  long int* njtsetups, long int* njtimes);
SUNDIALS_EXPORT int FIRKodeGetLastLinFlag(void* firkode_mem, long int* flag);
SUNDIALS_EXPORT char* FIRKodeGetLinReturnFlagName(long int flag);

/* Mass matrix optional outputs */

SUNDIALS_EXPORT int FIRKodeGetCurrentMassMatrix(void* firkode_mem,
                                                SUNMatrix* M); // nb::rv_policy::reference
SUNDIALS_EXPORT int FIRKodeGetNumMassSetups(void* firkode_mem,
                                            long int* nmsetups);
SUNDIALS_EXPORT int FIRKodeGetNumMassMultSetups(void* firkode_mem,
                                                long int* nmvsetups);
SUNDIALS_EXPORT int FIRKodeGetNumMassMult(void* firkode_mem, long int* nmvevals);
SUNDIALS_EXPORT int FIRKodeGetNumMassSolves(void* firkode_mem,
                                            long int* nmsolves);
SUNDIALS_EXPORT int FIRKodeGetNumMassPrecEvals(void* firkode_mem,
                                               long int* nmpevals);
SUNDIALS_EXPORT int FIRKodeGetNumMassPrecSolves(void* firkode_mem,
                                                long int* nmpsolves);
SUNDIALS_EXPORT int FIRKodeGetNumMassIters(void* firkode_mem, long int* nmiters);
SUNDIALS_EXPORT int FIRKodeGetNumMassConvFails(void* firkode_mem,
                                               long int* nmcfails);
SUNDIALS_EXPORT int FIRKodeGetLastMassFlag(void* firkode_mem, long int* flag);

#ifdef __cplusplus
}
#endif

#endif
