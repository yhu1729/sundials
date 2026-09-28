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
 * This is the header file for the main FIRKODE integrator: fully
 * implicit Radau IIA Runge-Kutta methods for stiff ODE systems
 *   M y' = f(t,y),  y(t0) = y0.
 * -----------------------------------------------------------------*/

#ifndef _FIRKODE_H
#define _FIRKODE_H

#include <stdio.h>

#include <firkode/firkode_tables.h>
#include <sundials/sundials_adaptcontroller.h>
#include <sundials/sundials_core.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -------------------
 * FIRKODE Constants
 * ------------------- */

/* itask */
#define FIRK_NORMAL   1
#define FIRK_ONE_STEP 2

/* stage predictor */
#define FIRK_PREDICT_EXTRAPOLATE 0
#define FIRK_PREDICT_TRIVIAL     1

/* return values */

#define FIRK_SUCCESS      0
#define FIRK_TSTOP_RETURN 1
#define FIRK_ROOT_RETURN  2

#define FIRK_WARNING 99

#define FIRK_TOO_MUCH_WORK -1
#define FIRK_TOO_MUCH_ACC  -2
#define FIRK_ERR_FAILURE   -3
#define FIRK_CONV_FAILURE  -4

#define FIRK_LINIT_FAIL        -5
#define FIRK_LSETUP_FAIL       -6
#define FIRK_LSOLVE_FAIL       -7
#define FIRK_RHSFUNC_FAIL      -8
#define FIRK_FIRST_RHSFUNC_ERR -9
#define FIRK_REPTD_RHSFUNC_ERR -10
#define FIRK_UNREC_RHSFUNC_ERR -11
#define FIRK_RTFUNC_FAIL       -12
#define FIRK_NLS_INIT_FAIL     -13
#define FIRK_NLS_SETUP_FAIL    -14
#define FIRK_NLS_FAIL          -15
#define FIRK_MASSINIT_FAIL     -16
#define FIRK_MASSSETUP_FAIL    -17
#define FIRK_MASSSOLVE_FAIL    -18
#define FIRK_MASSMULT_FAIL     -19

#define FIRK_MEM_FAIL       -20
#define FIRK_MEM_NULL       -21
#define FIRK_ILL_INPUT      -22
#define FIRK_NO_MALLOC      -23
#define FIRK_BAD_K          -24
#define FIRK_BAD_T          -25
#define FIRK_BAD_DKY        -26
#define FIRK_TOO_CLOSE      -27
#define FIRK_VECTOROP_ERR   -28
#define FIRK_CONTROLLER_ERR -29
#define FIRK_TABLE_FAIL     -30

#define FIRK_CONTEXT_ERR -32

#define FIRK_UNRECOGNIZED_ERR -99

/* ------------------------------
 * User-Supplied Function Types
 * ------------------------------ */

typedef int (*FIRKRhsFn)(sunrealtype t, N_Vector y, N_Vector ydot,
                         void* user_data);

typedef int (*FIRKRootFn)(sunrealtype t, N_Vector y, sunrealtype* gout_1d,
                          void* user_data);

typedef int (*FIRKEwtFn)(N_Vector y, N_Vector ewt, void* user_data);

/* -------------------
 * Exported Functions
 * ------------------- */

/* Initialization functions */

SUNDIALS_EXPORT void* FIRKodeCreate(SUNContext sunctx);

SUNDIALS_EXPORT int FIRKodeInit(void* firkode_mem, FIRKRhsFn f, sunrealtype t0,
                                N_Vector y0);
SUNDIALS_EXPORT int FIRKodeReInit(void* firkode_mem, sunrealtype t0, N_Vector y0);
SUNDIALS_EXPORT int FIRKodeReset(void* firkode_mem, sunrealtype tR, N_Vector yR);

/* Tolerance input functions */

SUNDIALS_EXPORT int FIRKodeSStolerances(void* firkode_mem, sunrealtype reltol,
                                        sunrealtype abstol);
SUNDIALS_EXPORT int FIRKodeSVtolerances(void* firkode_mem, sunrealtype reltol,
                                        N_Vector abstol);
SUNDIALS_EXPORT int FIRKodeWFtolerances(void* firkode_mem, FIRKEwtFn efun);

/* Method selection */

SUNDIALS_EXPORT int FIRKodeSetNumStages(void* firkode_mem, int s);
SUNDIALS_EXPORT int FIRKodeSetOrder(void* firkode_mem, int ord);
SUNDIALS_EXPORT int FIRKodeGetNumStages(void* firkode_mem, int* s);
SUNDIALS_EXPORT int FIRKodeGetCurrentTable(void* firkode_mem,
                                           FIRKodeTable* T); // nb::rv_policy::reference

/* Optional input functions */

SUNDIALS_EXPORT int FIRKodeSetUserData(void* firkode_mem, void* user_data);
SUNDIALS_EXPORT int FIRKodeSetMaxNumSteps(void* firkode_mem, long int mxsteps);
SUNDIALS_EXPORT int FIRKodeSetMaxHnilWarns(void* firkode_mem, int mxhnil);
SUNDIALS_EXPORT int FIRKodeSetInitStep(void* firkode_mem, sunrealtype hin);
SUNDIALS_EXPORT int FIRKodeSetMinStep(void* firkode_mem, sunrealtype hmin);
SUNDIALS_EXPORT int FIRKodeSetMaxStep(void* firkode_mem, sunrealtype hmax);
SUNDIALS_EXPORT int FIRKodeSetFixedStep(void* firkode_mem, sunrealtype hfixed);
SUNDIALS_EXPORT int FIRKodeSetStopTime(void* firkode_mem, sunrealtype tstop);
SUNDIALS_EXPORT int FIRKodeSetInterpolateStopTime(void* firkode_mem,
                                                  sunbooleantype interp);
SUNDIALS_EXPORT int FIRKodeClearStopTime(void* firkode_mem);
SUNDIALS_EXPORT int FIRKodeSetMaxErrTestFails(void* firkode_mem, int maxnef);
SUNDIALS_EXPORT int FIRKodeSetMaxConvFails(void* firkode_mem, int maxncf);
SUNDIALS_EXPORT int FIRKodeSetLSetupFrequency(void* firkode_mem, long int msbp);
SUNDIALS_EXPORT int FIRKodeSetDeltaGammaMaxLSetup(void* firkode_mem,
                                                  sunrealtype dgmax_lsetup);
SUNDIALS_EXPORT int FIRKodeSetPredictorMethod(void* firkode_mem, int method);
SUNDIALS_EXPORT int FIRKodeSetErrorRefilter(void* firkode_mem,
                                            sunbooleantype onoff);

/* Optional step adaptivity input functions */

SUNDIALS_EXPORT int FIRKodeSetAdaptController(void* firkode_mem,
                                              SUNAdaptController C);
SUNDIALS_EXPORT int FIRKodeSetSafetyFactor(void* firkode_mem, sunrealtype safety);
SUNDIALS_EXPORT int FIRKodeSetErrorBias(void* firkode_mem, sunrealtype bias);
SUNDIALS_EXPORT int FIRKodeSetMaxGrowth(void* firkode_mem, sunrealtype mx_growth);
SUNDIALS_EXPORT int FIRKodeSetMinReduction(void* firkode_mem,
                                           sunrealtype eta_min);
SUNDIALS_EXPORT int FIRKodeSetFixedStepBounds(void* firkode_mem, sunrealtype lb,
                                              sunrealtype ub);
SUNDIALS_EXPORT int FIRKodeSetMaxFirstGrowth(void* firkode_mem,
                                             sunrealtype etamx1);
SUNDIALS_EXPORT int FIRKodeSetMaxEFailGrowth(void* firkode_mem,
                                             sunrealtype etamxf);
SUNDIALS_EXPORT int FIRKodeSetSmallNumEFails(void* firkode_mem, int small_nef);
SUNDIALS_EXPORT int FIRKodeSetMaxCFailGrowth(void* firkode_mem,
                                             sunrealtype etacf);
SUNDIALS_EXPORT int FIRKodeSetNewtonCountSafety(void* firkode_mem,
                                                sunbooleantype onoff);

/* Optional nonlinear solver input functions */

SUNDIALS_EXPORT int FIRKodeSetMaxNonlinIters(void* firkode_mem, int maxcor);
SUNDIALS_EXPORT int FIRKodeSetNonlinConvCoef(void* firkode_mem,
                                             sunrealtype nlscoef);
SUNDIALS_EXPORT int FIRKodeSetNonlinDivergenceRate(void* firkode_mem,
                                                   sunrealtype theta_max);
SUNDIALS_EXPORT int FIRKodeSetJacBadConvRate(void* firkode_mem,
                                             sunrealtype theta_jbad);

/* Rootfinding initialization function */

SUNDIALS_EXPORT int FIRKodeRootInit(void* firkode_mem, int nrtfn, FIRKRootFn g);

/* Rootfinding optional input functions */

SUNDIALS_EXPORT int FIRKodeSetRootDirection(void* firkode_mem, int* rootdir_1d);
SUNDIALS_EXPORT int FIRKodeSetNoInactiveRootWarn(void* firkode_mem);

/* Solver function */

SUNDIALS_EXPORT int FIRKodeEvolve(void* firkode_mem, sunrealtype tout,
                                  N_Vector yout, sunrealtype* tret, int itask);

/* Dense output function */

SUNDIALS_EXPORT int FIRKodeGetDky(void* firkode_mem, sunrealtype t, int k,
                                  N_Vector dky);

/* Optional output functions */

SUNDIALS_EXPORT int FIRKodeGetNumSteps(void* firkode_mem, long int* nsteps);
SUNDIALS_EXPORT int FIRKodeGetNumStepAttempts(void* firkode_mem,
                                              long int* step_attempts);
SUNDIALS_EXPORT int FIRKodeGetNumRhsEvals(void* firkode_mem, long int* nfevals);
SUNDIALS_EXPORT int FIRKodeGetNumLinSolvSetups(void* firkode_mem,
                                               long int* nlinsetups);
SUNDIALS_EXPORT int FIRKodeGetNumErrTestFails(void* firkode_mem,
                                              long int* netfails);
SUNDIALS_EXPORT int FIRKodeGetNumStepSolveFails(void* firkode_mem,
                                                long int* nncfails);
SUNDIALS_EXPORT int FIRKodeGetNumNonlinSolvIters(void* firkode_mem,
                                                 long int* nniters);
SUNDIALS_EXPORT int FIRKodeGetNumNonlinSolvConvFails(void* firkode_mem,
                                                     long int* nnfails);
SUNDIALS_EXPORT int FIRKodeGetNonlinSolvStats(void* firkode_mem,
                                              long int* nniters,
                                              long int* nnfails);
SUNDIALS_EXPORT int FIRKodeGetActualInitStep(void* firkode_mem,
                                             sunrealtype* hinused);
SUNDIALS_EXPORT int FIRKodeGetLastStep(void* firkode_mem, sunrealtype* hlast);
SUNDIALS_EXPORT int FIRKodeGetCurrentStep(void* firkode_mem, sunrealtype* hcur);
SUNDIALS_EXPORT int FIRKodeGetCurrentTime(void* firkode_mem, sunrealtype* tcur);
SUNDIALS_EXPORT int FIRKodeGetCurrentState(void* firkode_mem,
                                           N_Vector* y); // nb::rv_policy::reference
SUNDIALS_EXPORT int FIRKodeGetCurrentGamma(void* firkode_mem, sunrealtype* gamma);
SUNDIALS_EXPORT int FIRKodeGetTolScaleFactor(void* firkode_mem,
                                             sunrealtype* tolsfac);
SUNDIALS_EXPORT int FIRKodeGetErrWeights(void* firkode_mem, N_Vector eweight);
SUNDIALS_EXPORT int FIRKodeGetEstLocalErrors(void* firkode_mem, N_Vector ele);
SUNDIALS_EXPORT int FIRKodeGetIntegratorStats(
  void* firkode_mem, long int* nsteps, long int* nfevals, long int* nlinsetups,
  long int* netfails, sunrealtype* hinused, sunrealtype* hlast,
  sunrealtype* hcur, sunrealtype* tcur);
SUNDIALS_EXPORT int FIRKodeGetNumGEvals(void* firkode_mem, long int* ngevals);
SUNDIALS_EXPORT int FIRKodeGetRootInfo(void* firkode_mem, int* rootsfound_1d);
SUNDIALS_EXPORT int FIRKodeGetUserData(void* firkode_mem, void** user_data);
SUNDIALS_EXPORT int FIRKodePrintAllStats(void* firkode_mem, FILE* outfile,
                                         SUNOutputFormat fmt);
SUNDIALS_EXPORT int FIRKodeWriteParameters(void* firkode_mem, FILE* fp);
SUNDIALS_EXPORT char* FIRKodeGetReturnFlagName(long int flag);

/* Free function */

SUNDIALS_EXPORT void FIRKodeFree(void** firkode_mem);

/* FIRKLS interface function that depends on FIRKRhsFn */

SUNDIALS_EXPORT int FIRKodeSetJacTimesRhsFn(void* firkode_mem,
                                            FIRKRhsFn jtimesRhsFn);

#ifdef __cplusplus
}
#endif

#endif
