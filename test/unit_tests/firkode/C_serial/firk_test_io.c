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
 * Unit test of the optional inputs and outputs that no other test calls, on
 * y' = lambda y with lambda = (-1, -10):
 *
 * - before the first step: no initial step, unit tolerance scale factor,
 *   zero gamma, the current state is y0, no mass matrix interface;
 * - after a step: the actual initial step is the first step, gamma is
 *   h gamma0 of the accepted attempt, the Jacobian time is that of its
 *   evaluation, the current state is the returned solution, the integrator
 *   statistics equal the individual optional outputs, the last linear
 *   solver flag is 0;
 * - FIRKodeSetMinStep rejects a negative value and, set before the first
 *   step, bounds the initial step and every reduction below;
 * - a stop time is landed on exactly, or passed and interpolated after
 *   FIRKodeSetInterpolateStopTime, and FIRKodeClearStopTime removes it;
 * - FIRKodeSetJacTimesRhsFn routes the difference-quotient products through
 *   the alternative right-hand side, NULL restores f, and a user jtimes
 *   disables it;
 * - with a mass matrix: FIRKodeSetMassLSNormFactor sets, computes or
 *   defaults the factor, FIRKodeGetCurrentMassMatrix returns the template
 *   and the last mass solver flag is 0 after a step.
 *
 * Pass criterion: every return value and value matches.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>

#include "firkode/firkode.h"
#include "firkode/firkode_impl.h"
#include "firkode/firkode_ls.h"
#include "firkode/firkode_ls_impl.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define NEQ  2
#define RTOL SUN_RCONST(1.0e-6)

static const sunrealtype lambda[NEQ] = {SUN_RCONST(-1.0), SUN_RCONST(-10.0)};

static long int nf = 0, nf_alt = 0;
static sunrealtype t_jac = ZERO;

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  int i;
  (void)t;
  (void)user_data;
  nf++;
  for (i = 0; i < NEQ; i++) { NV_Ith_S(ydot, i) = lambda[i] * NV_Ith_S(y, i); }
  return 0;
}

static int f_alt(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  nf_alt++;
  return f(t, y, ydot, user_data);
}

static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int i;
  (void)y;
  (void)fy;
  (void)user_data;
  (void)tmp1;
  (void)tmp2;
  (void)tmp3;
  t_jac = t;
  SUNMatZero(J);
  for (i = 0; i < NEQ; i++) { SM_ELEMENT_D(J, i, i) = lambda[i]; }
  return 0;
}

static int jtimes(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y,
                  N_Vector fy, void* user_data, N_Vector tmp)
{
  int i;
  (void)t;
  (void)y;
  (void)fy;
  (void)user_data;
  (void)tmp;
  for (i = 0; i < NEQ; i++) { NV_Ith_S(Jv, i) = lambda[i] * NV_Ith_S(v, i); }
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

static int check(int ok, const char* what)
{
  if (!ok)
  {
    printf("  FAIL: %s\n", what);
    return 1;
  }
  return 0;
}

static int same(N_Vector a, N_Vector b)
{
  int i;
  for (i = 0; i < NEQ; i++)
  {
    if (NV_Ith_S(a, i) != NV_Ith_S(b, i)) { return 0; }
  }
  return 1;
}

struct run
{
  N_Vector y;
  SUNMatrix A;
  SUNLinearSolver LS;
  void* mem;
};

/* a dense (or matrix-free SPGMR) run with the analytic Jacobian */
static int setup(SUNContext sunctx, struct run* R, sunbooleantype dense)
{
  R->y = N_VNew_Serial(NEQ, sunctx);
  N_VConst(ONE, R->y);
  R->A   = dense ? SUNDenseMatrix(NEQ, NEQ, sunctx) : NULL;
  R->LS  = dense ? SUNLinSol_Dense(R->y, R->A, sunctx)
                 : SUNLinSol_SPGMR(R->y, SUN_PREC_NONE, 0, sunctx);
  R->mem = FIRKodeCreate(sunctx);
  return FIRKodeInit(R->mem, f, ZERO, R->y) ||
         FIRKodeSStolerances(R->mem, RTOL, RTOL) ||
         FIRKodeSetLinearSolver(R->mem, R->LS, R->A) ||
         (dense && FIRKodeSetJacFn(R->mem, jac));
}

static void teardown(struct run* R)
{
  FIRKodeFree(&R->mem);
  SUNLinSolFree(R->LS);
  if (R->A != NULL) { SUNMatDestroy(R->A); }
  N_VDestroy(R->y);
}

static int check_outputs(SUNContext sunctx)
{
  int fails = 0, k;
  long int nst, nfe, nsetups, netf, s_nst, s_nfe, s_nsetups, s_netf, lflag;
  sunrealtype tret, h0u, hlast, hcur, tcur, s_h0u, s_hlast, s_hcur, s_tcur;
  sunrealtype gamma, tolsf, t_J, hmin_steps;
  FIRKodeTable tbl = NULL;
  N_Vector ycur    = NULL;
  SUNMatrix M      = NULL;
  struct run R;
  char what[200];

  if (setup(sunctx, &R, SUNTRUE)) { return 1; }

  /* before the first step */
  fails += check(FIRKodeGetActualInitStep(R.mem, &h0u) == FIRK_SUCCESS &&
                   h0u == ZERO,
                 "no actual initial step before the first step");
  fails += check(FIRKodeGetTolScaleFactor(R.mem, &tolsf) == FIRK_SUCCESS &&
                   tolsf == ONE,
                 "unit tolerance scale factor before the first step");
  fails += check(FIRKodeGetCurrentGamma(R.mem, &gamma) == FIRK_SUCCESS &&
                   gamma == ZERO,
                 "zero gamma before the first step");
  fails += check(FIRKodeGetCurrentState(R.mem, &ycur) == FIRK_SUCCESS &&
                   ycur != NULL && same(ycur, R.y),
                 "the current state is y0 before the first step");
  fails += check(FIRKodeGetCurrentMassMatrix(R.mem, &M) == FIRKLS_MASSMEM_NULL &&
                   FIRKodeGetLastMassFlag(R.mem, &lflag) == FIRKLS_MASSMEM_NULL,
                 "mass outputs without a mass solver are FIRKLS_MASSMEM_NULL");

  /* one step */
  fails +=
    check(FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP) == FIRK_SUCCESS,
          "the first step succeeds");
  FIRKodeGetLastStep(R.mem, &hlast);
  FIRKodeGetCurrentTable(R.mem, &tbl);
  fails += check(FIRKodeGetActualInitStep(R.mem, &h0u) == FIRK_SUCCESS &&
                   h0u == hlast && h0u > ZERO,
                 "the actual initial step is the first step");
  fails += check(FIRKodeGetCurrentGamma(R.mem, &gamma) == FIRK_SUCCESS &&
                   tbl != NULL && gamma == hlast * tbl->gamma0,
                 "gamma is h gamma0 of the accepted attempt");
  fails += check(FIRKodeGetJacTime(R.mem, &t_J) == FIRKLS_SUCCESS &&
                   t_J == t_jac && t_J == ZERO,
                 "the Jacobian time is that of its evaluation");
  fails += check(FIRKodeGetCurrentState(R.mem, &ycur) == FIRK_SUCCESS &&
                   same(ycur, R.y),
                 "the current state is the returned solution");
  fails += check(FIRKodeGetLastLinFlag(R.mem, &lflag) == FIRKLS_SUCCESS &&
                   lflag == 0,
                 "the last linear solver flag is 0 after a step");
  fails += check(FIRKodeGetTolScaleFactor(R.mem, &tolsf) == FIRK_SUCCESS &&
                   tolsf == ONE,
                 "unit tolerance scale factor after a step");

  /* the integrator statistics equal the individual outputs */
  fails += check(FIRKodeGetIntegratorStats(R.mem, &s_nst, &s_nfe, &s_nsetups,
                                           &s_netf, &s_h0u, &s_hlast, &s_hcur,
                                           &s_tcur) == FIRK_SUCCESS,
                 "GetIntegratorStats succeeds");
  FIRKodeGetNumSteps(R.mem, &nst);
  FIRKodeGetNumRhsEvals(R.mem, &nfe);
  FIRKodeGetNumLinSolvSetups(R.mem, &nsetups);
  FIRKodeGetNumErrTestFails(R.mem, &netf);
  FIRKodeGetCurrentStep(R.mem, &hcur);
  FIRKodeGetCurrentTime(R.mem, &tcur);
  fails += check(s_nst == nst && s_nfe == nfe && s_nsetups == nsetups &&
                   s_netf == netf && s_h0u == h0u && s_hlast == hlast &&
                   s_hcur == hcur && s_tcur == tcur && nst == 1 && nfe == nf,
                 "the integrator statistics equal the individual outputs");

  fails += check(FIRKodeSetMinStep(R.mem, -ONE) == FIRK_ILL_INPUT,
                 "a negative minimum step is FIRK_ILL_INPUT");
  teardown(&R);

  /* a minimum step set before the first step bounds the steps below: the
     initial step estimate is raised to it and no reduction goes under it */
  if (setup(sunctx, &R, SUNTRUE)) { return fails + 1; }
  hmin_steps = SUN_RCONST(0.01);
  fails += check(FIRKodeSetMinStep(R.mem, hmin_steps) == FIRK_SUCCESS,
                 "SetMinStep succeeds");
  for (k = 0; k < 4; k++)
  {
    if (FIRKodeEvolve(R.mem, SUN_RCONST(10.0), R.y, &tret, FIRK_ONE_STEP) !=
        FIRK_SUCCESS)
    {
      break;
    }
    FIRKodeGetLastStep(R.mem, &hlast);
    snprintf(what, sizeof(what), "step %d with SetMinStep: h = " SUN_FORMAT_G,
             k + 1, hlast);
    fails += check(hlast >= hmin_steps, what);
  }
  fails += check(k == 4, "the steps with SetMinStep succeed");
  teardown(&R);
  return fails;
}

static int check_stop_time(SUNContext sunctx)
{
  int fails = 0, flag;
  sunrealtype tret, tcur;
  const sunrealtype tstop1 = SUN_RCONST(0.3), tstop2 = SUN_RCONST(0.5),
                    tstop3 = SUN_RCONST(0.7), tout = ONE;
  struct run R;

  if (setup(sunctx, &R, SUNTRUE) || FIRKodeSetStopTime(R.mem, tstop1))
  {
    return 1;
  }
  flag = FIRKodeEvolve(R.mem, tout, R.y, &tret, FIRK_NORMAL);
  FIRKodeGetCurrentTime(R.mem, &tcur);
  fails += check(flag == FIRK_TSTOP_RETURN && tret == tstop1 && tcur == tstop1,
                 "a stop time is landed on exactly");

  fails += check(FIRKodeSetStopTime(R.mem, tstop2) == FIRK_SUCCESS &&
                   FIRKodeSetInterpolateStopTime(R.mem, SUNTRUE) == FIRK_SUCCESS,
                 "SetInterpolateStopTime succeeds");
  flag = FIRKodeEvolve(R.mem, tout, R.y, &tret, FIRK_NORMAL);
  FIRKodeGetCurrentTime(R.mem, &tcur);
  fails += check(flag == FIRK_TSTOP_RETURN && tret == tstop2 && tcur >= tstop2 &&
                   SUNRabs(NV_Ith_S(R.y, 0) - SUNRexp(lambda[0] * tstop2)) <=
                     SUN_RCONST(100.0) * RTOL,
                 "an interpolated stop time is returned from a step past it");

  fails += check(FIRKodeSetStopTime(R.mem, tstop3) == FIRK_SUCCESS &&
                   FIRKodeClearStopTime(R.mem) == FIRK_SUCCESS,
                 "ClearStopTime succeeds");
  flag = FIRKodeEvolve(R.mem, tout, R.y, &tret, FIRK_NORMAL);
  fails += check(flag == FIRK_SUCCESS && tret == tout,
                 "a cleared stop time is passed");
  teardown(&R);
  return fails;
}

static int check_jtimes_rhs(SUNContext sunctx)
{
  int fails      = 0;
  long int nfeDQ = 0, nf_alt_before;
  sunrealtype tret;
  struct run R;

  /* matrix-free, no user jtimes: the difference quotient uses the
     alternative right-hand side */
  if (setup(sunctx, &R, SUNFALSE)) { return 1; }
  nf = nf_alt = 0;
  fails += check(FIRKodeSetJacTimesRhsFn(R.mem, f_alt) == FIRKLS_SUCCESS,
                 "SetJacTimesRhsFn succeeds for the difference quotient");
  fails +=
    check(FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP) == FIRK_SUCCESS,
          "the step with the alternative right-hand side succeeds");
  FIRKodeGetNumLinRhsEvals(R.mem, &nfeDQ);
  fails += check(nfeDQ > 0 && nf_alt == nfeDQ && nf > nf_alt,
                 "the difference-quotient products call the alternative "
                 "right-hand side only");

  /* NULL restores f */
  fails += check(FIRKodeSetJacTimesRhsFn(R.mem, NULL) == FIRKLS_SUCCESS,
                 "SetJacTimesRhsFn(NULL) succeeds");
  nf_alt_before = nf_alt;
  fails += check(FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP) ==
                     FIRK_SUCCESS &&
                   nf_alt == nf_alt_before,
                 "f computes the products again after SetJacTimesRhsFn(NULL)");

  /* a user jtimes disables the difference quotient */
  fails += check(FIRKodeSetJacTimes(R.mem, NULL, jtimes) == FIRKLS_SUCCESS &&
                   FIRKodeSetJacTimesRhsFn(R.mem, f_alt) == FIRKLS_ILL_INPUT,
                 "SetJacTimesRhsFn with a user jtimes is FIRKLS_ILL_INPUT");
  teardown(&R);
  return fails;
}

static int check_mass_io(SUNContext sunctx)
{
  int fails      = 0;
  long int lflag = -1;
  sunrealtype tret;
  SUNMatrix M, Mcur = NULL;
  SUNLinearSolver MLS;
  FIRKLsMassMem mmem;
  struct run R;

  if (setup(sunctx, &R, SUNTRUE)) { return 1; }
  M   = SUNDenseMatrix(NEQ, NEQ, sunctx);
  MLS = SUNLinSol_Dense(R.y, M, sunctx);
  if (FIRKodeSetMassLinearSolver(R.mem, MLS, M, SUNFALSE) ||
      FIRKodeSetMassFn(R.mem, mass))
  {
    return 1;
  }
  mmem = (FIRKLsMassMem)((FIRKodeMem)R.mem)->mass_mem;

  fails += check(FIRKodeSetMassLSNormFactor(R.mem, SUN_RCONST(2.5)) ==
                     FIRKLS_SUCCESS &&
                   mmem->nrmfac == SUN_RCONST(2.5),
                 "a positive mass norm factor is set");
  fails += check(FIRKodeSetMassLSNormFactor(R.mem, -ONE) == FIRKLS_SUCCESS &&
                   mmem->nrmfac == SUNRsqrt((sunrealtype)NEQ),
                 "a negative mass norm factor is computed from a dot product");
  fails += check(FIRKodeSetMassLSNormFactor(R.mem, SUN_RCONST(2.5)) ==
                     FIRKLS_SUCCESS &&
                   FIRKodeSetMassLSNormFactor(R.mem, ZERO) == FIRKLS_SUCCESS &&
                   mmem->nrmfac == SUNRsqrt((sunrealtype)NEQ),
                 "a zero mass norm factor is the default");
  fails += check(FIRKodeGetCurrentMassMatrix(R.mem, &Mcur) == FIRKLS_SUCCESS &&
                   Mcur == M,
                 "GetCurrentMassMatrix returns the template");
  fails +=
    check(FIRKodeEvolve(R.mem, ONE, R.y, &tret, FIRK_ONE_STEP) == FIRK_SUCCESS,
          "the step with the mass matrix succeeds");
  fails += check(FIRKodeGetLastMassFlag(R.mem, &lflag) == FIRKLS_SUCCESS &&
                   lflag == 0,
                 "the last mass solver flag is 0 after a step");
  teardown(&R);
  SUNLinSolFree(MLS);
  SUNMatDestroy(M);
  return fails;
}

int main(void)
{
  SUNContext sunctx = NULL;
  int fails         = 0;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  fails += check_outputs(sunctx);
  fails += check_stop_time(sunctx);
  fails += check_jtimes_rhs(sunctx);
  fails += check_mass_io(sunctx);

  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  printf("All tests passed\n");
  return 0;
}
