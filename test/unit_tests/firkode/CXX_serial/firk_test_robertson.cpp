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
 * Adaptive-step test on the stiff Robertson kinetics problem
 *    u' = -0.04 u + 1e4 v w
 *    v' =  0.04 u - 1e4 v w - 3e7 v^2
 *    w' =  3e7 v^2
 * integrated to t = 1e11 and compared with the reference solution used by the
 * ARKODE examples (rtol = 1e-8, atol = 1e-14). Conservation u + v + w = 1 and
 * a bound on the number of steps are checked as well.
 * ---------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define ESYM "Le"
#else
#define ESYM "e"
#endif

namespace {

constexpr sunrealtype zero = SUN_RCONST(0.0);
constexpr sunrealtype one  = SUN_RCONST(1.0);

int f(sunrealtype, N_Vector y, N_Vector ydot, void*)
{
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);
  const sunrealtype w = NV_Ith_S(y, 2);
  NV_Ith_S(ydot, 0)   = -SUN_RCONST(0.04) * u + SUN_RCONST(1.0e4) * v * w;
  NV_Ith_S(ydot, 1)   = SUN_RCONST(0.04) * u - SUN_RCONST(1.0e4) * v * w -
                      SUN_RCONST(3.0e7) * v * v;
  NV_Ith_S(ydot, 2) = SUN_RCONST(3.0e7) * v * v;
  return 0;
}

int Jac(sunrealtype, N_Vector y, N_Vector, SUNMatrix J, void*, N_Vector,
        N_Vector, N_Vector)
{
  const sunrealtype v = NV_Ith_S(y, 1);
  const sunrealtype w = NV_Ith_S(y, 2);
  SUNMatZero(J);
  SM_ELEMENT_D(J, 0, 0) = -SUN_RCONST(0.04);
  SM_ELEMENT_D(J, 0, 1) = SUN_RCONST(1.0e4) * w;
  SM_ELEMENT_D(J, 0, 2) = SUN_RCONST(1.0e4) * v;
  SM_ELEMENT_D(J, 1, 0) = SUN_RCONST(0.04);
  SM_ELEMENT_D(J, 1, 1) = -SUN_RCONST(1.0e4) * w - SUN_RCONST(6.0e7) * v;
  SM_ELEMENT_D(J, 1, 2) = -SUN_RCONST(1.0e4) * v;
  SM_ELEMENT_D(J, 2, 1) = SUN_RCONST(6.0e7) * v;
  return 0;
}

} // namespace

int main(int argc, char* argv[])
{
  const int s              = (argc > 1) ? std::atoi(argv[1]) : 3;
  const sunrealtype rtol   = SUN_RCONST(1.0e-6);
  const sunrealtype tf     = SUN_RCONST(1.0e11);
  const sunrealtype ref[3] = {SUN_RCONST(2.0833403356917897e-08),
                              SUN_RCONST(8.1470714598028223e-14),
                              SUN_RCONST(9.9999997916651040e-01)};
  int fails                = 0;

  SUNContext sunctx;
  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  N_Vector y        = N_VNew_Serial(3, sunctx);
  N_Vector atol     = N_VClone(y);
  NV_Ith_S(y, 0)    = one;
  NV_Ith_S(y, 1)    = zero;
  NV_Ith_S(y, 2)    = zero;
  NV_Ith_S(atol, 0) = SUN_RCONST(1.0e-8);
  NV_Ith_S(atol, 1) = SUN_RCONST(1.0e-14);
  NV_Ith_S(atol, 2) = SUN_RCONST(1.0e-8);

  void* mem = FIRKodeCreate(sunctx);
  int flag  = FIRKodeInit(mem, f, zero, y);
  flag += FIRKodeSetNumStages(mem, s);
  flag += FIRKodeSVtolerances(mem, rtol, atol);
  flag += FIRKodeSetMaxNumSteps(mem, 100000);
  SUNMatrix A        = SUNDenseMatrix(3, 3, sunctx);
  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
  flag += FIRKodeSetLinearSolver(mem, LS, A);
  flag += FIRKodeSetJacFn(mem, Jac);
  if (flag)
  {
    std::printf("  FAIL: setup returned %d\n", flag);
    return 1;
  }

  sunrealtype t;
  flag = FIRKodeEvolve(mem, tf, y, &t, FIRK_NORMAL);
  if (flag != FIRK_SUCCESS)
  {
    std::printf("  FAIL: FIRKodeEvolve returned %d\n", flag);
    return 1;
  }

  long int nst, netf, nni, ncfn;
  FIRKodeGetNumSteps(mem, &nst);
  FIRKodeGetNumErrTestFails(mem, &netf);
  FIRKodeGetNumNonlinSolvIters(mem, &nni);
  FIRKodeGetNumStepSolveFails(mem, &ncfn);

  std::printf("Robertson, Radau IIA with s = %d stages, rtol = %.1" ESYM "\n",
              s, rtol);
  std::printf("  steps = %li, error test fails = %li, Newton iters = %li, "
              "solve fails = %li\n",
              nst, netf, nni, ncfn);

  sunrealtype sum = zero;
  for (int i = 0; i < 3; i++)
  {
    const sunrealtype err = std::abs(NV_Ith_S(y, i) - ref[i]);
    const sunrealtype tol = SUN_RCONST(50.0) * rtol * std::abs(ref[i]) +
                            SUN_RCONST(10.0) * NV_Ith_S(atol, i);
    sum += NV_Ith_S(y, i);
    std::printf("  y[%d] = %.10" ESYM "  error = %.3" ESYM "  tol = %.3" ESYM
                "\n",
                i, NV_Ith_S(y, i), err, tol);
    if (err > tol)
    {
      std::printf("  FAIL: component %d exceeds tolerance\n", i);
      fails++;
    }
  }
  if (std::abs(sum - one) > SUN_RCONST(1.0e-8))
  {
    std::printf("  FAIL: u + v + w - 1 = %.3" ESYM "\n", sum - one);
    fails++;
  }
  if (nst > 2000)
  {
    std::printf("  FAIL: too many steps (%li)\n", nst);
    fails++;
  }
  if (ncfn > nst / 10)
  {
    std::printf("  FAIL: too many nonlinear solver failures (%li)\n", ncfn);
    fails++;
  }

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  N_VDestroy(atol);
  SUNContext_Free(&sunctx);

  if (fails)
  {
    std::printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  std::printf("\nAll tests passed\n");
  return 0;
}
