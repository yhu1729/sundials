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
 * Adaptive-step test on the stiff Van der Pol oscillator in the form used by
 * Hairer and Wanner,
 *    y1' = y2,   y2' = ((1 - y1^2) y2 - y1) / eps,   eps = 1e-6,
 * with y(0) = (2, 0) integrated to t = 2. The solution at rtol = 1e-6 is
 * compared with a tight-tolerance reference computed in the same run.
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
constexpr sunrealtype eps  = SUN_RCONST(1.0e-6);

int f(sunrealtype, N_Vector y, N_Vector ydot, void*)
{
  const sunrealtype y1 = NV_Ith_S(y, 0);
  const sunrealtype y2 = NV_Ith_S(y, 1);
  NV_Ith_S(ydot, 0)    = y2;
  NV_Ith_S(ydot, 1)    = ((one - y1 * y1) * y2 - y1) / eps;
  return 0;
}

int Jac(sunrealtype, N_Vector y, N_Vector, SUNMatrix J, void*, N_Vector,
        N_Vector, N_Vector)
{
  const sunrealtype y1  = NV_Ith_S(y, 0);
  const sunrealtype y2  = NV_Ith_S(y, 1);
  SM_ELEMENT_D(J, 0, 0) = zero;
  SM_ELEMENT_D(J, 0, 1) = one;
  SM_ELEMENT_D(J, 1, 0) = (-SUN_RCONST(2.0) * y1 * y2 - one) / eps;
  SM_ELEMENT_D(J, 1, 1) = (one - y1 * y1) / eps;
  return 0;
}

// integrate to t = 2 with the given tolerances; returns the number of steps
long int run(SUNContext sunctx, int s, sunrealtype rtol, sunrealtype atol,
             sunrealtype* yout, int* err)
{
  N_Vector y     = N_VNew_Serial(2, sunctx);
  NV_Ith_S(y, 0) = SUN_RCONST(2.0);
  NV_Ith_S(y, 1) = zero;

  void* mem = FIRKodeCreate(sunctx);
  int flag  = FIRKodeInit(mem, f, zero, y);
  flag += FIRKodeSetNumStages(mem, s);
  flag += FIRKodeSStolerances(mem, rtol, atol);
  flag += FIRKodeSetMaxNumSteps(mem, 1000000);
  SUNMatrix A        = SUNDenseMatrix(2, 2, sunctx);
  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
  flag += FIRKodeSetLinearSolver(mem, LS, A);
  flag += FIRKodeSetJacFn(mem, Jac);
  if (flag)
  {
    *err = flag;
    return 0;
  }

  sunrealtype t;
  flag = FIRKodeEvolve(mem, SUN_RCONST(2.0), y, &t, FIRK_NORMAL);
  *err = (flag == FIRK_SUCCESS) ? 0 : flag;

  yout[0] = NV_Ith_S(y, 0);
  yout[1] = NV_Ith_S(y, 1);

  long int nst;
  FIRKodeGetNumSteps(mem, &nst);

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  return nst;
}

} // namespace

int main(int argc, char* argv[])
{
  const int s = (argc > 1) ? std::atoi(argv[1]) : 3;
  sunrealtype yref[2], y[2];
  int err, fails = 0;

  SUNContext sunctx;
  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  std::printf("Van der Pol (eps = 1e-6), Radau IIA with s = %d stages\n", s);

  const long int nst_ref = run(sunctx, s, SUN_RCONST(1.0e-10),
                               SUN_RCONST(1.0e-12), yref, &err);
  if (err)
  {
    std::printf("  FAIL: reference run failed with flag %d\n", err);
    return 1;
  }
  std::printf("  reference (rtol = 1e-10): y1(2) = %.10" ESYM
              ", y2(2) = %.10" ESYM ", steps = %li\n",
              yref[0], yref[1], nst_ref);

  const long int nst = run(sunctx, s, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-8), y,
                           &err);
  if (err)
  {
    std::printf("  FAIL: test run failed with flag %d\n", err);
    return 1;
  }
  const sunrealtype e1 = std::abs(y[0] - yref[0]);
  const sunrealtype e2 = std::abs(y[1] - yref[1]);
  std::printf("  test      (rtol = 1e-6):  y1(2) = %.10" ESYM
              ", y2(2) = %.10" ESYM ", steps = %li\n",
              y[0], y[1], nst);
  std::printf("  errors: %.3" ESYM " %.3" ESYM "\n", e1, e2);

  /* y1(2) = 1.706167732170469 (Hairer & Wanner, Solving ODEs II, VDPOL) */
  if (std::abs(yref[0] - SUN_RCONST(1.706167732170469)) > SUN_RCONST(1.0e-7))
  {
    std::printf("  FAIL: reference y1(2) disagrees with the literature\n");
    fails++;
  }
  if (e1 > SUN_RCONST(1.0e-4) || e2 > SUN_RCONST(1.0e-4) * std::abs(yref[1]))
  {
    std::printf("  FAIL: solution error exceeds tolerance\n");
    fails++;
  }
  if (nst > 2000)
  {
    std::printf("  FAIL: too many steps (%li)\n", nst);
    fails++;
  }

  SUNContext_Free(&sunctx);

  if (fails)
  {
    std::printf("\n%d test(s) FAILED\n", fails);
    return 1;
  }
  std::printf("\nAll tests passed\n");
  return 0;
}
