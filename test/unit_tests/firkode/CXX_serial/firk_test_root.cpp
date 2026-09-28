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
 * Rootfinding test on the Robertson kinetics problem: the roots of
 * g1 = u - 1e-4 and g2 = w - 1e-2 are located at t = 2.0792e7 (u decreasing)
 * and t = 2.640e-1 (w increasing) respectively (reference values from the
 * CVODE and ARKODE examples), and the returned states must satisfy the root
 * conditions to the dense-output accuracy.
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

int g(sunrealtype, N_Vector y, sunrealtype* gout, void*)
{
  gout[0] = NV_Ith_S(y, 0) - SUN_RCONST(1.0e-4);
  gout[1] = NV_Ith_S(y, 2) - SUN_RCONST(1.0e-2);
  return 0;
}

} // namespace

int main(int argc, char* argv[])
{
  const int s = (argc > 1) ? std::atoi(argv[1]) : 3;
  int fails   = 0;

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
  flag += FIRKodeSVtolerances(mem, SUN_RCONST(1.0e-6), atol);
  flag += FIRKodeSetMaxNumSteps(mem, 100000);
  flag += FIRKodeRootInit(mem, 2, g);
  SUNMatrix A        = SUNDenseMatrix(3, 3, sunctx);
  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
  flag += FIRKodeSetLinearSolver(mem, LS, A);
  flag += FIRKodeSetJacFn(mem, Jac);
  if (flag)
  {
    std::printf("  FAIL: setup returned %d\n", flag);
    return 1;
  }

  std::printf("Robertson rootfinding, Radau IIA with s = %d stages\n", s);

  /* expected roots: (time, tolerance, rootsfound[0], rootsfound[1]) */
  const sunrealtype troot[2] = {SUN_RCONST(2.640e-1), SUN_RCONST(2.0792e7)};
  const sunrealtype ttol[2]  = {SUN_RCONST(2.0e-3), SUN_RCONST(3.0e4)};
  const int info[2][2]       = {{0, 1}, {-1, 0}};

  sunrealtype t, gout[2];
  int rootsfound[2], nroots = 0;
  const sunrealtype tout = SUN_RCONST(1.0e8);
  while (nroots < 2)
  {
    flag = FIRKodeEvolve(mem, tout, y, &t, FIRK_NORMAL);
    if (flag == FIRK_ROOT_RETURN)
    {
      FIRKodeGetRootInfo(mem, rootsfound);
      g(t, y, gout, nullptr);
      std::printf("  root %d at t = %.6" ESYM
                  ": rootsfound = %d %d, g = %.3" ESYM " %.3" ESYM "\n",
                  nroots + 1, t, rootsfound[0], rootsfound[1], gout[0], gout[1]);
      if (std::abs(t - troot[nroots]) > ttol[nroots])
      {
        std::printf("  FAIL: root time differs from reference %.4" ESYM "\n",
                    troot[nroots]);
        fails++;
      }
      if (rootsfound[0] != info[nroots][0] || rootsfound[1] != info[nroots][1])
      {
        std::printf("  FAIL: unexpected root information\n");
        fails++;
      }
      const sunrealtype gval = (nroots == 0) ? gout[1] : gout[0];
      if (std::abs(gval) > SUN_RCONST(1.0e-8))
      {
        std::printf("  FAIL: root condition not satisfied at the returned "
                    "state\n");
        fails++;
      }
      nroots++;
    }
    else if (flag == FIRK_SUCCESS)
    {
      std::printf("  FAIL: reached tout = %.3" ESYM " with %d root(s) found\n",
                  tout, nroots);
      fails++;
      break;
    }
    else
    {
      std::printf("  FAIL: FIRKodeEvolve returned %d\n", flag);
      fails++;
      break;
    }
  }

  long int nge;
  FIRKodeGetNumGEvals(mem, &nge);
  std::printf("  root function evaluations = %li\n", nge);

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
