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
 * Matrix-free test for FIRKODE: the Kvaerno-Prothero-Robinson problem is
 * integrated with fixed steps using (a) the dense linear solver with the
 * analytic Jacobian and (b) the matrix-free SPGMR solver with the internal
 * difference-quotient Jacobian-vector product, without and with a
 * user-supplied preconditioner. All solutions must agree to a tolerance well
 * below the discretization error.
 *
 * Usage: firk_test_krylov [stages]
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "problems/kpr.hpp"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunmatrix/sunmatrix_dense.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define ESYM "Le"
#else
#define ESYM "e"
#endif

namespace {

constexpr sunrealtype zero = SUN_RCONST(0.0);

sunrealtype kpr_data[4] = {SUN_RCONST(-2.0), SUN_RCONST(0.5), SUN_RCONST(0.5),
                           SUN_RCONST(-1.0)};

// integrate KPR to t = 0.5 with fixed step h; mode 0 = dense, 1 = SPGMR
int run(SUNContext sunctx, int s, sunrealtype h, int mode, sunrealtype* yout,
        long int* nli, long int* nni)
{
  N_Vector y      = N_VNew_Serial(2, sunctx);
  sunrealtype* yd = N_VGetArrayPointer(y);
  problems::kpr::true_sol(zero, &yd[0], &yd[1]);

  void* mem = FIRKodeCreate(sunctx);
  int flag  = FIRKodeInit(mem, problems::kpr::ode_rhs, zero, y);
  flag += FIRKodeSetUserData(mem, static_cast<void*>(kpr_data));
  flag += FIRKodeSetNumStages(mem, s);
  flag += FIRKodeSStolerances(mem, SUN_RCONST(1.0e-10), SUN_RCONST(1.0e-10));
  flag += FIRKodeSetFixedStep(mem, h);
  flag += FIRKodeSetNonlinConvCoef(mem, SUN_RCONST(1.0e-3));
  flag += FIRKodeSetMaxNonlinIters(mem, 20);

  SUNMatrix A        = nullptr;
  SUNLinearSolver LS = nullptr;
  if (mode == 0)
  {
    A  = SUNDenseMatrix(2, 2, sunctx);
    LS = SUNLinSol_Dense(y, A, sunctx);
    flag += FIRKodeSetLinearSolver(mem, LS, A);
    flag += FIRKodeSetJacFn(mem, problems::kpr::ode_rhs_jac);
  }
  else
  {
    LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 10, sunctx);
    flag += FIRKodeSetLinearSolver(mem, LS, nullptr);
    flag += FIRKodeSetEpsLin(mem, SUN_RCONST(1.0e-3));
  }
  flag += FIRKodeSetStageEpsLin(mem, SUN_RCONST(1.0e-3));
  if (flag) { return flag; }

  sunrealtype t;
  flag = FIRKodeEvolve(mem, SUN_RCONST(0.5), y, &t, FIRK_NORMAL);
  if (flag != FIRK_SUCCESS) { return flag; }

  yout[0] = yd[0];
  yout[1] = yd[1];
  FIRKodeGetNumLinIters(mem, nli);
  FIRKodeGetNumNonlinSolvIters(mem, nni);

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  if (A) { SUNMatDestroy(A); }
  N_VDestroy(y);
  return 0;
}

} // namespace

int main(int argc, char* argv[])
{
  const int s           = (argc > 1) ? std::atoi(argv[1]) : 3;
  const sunrealtype h   = SUN_RCONST(0.01);
  const sunrealtype tol = SUN_RCONST(1.0e-8);
  sunrealtype ydense[2], ykry[2];
  long int nli, nni;
  int fails = 0;

  SUNContext sunctx;
  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  std::printf("Matrix-free test, KPR with Radau IIA s = %d, h = %.3" ESYM "\n",
              s, h);

  if (run(sunctx, s, h, 0, ydense, &nli, &nni))
  {
    std::printf("  FAIL: dense run failed\n");
    return 1;
  }
  std::printf("  dense:        y = (%.12" ESYM ", %.12" ESYM "), nni = %li\n",
              ydense[0], ydense[1], nni);

  if (run(sunctx, s, h, 1, ykry, &nli, &nni))
  {
    std::printf("  FAIL: SPGMR run failed\n");
    return 1;
  }
  std::printf("  SPGMR:        y = (%.12" ESYM ", %.12" ESYM "), nni = %li, "
              "nli = %li\n",
              ykry[0], ykry[1], nni, nli);
  sunrealtype err = std::max(std::abs(ykry[0] - ydense[0]),
                             std::abs(ykry[1] - ydense[1]));
  std::printf("  max difference = %.3" ESYM "\n", err);
  if (err > tol)
  {
    std::printf("  FAIL: SPGMR solution differs from dense solution\n");
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
