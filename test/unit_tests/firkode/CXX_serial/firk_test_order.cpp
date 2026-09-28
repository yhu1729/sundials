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
 * Fixed-step convergence test for the FIRKODE Radau IIA methods.
 *
 * Integrates a problem with a known solution over a fixed interval with a
 * sequence of halved step sizes and estimates the observed order of accuracy
 * from consecutive errors. The s-stage Radau IIA method has classical order
 * 2s-1 and stage order s; on the stiff Prothero-Robinson problem the observed
 * order reduces towards the stage order.
 *
 * Usage: firk_test_order <pr|kpr> <stages> [lambda]
 *   pr  : Prothero-Robinson, y' = lambda (y - phi(t)) + phi'(t), band Jacobian
 *   kpr : Kvaerno-Prothero-Robinson nonlinear 2x2 system, dense Jacobian
 * ---------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "firkode/firkode.h"
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "problems/kpr.hpp"
#include "problems/pr.hpp"
#include "sunlinsol/sunlinsol_band.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_band.h"
#include "sunmatrix/sunmatrix_dense.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define ESYM "Le"
#else
#define ESYM "e"
#endif

namespace {

constexpr int num_refinements = 4;
constexpr sunrealtype zero    = SUN_RCONST(0.0);
constexpr sunrealtype one     = SUN_RCONST(1.0);

struct Problem
{
  std::string name;
  sunindextype neq;
  sunrealtype t0, tf, h0;
  FIRKRhsFn rhs;
  FIRKLsJacFn jac;
  void* udata;
  bool band;
  int (*exact)(sunrealtype t, N_Vector y, void* udata);
};

int pr_exact(sunrealtype t, N_Vector y, void*)
{
  return problems::pr::true_solution(t, y);
}

int kpr_exact(sunrealtype t, N_Vector y, void*)
{
  sunrealtype* yd = N_VGetArrayPointer(y);
  return problems::kpr::true_sol(t, &yd[0], &yd[1]);
}

// integrate with fixed step h and return the max-norm error at tf
int run(SUNContext sunctx, const Problem& p, int s, sunrealtype h,
        sunrealtype& err)
{
  N_Vector y    = N_VNew_Serial(p.neq, sunctx);
  N_Vector yref = N_VClone(y);
  p.exact(p.t0, y, p.udata);

  void* mem = FIRKodeCreate(sunctx);
  int flag  = FIRKodeInit(mem, p.rhs, p.t0, y);
  flag += FIRKodeSetUserData(mem, p.udata);
  flag += FIRKodeSetNumStages(mem, s);
  flag += FIRKodeSStolerances(mem, SUN_RCONST(1.0e-12), SUN_RCONST(1.0e-12));
  flag += FIRKodeSetFixedStep(mem, h);
  flag += FIRKodeSetNonlinConvCoef(mem, SUN_RCONST(1.0e-3));
  flag += FIRKodeSetMaxNonlinIters(mem, 20);
  flag += FIRKodeSetMaxNumSteps(mem, 100000);

  SUNMatrix A;
  SUNLinearSolver LS;
  if (p.band)
  {
    A  = SUNBandMatrix(p.neq, 0, 0, sunctx);
    LS = SUNLinSol_Band(y, A, sunctx);
  }
  else
  {
    A  = SUNDenseMatrix(p.neq, p.neq, sunctx);
    LS = SUNLinSol_Dense(y, A, sunctx);
  }
  flag += FIRKodeSetLinearSolver(mem, LS, A);
  flag += FIRKodeSetJacFn(mem, p.jac);
  flag += FIRKodeSetStageEpsLin(mem, SUN_RCONST(1.0e-3));
  if (flag != 0)
  {
    std::fprintf(stderr, "  setup failed\n");
    return 1;
  }

  sunrealtype t = p.t0;
  flag          = FIRKodeEvolve(mem, p.tf, y, &t, FIRK_NORMAL);
  if (flag < 0)
  {
    std::fprintf(stderr, "  FIRKodeEvolve failed with flag %d\n", flag);
    return 1;
  }

  p.exact(t, yref, p.udata);
  N_VLinearSum(one, y, -one, yref, yref);
  err = N_VMaxNorm(yref);

  long int nst, nni;
  FIRKodeGetNumSteps(mem, &nst);
  FIRKodeGetNumNonlinSolvIters(mem, &nni);
  std::printf("  h = %8.4" ESYM
              "  steps = %6li  nni/step = %5.2f  error = %.4" ESYM "\n",
              h, nst, static_cast<double>(nni) / static_cast<double>(nst), err);

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  N_VDestroy(yref);
  return 0;
}

} // namespace

int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    std::fprintf(stderr, "usage: %s <pr|kpr> <stages> [lambda]\n", argv[0]);
    return 1;
  }
  const std::string which = argv[1];
  const int s             = std::atoi(argv[2]);

  sunrealtype lambda      = SUN_RCONST(-1.0);
  sunrealtype kpr_data[4] = {SUN_RCONST(-2.0), SUN_RCONST(0.5), SUN_RCONST(0.5),
                             SUN_RCONST(-1.0)};

  Problem p;
  sunrealtype expected;
  if (which == "pr")
  {
    if (argc > 3) { lambda = static_cast<sunrealtype>(std::atof(argv[3])); }
    p = {"Prothero-Robinson",
         4,
         zero,
         one,
         SUN_RCONST(0.2),
         problems::pr::ode_rhs,
         problems::pr::ode_rhs_jac,
         &lambda,
         true,
         pr_exact};
    // nonstiff: classical order 2s-1; stiff: order reduces towards the stage
    // order s
    expected = (lambda > SUN_RCONST(-100.0))
                 ? static_cast<sunrealtype>(2 * s - 1)
                 : static_cast<sunrealtype>(s);
  }
  else if (which == "kpr")
  {
    p        = {"Kvaerno-Prothero-Robinson",
                2,
                zero,
                SUN_RCONST(0.5),
                SUN_RCONST(0.02),
                problems::kpr::ode_rhs,
                problems::kpr::ode_rhs_jac,
                kpr_data,
                false,
                kpr_exact};
    expected = static_cast<sunrealtype>(2 * s - 1);
  }
  else
  {
    std::fprintf(stderr, "unknown problem %s\n", which.c_str());
    return 1;
  }

  SUNContext sunctx;
  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  std::printf("%s, Radau IIA with s = %d stages", p.name.c_str(), s);
  if (which == "pr")
  {
    std::printf(", lambda = %g", static_cast<double>(lambda));
  }
  std::printf("\n");

  std::vector<sunrealtype> errs;
  sunrealtype h = p.h0;
  for (int j = 0; j < num_refinements; j++)
  {
    sunrealtype err;
    if (run(sunctx, p, s, h, err)) { return 1; }
    errs.push_back(err);
    h /= SUN_RCONST(2.0);
  }

  // observed orders from consecutive halvings; the asymptotic order is judged
  // on the finest pair whose errors are above the roundoff level
  int fails             = 0;
  sunrealtype order_fin = SUN_RCONST(0.0);
  int num_orders        = 0;
  for (int j = 1; j < num_refinements; j++)
  {
    if (errs[j] < SUN_RCONST(1.0e4) * SUN_UNIT_ROUNDOFF) { break; }
    const sunrealtype order = std::log(errs[j - 1] / errs[j]) /
                              std::log(SUN_RCONST(2.0));
    std::printf("  observed order = %5.2f\n", static_cast<double>(order));
    order_fin = order;
    num_orders++;
  }

  if (num_orders == 0)
  {
    std::printf("  errors at roundoff level; order not measurable\n");
  }
  else if (order_fin < expected - SUN_RCONST(0.3))
  {
    std::printf("  FAIL: observed order %.2f < expected %.2f - 0.3\n",
                static_cast<double>(order_fin), static_cast<double>(expected));
    fails++;
  }
  else
  {
    std::printf("  observed order %.2f >= expected %.2f - 0.3\n",
                static_cast<double>(order_fin), static_cast<double>(expected));
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
