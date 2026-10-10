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
 * Mass matrix scaling test for FIRKODE: the heat equation
 * u_t = u_xx + (pi^2 - 1) e^{-t} sin(pi x) on [0, 1] with homogeneous
 * Dirichlet conditions, discretized with P1 finite elements on N nodes,
 * M u' = -K u + M f(t). The interior rows of the mass matrix are of order
 * dx = 1/(N - 1) and the boundary rows are those of the identity. The
 * semi-discrete solution is known exactly, u_i(t) = a(t) sin(pi x_i) with
 * a' = -mu a + (pi^2 - 1) e^{-t}, a(0) = 1, where mu is the discrete
 * eigenvalue of the mode sin(pi x_i), so the error of the time integration
 * alone is measured.
 *
 * The residual of the stacked stage system has the units of M u. If it were
 * compared with tolerances in the units of u, the Newton iteration would accept
 * the predictor on fine meshes and the step size would collapse. The test
 * checks that the number of steps does not grow with N and that the error
 * stays at the requested tolerance.
 *
 * Pass criterion (not run in single precision): for rtol = 1e-6 and 1e-8
 * (atol = 1e-4 rtol) the error is at most 10 rtol on N = 401 and N = 1601,
 * and the finer mesh takes at most twice the steps of the coarser.
 *
 * Usage: firk_test_mass_fem [stages]
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "firkode/firkode.h"
#include "firkode/firkode_impl.h" /* firkSetNumStagesAny */
#include "firkode/firkode_ls.h"
#include "nvector/nvector_serial.h"
#include "sunlinsol/sunlinsol_band.h"
#include "sunmatrix/sunmatrix_band.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define ESYM "Le"
#else
#define ESYM "e"
#endif

namespace {

constexpr sunrealtype zero = SUN_RCONST(0.0);
constexpr sunrealtype half = SUN_RCONST(0.5);
constexpr sunrealtype one  = SUN_RCONST(1.0);
constexpr sunrealtype two  = SUN_RCONST(2.0);
constexpr sunrealtype four = SUN_RCONST(4.0);
constexpr sunrealtype six  = SUN_RCONST(6.0);
constexpr sunrealtype pi   = SUN_RCONST(3.141592653589793238462643383279502884);
constexpr sunrealtype tf   = SUN_RCONST(1.0);

struct Mesh
{
  sunindextype N; // number of nodes
  sunrealtype dx; // mesh spacing
};

// f = -K u + M f_src(t) in the interior rows, zero in the boundary rows
int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  const Mesh* mesh     = static_cast<const Mesh*>(user_data);
  const sunindextype N = mesh->N;
  const sunrealtype dx = mesh->dx;
  const sunrealtype* Y = N_VGetArrayPointer(y);
  sunrealtype* Ydot    = N_VGetArrayPointer(ydot);
  const sunrealtype c  = (pi * pi - one) * std::exp(-t);

  Ydot[0]     = zero;
  Ydot[N - 1] = zero;
  for (sunindextype i = 1; i < N - 1; i++)
  {
    const sunrealtype fl = c *
                           std::sin(pi * dx * static_cast<sunrealtype>(i - 1));
    const sunrealtype fc = c * std::sin(pi * dx * static_cast<sunrealtype>(i));
    const sunrealtype fr = c *
                           std::sin(pi * dx * static_cast<sunrealtype>(i + 1));
    Ydot[i] = -(-Y[i - 1] + two * Y[i] - Y[i + 1]) / dx +
              dx * (fl + four * fc + fr) / six;
  }
  return 0;
}

// J = -K (band, mu = ml = 1)
int jac(sunrealtype, N_Vector, N_Vector, SUNMatrix J, void* user_data, N_Vector,
        N_Vector, N_Vector)
{
  const Mesh* mesh = static_cast<const Mesh*>(user_data);
  SUNMatZero(J);
  for (sunindextype i = 1; i < mesh->N - 1; i++)
  {
    SM_ELEMENT_B(J, i, i - 1) = one / mesh->dx;
    SM_ELEMENT_B(J, i, i)     = -two / mesh->dx;
    SM_ELEMENT_B(J, i, i + 1) = one / mesh->dx;
  }
  return 0;
}

// P1 mass matrix with identity boundary rows (band, mu = ml = 1)
int mass(sunrealtype, SUNMatrix M, void* user_data, N_Vector, N_Vector, N_Vector)
{
  const Mesh* mesh     = static_cast<const Mesh*>(user_data);
  const sunindextype N = mesh->N;
  SUNMatZero(M);
  SM_ELEMENT_B(M, 0, 0)         = one;
  SM_ELEMENT_B(M, N - 1, N - 1) = one;
  for (sunindextype i = 1; i < N - 1; i++)
  {
    SM_ELEMENT_B(M, i, i - 1) = mesh->dx / six;
    SM_ELEMENT_B(M, i, i)     = four * mesh->dx / six;
    SM_ELEMENT_B(M, i, i + 1) = mesh->dx / six;
  }
  return 0;
}

// amplitude a(t) of the exact semi-discrete solution: K v = mu M v for
// v_i = sin(pi x_i), mu = 12 sin^2(pi dx / 2) / (dx^2 (2 + cos(pi dx)))
sunrealtype amplitude(sunrealtype t, sunrealtype dx)
{
  const sunrealtype sn = std::sin(half * pi * dx);
  const sunrealtype mu = SUN_RCONST(12.0) * sn * sn /
                         (dx * dx * (two + std::cos(pi * dx)));
  const sunrealtype c = (pi * pi - one) / (mu - one);
  return c * std::exp(-t) + (one - c) * std::exp(-mu * t);
}

struct Result
{
  long int nst, nni, nsli;
  sunrealtype err; // max-norm error at tf
};

int run(SUNContext sunctx, int s, sunindextype N, sunrealtype rtol, Result* res)
{
  Mesh mesh{N, one / static_cast<sunrealtype>(N - 1)};

  N_Vector y      = N_VNew_Serial(N, sunctx);
  sunrealtype* yd = N_VGetArrayPointer(y);
  for (sunindextype i = 0; i < N; i++)
  {
    yd[i] = std::sin(pi * mesh.dx * static_cast<sunrealtype>(i));
  }

  SUNMatrix A         = SUNBandMatrix(N, 1, 1, sunctx);
  SUNLinearSolver LS  = SUNLinSol_Band(y, A, sunctx);
  SUNMatrix M         = SUNBandMatrix(N, 1, 1, sunctx);
  SUNLinearSolver MLS = SUNLinSol_Band(y, M, sunctx);

  void* mem = FIRKodeCreate(sunctx);
  int flag  = FIRKodeInit(mem, rhs, zero, y);
  flag += FIRKodeSetUserData(mem, &mesh);
  flag += firkSetNumStagesAny(mem, s);
  flag += FIRKodeSStolerances(mem, rtol, SUN_RCONST(1.0e-4) * rtol);
  flag += FIRKodeSetLinearSolver(mem, LS, A);
  flag += FIRKodeSetJacFn(mem, jac);
  flag += FIRKodeSetMassLinearSolver(mem, MLS, M, SUNFALSE);
  flag += FIRKodeSetMassFn(mem, mass);
  flag += FIRKodeSetStopTime(mem, tf);
  flag += FIRKodeSetMaxNumSteps(mem, 5000);

  sunrealtype t = zero;
  if (flag == 0) { flag = FIRKodeEvolve(mem, tf, y, &t, FIRK_NORMAL); }

  if (flag == FIRK_SUCCESS || flag == FIRK_TSTOP_RETURN)
  {
    flag                = 0;
    const sunrealtype a = amplitude(t, mesh.dx);
    res->err            = zero;
    for (sunindextype i = 0; i < N; i++)
    {
      const sunrealtype exact =
        a * std::sin(pi * mesh.dx * static_cast<sunrealtype>(i));
      res->err = std::max(res->err, std::abs(yd[i] - exact));
    }
    FIRKodeGetNumSteps(mem, &res->nst);
    FIRKodeGetNumNonlinSolvIters(mem, &res->nni);
    FIRKodeGetNumStageLinIters(mem, &res->nsli);
  }

  FIRKodeFree(&mem);
  SUNLinSolFree(MLS);
  SUNMatDestroy(M);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  return flag;
}

} // namespace

int main(int argc, char* argv[])
{
  const int s = (argc > 1) ? std::atoi(argv[1]) : 3;
  int fails   = 0;

  SUNContext sunctx;
  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  std::printf("Mass matrix scaling test, P1 FEM heat equation with Radau IIA "
              "s = %d\n",
              s);

  const sunindextype meshes[2] = {401, 1601};
  const sunrealtype rtols[2]   = {SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-8)};

  for (const sunrealtype rtol : rtols)
  {
    std::printf("rtol = %.0" ESYM ", atol = %.0" ESYM ":\n", rtol,
                SUN_RCONST(1.0e-4) * rtol);
    Result res[2];
    bool ran = true;
    for (int k = 0; k < 2; k++)
    {
      if (run(sunctx, s, meshes[k], rtol, &res[k]))
      {
        std::printf("  FAIL: integration failed for N = %li\n",
                    static_cast<long int>(meshes[k]));
        fails++;
        ran = false;
        continue;
      }
      // the stacked Krylov iteration count varies across platforms and is
      // not part of the compared output
      std::printf("  N = %4li: nst = %li, nni = %li, max error = %.1" ESYM "\n",
                  static_cast<long int>(meshes[k]), res[k].nst, res[k].nni,
                  res[k].err);

      // the error of a tolerance-controlled run is of the order of rtol |u|,
      // with |u| <= e^{-1} at tf
      if (res[k].err > SUN_RCONST(10.0) * rtol)
      {
        std::printf("  FAIL: error exceeds 10 rtol = %.1" ESYM "\n",
                    SUN_RCONST(10.0) * rtol);
        fails++;
      }
    }
    if (!ran) { continue; }

    // the time scales of the solution do not depend on N, so neither should
    // the number of steps
    if (res[1].nst > 2 * res[0].nst)
    {
      std::printf("  FAIL: %li steps for N = %li against %li for N = %li\n",
                  res[1].nst, static_cast<long int>(meshes[1]), res[0].nst,
                  static_cast<long int>(meshes[0]));
      fails++;
    }
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
