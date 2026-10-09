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
 * 1-D heat equation u_t = k u_xx + f(t,x) with P1 finite elements on
 * N nodes: M u' = -k K u + M f(t), with a constant tridiagonal mass
 * matrix (port of examples/firkode/C_serial/firk_heat1D_fem_mass.c
 * with a run-time mesh size).
 * -----------------------------------------------------------------*/

#include <cmath>

#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_band.h>

#include "problem.hpp"

namespace benchmark {

namespace {

constexpr sunrealtype ZERO = SUN_RCONST(0.0);
constexpr sunrealtype ONE  = SUN_RCONST(1.0);
constexpr sunrealtype TWO  = SUN_RCONST(2.0);
constexpr sunrealtype FOUR = SUN_RCONST(4.0);
constexpr sunrealtype SIX  = SUN_RCONST(6.0);
constexpr sunrealtype PI   = SUN_RCONST(3.14159265358979323846);

struct Heat1DData
{
  sunindextype N = 101;
  sunrealtype dx = ONE / SUN_RCONST(100.0);
  sunrealtype k  = ONE;
};

sunrealtype exact_pde(sunrealtype t, sunrealtype x)
{ return std::exp(-t) * std::sin(PI * x); }

int heat1d_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Heat1DData* d        = static_cast<Heat1DData*>(user_data);
  const sunindextype N = d->N;
  const sunrealtype k = d->k, dx = d->dx;
  const sunrealtype* Y = N_VGetArrayPointer(y);
  sunrealtype* Ydot    = N_VGetArrayPointer(ydot);
  const sunrealtype c  = (k * PI * PI - ONE) * std::exp(-t);
  Ydot[0]              = ZERO;
  Ydot[N - 1]          = ZERO;
  for (sunindextype i = 1; i < N - 1; i++)
  {
    const sunrealtype fl = c * std::sin(PI * dx * (sunrealtype)(i - 1));
    const sunrealtype fc = c * std::sin(PI * dx * (sunrealtype)i);
    const sunrealtype fr = c * std::sin(PI * dx * (sunrealtype)(i + 1));
    Ydot[i]              = -k * (-Y[i - 1] + TWO * Y[i] - Y[i + 1]) / dx +
                           dx * (fl + FOUR * fc + fr) / SIX;
  }
  return 0;
}

int heat1d_jac(sunrealtype, N_Vector, N_Vector, SUNMatrix J, void* user_data,
               N_Vector, N_Vector, N_Vector)
{
  Heat1DData* d        = static_cast<Heat1DData*>(user_data);
  const sunindextype N = d->N;
  const sunrealtype k = d->k, dx = d->dx;
  SUNMatZero(J);
  for (sunindextype i = 1; i < N - 1; i++)
  {
    SM_ELEMENT_B(J, i, i - 1) = k / dx;
    SM_ELEMENT_B(J, i, i)     = -TWO * k / dx;
    SM_ELEMENT_B(J, i, i + 1) = k / dx;
  }
  return 0;
}

int heat1d_mass(sunrealtype, SUNMatrix M, void* user_data, N_Vector, N_Vector,
                N_Vector)
{
  Heat1DData* d        = static_cast<Heat1DData*>(user_data);
  const sunindextype N = d->N;
  const sunrealtype dx = d->dx;
  SUNMatZero(M);
  SM_ELEMENT_B(M, 0, 0)         = ONE;
  SM_ELEMENT_B(M, N - 1, N - 1) = ONE;
  for (sunindextype i = 1; i < N - 1; i++)
  {
    SM_ELEMENT_B(M, i, i - 1) = dx / SIX;
    SM_ELEMENT_B(M, i, i)     = TWO * dx / SUN_RCONST(3.0);
    SM_ELEMENT_B(M, i, i + 1) = dx / SIX;
  }
  return 0;
}

class Heat1DFem : public Problem
{
public:
  Heat1DFem()
  {
    spec_.name         = "heat1dfem";
    spec_.t0           = ZERO;
    spec_.tf           = ONE;
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-10);
    spec_.ls_default   = LinSolKind::Band;
    spec_.mu = spec_.ml = 1;
    spec_.mass_mu = spec_.mass_ml = 1;
    spec_.has_jac                 = true;
    spec_.has_mass                = true;
    spec_.linear                  = true;
    spec_.symmetric               = true;
    set_size(101);
  }

  void configure(const Options& o) override
  {
    if (o.n > 1) { set_size((sunindextype)o.n); }
  }

  void init(N_Vector y) const override
  {
    sunrealtype* Y = N_VGetArrayPointer(y);
    for (sunindextype i = 0; i < data_.N; i++)
    {
      Y[i] = exact_pde(ZERO, data_.dx * (sunrealtype)i);
    }
  }

  void* user_data() override { return &data_; }

  FIRKRhsFn rhs() const override { return heat1d_rhs; }

  FIRKLsJacFn jac() const override { return heat1d_jac; }

  FIRKLsMassFn mass() const override { return heat1d_mass; }

private:
  void set_size(sunindextype N)
  {
    data_.N  = N;
    data_.dx = ONE / (sunrealtype)(N - 1);
    spec_.N  = N;
  }

  Heat1DData data_;
};

} // namespace

std::unique_ptr<Problem> make_heat1dfem()
{ return std::unique_ptr<Problem>(new Heat1DFem()); }

} // namespace benchmark
