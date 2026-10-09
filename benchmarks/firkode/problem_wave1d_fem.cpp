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
 * Damped 1-D elastodynamics with P1 finite elements (family B9):
 *
 *   M u'' + C u' + K u = g(t) b   on (0, 1), u = 0 at both ends,
 *   M = (dx/6) tridiag(1, 4, 1), K = (1/dx) tridiag(-1, 2, -1),
 *   n interior nodes, dx = 1/(n+1), first-order form with the
 *   interleaved state y = (u_1, v_1, ..., u_n, v_n):
 *     [I 0; 0 M] y' = [v; -K u - C v + g(t) b]
 *   (Jacobian band mu = ml = 3, mass band mu = ml = 2).
 *
 * Damping C: "mass" alpha M (high modes near the imaginary axis, the
 * test case), "kv" beta K (Kelvin-Voigt, real high-mode eigenvalues,
 * the control), "rayleigh" alpha M + beta K. With --mass lumped the
 * mass matrix is dx I and the system is an ODE (CVODE admissible).
 * Load b = M phi with phi = sin(pi x) + 0.5 sin(2 pi x), g(t) =
 * sin^2(pi t / 2) for t <= 1 and 1 afterwards; u = v = 0 at t = 0.
 *
 * Parameters: --n (interior nodes), --damping mass|kv|rayleigh,
 * --alpha, --beta or --beta-omega (beta times the largest undamped
 * frequency), --mass consistent|lumped, --urow-scale S|dx (scale of the
 * identity rows of the mass matrix, a diagnostic).
 * -----------------------------------------------------------------*/

#include <cmath>
#include <vector>

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

struct WaveData
{
  sunindextype n     = 100; // interior nodes
  sunrealtype dx     = ONE / SUN_RCONST(101.0);
  sunrealtype alpha  = ZERO; // mass-proportional damping
  sunrealtype beta   = ZERO; // stiffness-proportional damping
  bool lumped        = false;
  sunrealtype uscale = ONE;   // scale of the u-rows of M and f (diagnostic)
  std::vector<sunrealtype> b; // load vector
};

sunrealtype g_of_t(sunrealtype t)
{
  if (t >= ONE) { return ONE; }
  const sunrealtype s = std::sin(PI * t / TWO);
  return s * s;
}

// (K w)_i and (M w)_i on the interior nodes, w given with stride 2
sunrealtype apply_k(const WaveData& d, const sunrealtype* y, sunindextype i,
                    int offset)
{
  const sunrealtype wl = (i > 0) ? y[2 * (i - 1) + offset] : ZERO;
  const sunrealtype wc = y[2 * i + offset];
  const sunrealtype wr = (i + 1 < d.n) ? y[2 * (i + 1) + offset] : ZERO;
  return (-wl + TWO * wc - wr) / d.dx;
}

sunrealtype apply_m(const WaveData& d, const sunrealtype* y, sunindextype i,
                    int offset)
{
  if (d.lumped) { return d.dx * y[2 * i + offset]; }
  const sunrealtype wl = (i > 0) ? y[2 * (i - 1) + offset] : ZERO;
  const sunrealtype wc = y[2 * i + offset];
  const sunrealtype wr = (i + 1 < d.n) ? y[2 * (i + 1) + offset] : ZERO;
  return d.dx * (wl + FOUR * wc + wr) / SIX;
}

int wave_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  WaveData* d          = static_cast<WaveData*>(user_data);
  const sunrealtype* Y = N_VGetArrayPointer(y);
  sunrealtype* Ydot    = N_VGetArrayPointer(ydot);
  const sunrealtype g  = g_of_t(t);
  for (sunindextype i = 0; i < d->n; i++)
  {
    Ydot[2 * i]   = d->uscale * Y[2 * i + 1];
    sunrealtype r = -apply_k(*d, Y, i, 0) + g * d->b[(size_t)i];
    if (d->alpha != ZERO) { r -= d->alpha * apply_m(*d, Y, i, 1); }
    if (d->beta != ZERO) { r -= d->beta * apply_k(*d, Y, i, 1); }
    Ydot[2 * i + 1] = d->lumped ? r / d->dx : r;
  }
  return 0;
}

int wave_jac(sunrealtype, N_Vector, N_Vector, SUNMatrix J, void* user_data,
             N_Vector, N_Vector, N_Vector)
{
  WaveData* d          = static_cast<WaveData*>(user_data);
  const sunrealtype dx = d->dx;
  const sunrealtype sc = d->lumped ? ONE / dx : ONE;
  SUNMatZero(J);
  for (sunindextype i = 0; i < d->n; i++)
  {
    SM_ELEMENT_B(J, 2 * i, 2 * i + 1) = d->uscale;
    for (int s = -1; s <= 1; s++)
    {
      const sunindextype j = i + s;
      if (j < 0 || j >= d->n) { continue; }
      const sunrealtype k               = (s == 0 ? TWO : -ONE) / dx;
      const sunrealtype m               = d->lumped ? (s == 0 ? dx : ZERO)
                                                    : dx * (s == 0 ? FOUR : ONE) / SIX;
      SM_ELEMENT_B(J, 2 * i + 1, 2 * j) = -sc * k;
      SM_ELEMENT_B(J, 2 * i + 1, 2 * j + 1) = -sc * (d->alpha * m + d->beta * k);
    }
  }
  return 0;
}

int wave_mass(sunrealtype, SUNMatrix M, void* user_data, N_Vector, N_Vector,
              N_Vector)
{
  WaveData* d          = static_cast<WaveData*>(user_data);
  const sunrealtype dx = d->dx;
  SUNMatZero(M);
  for (sunindextype i = 0; i < d->n; i++)
  {
    SM_ELEMENT_B(M, 2 * i, 2 * i) = d->uscale;
    for (int s = -1; s <= 1; s++)
    {
      const sunindextype j = i + s;
      if (j < 0 || j >= d->n) { continue; }
      SM_ELEMENT_B(M, 2 * i + 1, 2 * j + 1) = dx * (s == 0 ? FOUR : ONE) / SIX;
    }
  }
  return 0;
}

class Wave1DFem : public Problem
{
public:
  Wave1DFem()
  {
    spec_.name         = "wave1dfem";
    spec_.t0           = ZERO;
    spec_.tf           = FOUR;
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-10);
    spec_.ls_default   = LinSolKind::Band;
    spec_.has_jac      = true;
    spec_.linear       = true;
  }

  std::vector<std::string> param_keys() const override
  {
    return {"damping", "alpha", "beta", "beta-omega", "mass", "urow-scale"};
  }

  void configure(const Options& o) override
  {
    d_.n                   = (o.n > 1) ? (sunindextype)o.n : 100;
    d_.dx                  = ONE / (sunrealtype)(d_.n + 1);
    const std::string mass = param(o, "mass", std::string("consistent"));
    if (mass != "consistent" && mass != "lumped")
    {
      throw std::runtime_error(
        "wave1dfem: --mass must be consistent or lumped");
    }
    d_.lumped = mass == "lumped";
    // the u-rows of the consistent mass matrix are the identity; a scale
    // of dx makes every mass entry O(dx) (diagnostic of the mixed-scale
    // behaviour of the mass-matrix integrators)
    const std::string us = param(o, "urow-scale", std::string("1"));
    d_.uscale            = (us == "dx") ? d_.dx : (sunrealtype)std::stod(us);
    if (d_.lumped) { d_.uscale = ONE; }
    const std::string damping = param(o, "damping", std::string("mass"));
    const double alpha        = param(o, "alpha", 0.1);
    // largest undamped frequency of the discretization
    const double omega_max = d_.lumped ? 2.0 / (double)d_.dx
                                       : 2.0 * std::sqrt(3.0) / (double)d_.dx;
    double beta            = param(o, "beta", 0.0);
    if (o.extra.count("beta-omega") > 0)
    {
      beta = param(o, "beta-omega", 0.0) / omega_max;
    }
    if (damping == "mass") { d_.alpha = (sunrealtype)alpha, d_.beta = ZERO; }
    else if (damping == "kv")
    {
      d_.alpha = ZERO;
      d_.beta  = (sunrealtype)(beta > 0.0 ? beta : 10.0 / omega_max);
    }
    else if (damping == "rayleigh")
    {
      d_.alpha = (sunrealtype)alpha;
      d_.beta  = (sunrealtype)(beta > 0.0 ? beta : 0.1 / omega_max);
    }
    else
    {
      throw std::runtime_error(
        "wave1dfem: --damping must be mass, kv or rayleigh");
    }
    // load b = M phi, phi = sin(pi x) + 0.5 sin(2 pi x) at the nodes
    std::vector<sunrealtype> phi((size_t)d_.n);
    for (sunindextype i = 0; i < d_.n; i++)
    {
      const sunrealtype x = d_.dx * (sunrealtype)(i + 1);
      phi[(size_t)i]      = std::sin(PI * x) +
                       SUN_RCONST(0.5) * std::sin(TWO * PI * x);
    }
    d_.b.assign((size_t)d_.n, ZERO);
    for (sunindextype i = 0; i < d_.n; i++)
    {
      if (d_.lumped) { d_.b[(size_t)i] = d_.dx * phi[(size_t)i]; }
      else
      {
        const sunrealtype pl = (i > 0) ? phi[(size_t)(i - 1)] : ZERO;
        const sunrealtype pr = (i + 1 < d_.n) ? phi[(size_t)(i + 1)] : ZERO;
        d_.b[(size_t)i]      = d_.dx * (pl + FOUR * phi[(size_t)i] + pr) / SIX;
      }
    }
    spec_.N  = 2 * d_.n;
    spec_.mu = spec_.ml = 3;
    spec_.mass_mu = spec_.mass_ml = 2;
    spec_.has_mass                = !d_.lumped;
  }

  void init(N_Vector y0) const override { N_VConst(ZERO, y0); }

  void* user_data() override { return &d_; }

  FIRKRhsFn rhs() const override { return wave_rhs; }

  FIRKLsJacFn jac() const override { return wave_jac; }

  FIRKLsMassFn mass() const override { return d_.lumped ? nullptr : wave_mass; }

private:
  WaveData d_;
};

} // namespace

std::unique_ptr<Problem> make_wave1dfem()
{
  return std::unique_ptr<Problem>(new Wave1DFem());
}

} // namespace benchmark
