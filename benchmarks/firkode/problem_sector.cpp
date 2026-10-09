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
 * SECTOR: a linear test problem whose stiff eigenvalues sit at a
 * prescribed angle theta from the negative real axis, with an exact
 * solution (benchmark family B8, and B11's rc_pwm variant).
 *
 *   y' = A (y - phi(t)) + phi'(t) + nonlin (y - phi) o (y - phi),
 *   y(0) = phi(0) + delta w,          t in [0, tf]
 *   A = Q blockdiag(B_1, ..., B_m) Q,  N = 2 m,
 *   B_k = r_k [[-cos theta, sin theta], [-sin theta, -cos theta]],
 *   eigenvalues r_k (-cos theta +- i sin theta), r_k log-spaced in
 *   [rmin, rmax], Q = I (decoupled, band Jacobian with mu = ml = 1) or a
 *   Householder reflection I - 2 q q^T (coupled, dense Jacobian).
 *
 * Forcing: "sine", phi_j = sin(omega_j t + p_j) with omega_j in
 * [0.5, 3) and fixed phases, or "square", phi = u(t) c with u a +-1
 * square wave of the given period, left-continuous (u(t) = u_k on
 * (k T, (k+1) T]) so that a stage at the end of an interval sees that
 * interval's value, and c a fixed vector in [0.5, 1.5). At the time of a
 * re-initialization (--restart-every) the forcing takes its right limit,
 * so that the derivative an integrator evaluates at its new initial time
 * (IDA's consistent y', CVODE's first step) belongs to the interval that
 * starts there.
 *
 * Exact solution: with the sine forcing y = phi + Q E(t) Q delta w,
 * E = blockdiag(exp(B_k t)), exact when nonlin = 0 or delta = 0 (then
 * y = phi); with the square forcing the modal recursion over the
 * switching intervals, exact when nonlin = 0.
 *
 * Parameters: --m, --theta (degrees), --rmin, --rmax, --coupled 0|1,
 * --delta, --nonlin, --forcing sine|square, --period.
 * -----------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>

#include "problem.hpp"

namespace benchmark {

namespace {

constexpr sunrealtype ZERO = SUN_RCONST(0.0);
constexpr sunrealtype ONE  = SUN_RCONST(1.0);
constexpr sunrealtype TWO  = SUN_RCONST(2.0);
constexpr sunrealtype PI   = SUN_RCONST(3.14159265358979323846);

// fractional part of j * golden ratio: a fixed, well-spread sequence
sunrealtype spread(sunindextype j)
{
  const double x = (double)(j + 1) * 0.6180339887498949;
  return (sunrealtype)(x - std::floor(x));
}

struct SectorData
{
  sunindextype N = 0;
  long m         = 100;
  sunrealtype c  = ZERO; // cos theta
  sunrealtype s  = ONE;  // sin theta
  std::vector<sunrealtype> r;
  bool coupled = false;
  std::vector<sunrealtype> q;      // Householder unit vector
  std::vector<sunrealtype> Adense; // row-major A = Q D Q (coupled)
  sunrealtype delta  = ZERO;
  sunrealtype nonlin = ZERO;
  std::vector<sunrealtype> w; // unit perturbation of the initial state
  bool square        = false;
  sunrealtype period = SUN_RCONST(0.01);
  // time of the last re-initialization, negative before the first one
  sunrealtype t_restart = -ONE;
  std::vector<sunrealtype> omega, phase, cvec;
  // work vectors
  std::vector<sunrealtype> phi, dphi, u, v1, v2;
};

// u(t) of the square forcing: +1 on (0, T], -1 on (T, 2T], ..., and the
// right limit at the time of the last re-initialization. A time within a
// few ulps of a switch counts as the switch, so that restart times
// computed as t0 + (tf - t0) k / K fall on k T.
sunrealtype square_wave(SectorData& d, sunrealtype t)
{
  // no run goes back before its last re-initialization: an earlier time
  // belongs to a new run (a repetition or the reference run)
  if (t < d.t_restart) { d.t_restart = -ONE; }
  const double x   = (double)(t / d.period);
  const double tol = 64.0 * std::numeric_limits<double>::epsilon() *
                     std::max(1.0, std::fabs(x));
  long k = (t == d.t_restart) ? (long)std::floor(x + tol)
                              : (long)std::ceil(x - tol) - 1;
  if (k < 0) { k = 0; }
  return (k % 2 == 0) ? ONE : -ONE;
}

void forcing(SectorData& d, sunrealtype t, sunrealtype* phi, sunrealtype* dphi)
{
  if (d.square)
  {
    const sunrealtype u = square_wave(d, t);
    for (sunindextype j = 0; j < d.N; j++)
    {
      phi[j]  = u * d.cvec[(size_t)j];
      dphi[j] = ZERO;
    }
  }
  else
  {
    for (sunindextype j = 0; j < d.N; j++)
    {
      const sunrealtype a = d.omega[(size_t)j] * t + d.phase[(size_t)j];
      phi[j]              = std::sin(a);
      dphi[j]             = d.omega[(size_t)j] * std::cos(a);
    }
  }
}

// out = Q v = v - 2 q (q . v)
void apply_q(const SectorData& d, const sunrealtype* v, sunrealtype* out)
{
  sunrealtype dot = ZERO;
  for (sunindextype j = 0; j < d.N; j++) { dot += d.q[(size_t)j] * v[j]; }
  for (sunindextype j = 0; j < d.N; j++)
  {
    out[j] = v[j] - TWO * d.q[(size_t)j] * dot;
  }
}

// out = blockdiag(B_k) v
void apply_d(const SectorData& d, const sunrealtype* v, sunrealtype* out)
{
  for (long k = 0; k < d.m; k++)
  {
    const sunrealtype rk = d.r[(size_t)k];
    const sunrealtype x = v[2 * k], y = v[2 * k + 1];
    out[2 * k]     = rk * (-d.c * x + d.s * y);
    out[2 * k + 1] = rk * (-d.s * x - d.c * y);
  }
}

// out = A v
void apply_a(SectorData& d, const sunrealtype* v, sunrealtype* out)
{
  if (d.coupled)
  {
    apply_q(d, v, d.v1.data());
    apply_d(d, d.v1.data(), d.v2.data());
    apply_q(d, d.v2.data(), out);
  }
  else { apply_d(d, v, out); }
}

// z <- exp(blockdiag(B_k) dt) z
void apply_exp(const SectorData& d, sunrealtype dt, sunrealtype* z)
{
  for (long k = 0; k < d.m; k++)
  {
    const sunrealtype rk = d.r[(size_t)k];
    const sunrealtype e  = std::exp(-rk * d.c * dt);
    const sunrealtype cs = std::cos(rk * d.s * dt), sn = std::sin(rk * d.s * dt);
    const sunrealtype x = z[2 * k], y = z[2 * k + 1];
    z[2 * k]     = e * (cs * x + sn * y);
    z[2 * k + 1] = e * (-sn * x + cs * y);
  }
}

int sector_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  SectorData* d        = static_cast<SectorData*>(user_data);
  const sunrealtype* Y = N_VGetArrayPointer(y);
  sunrealtype* Ydot    = N_VGetArrayPointer(ydot);
  forcing(*d, t, d->phi.data(), d->dphi.data());
  for (sunindextype j = 0; j < d->N; j++)
  {
    d->u[(size_t)j] = Y[j] - d->phi[(size_t)j];
  }
  apply_a(*d, d->u.data(), Ydot);
  for (sunindextype j = 0; j < d->N; j++)
  {
    Ydot[j] += d->dphi[(size_t)j];
    if (d->nonlin != ZERO)
    {
      Ydot[j] += d->nonlin * d->u[(size_t)j] * d->u[(size_t)j];
    }
  }
  return 0;
}

void matrix_set(SUNMatrix J, sunindextype i, sunindextype j, sunrealtype v)
{
  if (SUNMatGetID(J) == SUNMATRIX_DENSE) { SM_ELEMENT_D(J, i, j) = v; }
  else { SM_ELEMENT_B(J, i, j) = v; }
}

int sector_jac(sunrealtype t, N_Vector y, N_Vector, SUNMatrix J,
               void* user_data, N_Vector, N_Vector, N_Vector)
{
  SectorData* d        = static_cast<SectorData*>(user_data);
  const sunrealtype* Y = N_VGetArrayPointer(y);
  const sunindextype N = d->N;
  SUNMatZero(J);
  if (d->coupled)
  {
    for (sunindextype i = 0; i < N; i++)
    {
      for (sunindextype j = 0; j < N; j++)
      {
        matrix_set(J, i, j, d->Adense[(size_t)(i * N + j)]);
      }
    }
  }
  else
  {
    for (long k = 0; k < d->m; k++)
    {
      const sunrealtype rk = d->r[(size_t)k];
      matrix_set(J, 2 * k, 2 * k, -rk * d->c);
      matrix_set(J, 2 * k, 2 * k + 1, rk * d->s);
      matrix_set(J, 2 * k + 1, 2 * k, -rk * d->s);
      matrix_set(J, 2 * k + 1, 2 * k + 1, -rk * d->c);
    }
  }
  if (d->nonlin != ZERO)
  {
    forcing(*d, t, d->phi.data(), d->dphi.data());
    for (sunindextype j = 0; j < N; j++)
    {
      const sunrealtype v = (SUNMatGetID(J) == SUNMATRIX_DENSE)
                              ? SM_ELEMENT_D(J, j, j)
                              : SM_ELEMENT_B(J, j, j);
      matrix_set(J, j, j, v + TWO * d->nonlin * (Y[j] - d->phi[(size_t)j]));
    }
  }
  return 0;
}

class Sector : public Problem
{
public:
  Sector()
  {
    spec_.name         = "sector";
    spec_.t0           = ZERO;
    spec_.tf           = SUN_RCONST(10.0);
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-6);
    spec_.has_jac      = true;
    spec_.linear       = true;
  }

  std::vector<std::string> param_keys() const override
  {
    return {"m",     "theta",  "rmin",    "rmax",  "coupled",
            "delta", "nonlin", "forcing", "period"};
  }

  void configure(const Options& o) override
  {
    d_.m = param(o, "m", 100L);
    if (d_.m < 1) { throw std::runtime_error("sector: --m must be >= 1"); }
    const double theta = param(o, "theta", 80.0) * (double)PI / 180.0;
    d_.c               = (sunrealtype)std::cos(theta);
    d_.s               = (sunrealtype)std::sin(theta);
    const double rmin = param(o, "rmin", 1.0), rmax = param(o, "rmax", 1.0e3);
    d_.coupled                     = param(o, "coupled", 0L) != 0;
    d_.delta                       = (sunrealtype)param(o, "delta", 0.0);
    d_.nonlin                      = (sunrealtype)param(o, "nonlin", 0.0);
    const std::string forcing_name = param(o, "forcing", std::string("sine"));
    if (forcing_name != "sine" && forcing_name != "square")
    {
      throw std::runtime_error("sector: --forcing must be sine or square");
    }
    d_.square = forcing_name == "square";
    d_.period = (sunrealtype)param(o, "period", 0.01);
    if (d_.period <= ZERO)
    {
      throw std::runtime_error("sector: --period must be positive");
    }

    const sunindextype N = (sunindextype)(2 * d_.m);
    d_.N                 = N;
    d_.r.resize((size_t)d_.m);
    for (long k = 0; k < d_.m; k++)
    {
      const double f  = (d_.m > 1) ? (double)k / (double)(d_.m - 1) : 0.0;
      d_.r[(size_t)k] = (sunrealtype)(rmin * std::pow(rmax / rmin, f));
    }
    d_.omega.resize((size_t)N);
    d_.phase.resize((size_t)N);
    d_.cvec.resize((size_t)N);
    d_.w.resize((size_t)N);
    d_.q.resize((size_t)N);
    sunrealtype wn = ZERO, qn = ZERO;
    for (sunindextype j = 0; j < N; j++)
    {
      d_.omega[(size_t)j] = SUN_RCONST(0.5) + SUN_RCONST(2.5) * spread(j);
      d_.phase[(size_t)j] = TWO * PI * spread(j + N);
      d_.cvec[(size_t)j]  = SUN_RCONST(0.5) + spread(j + 2 * N);
      d_.w[(size_t)j]     = std::cos(ONE + (sunrealtype)j);
      d_.q[(size_t)j]     = std::sin(ONE + (sunrealtype)j);
      wn += d_.w[(size_t)j] * d_.w[(size_t)j];
      qn += d_.q[(size_t)j] * d_.q[(size_t)j];
    }
    for (sunindextype j = 0; j < N; j++)
    {
      d_.w[(size_t)j] /= std::sqrt(wn);
      d_.q[(size_t)j] /= std::sqrt(qn);
    }
    d_.phi.assign((size_t)N, ZERO);
    d_.dphi.assign((size_t)N, ZERO);
    d_.u.assign((size_t)N, ZERO);
    d_.v1.assign((size_t)N, ZERO);
    d_.v2.assign((size_t)N, ZERO);
    if (d_.coupled)
    {
      // dense A = Q D Q, column by column
      d_.Adense.assign((size_t)(N * N), ZERO);
      std::vector<sunrealtype> e((size_t)N, ZERO), col((size_t)N, ZERO);
      for (sunindextype j = 0; j < N; j++)
      {
        std::fill(e.begin(), e.end(), ZERO);
        e[(size_t)j] = ONE;
        apply_a(d_, e.data(), col.data());
        for (sunindextype i = 0; i < N; i++)
        {
          d_.Adense[(size_t)(i * N + j)] = col[(size_t)i];
        }
      }
    }

    spec_.N          = N;
    spec_.mu         = d_.coupled ? N - 1 : 1;
    spec_.ml         = spec_.mu;
    spec_.ls_default = d_.coupled ? LinSolKind::Dense : LinSolKind::Band;
    spec_.linear     = d_.nonlin == ZERO;
    spec_.has_exact  = d_.square ? (d_.nonlin == ZERO)
                                 : (d_.nonlin == ZERO || d_.delta == ZERO);
  }

  void init(N_Vector y0) const override
  {
    sunrealtype* Y = N_VGetArrayPointer(y0);
    std::vector<sunrealtype> phi((size_t)d_.N), dphi((size_t)d_.N);
    forcing(d_, ZERO, phi.data(), dphi.data());
    for (sunindextype j = 0; j < d_.N; j++)
    {
      Y[j] = phi[(size_t)j] + d_.delta * d_.w[(size_t)j];
    }
  }

  void* user_data() override { return &d_; }

  void restart_update(sunrealtype t, N_Vector) const override
  {
    d_.t_restart = t;
  }

  FIRKRhsFn rhs() const override { return sector_rhs; }

  FIRKLsJacFn jac() const override { return sector_jac; }

  int exact(sunrealtype t, N_Vector y) const override
  {
    const sunindextype N = d_.N;
    sunrealtype* Y       = N_VGetArrayPointer(y);
    std::vector<sunrealtype> z((size_t)N), tmp((size_t)N), phi((size_t)N),
      dphi((size_t)N);
    if (!d_.square)
    {
      forcing(d_, t, phi.data(), dphi.data());
      for (sunindextype j = 0; j < N; j++) { Y[j] = phi[(size_t)j]; }
      if (d_.delta != ZERO)
      {
        for (sunindextype j = 0; j < N; j++)
        {
          tmp[(size_t)j] = d_.delta * d_.w[(size_t)j];
        }
        if (d_.coupled) { apply_q(d_, tmp.data(), z.data()); }
        else { z = tmp; }
        apply_exp(d_, t, z.data());
        if (d_.coupled) { apply_q(d_, z.data(), tmp.data()); }
        else { tmp = z; }
        for (sunindextype j = 0; j < N; j++) { Y[j] += tmp[(size_t)j]; }
      }
      return 0;
    }
    // square forcing: modal recursion over the switching intervals
    std::vector<sunrealtype> cq((size_t)N);
    if (d_.coupled) { apply_q(d_, d_.cvec.data(), cq.data()); }
    else { cq = d_.cvec; }
    for (sunindextype j = 0; j < N; j++)
    {
      tmp[(size_t)j] = d_.cvec[(size_t)j] + d_.delta * d_.w[(size_t)j];
    }
    if (d_.coupled) { apply_q(d_, tmp.data(), z.data()); }
    else { z = tmp; }
    sunrealtype tk = ZERO;
    for (long k = 0; tk < t; k++)
    {
      const sunrealtype tnext = std::min(t, (sunrealtype)(k + 1) * d_.period);
      const sunrealtype uk    = (k % 2 == 0) ? ONE : -ONE;
      for (sunindextype j = 0; j < N; j++)
      {
        z[(size_t)j] -= uk * cq[(size_t)j];
      }
      apply_exp(d_, tnext - tk, z.data());
      for (sunindextype j = 0; j < N; j++)
      {
        z[(size_t)j] += uk * cq[(size_t)j];
      }
      tk = tnext;
    }
    if (d_.coupled) { apply_q(d_, z.data(), Y); }
    else
    {
      for (sunindextype j = 0; j < N; j++) { Y[j] = z[(size_t)j]; }
    }
    return 0;
  }

private:
  // mutable: restart_update records the re-initialization time
  mutable SectorData d_;
};

} // namespace

std::unique_ptr<Problem> make_sector()
{
  return std::unique_ptr<Problem>(new Sector());
}

} // namespace benchmark
