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
 * 2-D advection-diffusion on the unit square,
 *   u_t = u_xx + u_yy - v u_x + b(t,x,y),  u = 0 on the boundary,
 * second-order centered differences on an mx x my grid of interior
 * points (x fastest), N = mx*my. The forcing b is chosen so that
 * u*(t,x,y) = exp(-t) sin(pi x) sin(pi y) solves the PDE; u* is used
 * for the initial condition only. The semi-discrete system is linear
 * with a constant banded Jacobian (mu = ml = mx), symmetric when
 * v = 0. Matrix-free variants: Jacobian-vector product plus a Jacobi
 * or inner-PCG preconditioner solve on I - gamma J.
 * -----------------------------------------------------------------*/

#include <cmath>
#include <stdexcept>
#include <vector>

#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_band.h>

#include "problem.hpp"

namespace benchmark {

namespace {

constexpr sunrealtype ZERO = SUN_RCONST(0.0);
constexpr sunrealtype ONE  = SUN_RCONST(1.0);
constexpr sunrealtype TWO  = SUN_RCONST(2.0);
constexpr sunrealtype PI   = SUN_RCONST(3.14159265358979323846);

struct Heat2DData
{
  int mx = 32, my = 32;
  sunrealtype dx = ZERO, dy = ZERO, v = ZERO;
  sunrealtype cx = ZERO, cy = ZERO, ca = ZERO; // 1/dx^2, 1/dy^2, v/(2dx)
  // matrix-free preconditioner
  bool inner_pcg       = false;
  sunrealtype pcg_rtol = SUN_RCONST(1.0e-2);
  int pcg_maxit        = 200;
  std::vector<sunrealtype> p, q, r; // PCG work arrays
  long inner_iters = 0;

  sunindextype idx(int i, int j) const { return i + (sunindextype)j * mx; }

  void set(int mx_, int my_, sunrealtype v_)
  {
    mx = mx_;
    my = my_;
    v  = v_;
    dx = ONE / (sunrealtype)(mx + 1);
    dy = ONE / (sunrealtype)(my + 1);
    cx = ONE / (dx * dx);
    cy = ONE / (dy * dy);
    ca = v / (TWO * dx);
  }

  // y = J u (stencil without forcing)
  void apply_J(const sunrealtype* u, sunrealtype* out) const
  {
    for (int j = 0; j < my; j++)
    {
      for (int i = 0; i < mx; i++)
      {
        const sunrealtype uc = u[idx(i, j)];
        const sunrealtype ul = (i > 0) ? u[idx(i - 1, j)] : ZERO;
        const sunrealtype ur = (i < mx - 1) ? u[idx(i + 1, j)] : ZERO;
        const sunrealtype ud = (j > 0) ? u[idx(i, j - 1)] : ZERO;
        const sunrealtype uu = (j < my - 1) ? u[idx(i, j + 1)] : ZERO;
        out[idx(i, j)] = cx * (ul - TWO * uc + ur) + cy * (ud - TWO * uc + uu) -
                         ca * (ur - ul);
      }
    }
  }

  sunrealtype diag() const { return -TWO * cx - TWO * cy; }
};

sunrealtype exact_pde(sunrealtype t, sunrealtype x, sunrealtype y)
{ return std::exp(-t) * std::sin(PI * x) * std::sin(PI * y); }

int heat2d_rhs(sunrealtype t, N_Vector uv, N_Vector udot, void* user_data)
{
  Heat2DData* d        = static_cast<Heat2DData*>(user_data);
  const sunrealtype* u = N_VGetArrayPointer(uv);
  sunrealtype* f       = N_VGetArrayPointer(udot);
  d->apply_J(u, f);
  // forcing b = u*_t - (Lu* - v u*_x) with u* = e^{-t} sin(pi x) sin(pi y)
  const sunrealtype et = std::exp(-t);
  for (int j = 0; j < d->my; j++)
  {
    const sunrealtype y  = (j + 1) * d->dy;
    const sunrealtype sy = std::sin(PI * y);
    for (int i = 0; i < d->mx; i++)
    {
      const sunrealtype x  = (i + 1) * d->dx;
      const sunrealtype sx = std::sin(PI * x);
      const sunrealtype us = et * sx * sy;
      const sunrealtype b  = -us + TWO * PI * PI * us +
                             d->v * PI * et * std::cos(PI * x) * sy;
      f[d->idx(i, j)] += b;
    }
  }
  return 0;
}

int heat2d_jac(sunrealtype, N_Vector, N_Vector, SUNMatrix J, void* user_data,
               N_Vector, N_Vector, N_Vector)
{
  Heat2DData* d = static_cast<Heat2DData*>(user_data);
  SUNMatZero(J);
  for (int j = 0; j < d->my; j++)
  {
    for (int i = 0; i < d->mx; i++)
    {
      const sunindextype k  = d->idx(i, j);
      SM_ELEMENT_B(J, k, k) = d->diag();
      if (i > 0) { SM_ELEMENT_B(J, k, k - 1) = d->cx + d->ca; }
      if (i < d->mx - 1) { SM_ELEMENT_B(J, k, k + 1) = d->cx - d->ca; }
      if (j > 0) { SM_ELEMENT_B(J, k, k - d->mx) = d->cy; }
      if (j < d->my - 1) { SM_ELEMENT_B(J, k, k + d->mx) = d->cy; }
    }
  }
  return 0;
}

int heat2d_jtv(N_Vector v, N_Vector Jv, sunrealtype, N_Vector, N_Vector,
               void* user_data, N_Vector)
{
  Heat2DData* d = static_cast<Heat2DData*>(user_data);
  d->apply_J(N_VGetArrayPointer(v), N_VGetArrayPointer(Jv));
  return 0;
}

int heat2d_psetup(sunrealtype, N_Vector, N_Vector, sunbooleantype,
                  sunbooleantype* jcurPtr, sunrealtype, void*)
{
  *jcurPtr = SUNFALSE; // the Jacobian is constant; nothing to set up
  return 0;
}

// Solve (I - gamma J) z = r approximately: Jacobi, or Jacobi-preconditioned
// CG to a relative residual pcg_rtol (symmetric case only).
int heat2d_psolve(sunrealtype, N_Vector, N_Vector, N_Vector rv, N_Vector zv,
                  sunrealtype gamma, sunrealtype, int, void* user_data)
{
  Heat2DData* d          = static_cast<Heat2DData*>(user_data);
  const sunrealtype* r   = N_VGetArrayPointer(rv);
  sunrealtype* z         = N_VGetArrayPointer(zv);
  const sunindextype N   = (sunindextype)d->mx * d->my;
  const sunrealtype dinv = ONE / (ONE - gamma * d->diag());

  if (!d->inner_pcg)
  {
    for (sunindextype k = 0; k < N; k++) { z[k] = dinv * r[k]; }
    return 0;
  }

  d->p.assign((size_t)N, ZERO);
  d->q.assign((size_t)N, ZERO);
  d->r.assign((size_t)N, ZERO);
  sunrealtype* p     = d->p.data();
  sunrealtype* q     = d->q.data();
  sunrealtype* rr    = d->r.data();
  sunrealtype rnorm0 = ZERO, rz = ZERO;
  for (sunindextype k = 0; k < N; k++)
  {
    z[k]  = ZERO;
    rr[k] = r[k];
    p[k]  = dinv * rr[k];
    rz += rr[k] * p[k];
    rnorm0 += rr[k] * rr[k];
  }
  rnorm0 = std::sqrt(rnorm0);
  if (rnorm0 == ZERO) { return 0; }
  for (int it = 0; it < d->pcg_maxit; it++)
  {
    d->apply_J(p, q); // q = J p
    sunrealtype pq = ZERO;
    for (sunindextype k = 0; k < N; k++)
    {
      q[k] = p[k] - gamma * q[k]; // (I - gamma J) p
      pq += p[k] * q[k];
    }
    const sunrealtype alpha = rz / pq;
    sunrealtype rnorm = ZERO, rz_new = ZERO;
    for (sunindextype k = 0; k < N; k++)
    {
      z[k] += alpha * p[k];
      rr[k] -= alpha * q[k];
      rnorm += rr[k] * rr[k];
    }
    d->inner_iters++;
    if (std::sqrt(rnorm) <= d->pcg_rtol * rnorm0) { break; }
    for (sunindextype k = 0; k < N; k++) { rz_new += rr[k] * dinv * rr[k]; }
    const sunrealtype beta = rz_new / rz;
    rz                     = rz_new;
    for (sunindextype k = 0; k < N; k++) { p[k] = dinv * rr[k] + beta * p[k]; }
  }
  return 0;
}

class Heat2D : public Problem
{
public:
  Heat2D()
  {
    spec_.name         = "heat2d";
    spec_.t0           = ZERO;
    spec_.tf           = ONE;
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-10);
    spec_.ls_default   = LinSolKind::Band;
    spec_.has_jac      = true;
    spec_.has_jtimes   = true;
    spec_.has_prec     = true;
    spec_.linear       = true;
    set_grid(32, 32, ZERO);
  }

  void configure(const Options& o) override
  {
    const int mx = o.mx > 0 ? o.mx : 32;
    const int my = o.my > 0 ? o.my : mx;
    set_grid(mx, my, (sunrealtype)o.velocity);
  }

  void prepare(const Options& o) override
  {
    data_.inner_pcg = (o.psolve == "innerpcg");
    if (data_.inner_pcg && !spec_.symmetric)
    {
      throw std::runtime_error("innerpcg psolve requires velocity 0");
    }
    data_.pcg_rtol = (sunrealtype)o.psolve_rtol;
  }

  void init(N_Vector u) const override
  {
    sunrealtype* ud = N_VGetArrayPointer(u);
    for (int j = 0; j < data_.my; j++)
    {
      for (int i = 0; i < data_.mx; i++)
      {
        ud[data_.idx(i, j)] = exact_pde(ZERO, (i + 1) * data_.dx,
                                        (j + 1) * data_.dy);
      }
    }
  }

  void* user_data() override { return &data_; }

  FIRKRhsFn rhs() const override { return heat2d_rhs; }

  FIRKLsJacFn jac() const override { return heat2d_jac; }

  FIRKLsJacTimesVecFn jtimes() const override { return heat2d_jtv; }

  FIRKLsPrecSetupFn psetup() const override { return heat2d_psetup; }

  FIRKLsPrecSolveFn psolve() const override { return heat2d_psolve; }

private:
  void set_grid(int mx, int my, sunrealtype v)
  {
    data_.set(mx, my, v);
    spec_.N         = (sunindextype)mx * my;
    spec_.mu        = mx;
    spec_.ml        = mx;
    spec_.symmetric = (v == ZERO);
  }

  Heat2DData data_;
};

} // namespace

std::unique_ptr<Problem> make_heat2d()
{ return std::unique_ptr<Problem>(new Heat2D()); }

} // namespace benchmark
