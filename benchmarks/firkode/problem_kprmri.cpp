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
 * The multirate Kvaerno-Prothero-Robinson problem of
 * examples/arkode/C_serial/ark_kpr_mri.c (ported unchanged):
 *
 *   [u]' = [ G  e ] [(-1 + u^2 - r(t)) / (2u)] + [      r'(t) / (2u)        ]
 *   [v]    [ e -1 ] [(-2 + v^2 - s(t)) / (2v)]   [ s'(t) / (2 sqrt(2 + s(t))) ]
 *
 * r(t) = 0.5 cos(t), s(t) = cos(w t), t in [0, 5], exact solution
 * u = sqrt(1 + r(t)), v = sqrt(2 + s(t)). G (default -100) sets the
 * stiffness of the slow scale, w (default 100) the time-scale separation
 * and e (default 0.5) the coupling. The first row is slow, the second
 * fast. Splittings for MRIStep:
 *   mri:          fi = fs (slow implicit, Jacobian Js), ff
 *   mri-imex:     fe = r'/(2u) (slow explicit), fi = the rest of fs
 *                 (Jacobian Jsi), ff
 *   mri-explicit: fe = fs (slow explicit), ff
 * The harness's kpr problem (test/unit_tests/problems/kpr.hpp) is a
 * different system with fixed parameters.
 * -----------------------------------------------------------------*/

#include <cmath>

#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>

#include "problem.hpp"

namespace benchmark {

namespace {

constexpr sunrealtype ZERO = SUN_RCONST(0.0);
constexpr sunrealtype ONE  = SUN_RCONST(1.0);
constexpr sunrealtype TWO  = SUN_RCONST(2.0);

// rpar of the example: G, w, e
struct KprMriData
{
  sunrealtype G = SUN_RCONST(-100.0);
  sunrealtype w = SUN_RCONST(100.0);
  sunrealtype e = SUN_RCONST(0.5);
};

sunrealtype r(sunrealtype t, const KprMriData&)
{
  return SUN_RCONST(0.5) * std::cos(t);
}

sunrealtype s(sunrealtype t, const KprMriData& d) { return std::cos(d.w * t); }

sunrealtype rdot(sunrealtype t, const KprMriData&)
{
  return -SUN_RCONST(0.5) * std::sin(t);
}

sunrealtype sdot(sunrealtype t, const KprMriData& d)
{
  return -d.w * std::sin(d.w * t);
}

sunrealtype utrue(sunrealtype t, const KprMriData& d)
{
  return std::sqrt(ONE + r(t, d));
}

sunrealtype vtrue(sunrealtype t, const KprMriData& d)
{
  return std::sqrt(TWO + s(t, d));
}

// fn of the example: the full right-hand side
int kprmri_fn(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  const KprMriData& d    = *static_cast<KprMriData*>(user_data);
  const sunrealtype u    = NV_Ith_S(y, 0);
  const sunrealtype v    = NV_Ith_S(y, 1);
  const sunrealtype tmp1 = (-ONE + u * u - r(t, d)) / (TWO * u);
  const sunrealtype tmp2 = (-TWO + v * v - s(t, d)) / (TWO * v);
  NV_Ith_S(ydot, 0)      = d.G * tmp1 + d.e * tmp2 + rdot(t, d) / (TWO * u);
  NV_Ith_S(ydot, 1)      = d.e * tmp1 - tmp2 + sdot(t, d) / (TWO * vtrue(t, d));
  return 0;
}

// fs: the slow row
int kprmri_fs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  const KprMriData& d    = *static_cast<KprMriData*>(user_data);
  const sunrealtype u    = NV_Ith_S(y, 0);
  const sunrealtype v    = NV_Ith_S(y, 1);
  const sunrealtype tmp1 = (-ONE + u * u - r(t, d)) / (TWO * u);
  const sunrealtype tmp2 = (-TWO + v * v - s(t, d)) / (TWO * v);
  NV_Ith_S(ydot, 0)      = d.G * tmp1 + d.e * tmp2 + rdot(t, d) / (TWO * u);
  NV_Ith_S(ydot, 1)      = ZERO;
  return 0;
}

// fse: the explicit part of the slow row
int kprmri_fse(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  const KprMriData& d = *static_cast<KprMriData*>(user_data);
  const sunrealtype u = NV_Ith_S(y, 0);
  NV_Ith_S(ydot, 0)   = rdot(t, d) / (TWO * u);
  NV_Ith_S(ydot, 1)   = ZERO;
  return 0;
}

// fsi: the implicit part of the slow row
int kprmri_fsi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  const KprMriData& d    = *static_cast<KprMriData*>(user_data);
  const sunrealtype u    = NV_Ith_S(y, 0);
  const sunrealtype v    = NV_Ith_S(y, 1);
  const sunrealtype tmp1 = (-ONE + u * u - r(t, d)) / (TWO * u);
  const sunrealtype tmp2 = (-TWO + v * v - s(t, d)) / (TWO * v);
  NV_Ith_S(ydot, 0)      = d.G * tmp1 + d.e * tmp2;
  NV_Ith_S(ydot, 1)      = ZERO;
  return 0;
}

// ff: the fast row
int kprmri_ff(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  const KprMriData& d    = *static_cast<KprMriData*>(user_data);
  const sunrealtype u    = NV_Ith_S(y, 0);
  const sunrealtype v    = NV_Ith_S(y, 1);
  const sunrealtype tmp1 = (-ONE + u * u - r(t, d)) / (TWO * u);
  const sunrealtype tmp2 = (-TWO + v * v - s(t, d)) / (TWO * v);
  NV_Ith_S(ydot, 0)      = ZERO;
  NV_Ith_S(ydot, 1)      = d.e * tmp1 - tmp2 + sdot(t, d) / (TWO * vtrue(t, d));
  return 0;
}

// Jn: the Jacobian of fn
int kprmri_jn(sunrealtype t, N_Vector y, N_Vector, SUNMatrix J, void* user_data,
              N_Vector, N_Vector, N_Vector)
{
  const KprMriData& d   = *static_cast<KprMriData*>(user_data);
  const sunrealtype u   = NV_Ith_S(y, 0);
  const sunrealtype v   = NV_Ith_S(y, 1);
  SM_ELEMENT_D(J, 0, 0) = d.G / TWO +
                          (d.G * (ONE + r(t, d)) - rdot(t, d)) / (2 * u * u);
  SM_ELEMENT_D(J, 0, 1) = d.e / TWO + d.e * (TWO + s(t, d)) / (TWO * v * v);
  SM_ELEMENT_D(J, 1, 0) = d.e / TWO + d.e * (ONE + r(t, d)) / (TWO * u * u);
  SM_ELEMENT_D(J, 1, 1) = -ONE / TWO - (TWO + s(t, d)) / (TWO * v * v);
  return 0;
}

// Js: the Jacobian of fs
int kprmri_js(sunrealtype t, N_Vector y, N_Vector, SUNMatrix J, void* user_data,
              N_Vector, N_Vector, N_Vector)
{
  const KprMriData& d   = *static_cast<KprMriData*>(user_data);
  const sunrealtype u   = NV_Ith_S(y, 0);
  const sunrealtype v   = NV_Ith_S(y, 1);
  SM_ELEMENT_D(J, 0, 0) = d.G / TWO +
                          (d.G * (ONE + r(t, d)) - rdot(t, d)) / (2 * u * u);
  SM_ELEMENT_D(J, 0, 1) = d.e / TWO + d.e * (TWO + s(t, d)) / (TWO * v * v);
  SM_ELEMENT_D(J, 1, 0) = ZERO;
  SM_ELEMENT_D(J, 1, 1) = ZERO;
  return 0;
}

// Jsi: the Jacobian of fsi
int kprmri_jsi(sunrealtype t, N_Vector y, N_Vector, SUNMatrix J,
               void* user_data, N_Vector, N_Vector, N_Vector)
{
  const KprMriData& d   = *static_cast<KprMriData*>(user_data);
  const sunrealtype u   = NV_Ith_S(y, 0);
  const sunrealtype v   = NV_Ith_S(y, 1);
  SM_ELEMENT_D(J, 0, 0) = d.G / TWO + (d.G * (ONE + r(t, d))) / (2 * u * u);
  SM_ELEMENT_D(J, 0, 1) = d.e / TWO + d.e * (TWO + s(t, d)) / (TWO * v * v);
  SM_ELEMENT_D(J, 1, 0) = ZERO;
  SM_ELEMENT_D(J, 1, 1) = ZERO;
  return 0;
}

class KprMri : public Problem
{
public:
  KprMri()
  {
    spec_.name         = "kprmri";
    spec_.N            = 2;
    spec_.t0           = ZERO;
    spec_.tf           = SUN_RCONST(5.0);
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-11);
    spec_.ls_default   = LinSolKind::Dense;
    spec_.has_jac      = true;
    spec_.has_exact    = true;
  }

  std::vector<std::string> param_keys() const override
  {
    return {"G", "w", "e"};
  }

  void configure(const Options& o) override
  {
    d_.G = (sunrealtype)param(o, "G", -100.0);
    d_.w = (sunrealtype)param(o, "w", 100.0);
    d_.e = (sunrealtype)param(o, "e", 0.5);
  }

  void init(N_Vector y) const override { exact(spec_.t0, y); }

  void* user_data() override { return &d_; }

  FIRKRhsFn rhs() const override { return kprmri_fn; }

  FIRKLsJacFn jac() const override { return kprmri_jn; }

  int exact(sunrealtype t, N_Vector y) const override
  {
    NV_Ith_S(y, 0) = utrue(t, d_);
    NV_Ith_S(y, 1) = vtrue(t, d_);
    return 0;
  }

  std::vector<RhsSplit> splits() const override
  {
    RhsSplit mri, mri_imex, mri_explicit;
    mri.name          = "mri";
    mri.fi            = kprmri_fs;
    mri.ff            = kprmri_ff;
    mri.jac_i         = kprmri_js;
    mri_imex.name     = "mri-imex";
    mri_imex.fe       = kprmri_fse;
    mri_imex.fi       = kprmri_fsi;
    mri_imex.ff       = kprmri_ff;
    mri_imex.jac_i    = kprmri_jsi;
    mri_explicit.name = "mri-explicit";
    mri_explicit.fe   = kprmri_fs;
    mri_explicit.ff   = kprmri_ff;
    return {mri, mri_imex, mri_explicit};
  }

private:
  KprMriData d_;
};

} // namespace

std::unique_ptr<Problem> make_kprmri()
{
  return std::unique_ptr<Problem>(new KprMri());
}

} // namespace benchmark
