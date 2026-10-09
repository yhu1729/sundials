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
 * Small stiff test problems of the FIRKODE benchmark harness:
 *   rober  Robertson kinetics (Hairer-Wanner IV.10), N = 3, tf = 1e11
 *   vdpol  Van der Pol in the Hairer-Wanner form, eps = 1e-6, tf = 2
 *   orego  Oregonator, N = 3, tf = 360
 *   hires  HIRES (Schaefer 1975), N = 8, tf = 321.8122
 *   analytic  y' = lambda*y + 1/(1+t^2) - lambda*atan(t), N = 1
 *   pr     Prothero-Robinson, test/unit_tests/problems/pr.hpp
 *   kpr    Kvaerno-Prothero-Robinson, test/unit_tests/problems/kpr.hpp
 *
 * This is the only translation unit that includes kpr.hpp, which
 * defines a non-inline global.
 * -----------------------------------------------------------------*/

#include <cmath>
#include <stdexcept>

#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>

#include "problem.hpp"
#include "problems/kpr.hpp"
#include "problems/pr.hpp"

namespace benchmark {

namespace {

constexpr sunrealtype ZERO = SUN_RCONST(0.0);
constexpr sunrealtype ONE  = SUN_RCONST(1.0);
constexpr sunrealtype TWO  = SUN_RCONST(2.0);

/* ------------------------------------------------------------------
 * Robertson
 * ------------------------------------------------------------------*/

int rober_rhs(sunrealtype, N_Vector y, N_Vector ydot, void*)
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

int rober_jac(sunrealtype, N_Vector y, N_Vector, SUNMatrix J, void*, N_Vector,
              N_Vector, N_Vector)
{
  const sunrealtype v = NV_Ith_S(y, 1);
  const sunrealtype w = NV_Ith_S(y, 2);
  // dense, band or sparse (KLU) matrix: the full 3 x 3 block is written
  const sunrealtype a[9] = {-SUN_RCONST(0.04),
                            SUN_RCONST(1.0e4) * w,
                            SUN_RCONST(1.0e4) * v,
                            SUN_RCONST(0.04),
                            -SUN_RCONST(1.0e4) * w - SUN_RCONST(6.0e7) * v,
                            -SUN_RCONST(1.0e4) * v,
                            ZERO,
                            SUN_RCONST(6.0e7) * v,
                            ZERO};
  return write_full_jacobian(J, a, 3);
}

int rober_roots(sunrealtype, N_Vector y, sunrealtype* gout, void*)
{
  gout[0] = NV_Ith_S(y, 0) - SUN_RCONST(1.0e-4);
  gout[1] = NV_Ith_S(y, 2) - SUN_RCONST(1.0e-2);
  return 0;
}

class Rober : public Problem
{
public:
  Rober()
  {
    spec_.name         = "rober";
    spec_.N            = 3;
    spec_.t0           = ZERO;
    spec_.tf           = SUN_RCONST(1.0e11);
    spec_.rtol_default = SUN_RCONST(1.0e-4);
    spec_.atol_default = SUN_RCONST(1.0e-10);
    spec_.ls_default   = LinSolKind::Dense;
    spec_.has_jac      = true;
    spec_.nnz          = 9; // full block for --ls klu
    spec_.nroots       = 2;
  }

  void init(N_Vector y) const override
  {
    NV_Ith_S(y, 0) = ONE;
    NV_Ith_S(y, 1) = ZERO;
    NV_Ith_S(y, 2) = ZERO;
  }

  void* user_data() override { return nullptr; }

  FIRKRhsFn rhs() const override { return rober_rhs; }

  FIRKLsJacFn jac() const override { return rober_jac; }

  FIRKRootFn roots() const override { return rober_roots; }
};

/* ------------------------------------------------------------------
 * Van der Pol (Hairer-Wanner form)
 * ------------------------------------------------------------------*/

int vdpol_rhs(sunrealtype, N_Vector y, N_Vector ydot, void* user_data)
{
  const sunrealtype eps = *static_cast<sunrealtype*>(user_data);
  const sunrealtype y1  = NV_Ith_S(y, 0);
  const sunrealtype y2  = NV_Ith_S(y, 1);
  NV_Ith_S(ydot, 0)     = y2;
  NV_Ith_S(ydot, 1)     = ((ONE - y1 * y1) * y2 - y1) / eps;
  return 0;
}

int vdpol_jac(sunrealtype, N_Vector y, N_Vector, SUNMatrix J, void* user_data,
              N_Vector, N_Vector, N_Vector)
{
  const sunrealtype eps = *static_cast<sunrealtype*>(user_data);
  const sunrealtype y1  = NV_Ith_S(y, 0);
  const sunrealtype y2  = NV_Ith_S(y, 1);
  SM_ELEMENT_D(J, 0, 0) = ZERO;
  SM_ELEMENT_D(J, 0, 1) = ONE;
  SM_ELEMENT_D(J, 1, 0) = (-TWO * y1 * y2 - ONE) / eps;
  SM_ELEMENT_D(J, 1, 1) = (ONE - y1 * y1) / eps;
  return 0;
}

class Vdpol : public Problem
{
public:
  Vdpol()
  {
    spec_.name         = "vdpol";
    spec_.N            = 2;
    spec_.t0           = ZERO;
    spec_.tf           = TWO;
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-6);
    spec_.ls_default   = LinSolKind::Dense;
    spec_.has_jac      = true;
  }

  void configure(const Options& o) override { eps_ = (sunrealtype)o.eps; }

  void init(N_Vector y) const override
  {
    NV_Ith_S(y, 0) = TWO;
    NV_Ith_S(y, 1) = ZERO;
  }

  void* user_data() override { return &eps_; }

  FIRKRhsFn rhs() const override { return vdpol_rhs; }

  FIRKLsJacFn jac() const override { return vdpol_jac; }

private:
  sunrealtype eps_ = SUN_RCONST(1.0e-6);
};

/* ------------------------------------------------------------------
 * Oregonator
 * ------------------------------------------------------------------*/

constexpr sunrealtype OREGO_A = SUN_RCONST(77.27);
constexpr sunrealtype OREGO_B = SUN_RCONST(8.375e-6);
constexpr sunrealtype OREGO_C = SUN_RCONST(0.161);

int orego_rhs(sunrealtype, N_Vector y, N_Vector ydot, void*)
{
  const sunrealtype y1 = NV_Ith_S(y, 0);
  const sunrealtype y2 = NV_Ith_S(y, 1);
  const sunrealtype y3 = NV_Ith_S(y, 2);
  NV_Ith_S(ydot, 0)    = OREGO_A * (y2 + y1 * (ONE - OREGO_B * y1 - y2));
  NV_Ith_S(ydot, 1)    = (y3 - (ONE + y1) * y2) / OREGO_A;
  NV_Ith_S(ydot, 2)    = OREGO_C * (y1 - y3);
  return 0;
}

int orego_jac(sunrealtype, N_Vector y, N_Vector, SUNMatrix J, void*, N_Vector,
              N_Vector, N_Vector)
{
  const sunrealtype y1 = NV_Ith_S(y, 0);
  const sunrealtype y2 = NV_Ith_S(y, 1);
  SUNMatZero(J);
  SM_ELEMENT_D(J, 0, 0) = OREGO_A * (ONE - TWO * OREGO_B * y1 - y2);
  SM_ELEMENT_D(J, 0, 1) = OREGO_A * (ONE - y1);
  SM_ELEMENT_D(J, 1, 0) = -y2 / OREGO_A;
  SM_ELEMENT_D(J, 1, 1) = -(ONE + y1) / OREGO_A;
  SM_ELEMENT_D(J, 1, 2) = ONE / OREGO_A;
  SM_ELEMENT_D(J, 2, 0) = OREGO_C;
  SM_ELEMENT_D(J, 2, 2) = -OREGO_C;
  return 0;
}

class Orego : public Problem
{
public:
  Orego()
  {
    spec_.name         = "orego";
    spec_.N            = 3;
    spec_.t0           = ZERO;
    spec_.tf           = SUN_RCONST(360.0);
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-6);
    spec_.ls_default   = LinSolKind::Dense;
    spec_.has_jac      = true;
  }

  void init(N_Vector y) const override
  {
    NV_Ith_S(y, 0) = ONE;
    NV_Ith_S(y, 1) = TWO;
    NV_Ith_S(y, 2) = SUN_RCONST(3.0);
  }

  void* user_data() override { return nullptr; }

  FIRKRhsFn rhs() const override { return orego_rhs; }

  FIRKLsJacFn jac() const override { return orego_jac; }
};

/* ------------------------------------------------------------------
 * HIRES
 * ------------------------------------------------------------------*/

int hires_rhs(sunrealtype, N_Vector yv, N_Vector ydot, void*)
{
  const sunrealtype* y = N_VGetArrayPointer(yv);
  sunrealtype* f       = N_VGetArrayPointer(ydot);
  f[0]                 = -SUN_RCONST(1.71) * y[0] + SUN_RCONST(0.43) * y[1] +
         SUN_RCONST(8.32) * y[2] + SUN_RCONST(0.0007);
  f[1] = SUN_RCONST(1.71) * y[0] - SUN_RCONST(8.75) * y[1];
  f[2] = -SUN_RCONST(10.03) * y[2] + SUN_RCONST(0.43) * y[3] +
         SUN_RCONST(0.035) * y[4];
  f[3] = SUN_RCONST(8.32) * y[1] + SUN_RCONST(1.71) * y[2] -
         SUN_RCONST(1.12) * y[3];
  f[4] = -SUN_RCONST(1.745) * y[4] + SUN_RCONST(0.43) * y[5] +
         SUN_RCONST(0.43) * y[6];
  f[5] = -SUN_RCONST(280.0) * y[5] * y[7] + SUN_RCONST(0.69) * y[3] +
         SUN_RCONST(1.71) * y[4] - SUN_RCONST(0.43) * y[5] +
         SUN_RCONST(0.69) * y[6];
  f[6] = SUN_RCONST(280.0) * y[5] * y[7] - SUN_RCONST(1.81) * y[6];
  f[7] = -f[6];
  return 0;
}

int hires_jac(sunrealtype, N_Vector yv, N_Vector, SUNMatrix J, void*, N_Vector,
              N_Vector, N_Vector)
{
  const sunrealtype* y = N_VGetArrayPointer(yv);
  SUNMatZero(J);
  SM_ELEMENT_D(J, 0, 0) = -SUN_RCONST(1.71);
  SM_ELEMENT_D(J, 0, 1) = SUN_RCONST(0.43);
  SM_ELEMENT_D(J, 0, 2) = SUN_RCONST(8.32);
  SM_ELEMENT_D(J, 1, 0) = SUN_RCONST(1.71);
  SM_ELEMENT_D(J, 1, 1) = -SUN_RCONST(8.75);
  SM_ELEMENT_D(J, 2, 2) = -SUN_RCONST(10.03);
  SM_ELEMENT_D(J, 2, 3) = SUN_RCONST(0.43);
  SM_ELEMENT_D(J, 2, 4) = SUN_RCONST(0.035);
  SM_ELEMENT_D(J, 3, 1) = SUN_RCONST(8.32);
  SM_ELEMENT_D(J, 3, 2) = SUN_RCONST(1.71);
  SM_ELEMENT_D(J, 3, 3) = -SUN_RCONST(1.12);
  SM_ELEMENT_D(J, 4, 4) = -SUN_RCONST(1.745);
  SM_ELEMENT_D(J, 4, 5) = SUN_RCONST(0.43);
  SM_ELEMENT_D(J, 4, 6) = SUN_RCONST(0.43);
  SM_ELEMENT_D(J, 5, 3) = SUN_RCONST(0.69);
  SM_ELEMENT_D(J, 5, 4) = SUN_RCONST(1.71);
  SM_ELEMENT_D(J, 5, 5) = -SUN_RCONST(280.0) * y[7] - SUN_RCONST(0.43);
  SM_ELEMENT_D(J, 5, 6) = SUN_RCONST(0.69);
  SM_ELEMENT_D(J, 5, 7) = -SUN_RCONST(280.0) * y[5];
  SM_ELEMENT_D(J, 6, 5) = SUN_RCONST(280.0) * y[7];
  SM_ELEMENT_D(J, 6, 6) = -SUN_RCONST(1.81);
  SM_ELEMENT_D(J, 6, 7) = SUN_RCONST(280.0) * y[5];
  SM_ELEMENT_D(J, 7, 5) = -SUN_RCONST(280.0) * y[7];
  SM_ELEMENT_D(J, 7, 6) = SUN_RCONST(1.81);
  SM_ELEMENT_D(J, 7, 7) = -SUN_RCONST(280.0) * y[5];
  return 0;
}

class Hires : public Problem
{
public:
  Hires()
  {
    spec_.name         = "hires";
    spec_.N            = 8;
    spec_.t0           = ZERO;
    spec_.tf           = SUN_RCONST(321.8122);
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-6);
    spec_.ls_default   = LinSolKind::Dense;
    spec_.has_jac      = true;
  }

  void init(N_Vector y) const override
  {
    N_VConst(ZERO, y);
    NV_Ith_S(y, 0) = ONE;
    NV_Ith_S(y, 7) = SUN_RCONST(0.0057);
  }

  void* user_data() override { return nullptr; }

  FIRKRhsFn rhs() const override { return hires_rhs; }

  FIRKLsJacFn jac() const override { return hires_jac; }
};

/* ------------------------------------------------------------------
 * Analytic stiff scalar problem of examples/firkode/C_serial/
 * firk_analytic.c, y' = lambda*y + 1/(1+t^2) - lambda*atan(t), with
 * the example's floating-point expression (the PR problem below is
 * mathematically identical but rounds differently).
 * ------------------------------------------------------------------*/

int analytic_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  const sunrealtype lambda = *static_cast<sunrealtype*>(user_data);
  const sunrealtype u      = NV_Ith_S(y, 0);
  NV_Ith_S(ydot, 0) = lambda * u + ONE / (ONE + t * t) - lambda * std::atan(t);
  return 0;
}

int analytic_jac(sunrealtype, N_Vector, N_Vector, SUNMatrix J, void* user_data,
                 N_Vector, N_Vector, N_Vector)
{
  SM_ELEMENT_D(J, 0, 0) = *static_cast<sunrealtype*>(user_data);
  return 0;
}

int analytic_exact(sunrealtype t, N_Vector y)
{
  NV_Ith_S(y, 0) = std::atan(t);
  return 0;
}

class Analytic : public Problem
{
public:
  Analytic()
  {
    spec_.name         = "analytic";
    spec_.N            = 1;
    spec_.t0           = ZERO;
    spec_.tf           = SUN_RCONST(10.0);
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-10);
    spec_.ls_default   = LinSolKind::Dense;
    spec_.has_jac      = true;
    spec_.has_exact    = true;
    spec_.linear       = true;
    spec_.symmetric    = true;
  }

  void configure(const Options& o) override { lambda_ = (sunrealtype)o.lambda; }

  void init(N_Vector y) const override { analytic_exact(spec_.t0, y); }

  void* user_data() override { return &lambda_; }

  FIRKRhsFn rhs() const override { return analytic_rhs; }

  FIRKLsJacFn jac() const override { return analytic_jac; }

  int exact(sunrealtype t, N_Vector y) const override
  {
    return analytic_exact(t, y);
  }

private:
  sunrealtype lambda_ = SUN_RCONST(-100.0);
};

/* ------------------------------------------------------------------
 * Prothero-Robinson (band Jacobian with mu = ml = 0)
 * ------------------------------------------------------------------*/

class PR : public Problem
{
public:
  PR()
  {
    spec_.name         = "pr";
    spec_.N            = 4;
    spec_.t0           = ZERO;
    spec_.tf           = ONE;
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-10);
    spec_.ls_default   = LinSolKind::Band;
    spec_.mu = spec_.ml = 0;
    spec_.has_jac       = true;
    spec_.has_exact     = true;
    spec_.linear        = true;
    spec_.symmetric     = true;
  }

  void configure(const Options& o) override
  {
    if (o.n > 0) { spec_.N = (sunindextype)o.n; }
    lambda_ = (sunrealtype)o.lambda;
  }

  void init(N_Vector y) const override
  {
    problems::pr::true_solution(spec_.t0, y);
  }

  void* user_data() override { return &lambda_; }

  FIRKRhsFn rhs() const override { return problems::pr::ode_rhs; }

  FIRKLsJacFn jac() const override { return problems::pr::ode_rhs_jac; }

  int exact(sunrealtype t, N_Vector y) const override
  {
    return problems::pr::true_solution(t, y);
  }

private:
  sunrealtype lambda_ = SUN_RCONST(-1.0e4);
};

/* ------------------------------------------------------------------
 * Kvaerno-Prothero-Robinson
 * ------------------------------------------------------------------*/

int kpr_exact(sunrealtype t, N_Vector y)
{
  sunrealtype u, v;
  problems::kpr::true_sol(t, &u, &v);
  NV_Ith_S(y, 0) = u;
  NV_Ith_S(y, 1) = v;
  return 0;
}

class KPR : public Problem
{
public:
  KPR()
  {
    spec_.name         = "kpr";
    spec_.N            = 2;
    spec_.t0           = ZERO;
    spec_.tf           = ONE;
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-10);
    spec_.ls_default   = LinSolKind::Dense;
    spec_.has_jac      = true;
    spec_.has_exact    = true;
  }

  void init(N_Vector y) const override { kpr_exact(spec_.t0, y); }

  void* user_data() override { return problems::kpr::problem_data; }

  FIRKRhsFn rhs() const override { return problems::kpr::ode_rhs; }

  FIRKLsJacFn jac() const override { return problems::kpr::ode_rhs_jac; }

  int exact(sunrealtype t, N_Vector y) const override
  {
    return kpr_exact(t, y);
  }
};

} // namespace

std::unique_ptr<Problem> make_diurnal();
std::unique_ptr<Problem> make_heat2d();
std::unique_ptr<Problem> make_heat1dfem();
std::unique_ptr<Problem> make_sector();
std::unique_ptr<Problem> make_wave1dfem();
std::unique_ptr<Problem> make_pollu();
std::unique_ptr<Problem> make_ringmod();
std::unique_ptr<Problem> make_beam();
std::unique_ptr<Problem> make_plate();
std::unique_ptr<Problem> make_detest_b5();
std::unique_ptr<Problem> make_bruss3d();

std::unique_ptr<Problem> make_problem(const std::string& name)
{
  if (name == "rober") { return std::unique_ptr<Problem>(new Rober()); }
  if (name == "vdpol") { return std::unique_ptr<Problem>(new Vdpol()); }
  if (name == "orego") { return std::unique_ptr<Problem>(new Orego()); }
  if (name == "hires") { return std::unique_ptr<Problem>(new Hires()); }
  if (name == "analytic") { return std::unique_ptr<Problem>(new Analytic()); }
  if (name == "pr") { return std::unique_ptr<Problem>(new PR()); }
  if (name == "kpr") { return std::unique_ptr<Problem>(new KPR()); }
  if (name == "diurnal") { return make_diurnal(); }
  if (name == "heat2d") { return make_heat2d(); }
  if (name == "heat1dfem") { return make_heat1dfem(); }
  if (name == "sector") { return make_sector(); }
  if (name == "wave1dfem") { return make_wave1dfem(); }
  if (name == "pollu") { return make_pollu(); }
  if (name == "ringmod") { return make_ringmod(); }
  if (name == "beam") { return make_beam(); }
  if (name == "plate") { return make_plate(); }
  if (name == "detest_b5") { return make_detest_b5(); }
  if (name == "bruss3d") { return make_bruss3d(); }
  throw std::runtime_error("unknown problem: " + name);
}

} // namespace benchmark
