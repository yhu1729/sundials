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
 * Profiling shim of the benchmark harness (--profile): a forwarding
 * SUNLinearSolver that times the setup (factorization) and solve calls
 * of the solver it wraps, and a Problem decorator that times the
 * right-hand-side and Jacobian callbacks. The accumulated times live in
 * a process-wide ProfileStats, reset before every repetition. Profiled
 * runs pay for the timers, so their wall times are never compared with
 * unprofiled ones.
 * -----------------------------------------------------------------*/

#ifndef FIRKODE_BENCHMARK_TIMED_LINSOL_HPP_
#define FIRKODE_BENCHMARK_TIMED_LINSOL_HPP_

#include <chrono>

#include <sundials/sundials_core.h>

#include "benchmark_common.hpp"
#include "problem.hpp"

namespace benchmark {

struct ProfileStats
{
  double t_setup = 0.0, t_solve = 0.0, t_rhs = 0.0, t_jac = 0.0;
  long n_setup = 0, n_solve = 0, n_rhs = 0, n_jac = 0;

  void reset() { *this = ProfileStats(); }
};

inline ProfileStats& profile_stats()
{
  static ProfileStats stats;
  return stats;
}

namespace timed {

using clock = std::chrono::steady_clock;

inline double seconds_since(clock::time_point t0)
{
  return std::chrono::duration<double>(clock::now() - t0).count();
}

// ---- forwarding linear solver ----

struct Content
{
  SUNLinearSolver inner;
  ProfileStats* stats;
};

inline Content* content(SUNLinearSolver S)
{
  return static_cast<Content*>(S->content);
}

inline SUNLinearSolver_Type gettype(SUNLinearSolver S)
{
  return SUNLinSolGetType(content(S)->inner);
}

inline SUNLinearSolver_ID getid(SUNLinearSolver S)
{
  return SUNLinSolGetID(content(S)->inner);
}

inline SUNErrCode initialize(SUNLinearSolver S)
{
  return SUNLinSolInitialize(content(S)->inner);
}

inline int setup(SUNLinearSolver S, SUNMatrix A)
{
  Content* c     = content(S);
  const auto t0  = clock::now();
  const int flag = SUNLinSolSetup(c->inner, A);
  c->stats->t_setup += seconds_since(t0);
  c->stats->n_setup++;
  return flag;
}

inline int solve(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b,
                 sunrealtype tol)
{
  Content* c     = content(S);
  const auto t0  = clock::now();
  const int flag = SUNLinSolSolve(c->inner, A, x, b, tol);
  c->stats->t_solve += seconds_since(t0);
  c->stats->n_solve++;
  return flag;
}

inline SUNErrCode setatimes(SUNLinearSolver S, void* A_data, SUNATimesFn ATimes)
{
  return SUNLinSolSetATimes(content(S)->inner, A_data, ATimes);
}

inline SUNErrCode setpreconditioner(SUNLinearSolver S, void* P_data,
                                    SUNPSetupFn Pset, SUNPSolveFn Psol)
{
  return SUNLinSolSetPreconditioner(content(S)->inner, P_data, Pset, Psol);
}

inline SUNErrCode setscalingvectors(SUNLinearSolver S, N_Vector s1, N_Vector s2)
{
  return SUNLinSolSetScalingVectors(content(S)->inner, s1, s2);
}

inline SUNErrCode setzeroguess(SUNLinearSolver S, sunbooleantype onoff)
{
  return SUNLinSolSetZeroGuess(content(S)->inner, onoff);
}

inline int numiters(SUNLinearSolver S)
{
  return SUNLinSolNumIters(content(S)->inner);
}

inline sunrealtype resnorm(SUNLinearSolver S)
{
  return SUNLinSolResNorm(content(S)->inner);
}

inline N_Vector resid(SUNLinearSolver S)
{
  return SUNLinSolResid(content(S)->inner);
}

inline sunindextype lastflag(SUNLinearSolver S)
{
  return SUNLinSolLastFlag(content(S)->inner);
}

inline SUNErrCode free_solver(SUNLinearSolver S)
{
  Content* c = content(S);
  if (c != nullptr)
  {
    if (c->inner != nullptr) { SUNLinSolFree(c->inner); }
    delete c;
    S->content = nullptr;
  }
  SUNLinSolFreeEmpty(S);
  return SUN_SUCCESS;
}

} // namespace timed

// Wrap a linear solver; the wrapper owns the inner solver and frees it.
// Only the operations the inner solver implements are forwarded, so the
// packages see the same capabilities.
inline SUNLinearSolver wrap_timed_linsol(SUNLinearSolver inner, SUNContext ctx)
{
  SUNLinearSolver S = SUNLinSolNewEmpty(ctx);
  check_ptr(S, "SUNLinSolNewEmpty");
  S->content                   = new timed::Content{inner, &profile_stats()};
  const SUNLinearSolver_Ops in = inner->ops; // a pointer typedef
  S->ops->gettype              = timed::gettype;
  S->ops->getid                = in->getid ? timed::getid : nullptr;
  S->ops->initialize           = in->initialize ? timed::initialize : nullptr;
  S->ops->setup                = in->setup ? timed::setup : nullptr;
  S->ops->solve                = timed::solve;
  S->ops->setatimes            = in->setatimes ? timed::setatimes : nullptr;
  S->ops->setpreconditioner = in->setpreconditioner ? timed::setpreconditioner
                                                    : nullptr;
  S->ops->setscalingvectors = in->setscalingvectors ? timed::setscalingvectors
                                                    : nullptr;
  S->ops->setzeroguess      = in->setzeroguess ? timed::setzeroguess : nullptr;
  S->ops->numiters          = in->numiters ? timed::numiters : nullptr;
  S->ops->resnorm           = in->resnorm ? timed::resnorm : nullptr;
  S->ops->resid             = in->resid ? timed::resid : nullptr;
  S->ops->lastflag          = in->lastflag ? timed::lastflag : nullptr;
  S->ops->free              = timed::free_solver;
  return S;
}

// ---- timed problem callbacks ----

struct TimedCallbacks
{
  Problem* inner = nullptr;
  void* udata    = nullptr;
  ProfileStats* stats;
};

namespace timed {

inline int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* ud)
{
  TimedCallbacks* c = static_cast<TimedCallbacks*>(ud);
  const auto t0     = clock::now();
  const int flag    = c->inner->rhs()(t, y, ydot, c->udata);
  c->stats->t_rhs += seconds_since(t0);
  c->stats->n_rhs++;
  return flag;
}

inline int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* ud,
               N_Vector t1, N_Vector t2, N_Vector t3)
{
  TimedCallbacks* c = static_cast<TimedCallbacks*>(ud);
  const auto t0     = clock::now();
  const int flag    = c->inner->jac()(t, y, fy, J, c->udata, t1, t2, t3);
  c->stats->t_jac += seconds_since(t0);
  c->stats->n_jac++;
  return flag;
}

inline int jtimes(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y,
                  N_Vector fy, void* ud, N_Vector tmp)
{
  TimedCallbacks* c = static_cast<TimedCallbacks*>(ud);
  return c->inner->jtimes()(v, Jv, t, y, fy, c->udata, tmp);
}

inline int psetup(sunrealtype t, N_Vector y, N_Vector fy, sunbooleantype jok,
                  sunbooleantype* jcur, sunrealtype gamma, void* ud)
{
  TimedCallbacks* c = static_cast<TimedCallbacks*>(ud);
  const auto t0     = clock::now();
  const int flag    = c->inner->psetup()(t, y, fy, jok, jcur, gamma, c->udata);
  c->stats->t_setup += seconds_since(t0);
  c->stats->n_setup++;
  return flag;
}

inline int psolve(sunrealtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z,
                  sunrealtype gamma, sunrealtype delta, int lr, void* ud)
{
  TimedCallbacks* c = static_cast<TimedCallbacks*>(ud);
  const auto t0     = clock::now();
  const int flag = c->inner->psolve()(t, y, fy, r, z, gamma, delta, lr, c->udata);
  c->stats->t_solve += seconds_since(t0);
  c->stats->n_solve++;
  return flag;
}

inline int mass(sunrealtype t, SUNMatrix M, void* ud, N_Vector t1, N_Vector t2,
                N_Vector t3)
{
  TimedCallbacks* c = static_cast<TimedCallbacks*>(ud);
  return c->inner->mass()(t, M, c->udata, t1, t2, t3);
}

inline int roots(sunrealtype t, N_Vector y, sunrealtype* gout, void* ud)
{
  TimedCallbacks* c = static_cast<TimedCallbacks*>(ud);
  return c->inner->roots()(t, y, gout, c->udata);
}

} // namespace timed

// Problem decorator whose callbacks time the inner problem's right-hand
// side and Jacobian (preconditioner setup and solve count as linear
// algebra); every other callback forwards with the inner user data.
class TimedProblem : public Problem
{
public:
  explicit TimedProblem(std::unique_ptr<Problem> inner)
    : inner_(std::move(inner))
  {
    spec_     = inner_->spec();
    cb_.inner = inner_.get();
    cb_.udata = inner_->user_data();
    cb_.stats = &profile_stats();
  }

  std::vector<std::string> param_keys() const override
  {
    return inner_->param_keys();
  }

  void configure(const Options& o) override
  {
    inner_->configure(o);
    spec_ = inner_->spec();
  }

  void prepare(const Options& o) override
  {
    inner_->prepare(o);
    spec_     = inner_->spec();
    cb_.udata = inner_->user_data();
  }

  void init(N_Vector y0) const override { inner_->init(y0); }

  void* user_data() override
  {
    cb_.udata = inner_->user_data();
    return &cb_;
  }

  FIRKRhsFn rhs() const override { return timed::rhs; }

  FIRKLsJacFn jac() const override
  {
    return inner_->jac() ? timed::jac : nullptr;
  }

  FIRKLsJacTimesVecFn jtimes() const override
  {
    return inner_->jtimes() ? timed::jtimes : nullptr;
  }

  FIRKLsPrecSetupFn psetup() const override
  {
    return inner_->psetup() ? timed::psetup : nullptr;
  }

  FIRKLsPrecSolveFn psolve() const override
  {
    return inner_->psolve() ? timed::psolve : nullptr;
  }

  FIRKLsMassFn mass() const override
  {
    return inner_->mass() ? timed::mass : nullptr;
  }

  FIRKRootFn roots() const override
  {
    return inner_->roots() ? timed::roots : nullptr;
  }

  int exact(sunrealtype t, N_Vector y) const override
  {
    return inner_->exact(t, y);
  }

  void restart_update(sunrealtype t, N_Vector y) const override
  {
    inner_->restart_update(t, y);
  }

private:
  std::unique_ptr<Problem> inner_;
  TimedCallbacks cb_;
};

} // namespace benchmark

#endif // FIRKODE_BENCHMARK_TIMED_LINSOL_HPP_
