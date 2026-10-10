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
 * Uniform integrator interface of the FIRKODE benchmark harness and the
 * shared linear-solver construction helper.
 * -----------------------------------------------------------------*/

#ifndef FIRKODE_BENCHMARK_INTEGRATOR_HPP_
#define FIRKODE_BENCHMARK_INTEGRATOR_HPP_

#include <memory>
#include <string>
#include <vector>

#include <sundials/sundials_core.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_pcg.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>
#ifdef BENCHMARK_HAVE_KLU
#include <sunlinsol/sunlinsol_klu.h>
#endif

#include "benchmark_common.hpp"
#include "problem.hpp"
#include "timed_linsol.hpp"

namespace benchmark {

// A matrix/linear-solver pair owned by an integrator
struct LinSys
{
  SUNMatrix A        = nullptr;
  SUNLinearSolver LS = nullptr;
  LinSolKind kind    = LinSolKind::Dense;
  std::string prec   = "none";

  void destroy()
  {
    if (LS != nullptr) { SUNLinSolFree(LS); }
    if (A != nullptr) { SUNMatDestroy(A); }
    LS = nullptr;
    A  = nullptr;
  }
};

inline LinSolKind parse_ls(const std::string& s, LinSolKind dflt)
{
  if (s.empty()) { return dflt; }
  if (s == "dense") { return LinSolKind::Dense; }
  if (s == "band") { return LinSolKind::Band; }
  if (s == "klu") { return LinSolKind::KLU; }
  if (s == "spgmr") { return LinSolKind::SPGMR; }
  if (s == "pcg") { return LinSolKind::PCG; }
  throw std::runtime_error("unknown linear solver: " + s);
}

inline int parse_prec(const std::string& s)
{
  if (s == "none") { return SUN_PREC_NONE; }
  if (s == "left") { return SUN_PREC_LEFT; }
  if (s == "right") { return SUN_PREC_RIGHT; }
  throw std::runtime_error("unknown preconditioning type: " + s);
}

// Build the Newton-matrix linear solver for a problem. extra_mu/extra_ml
// widen the band (IDA needs max of the Jacobian and mass bandwidths).
inline LinSys make_linsys(LinSolKind kind, const ProblemSpec& spec, N_Vector y,
                          const Options& opts, SUNContext ctx,
                          sunindextype extra_mu = 0, sunindextype extra_ml = 0)
{
  LinSys ls;
  ls.kind = kind;
  switch (kind)
  {
  case LinSolKind::Dense:
    ls.A = SUNDenseMatrix(spec.N, spec.N, ctx);
    check_ptr(ls.A, "SUNDenseMatrix");
    ls.LS = SUNLinSol_Dense(y, ls.A, ctx);
    check_ptr(ls.LS, "SUNLinSol_Dense");
    break;
  case LinSolKind::Band:
  {
    const sunindextype mu = std::max(spec.mu, extra_mu);
    const sunindextype ml = std::max(spec.ml, extra_ml);
    ls.A                  = SUNBandMatrix(spec.N, mu, ml, ctx);
    check_ptr(ls.A, "SUNBandMatrix");
    ls.LS = SUNLinSol_Band(y, ls.A, ctx);
    check_ptr(ls.LS, "SUNLinSol_Band");
    break;
  }
  case LinSolKind::KLU:
#ifdef BENCHMARK_HAVE_KLU
  {
    if (spec.nnz <= 0)
    {
      throw std::runtime_error(
        "unsupported: the problem has no sparse Jacobian");
    }
    ls.A = SUNSparseMatrix(spec.N, spec.N, spec.nnz, SUN_CSR_MAT, ctx);
    check_ptr(ls.A, "SUNSparseMatrix");
    ls.LS = SUNLinSol_KLU(y, ls.A, ctx);
    check_ptr(ls.LS, "SUNLinSol_KLU");
    break;
  }
#else
    throw std::runtime_error("unsupported: KLU is not built");
#endif
  case LinSolKind::SPGMR:
    ls.prec = spec.has_prec ? opts.prec : "none";
    ls.LS   = SUNLinSol_SPGMR(y, parse_prec(ls.prec), opts.maxl, ctx);
    check_ptr(ls.LS, "SUNLinSol_SPGMR");
    if (opts.maxrs >= 0)
    {
      check(SUNLinSol_SPGMRSetMaxRestarts(ls.LS, opts.maxrs),
            "SUNLinSol_SPGMRSetMaxRestarts");
    }
    break;
  case LinSolKind::PCG:
    if (!spec.symmetric)
    {
      throw std::runtime_error("PCG requires a symmetric Jacobian");
    }
    ls.prec = spec.has_prec ? opts.prec : "none";
    ls.LS   = SUNLinSol_PCG(y, parse_prec(ls.prec), opts.maxl, ctx);
    check_ptr(ls.LS, "SUNLinSol_PCG");
    break;
  }
  if (opts.profile) { ls.LS = wrap_timed_linsol(ls.LS, ctx); }
  return ls;
}

inline bool is_direct(LinSolKind k)
{
  return k == LinSolKind::Dense || k == LinSolKind::Band || k == LinSolKind::KLU;
}

inline LinSys make_mass_linsys(const ProblemSpec& spec, N_Vector y,
                               const Options& opts, SUNContext ctx)
{
  LinSys ls;
  ls.kind = LinSolKind::Band;
  ls.A    = SUNBandMatrix(spec.N, spec.mass_mu, spec.mass_ml, ctx);
  check_ptr(ls.A, "SUNBandMatrix(mass)");
  ls.LS = SUNLinSol_Band(y, ls.A, ctx);
  check_ptr(ls.LS, "SUNLinSol_Band(mass)");
  if (opts.profile) { ls.LS = wrap_timed_linsol(ls.LS, ctx); }
  return ls;
}

// ARKODE's name of a --controller value (ARKodeSetAdaptControllerByName)
inline std::string arkode_controller_name(const std::string& lower)
{
  if (lower == "impgus") { return "ImpGus"; }
  if (lower == "expgus") { return "ExpGus"; }
  if (lower == "imexgus") { return "ImExGus"; }
  if (lower == "i") { return "I"; }
  if (lower == "pi") { return "PI"; }
  if (lower == "pid") { return "PID"; }
  if (lower == "soderlind") { return "Soderlind"; }
  if (lower == "h0211") { return "H0211"; }
  if (lower == "h0321") { return "H0321"; }
  if (lower == "h211") { return "H211"; }
  if (lower == "h312") { return "H312"; }
  throw std::runtime_error("unknown controller: " + lower);
}

class Integrator
{
public:
  virtual ~Integrator() = default;

  // Create the package memory, attach the problem and apply the options.
  // y holds the initial condition and is advanced in place by evolve().
  virtual void setup(Problem& p, const Options& opts, N_Vector y,
                     SUNContext ctx) = 0;

  // Advance to tout (or one internal step). Returns the package flag.
  virtual int evolve(sunrealtype tout, N_Vector y, sunrealtype& tret,
                     bool one_step) = 0;

  virtual bool is_success(int flag) const        = 0; // SUCCESS or TSTOP_RETURN
  virtual bool is_root(int flag) const           = 0;
  virtual void root_info(std::vector<int>& dirs) = 0;
  virtual int get_dky(sunrealtype t, int k, N_Vector dky) = 0;
  virtual int est_local_errors(N_Vector ele)              = 0;
  virtual int err_weights(N_Vector ewt)                   = 0;
  virtual sunrealtype last_step()                         = 0;

  // Order of the last step (CVODE, IDA); -1 when the package has no such
  // getter (FIRKODE's order is fixed, ARKODE has none)
  virtual int last_order() { return -1; }

  // Re-initialize at (t, y) with the same options (counters restart from
  // zero in every package) and set a new stop time
  virtual void reinit(sunrealtype t, N_Vector y) = 0;
  virtual void set_stop_time(sunrealtype tstop)  = 0;

  // Copy counters into s (NA where the package lacks a counter)
  virtual void harvest(RunStats& s) = 0;

  virtual std::string name() const   = 0;
  virtual std::string method() const = 0;
  virtual int order() const          = 0;
  virtual int stages() const         = 0;

  // FIRKODE: "beginner" or "expert"; "NA" for the other packages
  virtual std::string mode_name() const { return "NA"; }

  virtual std::string ls_name() const   = 0;
  virtual std::string prec_name() const = 0;

  // Release all package memory (called before the next repetition)
  virtual void teardown() = 0;
};

// Returns nullptr with a message on stderr when the package is not built
std::unique_ptr<Integrator> make_integrator(const std::string& name);

} // namespace benchmark

#endif // FIRKODE_BENCHMARK_INTEGRATOR_HPP_
