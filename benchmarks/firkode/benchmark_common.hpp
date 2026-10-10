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
 * Shared utilities for the FIRKODE benchmark harness: error checking,
 * timing, run statistics, the CSV output format, error metrics and the
 * reference-solution file format.
 * -----------------------------------------------------------------*/

#ifndef FIRKODE_BENCHMARK_COMMON_HPP_
#define FIRKODE_BENCHMARK_COMMON_HPP_

#include <cfloat>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_core.h>

namespace benchmark {

// Sentinel for counters a package does not provide
constexpr long NA_LONG = -1;

inline double nan_value() { return std::numeric_limits<double>::quiet_NaN(); }

// Throwing error checks for SUNDIALS calls made by the harness (never used
// inside C callbacks)
inline void check(int flag, const char* fn)
{
  if (flag < 0)
  {
    throw std::runtime_error(std::string(fn) + " failed with flag " +
                             std::to_string(flag));
  }
}

inline void check_ptr(const void* p, const char* fn)
{
  if (p == nullptr)
  {
    throw std::runtime_error(std::string(fn) + " returned NULL");
  }
}

// Wall-clock timer
class Timer
{
public:
  using clock = std::chrono::steady_clock;

  void start() { t0_ = clock::now(); }

  double stop() const
  {
    return std::chrono::duration<double>(clock::now() - t0_).count();
  }

private:
  clock::time_point t0_{};
};

// Statistics of one benchmark run. Counters default to NA_LONG and reals
// to NaN so that columns a package lacks print as "NA".
struct RunStats
{
  std::string problem, integrator, method, ls, prec;
  int order        = 0;
  int stages       = 0;
  std::string mode = "NA"; // FIRKODE: beginner|expert
  double rtol = nan_value(), atol = nan_value();
  long N                   = NA_LONG;
  double fixed_h           = nan_value();
  long nst                 = NA_LONG;
  long nst_att             = NA_LONG;
  long nfe                 = NA_LONG;
  long nfe_ls              = NA_LONG;
  long nsetups             = NA_LONG;
  long nje                 = NA_LONG;
  long nni                 = NA_LONG;
  long nncf                = NA_LONG; // nonlinear convergence failures
  long nsf                 = NA_LONG; // step solve failures
  long netf                = NA_LONG;
  long nstab               = NA_LONG; // BDF order reductions by STALD (CVODE)
  long nli                 = NA_LONG; // user Krylov iterations
  long nlcf                = NA_LONG;
  long npe                 = NA_LONG;
  long nps                 = NA_LONG;
  long njtv                = NA_LONG;
  long nsli                = NA_LONG; // stacked FGMRES iterations (FIRKODE)
  long nscf                = NA_LONG; // stacked FGMRES failures (FIRKODE)
  long nbs                 = NA_LONG; // block solves (FIRKODE)
  long nmass_setups        = NA_LONG;
  long nmass_solves        = NA_LONG;
  long nmass_mult          = NA_LONG;
  long ngeval              = NA_LONG;
  long nroots              = NA_LONG;
  long nrestarts           = NA_LONG; // re-initializations (--restart-every)
  std::string root_times   = "NA";    // "t:d0/d1;t:d0/d1;..."
  double err_mixed         = nan_value();
  double err_wrms          = nan_value();
  double err_abs_max       = nan_value();
  double dense_err_mixed   = nan_value();
  double dense_err_abs_max = nan_value();
  double hlast             = nan_value();
  double ele_wrms          = nan_value();
  double ele_max           = nan_value();
  std::string q_hist       = "NA"; // steps per BDF order 1..5 (--order-hist)
  double q_mean_step       = nan_value(); // step-weighted mean order
  double q_mean_time       = nan_value(); // time-weighted mean order
  long nst_ref             = NA_LONG;
  double ref_disagreement  = nan_value();
  double time_min          = nan_value();
  double time_median       = nan_value();
  int reps                 = 0;
  double t_setup           = nan_value(); // --profile: linear solver setups
  double t_solve           = nan_value(); // --profile: linear solves
  double t_rhs             = nan_value(); // --profile: right-hand sides
  double t_jac             = nan_value(); // --profile: Jacobian evaluations
  double t_profiled        = nan_value(); // --profile: wall time of that run
  int profiled             = 0;
  int flag                 = 0;
  // ARKODE steppers other than implicit ARKStep (NA elsewhere)
  std::string split     = "NA"; // the splitting of f (arkode-imex, arkode-mri)
  long nfe_e            = NA_LONG; // explicit (IMEX) or slow explicit (MRI) f
  long nfe_f            = NA_LONG; // fast right-hand sides (MRI)
  long nst_f            = NA_LONG; // fast steps (MRI)
  long nst_f_att        = NA_LONG; // fast step attempts (MRI)
  long netf_f           = NA_LONG; // fast error-test failures (MRI)
  long ndomeig          = NA_LONG; // dominant-eigenvalue updates (LSRK)
  long nfe_domeig       = NA_LONG; // f evaluations of the power iteration
  long lsrk_max_stages  = NA_LONG; // most stages of a step (LSRK)
  long nst_expstab      = NA_LONG; // steps limited by stability (ARKODE)
  long nst_acc          = NA_LONG; // steps limited by accuracy (ARKODE)
  double split_residual = nan_value(); // --check-split
  double split_jac_err  = nan_value(); // --check-split
  double domeig_ratio   = nan_value(); // --check-domeig
};

// Counters that are summed over the intervals of a restarted run
inline const std::vector<long RunStats::*>& summed_counters()
{
  static const std::vector<long RunStats::*> fields =
    {&RunStats::nst,         &RunStats::nst_att,      &RunStats::nfe,
     &RunStats::nfe_ls,      &RunStats::nsetups,      &RunStats::nje,
     &RunStats::nni,         &RunStats::nncf,         &RunStats::nsf,
     &RunStats::netf,        &RunStats::nstab,        &RunStats::nli,
     &RunStats::nlcf,        &RunStats::npe,          &RunStats::nps,
     &RunStats::njtv,        &RunStats::nsli,         &RunStats::nscf,
     &RunStats::nbs,         &RunStats::nmass_setups, &RunStats::nmass_solves,
     &RunStats::nmass_mult,  &RunStats::ngeval,       &RunStats::nfe_e,
     &RunStats::nfe_f,       &RunStats::nst_f,        &RunStats::nst_f_att,
     &RunStats::netf_f,      &RunStats::ndomeig,      &RunStats::nfe_domeig,
     &RunStats::nst_expstab, &RunStats::nst_acc};
  return fields;
}

// acc += part for every summed counter that part provides
inline void add_counters(RunStats& acc, const RunStats& part)
{
  for (auto f : summed_counters())
  {
    if (part.*f != NA_LONG)
    {
      acc.*f = (acc.*f == NA_LONG ? 0 : acc.*f) + part.*f;
    }
  }
}

inline void copy_counters(RunStats& dst, const RunStats& src)
{
  for (auto f : summed_counters()) { dst.*f = src.*f; }
}

inline std::string csv_header()
{
  return "problem,integrator,method,order,stages,mode,ls,prec,rtol,atol,N,"
         "fixed_h,"
         "nst,nst_att,nfe,nfe_ls,nsetups,nje,nni,nncf,nsf,netf,nstab,nli,nlcf,"
         "npe,"
         "nps,njtv,nsli,nscf,nbs,nmass_setups,nmass_solves,nmass_mult,ngeval,"
         "nroots,nrestarts,root_times,err_mixed,err_wrms,err_abs_max,dense_err_"
         "mixed,"
         "dense_err_abs_max,hlast,ele_wrms,ele_max,q_hist,q_mean_step,"
         "q_mean_time,nst_ref,ref_disagreement,"
         "time_min,time_median,reps,t_setup,t_solve,t_rhs,t_jac,t_profiled,"
         "profiled,flag,"
         "split,nfe_e,nfe_f,nst_f,nst_f_att,netf_f,ndomeig,nfe_domeig,"
         "lsrk_max_stages,nst_expstab,nst_acc,split_residual,split_jac_err,"
         "domeig_ratio";
}

inline std::string fmt_long(long v)
{
  return v == NA_LONG ? std::string("NA") : std::to_string(v);
}

inline std::string fmt_real(double v, const char* fmt = "%.9e")
{
  if (std::isnan(v)) { return "NA"; }
  char buf[64];
  std::snprintf(buf, sizeof(buf), fmt, v);
  return buf;
}

inline std::string csv_row(const RunStats& s)
{
  std::string r;
  r += s.problem + "," + s.integrator + "," + s.method + ",";
  r += std::to_string(s.order) + "," + std::to_string(s.stages) + ",";
  r += s.mode + ",";
  r += s.ls + "," + s.prec + ",";
  r += fmt_real(s.rtol, "%.3e") + "," + fmt_real(s.atol, "%.3e") + ",";
  r += fmt_long(s.N) + "," + fmt_real(s.fixed_h, "%.9e") + ",";
  for (long v :
       {s.nst,          s.nst_att,    s.nfe,    s.nfe_ls, s.nsetups,
        s.nje,          s.nni,        s.nncf,   s.nsf,    s.netf,
        s.nstab,        s.nli,        s.nlcf,   s.npe,    s.nps,
        s.njtv,         s.nsli,       s.nscf,   s.nbs,    s.nmass_setups,
        s.nmass_solves, s.nmass_mult, s.ngeval, s.nroots, s.nrestarts})
  {
    r += fmt_long(v) + ",";
  }
  r += s.root_times + ",";
  for (double v : {s.err_mixed, s.err_wrms, s.err_abs_max, s.dense_err_mixed,
                   s.dense_err_abs_max, s.hlast, s.ele_wrms, s.ele_max})
  {
    r += fmt_real(v) + ",";
  }
  r += s.q_hist + "," + fmt_real(s.q_mean_step, "%.4f") + "," +
       fmt_real(s.q_mean_time, "%.4f") + ",";
  r += fmt_long(s.nst_ref) + "," + fmt_real(s.ref_disagreement) + ",";
  r += fmt_real(s.time_min, "%.6e") + "," + fmt_real(s.time_median, "%.6e") + ",";
  r += std::to_string(s.reps) + ",";
  for (double v : {s.t_setup, s.t_solve, s.t_rhs, s.t_jac, s.t_profiled})
  {
    r += fmt_real(v, "%.6e") + ",";
  }
  r += std::to_string(s.profiled) + "," + std::to_string(s.flag) + ",";
  r += s.split + ",";
  for (long v : {s.nfe_e, s.nfe_f, s.nst_f, s.nst_f_att, s.netf_f, s.ndomeig,
                 s.nfe_domeig, s.lsrk_max_stages, s.nst_expstab, s.nst_acc})
  {
    r += fmt_long(v) + ",";
  }
  r += fmt_real(s.split_residual, "%.3e") + "," +
       fmt_real(s.split_jac_err, "%.3e") + "," + fmt_real(s.domeig_ratio, "%.6f");
  return r;
}

// Error metrics between y and yref (both serial vectors of equal length):
//   mixed  = max_i |y_i - yref_i| / (atol/rtol + |yref_i|)
//   wrms   = sqrt( mean_i ((y_i - yref_i) / (rtol |yref_i| + atol))^2 )
//   absmax = max_i |y_i - yref_i|
struct ErrorMetrics
{
  double mixed = nan_value(), wrms = nan_value(), absmax = nan_value();
};

inline ErrorMetrics compute_errors(const sunrealtype* y, const sunrealtype* yref,
                                   long N, double rtol, double atol)
{
  ErrorMetrics e;
  e.mixed = e.wrms = e.absmax = 0.0;
  const double ratio          = atol / rtol;
  for (long i = 0; i < N; i++)
  {
    const double d = std::fabs((double)y[i] - (double)yref[i]);
    const double m = d / (ratio + std::fabs((double)yref[i]));
    const double w = d / (rtol * std::fabs((double)yref[i]) + atol);
    if (m > e.mixed) { e.mixed = m; }
    if (d > e.absmax) { e.absmax = d; }
    e.wrms += w * w;
  }
  e.wrms = std::sqrt(e.wrms / (double)N);
  return e;
}

// Full-precision text form of a sunrealtype (every digit needed to round
// trip: 17 for double, 21 for the 80-bit extended type)
inline std::string fmt_full(sunrealtype v)
{
  char buf[64];
#if defined(SUNDIALS_EXTENDED_PRECISION)
  std::snprintf(buf, sizeof(buf), "%.*Le",
                std::numeric_limits<long double>::max_digits10, (long double)v);
#else
  std::snprintf(buf, sizeof(buf), "%.*e",
                std::numeric_limits<double>::max_digits10, (double)v);
#endif
  return buf;
}

// Reference-solution file: "# key=value" comment lines, then a line
// "N tf", then N lines with one full-precision value each (a double
// build reads a 21-digit extended-precision reference rounded to double).
struct Reference
{
  long N    = 0;
  double tf = nan_value();
  std::vector<sunrealtype> y;
  std::vector<std::string> comments;
};

inline bool read_reference(const std::string& path, Reference& ref)
{
  std::ifstream in(path);
  if (!in) { return false; }
  std::string line;
  ref.comments.clear();
  while (std::getline(in, line))
  {
    if (!line.empty() && line[0] == '#')
    {
      ref.comments.push_back(line);
      continue;
    }
    break;
  }
  if (std::sscanf(line.c_str(), "%ld %lf", &ref.N, &ref.tf) != 2)
  {
    return false;
  }
  ref.y.assign((size_t)ref.N, SUN_RCONST(0.0));
  for (long i = 0; i < ref.N; i++)
  {
    if (!(in >> ref.y[(size_t)i])) { return false; }
  }
  return true;
}

inline bool write_reference(const std::string& path, const Reference& ref)
{
  std::FILE* fp = std::fopen(path.c_str(), "w");
  if (fp == nullptr) { return false; }
  for (const auto& c : ref.comments) { std::fprintf(fp, "%s\n", c.c_str()); }
  std::fprintf(fp, "%ld %.17e\n", ref.N, ref.tf);
  for (long i = 0; i < ref.N; i++)
  {
    std::fprintf(fp, "%s\n", fmt_full(ref.y[(size_t)i]).c_str());
  }
  std::fclose(fp);
  return true;
}

} // namespace benchmark

#endif // FIRKODE_BENCHMARK_COMMON_HPP_
