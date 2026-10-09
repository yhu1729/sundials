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
 * firkode_benchmark: run one benchmark cell (problem x integrator x
 * options) and print a CSV header and one data row to stdout.
 *
 * Exit codes: 0 success, 1 harness/usage error, 2 unsupported
 * combination, 3 the integrator returned a failure flag (a row was
 * still printed).
 * -----------------------------------------------------------------*/

#include <algorithm>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <firkode/firkode.h>
#include <sundials/sundials_core.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>

#include "benchmark_cli.hpp"
#include "benchmark_common.hpp"
#include "integrator.hpp"
#include "problem.hpp"
#include "timed_linsol.hpp"

namespace benchmark {

std::unique_ptr<Integrator> make_firkode_integrator();
#ifdef BENCHMARK_HAVE_CVODE
std::unique_ptr<Integrator> make_cvode_integrator();
#endif
#ifdef BENCHMARK_HAVE_ARKODE
std::unique_ptr<Integrator> make_arkode_integrator();
#endif
#ifdef BENCHMARK_HAVE_IDA
std::unique_ptr<Integrator> make_ida_integrator();
#endif

std::unique_ptr<Integrator> make_integrator(const std::string& name)
{
  if (name == "firkode") { return make_firkode_integrator(); }
#ifdef BENCHMARK_HAVE_CVODE
  if (name == "cvode") { return make_cvode_integrator(); }
#endif
#ifdef BENCHMARK_HAVE_ARKODE
  if (name == "arkode") { return make_arkode_integrator(); }
#endif
#ifdef BENCHMARK_HAVE_IDA
  if (name == "ida") { return make_ida_integrator(); }
#endif
  throw std::runtime_error("unsupported: integrator '" + name +
                           "' is unknown or not built");
}

namespace {

// Output times of the evolve loop
std::vector<double> output_times(const Options& o, double t0, double tf)
{
  std::vector<double> ts;
  if (o.grid_geom_t1 > 0.0 && o.grid_geom_mult > 0.0)
  {
    double t = o.grid_geom_t1;
    for (int k = 0; k < std::max(o.grid, 1); k++, t *= o.grid_geom_mult)
    {
      ts.push_back(t);
    }
  }
  else
  {
    const int K = std::max(o.grid, 1);
    for (int k = 1; k <= K; k++) { ts.push_back(t0 + (tf - t0) * k / K); }
    ts.back() = tf;
  }
  return ts;
}

// Options of an in-process reference run: same problem and linear
// solver, everything else reset to the integrator defaults.
Options reference_options(const Options& o, const std::string& integ, double rtol)
{
  Options r    = o;
  r.integrator = integ;
  r.mode       = Mode::Evolve;
  r.rtol       = rtol;
  r.atol       = rtol * (o.atol / o.rtol);
  r.atol_vec.clear();
  r.stages     = 0;
  r.order      = 0;
  r.tstop      = true;
  r.maxsteps   = 200000000;
  r.maxnef     = -1;
  r.h0         = 0.0;
  r.fixed_h    = 0.0;
  r.order_hist = false;
  r.profile    = false;
  r.dump_times.clear();
  r.h_one = 0.0;
  r.predictor.clear();
  r.refilter      = true;
  r.newton_safety = true;
  r.controller.clear();
  r.impgus_k1 = r.impgus_k2 = -1.0;
  r.mode_name               = "beginner";
  r.stald                   = false;
  r.arkode_linear           = false;
  r.ark_predictor           = -1;
  r.msbp = r.msbj = -1;
  r.dgmax = r.dgmax_jbad = -1.0;
  r.reuse_policy.clear();
  r.maxcor       = -1;
  r.nlscoef      = -1.0;
  r.epslin       = -1.0;
  r.stage_maxl   = 0;
  r.stage_maxrs  = -1;
  r.stage_epslin = -1.0;
  return r;
}

// Internal boundaries of a restarted run: K - 1 equally spaced times
std::vector<double> restart_boundaries(const Options& o, double t0, double tf)
{
  std::vector<double> b;
  for (int k = 1; k < o.restart_every; k++)
  {
    b.push_back(t0 + (tf - t0) * k / o.restart_every);
  }
  return b;
}

bool is_boundary(const std::vector<double>& b, double t)
{
  return std::find(b.begin(), b.end(), t) != b.end();
}

// Write a dense matrix (row-major) with a one-line header
void write_dense(const std::string& path, const std::string& header,
                 const std::vector<double>& a, long N)
{
  std::FILE* fp = std::fopen(path.c_str(), "w");
  if (fp == nullptr) { throw std::runtime_error("cannot write " + path); }
  std::fprintf(fp, "%s\n", header.c_str());
  for (long i = 0; i < N; i++)
  {
    for (long j = 0; j < N; j++)
    {
      std::fprintf(fp, "%.17e%c", a[(size_t)(i * N + j)], j + 1 < N ? ' ' : '\n');
    }
  }
  std::fclose(fp);
}

// Copy a dense or band SUNMatrix into a row-major dense array
void to_dense(SUNMatrix A, std::vector<double>& out, long N)
{
  out.assign((size_t)(N * N), 0.0);
  if (SUNMatGetID(A) == SUNMATRIX_DENSE)
  {
    for (long i = 0; i < N; i++)
    {
      for (long j = 0; j < N; j++)
      {
        out[(size_t)(i * N + j)] = (double)SM_ELEMENT_D(A, i, j);
      }
    }
  }
  else if (SUNMatGetID(A) == SUNMATRIX_SPARSE)
  {
    const sunindextype* ptrs = SUNSparseMatrix_IndexPointers(A);
    const sunindextype* vals = SUNSparseMatrix_IndexValues(A);
    const sunrealtype* data  = SUNSparseMatrix_Data(A);
    const bool csr           = SUNSparseMatrix_SparseType(A) == SUN_CSR_MAT;
    for (long r = 0; r < N; r++)
    {
      for (sunindextype k = ptrs[r]; k < ptrs[r + 1]; k++)
      {
        const long c = (long)vals[k];
        out[(size_t)(csr ? r * N + c : c * N + r)] += (double)data[k];
      }
    }
  }
  else
  {
    const long mu = (long)SM_UBAND_B(A), ml = (long)SM_LBAND_B(A);
    for (long j = 0; j < N; j++)
    {
      for (long i = std::max(0L, j - mu); i <= std::min(N - 1, j + ml); i++)
      {
        out[(size_t)(i * N + j)] = (double)SM_ELEMENT_B(A, i, j);
      }
    }
  }
}

// Jacobian dump for the spectrum analysis: the user's Jacobian through a
// matrix of the problem's default kind, or forward differences of the
// right-hand side (when no Jacobian exists or --dump-jac-dq is given).
// The mass matrix, if any, goes to mass.txt once.
void dump_jacobian(Problem& p, const Options& o, double t, N_Vector y,
                   SUNContext ctx, int index, bool& mass_written)
{
  const ProblemSpec& sp = p.spec();
  const long N          = (long)sp.N;
  if (N > 2000) { throw std::runtime_error("--dump-jac: N exceeds 2000"); }
  std::vector<double> J;
  std::string source;
  const bool direct = is_direct(sp.ls_default);
  if (!o.dump_dq && p.jac() != nullptr && direct)
  {
    SUNMatrix A = nullptr;
    if (sp.ls_default == LinSolKind::Dense)
    {
      A = SUNDenseMatrix(sp.N, sp.N, ctx);
    }
    else if (sp.ls_default == LinSolKind::Band)
    {
      A = SUNBandMatrix(sp.N, sp.mu, sp.ml, ctx);
    }
    else { A = SUNSparseMatrix(sp.N, sp.N, sp.nnz, SUN_CSR_MAT, ctx); }
    check_ptr(A, "SUNMatrix(dump)");
    N_Vector t1 = N_VClone(y), t2 = N_VClone(y), t3 = N_VClone(y);
    check(SUNMatZero(A), "SUNMatZero");
    check(p.jac()((sunrealtype)t, y, nullptr, A, p.user_data(), t1, t2, t3),
          "jac(dump)");
    to_dense(A, J, N);
    N_VDestroy(t1);
    N_VDestroy(t2);
    N_VDestroy(t3);
    SUNMatDestroy(A);
    source = "jac";
  }
  else
  {
    J.assign((size_t)(N * N), 0.0);
    N_Vector f0 = N_VClone(y), f1 = N_VClone(y), yp = N_VClone(y);
    check(p.rhs()((sunrealtype)t, y, f0, p.user_data()), "rhs(dump)");
    const double eps       = std::sqrt(std::numeric_limits<double>::epsilon());
    const sunrealtype* yd  = N_VGetArrayPointer(y);
    const sunrealtype* f0d = N_VGetArrayPointer(f0);
    sunrealtype* ypd       = N_VGetArrayPointer(yp);
    const sunrealtype* f1d = N_VGetArrayPointer(f1);
    for (long j = 0; j < N; j++)
    {
      N_VScale(SUN_RCONST(1.0), y, yp);
      const double del = eps * std::max(std::fabs((double)yd[j]), 1.0);
      ypd[j]           = (sunrealtype)((double)yd[j] + del);
      check(p.rhs()((sunrealtype)t, yp, f1, p.user_data()), "rhs(dump)");
      for (long i = 0; i < N; i++)
      {
        J[(size_t)(i * N + j)] = ((double)f1d[i] - (double)f0d[i]) / del;
      }
    }
    N_VDestroy(f0);
    N_VDestroy(f1);
    N_VDestroy(yp);
    source = "dq";
  }
  char hdr[128];
  std::snprintf(hdr, sizeof(hdr), "# t=%.17e N=%ld source=%s", t, N,
                source.c_str());
  write_dense(o.dump_dir + "/jac_" + std::to_string(index) + ".txt", hdr, J, N);
  if (sp.has_mass && !mass_written)
  {
    SUNMatrix M = SUNBandMatrix(sp.N, sp.mass_mu, sp.mass_ml, ctx);
    check_ptr(M, "SUNBandMatrix(mass dump)");
    N_Vector t1 = N_VClone(y), t2 = N_VClone(y), t3 = N_VClone(y);
    check(SUNMatZero(M), "SUNMatZero");
    check(p.mass()((sunrealtype)t, M, p.user_data(), t1, t2, t3), "mass(dump)");
    std::vector<double> Md;
    to_dense(M, Md, N);
    N_VDestroy(t1);
    N_VDestroy(t2);
    N_VDestroy(t3);
    SUNMatDestroy(M);
    std::snprintf(hdr, sizeof(hdr), "# mass N=%ld", N);
    write_dense(o.dump_dir + "/mass.txt", hdr, Md, N);
    mass_written = true;
  }
}

// Integrate with the given options and record y at each time in ts
long run_trajectory(Problem& p, const Options& o, N_Vector y0, SUNContext ctx,
                    const std::vector<double>& ts,
                    std::vector<std::vector<sunrealtype>>& out)
{
  auto integ = make_integrator(o.integrator);
  N_Vector y = N_VClone(y0);
  N_VScale(SUN_RCONST(1.0), y0, y);
  integ->setup(p, o, y, ctx);
  out.assign(ts.size(),
             std::vector<sunrealtype>((size_t)p.spec().N, SUN_RCONST(0.0)));
  // the reference follows the same restart schedule as the benchmark run
  const double t0                  = (double)p.spec().t0;
  const double tf                  = ts.back();
  const std::vector<double> bounds = restart_boundaries(o, t0, tf);
  std::vector<double> events       = ts;
  events.insert(events.end(), bounds.begin(), bounds.end());
  std::sort(events.begin(), events.end());
  events.erase(std::unique(events.begin(), events.end()), events.end());
  if (!bounds.empty() && o.tstop)
  {
    integ->set_stop_time((sunrealtype)events[0]);
  }
  RunStats acc;
  sunrealtype tret   = (sunrealtype)t0;
  size_t next_sample = 0;
  for (size_t k = 0; k < events.size(); k++)
  {
    const double t_event = events[k];
    int flag             = integ->evolve((sunrealtype)t_event, y, tret, false);
    if (!integ->is_success(flag))
    {
      N_VDestroy(y);
      throw std::runtime_error("reference run failed with flag " +
                               std::to_string(flag));
    }
    if (next_sample < ts.size() && ts[next_sample] == t_event)
    {
      const sunrealtype* yd = N_VGetArrayPointer(y);
      for (long i = 0; i < p.spec().N; i++)
      {
        out[next_sample][(size_t)i] = yd[i];
      }
      next_sample++;
    }
    if (is_boundary(bounds, t_event))
    {
      RunStats part;
      integ->harvest(part);
      add_counters(acc, part);
      p.restart_update((sunrealtype)t_event, y);
      integ->reinit((sunrealtype)t_event, y);
      if (o.tstop) { integ->set_stop_time((sunrealtype)events[k + 1]); }
    }
  }
  RunStats s;
  integ->harvest(s);
  add_counters(acc, s);
  integ->teardown();
  N_VDestroy(y);
  return acc.nst;
}

void exact_trajectory(Problem& p, N_Vector tmpl, const std::vector<double>& ts,
                      std::vector<std::vector<sunrealtype>>& out)
{
  N_Vector y = N_VClone(tmpl);
  out.assign(ts.size(),
             std::vector<sunrealtype>((size_t)p.spec().N, SUN_RCONST(0.0)));
  for (size_t k = 0; k < ts.size(); k++)
  {
    check(p.exact((sunrealtype)ts[k], y), "exact");
    const sunrealtype* yd = N_VGetArrayPointer(y);
    for (long i = 0; i < p.spec().N; i++) { out[k][(size_t)i] = yd[i]; }
  }
  N_VDestroy(y);
}

std::string format_real(double v)
{
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.17e", v);
  return buf;
}

double median(std::vector<double> v)
{
  std::sort(v.begin(), v.end());
  const size_t n = v.size();
  return (n % 2 == 1) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

int run(int argc, char* argv[])
{
  Options opts = parse_args(argc, argv);

  if (opts.version)
  {
    std::printf("SUNDIALS %s\n", SUNDIALS_VERSION);
    std::printf("firkode_benchmark packages: firkode");
#ifdef BENCHMARK_HAVE_CVODE
    std::printf(" cvode");
#endif
#ifdef BENCHMARK_HAVE_ARKODE
    std::printf(" arkode");
#endif
#ifdef BENCHMARK_HAVE_IDA
    std::printf(" ida");
#endif
    std::
      printf("\nprecision: %s, sizeof(sunrealtype) = %zu, unit roundoff = %s, "
             "LDBL_MANT_DIG = %d\n",
#if defined(SUNDIALS_EXTENDED_PRECISION)
             "extended",
#elif defined(SUNDIALS_SINGLE_PRECISION)
             "single",
#else
             "double",
#endif
             sizeof(sunrealtype), fmt_full(SUN_UNIT_ROUNDOFF).c_str(),
             (int)LDBL_MANT_DIG);
    std::printf("FIRKODE default parameters:\n");
    SUNContext ctx;
    check(SUNContext_Create(SUN_COMM_NULL, &ctx), "SUNContext_Create");
    void* mem = FIRKodeCreate(ctx);
    FIRKodeWriteParameters(mem, stdout);
    FIRKodeFree(&mem);
    SUNContext_Free(&ctx);
    return 0;
  }

  SUNContext ctx;
  check(SUNContext_Create(SUN_COMM_NULL, &ctx), "SUNContext_Create");

  std::unique_ptr<Problem> problem = make_problem(opts.problem);
  if (opts.profile) { problem.reset(new TimedProblem(std::move(problem))); }
  {
    const std::vector<std::string> keys = problem->param_keys();
    for (const auto& kv : opts.extra)
    {
      if (std::find(keys.begin(), keys.end(), kv.first) == keys.end())
      {
        throw std::runtime_error("unrecognized argument: --" + kv.first);
      }
    }
  }
  problem->configure(opts);
  const ProblemSpec& spec = problem->spec();
  if (opts.rtol <= 0.0) { opts.rtol = (double)spec.rtol_default; }
  if (opts.atol <= 0.0) { opts.atol = (double)spec.atol_default; }
  if (opts.tf <= 0.0) { opts.tf = (double)spec.tf; }
  const double t0 = (double)spec.t0;
  const double tf = opts.tf;
  if (opts.restart_every > 1 &&
      (opts.mode != Mode::Evolve || opts.grid != 1 || opts.grid_geom_t1 > 0.0 ||
       !opts.dump_times.empty()))
  {
    throw std::runtime_error("--restart-every needs the plain evolve mode with "
                             "a single output time and no Jacobian dump");
  }
  if (opts.order_hist && (opts.mode != Mode::Evolve || !opts.tstop ||
                          opts.grid != 1 || opts.grid_geom_t1 > 0.0))
  {
    throw std::runtime_error("--order-hist needs the plain evolve mode with "
                             "a single output time and a stop time");
  }
  problem->prepare(opts);

  N_Vector y0 = N_VNew_Serial(spec.N, ctx);
  check_ptr(y0, "N_VNew_Serial");
  problem->init(y0);

  RunStats st;
  st.problem    = spec.name;
  st.integrator = opts.integrator;
  st.rtol       = opts.rtol;
  st.atol       = opts.atol;
  st.N          = (long)spec.N;
  st.reps       = opts.reps;

  // ---- reference solution at tf (and at the dense-output samples) ----
  std::vector<double> sample_times;
  if (opts.mode == Mode::Dense)
  {
    for (int k = 1; k <= opts.dense_k; k++)
    {
      sample_times.push_back(t0 + (tf - t0) * k / opts.dense_k);
    }
    sample_times.back() = tf;
  }
  else { sample_times.push_back(tf); }

  std::vector<std::vector<sunrealtype>> ref; // ref[k] at sample_times[k]
  bool have_ref = false;
  if (spec.has_exact && opts.mode != Mode::OneStep)
  {
    exact_trajectory(*problem, y0, sample_times, ref);
    have_ref = true;
  }
  else if (opts.mode != Mode::OneStep)
  {
    Reference rf;
    bool from_file = false;
    if (!opts.ref_file.empty() && read_reference(opts.ref_file, rf))
    {
      if (rf.N != (long)spec.N || std::fabs(rf.tf - tf) > 1e-12 * std::fabs(tf))
      {
        throw std::runtime_error("reference file does not match N or tf");
      }
      from_file = true;
    }
    if (from_file && opts.mode != Mode::Dense)
    {
      ref.assign(1, rf.y);
      have_ref = true;
      for (const auto& c : rf.comments)
      {
        if (c.rfind("# nst_ref=", 0) == 0)
        {
          st.nst_ref = std::stol(c.substr(10));
        }
        if (c.rfind("# ref_disagreement=", 0) == 0)
        {
          st.ref_disagreement = std::stod(c.substr(19));
        }
      }
    }
    else if (opts.ref_inprocess || !opts.ref_file.empty() ||
             !opts.write_ref.empty() || opts.mode == Mode::Dense)
    {
      Options ro = reference_options(opts, opts.ref_integrator, opts.ref_rtol);
      std::string chk = "none";
      st.nst_ref = run_trajectory(*problem, ro, y0, ctx, sample_times, ref);
      have_ref   = true;
      if (opts.ref_check)
      {
        // second reference with a different integrator; report the mixed
        // error between the two
        std::vector<std::vector<sunrealtype>> ref2;
        chk = opts.ref_check_integrator;
        if (chk.empty()) { chk = spec.has_mass ? "ida" : "cvode"; }
        if (chk == opts.ref_integrator)
        {
          throw std::runtime_error(
            "--ref-check-integrator must differ from --ref-integrator");
        }
        Options ro2 = reference_options(opts, chk, opts.ref_rtol);
        run_trajectory(*problem, ro2, y0, ctx, sample_times, ref2);
        double dis = 0.0;
        for (size_t k = 0; k < ref.size(); k++)
        {
          ErrorMetrics e = compute_errors(ref[k].data(), ref2[k].data(),
                                          (long)spec.N, opts.rtol, opts.atol);
          dis            = std::max(dis, e.mixed);
        }
        st.ref_disagreement = dis;
      }
      const std::string path = !opts.write_ref.empty() ? opts.write_ref
                                                       : opts.ref_file;
      if (!path.empty() && opts.mode != Mode::Dense)
      {
        rf.N  = (long)spec.N;
        rf.tf = tf;
        rf.y  = ref.back();
        rf.comments.clear();
        rf.comments.push_back("# problem=" + spec.name);
        rf.comments.push_back("# ref_integrator=" + opts.ref_integrator);
        rf.comments.push_back("# ref_check_integrator=" + chk);
        rf.comments.push_back("# ref_rtol=" + format_real(opts.ref_rtol));
        rf.comments.push_back("# ref_atol=" + format_real(ro.atol));
        rf.comments.push_back("# nst_ref=" + std::to_string(st.nst_ref));
        if (!std::isnan(st.ref_disagreement))
        {
          rf.comments.push_back("# ref_disagreement=" +
                                format_real(st.ref_disagreement));
        }
        if (!write_reference(path, rf))
        {
          throw std::runtime_error("cannot write reference file " + path);
        }
      }
    }
  }

  // ---- the benchmark run (repeated opts.reps times) ----
  std::vector<double> times;
  std::unique_ptr<Integrator> integ;
  N_Vector y = nullptr;
  int flag   = 0;
  RunStats restart_acc; // counters of the completed intervals (--restart-every)
  int exit_code = 0;
  Timer timer;

  for (int rep = 0; rep < opts.reps; rep++)
  {
    if (integ) { integ->teardown(); }
    if (y != nullptr) { N_VDestroy(y); }
    integ = make_integrator(opts.integrator);
    y     = N_VClone(y0);
    N_VScale(SUN_RCONST(1.0), y0, y);

    profile_stats().reset();
    timer.start();
    integ->setup(*problem, opts, y, ctx);
    sunrealtype tret = (sunrealtype)t0;
    flag             = 0;

    if (opts.mode == Mode::Evolve || opts.mode == Mode::FixedStep ||
        opts.mode == Mode::Roots)
    {
      // Follow the output grid exactly as the examples do (the first output
      // time enters the initial-step heuristic); root returns resume the
      // integration towards the same output time.
      std::string rts;
      long nroots = 0;
      std::vector<int> dirs;
      constexpr int QMAX = 5;
      std::vector<long> qhist(QMAX + 1, 0);
      std::vector<double> qtime(QMAX + 1, 0.0);
      std::vector<double> touts = output_times(opts, t0, tf);
      for (double td : opts.dump_times)
      {
        if (td > t0 && td <= tf) { touts.push_back(td); }
      }
      std::sort(touts.begin(), touts.end());
      const std::vector<double> bounds = restart_boundaries(opts, t0, tf);
      touts.insert(touts.end(), bounds.begin(), bounds.end());
      std::sort(touts.begin(), touts.end());
      touts.erase(std::unique(touts.begin(), touts.end()), touts.end());
      if (touts.size() > 1 && opts.tstop)
      {
        integ->set_stop_time((sunrealtype)touts[0]);
      }
      RunStats acc;
      bool mass_written = false;
      int dump_index    = 0;
      for (size_t it = 0; it < touts.size(); it++)
      {
        const double tout = touts[it];
        // every output time is a stop time (not only the final one), so
        // that a grid of switching times or a restart schedule is honoured
        // exactly rather than by interpolation
        if (opts.tstop && touts.size() > 1 && it > 0)
        {
          integ->set_stop_time((sunrealtype)tout);
        }
        while (true)
        {
          flag = integ->evolve((sunrealtype)tout, y, tret, opts.order_hist);
          if (opts.order_hist && integ->is_success(flag))
          {
            const int q = integ->last_order();
            if (q >= 1 && q <= QMAX)
            {
              qhist[(size_t)q]++;
              qtime[(size_t)q] += std::fabs((double)integ->last_step());
            }
            if ((double)tret < tout) { continue; }
          }
          if (integ->is_root(flag))
          {
            integ->root_info(dirs);
            rts += (rts.empty() ? "" : ";") + format_real((double)tret) + ":";
            for (size_t i = 0; i < dirs.size(); i++)
            {
              rts += (i ? "/" : "") + std::to_string(dirs[i]);
            }
            nroots++;
            if ((double)tret < tout) { continue; }
          }
          break;
        }
        if (!integ->is_success(flag) && !integ->is_root(flag)) { break; }
        if (rep == opts.reps - 1 &&
            std::find(opts.dump_times.begin(), opts.dump_times.end(), tout) !=
              opts.dump_times.end())
        {
          dump_jacobian(*problem, opts, tout, y, ctx, dump_index++, mass_written);
        }
        if (is_boundary(bounds, tout))
        {
          RunStats part;
          integ->harvest(part);
          add_counters(acc, part);
          problem->restart_update((sunrealtype)tout, y);
          integ->reinit((sunrealtype)tout, y);
        }
      }
      if (!bounds.empty() && rep == opts.reps - 1)
      {
        restart_acc  = acc;
        st.nrestarts = (long)bounds.size();
      }
      if (opts.mode == Mode::Roots)
      {
        st.root_times = rts.empty() ? "NA" : rts;
        st.nroots     = nroots;
      }
      if (opts.order_hist && integ->last_order() >= 0)
      {
        long nq   = 0;
        double sq = 0.0, st_sum = 0.0, sqt = 0.0;
        st.q_hist.clear();
        for (int q = 1; q <= QMAX; q++)
        {
          st.q_hist += (q > 1 ? ";" : "") + std::to_string(qhist[(size_t)q]);
          nq += qhist[(size_t)q];
          sq += q * (double)qhist[(size_t)q];
          st_sum += qtime[(size_t)q];
          sqt += q * qtime[(size_t)q];
        }
        if (nq > 0) { st.q_mean_step = sq / (double)nq; }
        if (st_sum > 0.0) { st.q_mean_time = sqt / st_sum; }
      }
    }
    else if (opts.mode == Mode::OneStep)
    {
      flag = integ->evolve((sunrealtype)tf, y, tret, true);
    }
    else if (opts.mode == Mode::Dense)
    {
      N_Vector dky = N_VClone(y);
      size_t next  = 0;
      double dmix = 0.0, dabs = 0.0;
      while ((double)tret < tf)
      {
        flag = integ->evolve((sunrealtype)tf, y, tret, true);
        if (!integ->is_success(flag)) { break; }
        while (next < sample_times.size() && sample_times[next] <= (double)tret)
        {
          check(integ->get_dky((sunrealtype)sample_times[next], 0, dky),
                "GetDky");
          ErrorMetrics e = compute_errors(N_VGetArrayPointer(dky),
                                          ref[next].data(), (long)spec.N,
                                          opts.rtol, opts.atol);
          dmix           = std::max(dmix, e.mixed);
          dabs           = std::max(dabs, e.absmax);
          next++;
        }
      }
      st.dense_err_mixed   = dmix;
      st.dense_err_abs_max = dabs;
      N_VDestroy(dky);
    }
    times.push_back(timer.stop());

    if (rep == opts.reps - 1)
    {
      integ->harvest(st);
      if (st.nrestarts > 0)
      {
        add_counters(restart_acc, st);
        copy_counters(st, restart_acc);
      }
      st.method = integ->method();
      st.order  = integ->order();
      st.stages = integ->stages();
      st.mode   = integ->mode_name();
      st.ls     = integ->ls_name();
      st.prec   = integ->prec_name();
      st.flag   = flag;
      st.hlast  = (double)integ->last_step();
      if (opts.mode == Mode::FixedStep) { st.fixed_h = opts.fixed_h; }
      if (opts.profile)
      {
        const ProfileStats& ps = profile_stats();
        st.t_setup             = ps.t_setup;
        st.t_solve             = ps.t_solve;
        st.t_rhs               = ps.t_rhs;
        st.t_jac               = ps.t_jac;
        st.t_profiled          = times.back();
        st.profiled            = 1;
      }

      if (opts.mode == Mode::OneStep)
      {
        N_Vector ele = N_VClone(y), ewt = N_VClone(y);
        if (integ->est_local_errors(ele) == 0 && integ->err_weights(ewt) == 0)
        {
          st.ele_wrms = (double)N_VWrmsNorm(ele, ewt);
          st.ele_max  = (double)N_VMaxNorm(ele);
        }
        if (spec.has_exact)
        {
          N_Vector yex = N_VClone(y);
          check(problem->exact(tret, yex), "exact");
          N_VLinearSum(SUN_RCONST(1.0), y, SUN_RCONST(-1.0), yex, yex);
          st.err_abs_max = (double)N_VMaxNorm(yex);
          st.err_wrms    = (double)N_VWrmsNorm(yex, ewt);
          N_VDestroy(yex);
        }
        N_VDestroy(ele);
        N_VDestroy(ewt);
      }
      else if (have_ref && integ->is_success(flag))
      {
        ErrorMetrics e = compute_errors(N_VGetArrayPointer(y), ref.back().data(),
                                        (long)spec.N, opts.rtol, opts.atol);
        st.err_mixed   = e.mixed;
        st.err_wrms    = e.wrms;
        st.err_abs_max = e.absmax;
      }
      if (!integ->is_success(flag) && !integ->is_root(flag)) { exit_code = 3; }
      if (opts.verbose)
      {
        std::fprintf(stderr, "flag=%d tret=%.17e nst=%ld nfe=%ld nni=%ld\n",
                     flag, (double)tret, st.nst, st.nfe, st.nni);
      }
    }
    if (!integ->is_success(flag) && !integ->is_root(flag)) { break; }
  }

  st.time_min    = *std::min_element(times.begin(), times.end());
  st.time_median = median(times);
  st.reps        = (int)times.size();

  if (!opts.no_header) { std::printf("%s\n", csv_header().c_str()); }
  std::printf("%s\n", csv_row(st).c_str());

  integ->teardown();
  integ.reset();
  if (y != nullptr) { N_VDestroy(y); }
  N_VDestroy(y0);
  problem.reset();
  SUNContext_Free(&ctx);
  return exit_code;
}

} // namespace
} // namespace benchmark

int main(int argc, char* argv[])
{
  try
  {
    return benchmark::run(argc, argv);
  }
  catch (const std::exception& e)
  {
    std::fprintf(stderr, "firkode_benchmark: %s\n", e.what());
    const std::string msg = e.what();
    return msg.rfind("unsupported", 0) == 0 ? 2 : 1;
  }
}
