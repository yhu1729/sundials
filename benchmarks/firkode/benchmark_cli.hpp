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
 * Command-line options of firkode_benchmark.
 * -----------------------------------------------------------------*/

#ifndef FIRKODE_BENCHMARK_CLI_HPP_
#define FIRKODE_BENCHMARK_CLI_HPP_

#include <algorithm>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <sundials/sundials_types.h>

namespace benchmark {

enum class Mode
{
  Evolve,
  FixedStep,
  OneStep,
  Roots,
  Dense
};

struct Options
{
  // what to run
  std::string problem;
  std::string integrator;
  Mode mode = Mode::Evolve;

  // tolerances (negative = problem default)
  double rtol = -1.0;
  double atol = -1.0;
  std::vector<double> atol_vec; // optional vector absolute tolerance

  // method
  int stages = 0; // FIRKODE (0 = default 3)
  int order  = 0; // ARKODE DIRK order, CVODE/IDA max order (0 = default)

  // linear solver
  std::string ls;            // dense|band|spgmr|pcg ("" = problem default)
  std::string prec = "left"; // none|left|right
  int maxl         = 0;
  int maxrs        = -1;
  double epslin    = -1.0;
  long msbp        = -1;
  long msbj        = -1;
  double dgmax =
    -1.0; // gamma change forcing a setup (FIRKODE periodic, CVODE, ARKODE)
  double dgmax_jbad =
    -1.0; // gamma change forcing a Jacobian (FIRKODE periodic, CVODE)
  std::string reuse_policy; // FIRKODE: radau5|periodic ("" = default radau5)
  int stage_maxl      = 0;
  int stage_maxrs     = -1;
  double stage_epslin = -1.0;
  std::string psolve  = "jacobi"; // matrix-free heat2d: jacobi|innerpcg
  double psolve_rtol  = 1.0e-2;

  // nonlinear solver
  int maxcor     = -1;
  double nlscoef = -1.0;

  // FIRKODE extras
  std::string predictor; // extrapolate|trivial ("" = default)
  bool refilter      = true;
  bool newton_safety = true;
  std::string controller; // "" = default
  double impgus_k1 = -1.0, impgus_k2 = -1.0;
  std::string mode_name = "beginner"; // beginner|expert

  // CVODE extras
  bool stald = false; // BDF stability-limit detection

  // ARKODE extras
  bool arkode_linear = false;
  int ark_predictor  = -1;

  // limits and interval
  long maxsteps = 10000000;
  int maxnef    = -1;
  double h0     = 0.0;
  double tf     = -1.0; // override final time (negative = problem default)
  int grid      = 1;    // number of output times
  double grid_geom_t1   = 0.0; // geometric output grid: t1 * mult^k
  double grid_geom_mult = 0.0;
  bool tstop            = true;

  // modes
  int restart_every = 1; // split [t0, tf] into K intervals, re-init in between
  bool order_hist   = false; // step-by-step order histogram (CVODE, IDA)
  double fixed_h    = 0.0;
  double h_one      = 0.0;
  int dense_k       = 0;

  // Jacobian dump (spectrum analysis with scripts/spectrum.py)
  std::vector<double> dump_times; // --dump-jac t1,t2,...
  std::string dump_dir = ".";
  bool dump_dq         = false; // forward differences even if jac() exists

  // reference
  std::string ref_file;
  std::string write_ref;
  bool ref_inprocess         = false;
  bool ref_check             = false;
  std::string ref_integrator = "firkode";
  std::string ref_check_integrator; // "" = cvode (ida with a mass matrix)
  double ref_rtol = 1.0e-13;

  // repetitions and output
  bool profile   = false; // time setup/solve/rhs/jac (timings not comparable)
  int reps       = 1;
  bool no_header = false;
  bool verbose   = false;
  bool version   = false;

  // problem parameters
  int mx          = 0; // diurnal / heat2d grid
  int my          = 0;
  long n          = 0;      // pr / heat1dfem size
  double velocity = 0.0;    // heat2d advection velocity
  double lambda   = -1.0e4; // pr
  double eps      = 1.0e-6; // vdpol
  // every other "--key value" pair; the problem declares the keys it reads
  std::map<std::string, std::string> extra;
};

// Problem parameters from the generic pairs, with defaults
inline double param(const Options& o, const std::string& key, double dflt)
{
  auto it = o.extra.find(key);
  return it == o.extra.end() ? dflt : std::stod(it->second);
}

inline long param(const Options& o, const std::string& key, long dflt)
{
  auto it = o.extra.find(key);
  return it == o.extra.end() ? dflt : std::stol(it->second);
}

inline std::string param(const Options& o, const std::string& key,
                         const std::string& dflt)
{
  auto it = o.extra.find(key);
  return it == o.extra.end() ? dflt : it->second;
}

namespace detail {

inline bool take(std::vector<std::string>& args, const std::string& key,
                 std::string& dest)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it == args.end()) { return false; }
  if (it + 1 == args.end())
  {
    throw std::runtime_error("missing value for " + key);
  }
  dest = *(it + 1);
  args.erase(it, it + 2);
  return true;
}

inline void take(std::vector<std::string>& args, const std::string& key,
                 double& dest)
{
  std::string v;
  if (take(args, key, v)) { dest = std::stod(v); }
}

inline void take(std::vector<std::string>& args, const std::string& key, int& dest)
{
  std::string v;
  if (take(args, key, v)) { dest = std::stoi(v); }
}

inline void take(std::vector<std::string>& args, const std::string& key,
                 long& dest)
{
  std::string v;
  if (take(args, key, v)) { dest = std::stol(v); }
}

// Two values following the key (e.g. --impgus-k K1 K2)
inline bool take2(std::vector<std::string>& args, const std::string& key,
                  double& a, double& b)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it == args.end()) { return false; }
  if (it + 2 >= args.end())
  {
    throw std::runtime_error("missing values for " + key);
  }
  a = std::stod(*(it + 1));
  b = std::stod(*(it + 2));
  args.erase(it, it + 3);
  return true;
}

inline bool flag(std::vector<std::string>& args, const std::string& key)
{
  auto it = std::find(args.begin(), args.end(), key);
  if (it == args.end()) { return false; }
  args.erase(it);
  return true;
}

} // namespace detail

inline void usage(std::FILE* fp)
{
  std::
    fprintf(fp,
            "Usage: firkode_benchmark --problem P --integrator I [options]\n"
            "  P: "
            "rober|vdpol|orego|hires|analytic|pr|kpr|diurnal|heat2d|heat1dfem|"
            "sector|wave1dfem|pollu|ringmod|beam|plate|bruss3d\n"
            "  I: firkode|cvode|arkode|ida\n"
            "Tolerances: --rtol R --atol A --atol-vec a1,a2,...\n"
            "Method:     --stages S (firkode, 3 only) --order Q (arkode DIRK "
            "order; cvode/ida max order)\n"
            "Linear:     --ls dense|band|spgmr|pcg --prec none|left|right "
            "--maxl M\n"
            "            --maxrs R --epslin E --msbp N --msbj N --dgmax D "
            "--dgmax-jbad D\n"
            "            --reuse-policy radau5|periodic (firkode; msbp, msbj, "
            "dgmax act under periodic only)\n"
            "            --stage-maxl M --stage-maxrs R --stage-epslin E "
            "(firkode)\n"
            "            --psolve jacobi|innerpcg --psolve-rtol R (matrix-free "
            "heat2d)\n"
            "Nonlinear:  --maxcor N --nlscoef C\n"
            "FIRKODE:    --predictor extrapolate|trivial --no-refilter "
            "--no-newton-safety\n"
            "            --controller NAME --impgus-k K1 K2 --mode "
            "beginner|expert\n"
            "CVODE:      --stald (stability-limit detection)\n"
            "ARKODE:     --arkode-linear --ark-predictor P\n"
            "Limits:     --maxsteps N --maxnef N --h0 H --tf T --grid K\n"
            "            --grid-geom T1 MULT --no-tstop\n"
            "Modes:      --fixed-h H | --one-step --h H | --roots | --dense-k "
            "K\n"
            "            --restart-every K (K equal intervals with a re-init "
            "between them)\n"
            "            --order-hist (BDF order histogram; untimed, "
            "cvode/ida)\n"
            "Dump:       --dump-jac t1,t2,... [--dump-dir DIR] "
            "[--dump-jac-dq]\n"
            "Reference:  --ref FILE | --write-ref FILE | --ref-inprocess "
            "[--ref-check]\n"
            "            --ref-integrator firkode|cvode|ida "
            "--ref-check-integrator I --ref-rtol R\n"
            "Problem:    --mx MX --my MY --n N --velocity V --lambda L --eps "
            "E\n"
            "            --KEY VALUE (further parameters a problem declares)\n"
            "Output:     --reps R --no-header --verbose --version\n"
            "            --profile (time factorizations, solves, rhs and "
            "Jacobians)\n");
}

inline Options parse_args(int argc, char* argv[])
{
  Options o;
  std::vector<std::string> args(argv + 1, argv + argc);
  using namespace detail;

  if (flag(args, "--help") || flag(args, "-h"))
  {
    usage(stdout);
    std::exit(0);
  }
  o.version = flag(args, "--version");

  take(args, "--problem", o.problem);
  take(args, "--integrator", o.integrator);
  take(args, "--rtol", o.rtol);
  take(args, "--atol", o.atol);
  std::string s;
  if (take(args, "--atol-vec", s))
  {
    size_t pos = 0;
    while (pos <= s.size())
    {
      size_t next = s.find(',', pos);
      if (next == std::string::npos) { next = s.size(); }
      o.atol_vec.push_back(std::stod(s.substr(pos, next - pos)));
      pos = next + 1;
    }
  }
  take(args, "--stages", o.stages);
  take(args, "--order", o.order);
  take(args, "--ls", o.ls);
  take(args, "--prec", o.prec);
  take(args, "--maxl", o.maxl);
  take(args, "--maxrs", o.maxrs);
  take(args, "--epslin", o.epslin);
  take(args, "--msbp", o.msbp);
  take(args, "--msbj", o.msbj);
  take(args, "--dgmax", o.dgmax);
  take(args, "--dgmax-jbad", o.dgmax_jbad);
  take(args, "--reuse-policy", o.reuse_policy);
  take(args, "--stage-maxl", o.stage_maxl);
  take(args, "--stage-maxrs", o.stage_maxrs);
  take(args, "--stage-epslin", o.stage_epslin);
  take(args, "--psolve", o.psolve);
  take(args, "--psolve-rtol", o.psolve_rtol);
  take(args, "--maxcor", o.maxcor);
  take(args, "--nlscoef", o.nlscoef);
  take(args, "--predictor", o.predictor);
  if (flag(args, "--no-refilter")) { o.refilter = false; }
  if (flag(args, "--no-newton-safety")) { o.newton_safety = false; }
  take(args, "--controller", o.controller);
  take2(args, "--impgus-k", o.impgus_k1, o.impgus_k2);
  take(args, "--mode", o.mode_name);
  o.stald         = flag(args, "--stald");
  o.arkode_linear = flag(args, "--arkode-linear");
  take(args, "--ark-predictor", o.ark_predictor);
  take(args, "--maxsteps", o.maxsteps);
  take(args, "--maxnef", o.maxnef);
  take(args, "--h0", o.h0);
  take(args, "--tf", o.tf);
  take(args, "--grid", o.grid);
  take2(args, "--grid-geom", o.grid_geom_t1, o.grid_geom_mult);
  if (flag(args, "--no-tstop")) { o.tstop = false; }

  o.order_hist = flag(args, "--order-hist");
  take(args, "--restart-every", o.restart_every);
  if (o.restart_every < 1)
  {
    throw std::runtime_error("--restart-every needs K >= 1");
  }
  take(args, "--fixed-h", o.fixed_h);
  bool one_step = flag(args, "--one-step");
  take(args, "--h", o.h_one);
  take(args, "--dense-k", o.dense_k);
  bool roots = flag(args, "--roots");
  int nmodes = (o.fixed_h > 0.0) + one_step + roots + (o.dense_k > 0);
  if (nmodes > 1) { throw std::runtime_error("modes are mutually exclusive"); }
  if (o.fixed_h > 0.0) { o.mode = Mode::FixedStep; }
  else if (one_step)
  {
    if (o.h_one <= 0.0) { throw std::runtime_error("--one-step needs --h H"); }
    o.mode = Mode::OneStep;
  }
  else if (roots) { o.mode = Mode::Roots; }
  else if (o.dense_k > 0) { o.mode = Mode::Dense; }

  if (take(args, "--dump-jac", s))
  {
    size_t pos = 0;
    while (pos <= s.size())
    {
      size_t next = s.find(',', pos);
      if (next == std::string::npos) { next = s.size(); }
      o.dump_times.push_back(std::stod(s.substr(pos, next - pos)));
      pos = next + 1;
    }
  }
  take(args, "--dump-dir", o.dump_dir);
  o.dump_dq = flag(args, "--dump-jac-dq");
  take(args, "--ref", o.ref_file);
  take(args, "--write-ref", o.write_ref);
  o.ref_inprocess = flag(args, "--ref-inprocess");
  o.ref_check     = flag(args, "--ref-check");
  take(args, "--ref-integrator", o.ref_integrator);
  take(args, "--ref-check-integrator", o.ref_check_integrator);
  take(args, "--ref-rtol", o.ref_rtol);

  o.profile = flag(args, "--profile");
  take(args, "--reps", o.reps);
  o.no_header = flag(args, "--no-header");
  o.verbose   = flag(args, "--verbose");

  take(args, "--mx", o.mx);
  take(args, "--my", o.my);
  take(args, "--n", o.n);
  take(args, "--velocity", o.velocity);
  take(args, "--lambda", o.lambda);
  take(args, "--eps", o.eps);

  // leftover "--key value" pairs are problem parameters, validated against
  // the chosen problem's keys by the harness
  while (args.size() >= 2 && args[0].rfind("--", 0) == 0 &&
         args[1].rfind("--", 0) != 0)
  {
    o.extra[args[0].substr(2)] = args[1];
    args.erase(args.begin(), args.begin() + 2);
  }
  if (!args.empty())
  {
    throw std::runtime_error("unrecognized argument: " + args.front());
  }
  if (!o.version && (o.problem.empty() || o.integrator.empty()))
  {
    throw std::runtime_error("--problem and --integrator are required");
  }
  if (o.reps < 1) { o.reps = 1; }
  return o;
}

} // namespace benchmark

#endif // FIRKODE_BENCHMARK_CLI_HPP_
