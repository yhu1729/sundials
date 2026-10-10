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
 * ARKODE LSRKStep adapter (arkode-lsrk): the explicit stabilized
 * methods RKC and RKL (order 2) on f, with the problem's bound of the
 * dominant eigenvalue (--dom-eig user) or the power iteration on
 * difference quotients of f (--dom-eig power). No linear algebra and no
 * mass matrix.
 * -----------------------------------------------------------------*/

#include <cmath>

#include <arkode/arkode.h>
#include <arkode/arkode_lsrkstep.h>
#include <sundomeigest/sundomeigest_power.h>

#include "integrator.hpp"

namespace benchmark {

namespace {

// Settings of the power iteration (those of
// examples/arkode/C_serial/ark_analytic_lsrk_domeigest.c)
constexpr long POWER_MAX_ITERS      = 100;
constexpr double POWER_REL_TOL      = 5.0e-3;
constexpr int POWER_INIT_PREPROCESS = 10;

class LsrkIntegrator : public Integrator
{
public:
  ~LsrkIntegrator() override { teardown(); }

  void setup(Problem& p, const Options& o, N_Vector y, SUNContext ctx) override
  {
    const ProblemSpec& sp = p.spec();
    p_                    = &p;
    h0_                   = (o.mode == Mode::Evolve) ? o.h0 : 0.0;
    if (sp.has_mass)
    {
      throw std::runtime_error("unsupported: LSRKStep needs M = I");
    }

    mem_ = LSRKStepCreateSTS(p.rhs(), sp.t0, y, ctx);
    check_ptr(mem_, "LSRKStepCreateSTS");
    check(ARKodeSetUserData(mem_, p.user_data()), "ARKodeSetUserData");

    if (!o.atol_vec.empty())
    {
      N_Vector avec = N_VClone(y);
      for (sunindextype i = 0; i < sp.N; i++)
      {
        NV_Ith_S(avec,
                 i) = (sunrealtype)o.atol_vec[(size_t)i % o.atol_vec.size()];
      }
      check(ARKodeSVtolerances(mem_, (sunrealtype)o.rtol, avec),
            "ARKodeSVtolerances");
      N_VDestroy(avec);
    }
    else
    {
      check(ARKodeSStolerances(mem_, (sunrealtype)o.rtol, (sunrealtype)o.atol),
            "ARKodeSStolerances");
    }

    if (o.lsrk_method == "rkc2")
    {
      check(LSRKStepSetSTSMethod(mem_, ARKODE_LSRK_RKC_2),
            "LSRKStepSetSTSMethod");
      method_ = "ARKODE_LSRK_RKC_2";
    }
    else if (o.lsrk_method == "rkl2")
    {
      check(LSRKStepSetSTSMethod(mem_, ARKODE_LSRK_RKL_2),
            "LSRKStepSetSTSMethod");
      method_ = "ARKODE_LSRK_RKL_2";
    }
    else
    {
      throw std::runtime_error("unknown --lsrk-method: " + o.lsrk_method);
    }

    // the spectrum is constant when the bound says so (user) or the
    // problem is linear with a constant Jacobian (power): one estimate
    bool constant = false;
    if (o.dom_eig == "user")
    {
      if (p.dom_eig() == nullptr)
      {
        throw std::runtime_error("unsupported: problem " + sp.name +
                                 " has no eigenvalue bound");
      }
      // DomEigFn repeats the signature of ARKDomEigFn
      check(LSRKStepSetDomEigFn(mem_, p.dom_eig()), "LSRKStepSetDomEigFn");
      constant = sp.dom_eig_const;
    }
    else if (o.dom_eig == "power")
    {
      // a deterministic start vector (the example draws it with rand())
      N_Vector q      = N_VClone(y);
      sunrealtype* qd = N_VGetArrayPointer(q);
      for (sunindextype i = 0; i < sp.N; i++)
      {
        qd[i] = SUN_RCONST(1.0) +
                SUN_RCONST(0.5) * std::sin((sunrealtype)(i + 1));
      }
      dee_ = SUNDomEigEstimator_Power(q, POWER_MAX_ITERS,
                                      (sunrealtype)POWER_REL_TOL, ctx);
      N_VDestroy(q);
      check_ptr(dee_, "SUNDomEigEstimator_Power");
      check(LSRKStepSetDomEigEstimator(mem_, dee_),
            "LSRKStepSetDomEigEstimator");
      check(LSRKStepSetNumDomEigEstInitPreprocessIters(mem_,
                                                       POWER_INIT_PREPROCESS),
            "LSRKStepSetNumDomEigEstInitPreprocessIters");
      constant = sp.linear;
    }
    else { throw std::runtime_error("unknown --dom-eig: " + o.dom_eig); }
    method_ += "+" + o.dom_eig;
    // never negative: LSRKStepSetDomEigFrequency stores a negative value
    // instead of resetting to the default (see BENCHMARK_ARKODE.md)
    const long freq = (o.dom_eig_freq >= 0) ? o.dom_eig_freq
                      : constant            ? 0L
                                            : 25L;
    check(LSRKStepSetDomEigFrequency(mem_, freq), "LSRKStepSetDomEigFrequency");
    if (o.lsrk_max_stages > 0)
    {
      check(LSRKStepSetMaxNumStages(mem_, o.lsrk_max_stages),
            "LSRKStepSetMaxNumStages");
    }
    if (o.dom_eig_safety > 0.0)
    {
      check(LSRKStepSetDomEigSafetyFactor(mem_, (sunrealtype)o.dom_eig_safety),
            "LSRKStepSetDomEigSafetyFactor");
    }

    check(ARKodeSetMaxNumSteps(mem_, o.maxsteps), "ARKodeSetMaxNumSteps");
    if (o.maxnef > 0)
    {
      check(ARKodeSetMaxErrTestFails(mem_, o.maxnef),
            "ARKodeSetMaxErrTestFails");
    }
    if (o.tstop)
    {
      check(ARKodeSetStopTime(mem_, (sunrealtype)o.tf), "ARKodeSetStopTime");
    }
    if (o.mode == Mode::FixedStep)
    {
      check(ARKodeSetFixedStep(mem_, (sunrealtype)o.fixed_h),
            "ARKodeSetFixedStep");
    }
    else if (o.mode == Mode::OneStep)
    {
      check(ARKodeSetInitStep(mem_, (sunrealtype)o.h_one), "ARKodeSetInitStep");
    }
    else if (o.h0 > 0.0)
    {
      check(ARKodeSetInitStep(mem_, (sunrealtype)o.h0), "ARKodeSetInitStep");
    }
    if (!o.controller.empty())
    {
      check(ARKodeSetAdaptControllerByName(mem_,
                                           arkode_controller_name(o.controller)
                                             .c_str()),
            "ARKodeSetAdaptControllerByName");
    }
    if (o.mode == Mode::Roots && sp.nroots > 0)
    {
      check(ARKodeRootInit(mem_, sp.nroots, p.roots()), "ARKodeRootInit");
      roots_ = true;
    }
  }

  void reinit(sunrealtype t, N_Vector y) override
  {
    check(LSRKStepReInitSTS(mem_, p_->rhs(), t, y), "LSRKStepReInitSTS");
    if (h0_ > 0.0)
    {
      check(ARKodeSetInitStep(mem_, (sunrealtype)h0_), "ARKodeSetInitStep");
    }
  }

  void set_stop_time(sunrealtype tstop) override
  {
    check(ARKodeSetStopTime(mem_, tstop), "ARKodeSetStopTime");
  }

  int evolve(sunrealtype tout, N_Vector y, sunrealtype& tret, bool one_step) override
  {
    return ARKodeEvolve(mem_, tout, y, &tret,
                        one_step ? ARK_ONE_STEP : ARK_NORMAL);
  }

  bool is_success(int flag) const override
  {
    return flag == ARK_SUCCESS || flag == ARK_TSTOP_RETURN;
  }

  bool is_root(int flag) const override { return flag == ARK_ROOT_RETURN; }

  void root_info(std::vector<int>& dirs) override
  {
    dirs.assign((size_t)p_->spec().nroots, 0);
    check(ARKodeGetRootInfo(mem_, dirs.data()), "ARKodeGetRootInfo");
  }

  int get_dky(sunrealtype t, int k, N_Vector dky) override
  {
    return ARKodeGetDky(mem_, t, k, dky);
  }

  int est_local_errors(N_Vector ele) override
  {
    return ARKodeGetEstLocalErrors(mem_, ele);
  }

  int err_weights(N_Vector ewt) override
  {
    return ARKodeGetErrWeights(mem_, ewt);
  }

  sunrealtype last_step() override
  {
    sunrealtype h = 0;
    ARKodeGetLastStep(mem_, &h);
    return h;
  }

  void harvest(RunStats& s) override
  {
    ARKodeGetNumSteps(mem_, &s.nst);
    ARKodeGetNumStepAttempts(mem_, &s.nst_att);
    ARKodeGetNumRhsEvals(mem_, 0, &s.nfe);
    ARKodeGetNumErrTestFails(mem_, &s.netf);
    ARKodeGetNumExpSteps(mem_, &s.nst_expstab);
    ARKodeGetNumAccSteps(mem_, &s.nst_acc);
    LSRKStepGetNumDomEigUpdates(mem_, &s.ndomeig);
    if (dee_ != nullptr)
    {
      LSRKStepGetNumDomEigEstRhsEvals(mem_, &s.nfe_domeig);
    }
    int smax = 0;
    if (LSRKStepGetMaxNumStages(mem_, &smax) == ARK_SUCCESS)
    {
      s.lsrk_max_stages = smax;
      stages_           = smax;
    }
    if (roots_) { ARKodeGetNumGEvals(mem_, &s.ngeval); }
  }

  std::string name() const override { return "arkode-lsrk"; }

  std::string method() const override { return method_; }

  int order() const override { return 2; }

  // the most stages of a step
  int stages() const override { return stages_; }

  std::string ls_name() const override { return "none"; }

  std::string prec_name() const override { return "none"; }

  void teardown() override
  {
    roots_ = false;
    if (mem_ != nullptr) { ARKodeFree(&mem_); }
    // the estimator belongs to the caller and outlives the integrator
    if (dee_ != nullptr) { SUNDomEigEstimator_Destroy(&dee_); }
  }

private:
  void* mem_              = nullptr;
  SUNDomEigEstimator dee_ = nullptr;
  double h0_              = 0.0;
  Problem* p_             = nullptr;
  bool roots_             = false;
  int stages_             = 0;
  std::string method_;
};

} // namespace

std::unique_ptr<Integrator> make_arkode_lsrk_integrator()
{
  return std::unique_ptr<Integrator>(new LsrkIntegrator());
}

} // namespace benchmark
