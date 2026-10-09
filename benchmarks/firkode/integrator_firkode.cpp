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
 * FIRKODE (Radau IIA) integrator wrapper of the benchmark harness.
 * -----------------------------------------------------------------*/

#include <firkode/firkode.h>
#include <firkode/firkode_ls.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>

#include "integrator.hpp"

namespace benchmark {

namespace {

SUNAdaptController make_controller(const std::string& name, SUNContext ctx)
{
  if (name == "impgus") { return SUNAdaptController_ImpGus(ctx); }
  if (name == "expgus") { return SUNAdaptController_ExpGus(ctx); }
  if (name == "i") { return SUNAdaptController_I(ctx); }
  if (name == "pi") { return SUNAdaptController_PI(ctx); }
  if (name == "pid") { return SUNAdaptController_PID(ctx); }
  if (name == "soderlind") { return SUNAdaptController_Soderlind(ctx); }
  if (name == "h0211") { return SUNAdaptController_H0211(ctx); }
  if (name == "h0321") { return SUNAdaptController_H0321(ctx); }
  if (name == "h211") { return SUNAdaptController_H211(ctx); }
  if (name == "h312") { return SUNAdaptController_H312(ctx); }
  throw std::runtime_error("unknown controller: " + name);
}

class FirkodeIntegrator : public Integrator
{
public:
  ~FirkodeIntegrator() override { teardown(); }

  void setup(Problem& p, const Options& o, N_Vector y, SUNContext ctx) override
  {
    const ProblemSpec& sp = p.spec();
    p_                    = &p;
    h0_                   = (o.mode == Mode::Evolve) ? o.h0 : 0.0;

    // method: this branch supports the three-stage method only
    if (o.stages > 0 && o.stages != FIRK_DEFAULT_STAGES)
    {
      throw std::runtime_error("unsupported: FIRKODE supports only s = 3 on "
                               "this branch (--stages " +
                               std::to_string(o.stages) + ")");
    }
    if (o.order > 0 && o.order != 5)
    {
      throw std::runtime_error("unsupported: FIRKODE order is 5 on this "
                               "branch (--order " +
                               std::to_string(o.order) + ")");
    }
    s_ = FIRK_DEFAULT_STAGES;

    mem_ = FIRKodeCreate(ctx);
    check_ptr(mem_, "FIRKodeCreate");
    check(FIRKodeInit(mem_, p.rhs(), sp.t0, y), "FIRKodeInit");
    check(FIRKodeSetUserData(mem_, p.user_data()), "FIRKodeSetUserData");

    if (!o.atol_vec.empty())
    {
      N_Vector avec = N_VClone(y);
      for (sunindextype i = 0; i < sp.N; i++)
      {
        NV_Ith_S(avec,
                 i) = (sunrealtype)o.atol_vec[(size_t)i % o.atol_vec.size()];
      }
      check(FIRKodeSVtolerances(mem_, (sunrealtype)o.rtol, avec),
            "FIRKodeSVtolerances");
      N_VDestroy(avec);
    }
    else
    {
      check(FIRKodeSStolerances(mem_, (sunrealtype)o.rtol, (sunrealtype)o.atol),
            "FIRKodeSStolerances");
    }

    // EXPERT mode: requested explicitly, or implied by any option that is
    // rejected in BEGINNER mode. EXPERT mode also changes the step-size
    // formula, so the mode is recorded in the CSV row.
    if (o.mode_name != "beginner" && o.mode_name != "expert")
    {
      throw std::runtime_error("unknown --mode: " + o.mode_name);
    }
    expert_ = o.mode_name == "expert" || o.maxnef > 0 || !o.predictor.empty() ||
              !o.refilter || !o.newton_safety || !o.controller.empty() ||
              o.impgus_k1 > 0.0 || o.msbp > 0 || o.msbj > 0 || o.maxcor > 0 ||
              o.nlscoef > 0.0 || o.stage_maxl > 0 || o.stage_maxrs >= 0 ||
              o.stage_epslin > 0.0 || !o.reuse_policy.empty() ||
              o.dgmax > 0.0 || o.dgmax_jbad > 0.0;
    if (expert_) { check(FIRKodeSetMode(mem_, FIRK_EXPERT), "FIRKodeSetMode"); }

    // limits (BEGINNER options)
    check(FIRKodeSetMaxNumSteps(mem_, o.maxsteps), "FIRKodeSetMaxNumSteps");
    if (o.tstop)
    {
      check(FIRKodeSetStopTime(mem_, (sunrealtype)o.tf), "FIRKodeSetStopTime");
    }
    if (o.mode == Mode::FixedStep)
    {
      check(FIRKodeSetFixedStep(mem_, (sunrealtype)o.fixed_h),
            "FIRKodeSetFixedStep");
    }
    else if (o.mode == Mode::OneStep)
    {
      check(FIRKodeSetInitStep(mem_, (sunrealtype)o.h_one),
            "FIRKodeSetInitStep");
    }
    else if (o.h0 > 0.0)
    {
      check(FIRKodeSetInitStep(mem_, (sunrealtype)o.h0), "FIRKodeSetInitStep");
    }

    // EXPERT integrator options
    if (o.maxnef > 0)
    {
      check(FIRKodeSetMaxErrTestFails(mem_, o.maxnef),
            "FIRKodeSetMaxErrTestFails");
    }
    if (!o.predictor.empty())
    {
      int method = o.predictor == "trivial" ? FIRK_PREDICT_TRIVIAL
                                            : FIRK_PREDICT_EXTRAPOLATE;
      check(FIRKodeSetPredictorMethod(mem_, method),
            "FIRKodeSetPredictorMethod");
    }
    if (!o.refilter)
    {
      check(FIRKodeSetErrorRefilter(mem_, SUNFALSE), "FIRKodeSetErrorRefilter");
    }
    if (!o.newton_safety)
    {
      check(FIRKodeSetNewtonCountSafety(mem_, SUNFALSE),
            "FIRKodeSetNewtonCountSafety");
    }
    if (!o.controller.empty() || o.impgus_k1 > 0.0)
    {
      // the default is ImpGus with both gains 1; --controller impgus gives
      // the stock SUNDIALS gains unless --impgus-k is given too
      ctrl_ = make_controller(o.controller.empty() ? "impgus" : o.controller,
                              ctx);
      check_ptr(ctrl_, "SUNAdaptController");
      if (o.impgus_k1 > 0.0)
      {
        check(SUNAdaptController_SetParams_ImpGus(ctrl_, (sunrealtype)o.impgus_k1,
                                                  (sunrealtype)o.impgus_k2),
              "SUNAdaptController_SetParams_ImpGus");
      }
      check(FIRKodeSetAdaptController(mem_, ctrl_),
            "FIRKodeSetAdaptController");
    }
    if (!o.reuse_policy.empty())
    {
      if (o.reuse_policy != "radau5" && o.reuse_policy != "periodic")
      {
        throw std::runtime_error("unknown --reuse-policy: " + o.reuse_policy);
      }
      periodic_ = o.reuse_policy == "periodic";
      check(FIRKodeSetReusePolicy(mem_, periodic_ ? FIRK_POLICY_PERIODIC
                                                  : FIRK_POLICY_RADAU5),
            "FIRKodeSetReusePolicy");
    }
    // msbp, msbj, dgmax and dgmax_jbad act under the periodic policy only
    if (o.msbp > 0)
    {
      check(FIRKodeSetLSetupFrequency(mem_, o.msbp),
            "FIRKodeSetLSetupFrequency");
    }
    if (o.dgmax > 0.0)
    {
      check(FIRKodeSetDeltaGammaMaxLSetup(mem_, (sunrealtype)o.dgmax),
            "FIRKodeSetDeltaGammaMaxLSetup");
    }
    if (o.maxcor > 0)
    {
      check(FIRKodeSetMaxNonlinIters(mem_, o.maxcor),
            "FIRKodeSetMaxNonlinIters");
    }
    if (o.nlscoef > 0.0)
    {
      check(FIRKodeSetNonlinConvCoef(mem_, (sunrealtype)o.nlscoef),
            "FIRKodeSetNonlinConvCoef");
    }

    // linear solver
    const LinSolKind kind = parse_ls(o.ls, sp.ls_default);
    ls_                   = make_linsys(kind, sp, y, o, ctx);
    check(FIRKodeSetLinearSolver(mem_, ls_.LS, ls_.A),
          "FIRKodeSetLinearSolver");
    const bool direct = is_direct(kind);
    if (kind == LinSolKind::KLU && p.jac() == nullptr)
    {
      throw std::runtime_error("unsupported: KLU needs a user Jacobian");
    }
    if (direct)
    {
      if (p.jac() != nullptr)
      {
        check(FIRKodeSetJacFn(mem_, p.jac()), "FIRKodeSetJacFn");
      }
    }
    else
    {
      if (p.jtimes() != nullptr)
      {
        check(FIRKodeSetJacTimes(mem_, nullptr, p.jtimes()),
              "FIRKodeSetJacTimes");
      }
      if (sp.has_prec && ls_.prec != "none")
      {
        check(FIRKodeSetPreconditioner(mem_, p.psetup(), p.psolve()),
              "FIRKodeSetPreconditioner");
      }
    }
    if (o.epslin > 0.0)
    {
      check(FIRKodeSetEpsLin(mem_, (sunrealtype)o.epslin), "FIRKodeSetEpsLin");
    }

    // EXPERT linear-solver options: FIRKodeSetLinearSolver resets them, so
    // they are set after the solver is attached
    if (o.msbj > 0)
    {
      check(FIRKodeSetJacEvalFrequency(mem_, o.msbj),
            "FIRKodeSetJacEvalFrequency");
    }
    if (o.dgmax_jbad > 0.0)
    {
      check(FIRKodeSetDeltaGammaMaxBadJac(mem_, (sunrealtype)o.dgmax_jbad),
            "FIRKodeSetDeltaGammaMaxBadJac");
    }
    if (o.stage_maxl > 0)
    {
      check(FIRKodeSetStageSolverMaxl(mem_, o.stage_maxl),
            "FIRKodeSetStageSolverMaxl");
    }
    if (o.stage_maxrs >= 0)
    {
      check(FIRKodeSetStageSolverMaxRestarts(mem_, o.stage_maxrs),
            "FIRKodeSetStageSolverMaxRestarts");
    }
    if (o.stage_epslin > 0.0)
    {
      check(FIRKodeSetStageEpsLin(mem_, (sunrealtype)o.stage_epslin),
            "FIRKodeSetStageEpsLin");
    }

    // mass matrix
    if (sp.has_mass)
    {
      mls_ = make_mass_linsys(sp, y, o, ctx);
      check(FIRKodeSetMassLinearSolver(mem_, mls_.LS, mls_.A, SUNFALSE),
            "FIRKodeSetMassLinearSolver");
      check(FIRKodeSetMassFn(mem_, p.mass()), "FIRKodeSetMassFn");
    }

    // rootfinding
    if (o.mode == Mode::Roots && sp.nroots > 0)
    {
      check(FIRKodeRootInit(mem_, sp.nroots, p.roots()), "FIRKodeRootInit");
      roots_ = true;
    }
  }

  void reinit(sunrealtype t, N_Vector y) override
  {
    check(FIRKodeReInit(mem_, t, y), "FIRKodeReInit");
    if (h0_ > 0.0)
    {
      check(FIRKodeSetInitStep(mem_, (sunrealtype)h0_), "FIRKodeSetInitStep");
    }
  }

  void set_stop_time(sunrealtype tstop) override
  {
    check(FIRKodeSetStopTime(mem_, tstop), "FIRKodeSetStopTime");
  }

  int evolve(sunrealtype tout, N_Vector y, sunrealtype& tret, bool one_step) override
  {
    return FIRKodeEvolve(mem_, tout, y, &tret,
                         one_step ? FIRK_ONE_STEP : FIRK_NORMAL);
  }

  bool is_success(int flag) const override
  {
    return flag == FIRK_SUCCESS || flag == FIRK_TSTOP_RETURN;
  }

  bool is_root(int flag) const override { return flag == FIRK_ROOT_RETURN; }

  void root_info(std::vector<int>& dirs) override
  {
    dirs.assign((size_t)p_->spec().nroots, 0);
    check(FIRKodeGetRootInfo(mem_, dirs.data()), "FIRKodeGetRootInfo");
  }

  int get_dky(sunrealtype t, int k, N_Vector dky) override
  {
    return FIRKodeGetDky(mem_, t, k, dky);
  }

  int est_local_errors(N_Vector ele) override
  {
    return FIRKodeGetEstLocalErrors(mem_, ele);
  }

  int err_weights(N_Vector ewt) override
  {
    return FIRKodeGetErrWeights(mem_, ewt);
  }

  sunrealtype last_step() override
  {
    sunrealtype h = 0;
    FIRKodeGetLastStep(mem_, &h);
    return h;
  }

  void harvest(RunStats& s) override
  {
    FIRKodeGetNumSteps(mem_, &s.nst);
    FIRKodeGetNumStepAttempts(mem_, &s.nst_att);
    FIRKodeGetNumRhsEvals(mem_, &s.nfe);
    FIRKodeGetNumLinSolvSetups(mem_, &s.nsetups);
    FIRKodeGetNumErrTestFails(mem_, &s.netf);
    FIRKodeGetNumStepSolveFails(mem_, &s.nsf);
    FIRKodeGetNumNonlinSolvIters(mem_, &s.nni);
    FIRKodeGetNumNonlinSolvConvFails(mem_, &s.nncf);
    FIRKodeGetNumJacEvals(mem_, &s.nje);
    FIRKodeGetNumLinRhsEvals(mem_, &s.nfe_ls);
    FIRKodeGetNumPrecEvals(mem_, &s.npe);
    FIRKodeGetNumPrecSolves(mem_, &s.nps);
    FIRKodeGetNumLinIters(mem_, &s.nli);
    FIRKodeGetNumLinConvFails(mem_, &s.nlcf);
    FIRKodeGetNumJtimesEvals(mem_, &s.njtv);
    FIRKodeGetNumBlockSolves(mem_, &s.nbs);
    FIRKodeGetNumStageLinIters(mem_, &s.nsli);
    FIRKodeGetNumStageLinConvFails(mem_, &s.nscf);
    if (p_->spec().has_mass)
    {
      FIRKodeGetNumMassSetups(mem_, &s.nmass_setups);
      FIRKodeGetNumMassSolves(mem_, &s.nmass_solves);
      FIRKodeGetNumMassMult(mem_, &s.nmass_mult);
    }
    if (roots_) { FIRKodeGetNumGEvals(mem_, &s.ngeval); }
  }

  std::string name() const override { return "firkode"; }

  std::string method() const override
  {
    return "RadauIIA-s" + std::to_string(s_) + (periodic_ ? "+periodic" : "");
  }

  int order() const override { return 2 * s_ - 1; }

  int stages() const override { return s_; }

  std::string mode_name() const override
  {
    return expert_ ? "expert" : "beginner";
  }

  std::string ls_name() const override { return to_string(ls_.kind); }

  std::string prec_name() const override { return ls_.prec; }

  void teardown() override
  {
    roots_ = false;
    if (mem_ != nullptr) { FIRKodeFree(&mem_); }
    ls_.destroy();
    mls_.destroy();
    if (ctrl_ != nullptr)
    {
      SUNAdaptController_Destroy(ctrl_);
      ctrl_ = nullptr;
    }
  }

private:
  void* mem_ = nullptr;
  LinSys ls_, mls_;
  SUNAdaptController ctrl_ = nullptr;
  Problem* p_              = nullptr;
  bool roots_              = false;
  bool expert_             = false;
  bool periodic_           = false;
  double h0_               = 0.0;
  int s_                   = FIRK_DEFAULT_STAGES;
};

} // namespace

std::unique_ptr<Integrator> make_firkode_integrator()
{
  return std::unique_ptr<Integrator>(new FirkodeIntegrator());
}

} // namespace benchmark
