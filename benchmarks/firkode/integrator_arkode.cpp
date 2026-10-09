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
 * ARKODE (ARKStep, implicit DIRK/ESDIRK) integrator wrapper of the
 * benchmark harness. The problem callbacks are passed unchanged: the
 * ARKODE function-pointer types are structurally identical to the
 * FIRKODE ones.
 * -----------------------------------------------------------------*/

#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_butcher.h>
#include <arkode/arkode_ls.h>

#include "integrator.hpp"

namespace benchmark {

namespace {

const char* default_dirk_name(int order)
{
  switch (order)
  {
  case 1: return "ARKODE_BACKWARD_EULER_1_1";
  case 2: return "ARKODE_GKC21_ESDIRK_3_1_2";
  case 3: return "ARKODE_ESDIRK325L2SA_5_2_3";
  case 4: return "ARKODE_ESDIRK436L2SA_6_3_4";
  case 5: return "ARKODE_ESDIRK547L2SA2_7_4_5";
  default: return "ARKODE_DIRK_unknown";
  }
}

std::string controller_name(const std::string& lower)
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

class ArkodeIntegrator : public Integrator
{
public:
  ~ArkodeIntegrator() override { teardown(); }

  void setup(Problem& p, const Options& o, N_Vector y, SUNContext ctx) override
  {
    const ProblemSpec& sp = p.spec();
    p_                    = &p;
    h0_                   = (o.mode == Mode::Evolve) ? o.h0 : 0.0;

    mem_ = ARKStepCreate(nullptr, p.rhs(), sp.t0, y, ctx);
    check_ptr(mem_, "ARKStepCreate");
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

    order_ = (o.order > 0) ? o.order : 4;
    if (o.order > 0) { check(ARKodeSetOrder(mem_, o.order), "ARKodeSetOrder"); }
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
    if (o.arkode_linear && sp.linear)
    {
      check(ARKodeSetLinear(mem_, 0), "ARKodeSetLinear");
      linear_ = true;
    }
    if (o.ark_predictor >= 0)
    {
      check(ARKodeSetPredictorMethod(mem_, o.ark_predictor),
            "ARKodeSetPredictorMethod");
    }
    if (!o.controller.empty())
    {
      check(ARKodeSetAdaptControllerByName(mem_,
                                           controller_name(o.controller).c_str()),
            "ARKodeSetAdaptControllerByName");
    }
    if (o.msbp > 0)
    {
      check(ARKodeSetLSetupFrequency(mem_, (int)o.msbp),
            "ARKodeSetLSetupFrequency");
    }
    if (o.dgmax > 0.0)
    {
      check(ARKodeSetDeltaGammaMax(mem_, (sunrealtype)o.dgmax),
            "ARKodeSetDeltaGammaMax");
    }
    if (o.maxcor > 0)
    {
      check(ARKodeSetMaxNonlinIters(mem_, o.maxcor), "ARKodeSetMaxNonlinIters");
    }
    if (o.nlscoef > 0.0)
    {
      check(ARKodeSetNonlinConvCoef(mem_, (sunrealtype)o.nlscoef),
            "ARKodeSetNonlinConvCoef");
    }

    const LinSolKind kind = parse_ls(o.ls, sp.ls_default);
    ls_                   = make_linsys(kind, sp, y, o, ctx);
    check(ARKodeSetLinearSolver(mem_, ls_.LS, ls_.A), "ARKodeSetLinearSolver");
    const bool direct = is_direct(kind);
    if (kind == LinSolKind::KLU && p.jac() == nullptr)
    {
      throw std::runtime_error("unsupported: KLU needs a user Jacobian");
    }
    if (direct)
    {
      if (p.jac() != nullptr)
      {
        check(ARKodeSetJacFn(mem_, p.jac()), "ARKodeSetJacFn");
      }
    }
    else
    {
      if (p.jtimes() != nullptr)
      {
        check(ARKodeSetJacTimes(mem_, nullptr, p.jtimes()),
              "ARKodeSetJacTimes");
      }
      if (sp.has_prec && ls_.prec != "none")
      {
        check(ARKodeSetPreconditioner(mem_, p.psetup(), p.psolve()),
              "ARKodeSetPreconditioner");
      }
    }
    if (o.epslin > 0.0)
    {
      check(ARKodeSetEpsLin(mem_, (sunrealtype)o.epslin), "ARKodeSetEpsLin");
    }
    if (o.msbj > 0)
    {
      check(ARKodeSetJacEvalFrequency(mem_, o.msbj),
            "ARKodeSetJacEvalFrequency");
    }

    if (sp.has_mass)
    {
      mls_ = make_mass_linsys(sp, y, o, ctx);
      check(ARKodeSetMassLinearSolver(mem_, mls_.LS, mls_.A, SUNFALSE),
            "ARKodeSetMassLinearSolver");
      check(ARKodeSetMassFn(mem_, p.mass()), "ARKodeSetMassFn");
    }

    if (o.mode == Mode::Roots && sp.nroots > 0)
    {
      check(ARKodeRootInit(mem_, sp.nroots, p.roots()), "ARKodeRootInit");
      roots_ = true;
    }
  }

  void reinit(sunrealtype t, N_Vector y) override
  {
    check(ARKStepReInit(mem_, nullptr, p_->rhs(), t, y), "ARKStepReInit");
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
    ARKodeGetNumRhsEvals(mem_, 1, &s.nfe); // implicit RHS evaluations
    ARKodeGetNumLinSolvSetups(mem_, &s.nsetups);
    ARKodeGetNumErrTestFails(mem_, &s.netf);
    ARKodeGetNumNonlinSolvIters(mem_, &s.nni);
    ARKodeGetNumNonlinSolvConvFails(mem_, &s.nncf);
    ARKodeGetNumStepSolveFails(mem_, &s.nsf);
    ARKodeGetNumJacEvals(mem_, &s.nje);
    ARKodeGetNumLinRhsEvals(mem_, &s.nfe_ls);
    ARKodeGetNumPrecEvals(mem_, &s.npe);
    ARKodeGetNumPrecSolves(mem_, &s.nps);
    ARKodeGetNumLinIters(mem_, &s.nli);
    ARKodeGetNumLinConvFails(mem_, &s.nlcf);
    ARKodeGetNumJtimesEvals(mem_, &s.njtv);
    if (p_->spec().has_mass)
    {
      ARKodeGetNumMassSetups(mem_, &s.nmass_setups);
      ARKodeGetNumMassSolves(mem_, &s.nmass_solves);
      ARKodeGetNumMassMult(mem_, &s.nmass_mult);
    }
    if (roots_) { ARKodeGetNumGEvals(mem_, &s.ngeval); }
    ARKodeButcherTable Bi = nullptr, Be = nullptr;
    if (ARKStepGetCurrentButcherTables(mem_, &Bi, &Be) == ARK_SUCCESS &&
        Bi != nullptr)
    {
      stages_ = Bi->stages;
      order_  = Bi->q;
    }
  }

  std::string name() const override { return "arkode"; }

  std::string method() const override
  {
    std::string m = default_dirk_name(order_);
    if (linear_) { m += "+linear"; }
    return m;
  }

  int order() const override { return order_; }

  int stages() const override { return stages_; }

  std::string ls_name() const override { return to_string(ls_.kind); }

  std::string prec_name() const override { return ls_.prec; }

  void teardown() override
  {
    roots_ = false;
    if (mem_ != nullptr) { ARKodeFree(&mem_); }
    ls_.destroy();
    mls_.destroy();
  }

private:
  void* mem_ = nullptr;
  LinSys ls_, mls_;
  double h0_   = 0.0;
  Problem* p_  = nullptr;
  bool roots_  = false;
  int order_   = 4;
  int stages_  = 0;
  bool linear_ = false;
};

} // namespace

std::unique_ptr<Integrator> make_arkode_integrator()
{
  return std::unique_ptr<Integrator>(new ArkodeIntegrator());
}

} // namespace benchmark
