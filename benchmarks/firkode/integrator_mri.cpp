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
 * ARKODE MRIStep adapter (arkode-mri) on a splitting f = fe + fi + ff of
 * the problem (--split NAME): the slow scale integrates fe (explicit) and
 * fi (implicit, with the splitting's Jacobian or product and
 * preconditioner) with an MRI coupling (--mri-table, default MRIStep's),
 * the fast scale ff with an explicit ARKStep inner stepper
 * (--mri-fast-table, default ARKStep's). The slow step is adaptive
 * (decoupled control: an I controller on the slow step, the fast stepper
 * adaptive with the same tolerances) or fixed (--fixed-h); --mri-controller
 * htol couples the fast tolerance to the slow step (MRIHTol with two I
 * controllers); --mri-fast-h fixes the fast step. M = I only.
 * -----------------------------------------------------------------*/

#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_ls.h>
#include <arkode/arkode_mristep.h>
#include <sunadaptcontroller/sunadaptcontroller_mrihtol.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>

#include "integrator.hpp"

namespace benchmark {

namespace {

class MriIntegrator : public Integrator
{
public:
  ~MriIntegrator() override { teardown(); }

  void setup(Problem& p, const Options& o, N_Vector y, SUNContext ctx) override
  {
    const ProblemSpec& sp = p.spec();
    p_                    = &p;
    if (sp.has_mass)
    {
      throw std::runtime_error("unsupported: MRIStep needs M = I");
    }
    if (o.restart_every > 1)
    {
      throw std::runtime_error("unsupported: arkode-mri with --restart-every");
    }
    split_ = find_split(p, o.split);
    if (split_.ff == nullptr || (split_.fe == nullptr && split_.fi == nullptr))
    {
      throw std::runtime_error("unsupported: splitting '" + split_.name +
                               "' has no fast and slow parts");
    }

    // ---- fast scale: explicit ARKStep on ff ----
    inner_mem_ = ARKStepCreate(split_.ff, nullptr, sp.t0, y, ctx);
    check_ptr(inner_mem_, "ARKStepCreate(fast)");
    check(ARKodeSetUserData(inner_mem_, p.user_data()), "ARKodeSetUserData");
    check(ARKodeSStolerances(inner_mem_, (sunrealtype)o.rtol, (sunrealtype)o.atol),
          "ARKodeSStolerances(fast)");
    if (!o.mri_fast_table.empty())
    {
      check(ARKStepSetTableName(inner_mem_, "ARKODE_DIRK_NONE",
                                o.mri_fast_table.c_str()),
            "ARKStepSetTableName(fast)");
    }
    if (o.mri_fast_h > 0.0)
    {
      check(ARKodeSetFixedStep(inner_mem_, (sunrealtype)o.mri_fast_h),
            "ARKodeSetFixedStep(fast)");
    }
    check(ARKodeSetMaxNumSteps(inner_mem_, o.maxsteps),
          "ARKodeSetMaxNumSteps(fast)");
    if (o.mri_controller == "htol")
    {
      // the fast error is accumulated over a slow step for the controller
      check(ARKodeSetAccumulatedErrorType(inner_mem_, ARK_ACCUMERROR_MAX),
            "ARKodeSetAccumulatedErrorType");
    }
    else if (o.mri_controller != "decoupled")
    {
      throw std::runtime_error("unknown --mri-controller: " + o.mri_controller);
    }
    check(ARKodeCreateMRIStepInnerStepper(inner_mem_, &inner_),
          "ARKodeCreateMRIStepInnerStepper");

    // ---- slow scale ----
    mem_ = MRIStepCreate(split_.fe, split_.fi, sp.t0, y, inner_, ctx);
    check_ptr(mem_, "MRIStepCreate");
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
    if (!o.mri_table.empty())
    {
      MRIStepCoupling C = MRIStepCoupling_LoadTableByName(o.mri_table.c_str());
      check_ptr(C, "MRIStepCoupling_LoadTableByName");
      const int flag = MRIStepSetCoupling(mem_, C);
      MRIStepCoupling_Free(C);
      check(flag, "MRIStepSetCoupling");
    }
    else if (o.order > 0)
    {
      check(ARKodeSetOrder(mem_, o.order), "ARKodeSetOrder");
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
    if (o.mri_controller == "htol" && o.mode != Mode::FixedStep)
    {
      c_slow_ = SUNAdaptController_I(ctx);
      c_tol_  = SUNAdaptController_I(ctx);
      check_ptr(c_slow_, "SUNAdaptController_I");
      check_ptr(c_tol_, "SUNAdaptController_I");
      c_htol_ = SUNAdaptController_MRIHTol(c_slow_, c_tol_, ctx);
      check_ptr(c_htol_, "SUNAdaptController_MRIHTol");
      check(ARKodeSetAdaptController(mem_, c_htol_),
            "ARKodeSetAdaptController");
      htol_ = true;
    }

    if (split_.fi != nullptr)
    {
      if (o.arkode_linear && split_.linear_i)
      {
        check(ARKodeSetLinear(mem_, 0), "ARKodeSetLinear");
        linear_ = true;
      }
      if (o.maxcor > 0)
      {
        check(ARKodeSetMaxNonlinIters(mem_, o.maxcor),
              "ARKodeSetMaxNonlinIters");
      }
      if (o.nlscoef > 0.0)
      {
        check(ARKodeSetNonlinConvCoef(mem_, (sunrealtype)o.nlscoef),
              "ARKodeSetNonlinConvCoef");
      }
      if (o.msbp > 0)
      {
        check(ARKodeSetLSetupFrequency(mem_, (int)o.msbp),
              "ARKodeSetLSetupFrequency");
      }
      const LinSolKind kind = parse_ls(o.ls, sp.ls_default);
      ls_                   = make_linsys(kind, sp, y, o, ctx);
      check(ARKodeSetLinearSolver(mem_, ls_.LS, ls_.A),
            "ARKodeSetLinearSolver");
      if (kind == LinSolKind::KLU && split_.jac_i == nullptr)
      {
        throw std::runtime_error("unsupported: KLU needs a user Jacobian");
      }
      if (is_direct(kind))
      {
        if (split_.jac_i != nullptr)
        {
          check(ARKodeSetJacFn(mem_, split_.jac_i), "ARKodeSetJacFn");
        }
      }
      else
      {
        if (split_.jtimes_i != nullptr)
        {
          check(ARKodeSetJacTimes(mem_, nullptr, split_.jtimes_i),
                "ARKodeSetJacTimes");
        }
        if (sp.has_prec && ls_.prec != "none")
        {
          check(ARKodeSetPreconditioner(mem_, split_.psetup_i, split_.psolve_i),
                "ARKodeSetPreconditioner");
        }
      }
      if (o.epslin > 0.0)
      {
        check(ARKodeSetEpsLin(mem_, (sunrealtype)o.epslin), "ARKodeSetEpsLin");
      }
    }

    if (o.mode == Mode::Roots && sp.nroots > 0)
    {
      check(ARKodeRootInit(mem_, sp.nroots, p.roots()), "ARKodeRootInit");
      roots_ = true;
    }
    table_      = o.mri_table.empty() ? "MRIStep_default" : o.mri_table;
    fast_table_ = o.mri_fast_table.empty() ? "ARKStep_ERK_default"
                                           : o.mri_fast_table;
  }

  void reinit(sunrealtype, N_Vector) override
  {
    throw std::runtime_error("unsupported: arkode-mri with --restart-every");
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
    // slow scale; nfe counts the slow implicit, nfe_e the slow explicit part
    ARKodeGetNumSteps(mem_, &s.nst);
    ARKodeGetNumStepAttempts(mem_, &s.nst_att);
    if (split_.fi != nullptr) { ARKodeGetNumRhsEvals(mem_, 1, &s.nfe); }
    if (split_.fe != nullptr) { ARKodeGetNumRhsEvals(mem_, 0, &s.nfe_e); }
    ARKodeGetNumErrTestFails(mem_, &s.netf);
    if (split_.fi != nullptr)
    {
      ARKodeGetNumLinSolvSetups(mem_, &s.nsetups);
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
    }
    if (roots_) { ARKodeGetNumGEvals(mem_, &s.ngeval); }
    // fast scale
    ARKodeGetNumSteps(inner_mem_, &s.nst_f);
    ARKodeGetNumStepAttempts(inner_mem_, &s.nst_f_att);
    ARKodeGetNumRhsEvals(inner_mem_, 0, &s.nfe_f);
    ARKodeGetNumErrTestFails(inner_mem_, &s.netf_f);
    s.split           = split_.name;
    MRIStepCoupling C = nullptr;
    if (MRIStepGetCurrentCoupling(mem_, &C) == ARK_SUCCESS && C != nullptr)
    {
      order_  = C->q;
      stages_ = C->stages;
    }
  }

  std::string name() const override { return "arkode-mri"; }

  std::string method() const override
  {
    std::string m = table_ + "/" + fast_table_;
    if (htol_) { m += "+htol"; }
    if (linear_) { m += "+linear"; }
    return m;
  }

  int order() const override { return order_; }

  int stages() const override { return stages_; }

  std::string ls_name() const override
  {
    return split_.fi != nullptr ? to_string(ls_.kind) : "none";
  }

  std::string prec_name() const override { return ls_.prec; }

  void teardown() override
  {
    roots_ = false;
    // the order of examples/arkode/C_serial/ark_kpr_mri.c
    if (inner_mem_ != nullptr) { ARKodeFree(&inner_mem_); }
    if (inner_ != nullptr) { MRIStepInnerStepper_Free(&inner_); }
    if (mem_ != nullptr) { ARKodeFree(&mem_); }
    // the controllers belong to the caller; MRIHTol does not own its parts
    for (SUNAdaptController* c : {&c_htol_, &c_slow_, &c_tol_})
    {
      if (*c != nullptr)
      {
        SUNAdaptController_Destroy(*c);
        *c = nullptr;
      }
    }
    ls_.destroy();
  }

private:
  void* mem_                 = nullptr;
  void* inner_mem_           = nullptr;
  MRIStepInnerStepper inner_ = nullptr;
  SUNAdaptController c_htol_ = nullptr;
  SUNAdaptController c_slow_ = nullptr;
  SUNAdaptController c_tol_  = nullptr;
  LinSys ls_;
  Problem* p_ = nullptr;
  RhsSplit split_;
  bool roots_  = false;
  bool htol_   = false;
  bool linear_ = false;
  int order_   = 0;
  int stages_  = 0;
  std::string table_, fast_table_;
};

} // namespace

std::unique_ptr<Integrator> make_arkode_mri_integrator()
{
  return std::unique_ptr<Integrator>(new MriIntegrator());
}

} // namespace benchmark
