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

// Name of the DIRK table ARKStep selects for an order (its defaults)
std::string default_dirk_name(int order)
{
  int id = -1;
  switch (order)
  {
  case 1: id = ARKSTEP_DEFAULT_DIRK_1; break;
  case 2: id = ARKSTEP_DEFAULT_DIRK_2; break;
  case 3: id = ARKSTEP_DEFAULT_DIRK_3; break;
  case 4: id = ARKSTEP_DEFAULT_DIRK_4; break;
  case 5: id = ARKSTEP_DEFAULT_DIRK_5; break;
  default: return "ARKODE_DIRK_unknown";
  }
  return ARKodeButcherTable_DIRKIDToName((ARKODE_DIRKTableID)id);
}

// Names of the IMEX pair ARKStep selects for an order (its defaults)
std::string default_ark_name(int order)
{
  int it = -1, et = -1;
  switch (order)
  {
  case 2:
    it = ARKSTEP_DEFAULT_ARK_ITABLE_2;
    et = ARKSTEP_DEFAULT_ARK_ETABLE_2;
    break;
  case 3:
    it = ARKSTEP_DEFAULT_ARK_ITABLE_3;
    et = ARKSTEP_DEFAULT_ARK_ETABLE_3;
    break;
  case 4:
    it = ARKSTEP_DEFAULT_ARK_ITABLE_4;
    et = ARKSTEP_DEFAULT_ARK_ETABLE_4;
    break;
  case 5:
    it = ARKSTEP_DEFAULT_ARK_ITABLE_5;
    et = ARKSTEP_DEFAULT_ARK_ETABLE_5;
    break;
  default: return "ARKODE_ARK_unknown";
  }
  return std::string(ARKodeButcherTable_DIRKIDToName((ARKODE_DIRKTableID)it)) +
         "/" + ARKodeButcherTable_ERKIDToName((ARKODE_ERKTableID)et);
}

// ARKStepSetOptimalParams is deprecated (its parameter sets were tuned for
// older default tables) but is the package's own a-priori tuning
int set_optimal_params(void* mem)
{
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#endif
  return ARKStepSetOptimalParams(mem);
#if defined(_MSC_VER)
#pragma warning(pop)
#endif
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
}

// ARKStep: implicit (fe = 0, fi = f) or IMEX with a splitting of the
// problem (fe + fi = f, --split NAME)
class ArkodeIntegrator : public Integrator
{
public:
  explicit ArkodeIntegrator(bool imex) : imex_(imex) {}

  ~ArkodeIntegrator() override { teardown(); }

  void setup(Problem& p, const Options& o, N_Vector y, SUNContext ctx) override
  {
    const ProblemSpec& sp = p.spec();
    p_                    = &p;
    h0_                   = (o.mode == Mode::Evolve) ? o.h0 : 0.0;

    if (imex_)
    {
      split_ = find_split(p, o.split);
      if (split_.fe == nullptr || split_.fi == nullptr || split_.ff != nullptr)
      {
        throw std::runtime_error("unsupported: splitting '" + split_.name +
                                 "' is not an IMEX splitting (fe + fi)");
      }
      fe_ = split_.fe;
      fi_ = split_.fi;
    }
    else
    {
      fe_ = nullptr;
      fi_ = p.rhs();
    }
    mem_ = ARKStepCreate(fe_, fi_, sp.t0, y, ctx);
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
    if (!o.ark_itable.empty() || !o.ark_etable.empty())
    {
      if (o.order > 0)
      {
        throw std::runtime_error("table names and --order are exclusive");
      }
      if (imex_ && (o.ark_itable.empty() || o.ark_etable.empty()))
      {
        throw std::runtime_error("arkode-imex needs both --ark-itable and "
                                 "--ark-etable");
      }
      if (!imex_ && !o.ark_etable.empty())
      {
        throw std::runtime_error("--ark-etable applies to arkode-imex");
      }
      const std::string et = imex_ ? o.ark_etable : "ARKODE_ERK_NONE";
      check(ARKStepSetTableName(mem_, o.ark_itable.c_str(), et.c_str()),
            "ARKStepSetTableName");
      table_ = imex_ ? o.ark_itable + "/" + o.ark_etable : o.ark_itable;
    }
    // before the options below, which then override the tuned values
    if (o.ark_optimal)
    {
      check(set_optimal_params(mem_), "ARKStepSetOptimalParams");
      optimal_ = true;
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
    if (o.arkode_linear && (imex_ ? split_.linear_i : sp.linear))
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
                                           arkode_controller_name(o.controller)
                                             .c_str()),
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
    // the implicit part's Jacobian: the problem's, or the splitting's
    const FIRKLsJacFn jac         = imex_ ? split_.jac_i : p.jac();
    const FIRKLsJacTimesVecFn jtv = imex_ ? split_.jtimes_i : p.jtimes();
    const FIRKLsPrecSetupFn pset  = imex_ ? split_.psetup_i : p.psetup();
    const FIRKLsPrecSolveFn psol  = imex_ ? split_.psolve_i : p.psolve();
    if (kind == LinSolKind::KLU && jac == nullptr)
    {
      throw std::runtime_error("unsupported: KLU needs a user Jacobian");
    }
    if (direct)
    {
      if (jac != nullptr)
      {
        check(ARKodeSetJacFn(mem_, jac), "ARKodeSetJacFn");
      }
    }
    else
    {
      if (jtv != nullptr)
      {
        check(ARKodeSetJacTimes(mem_, nullptr, jtv), "ARKodeSetJacTimes");
      }
      if (sp.has_prec && ls_.prec != "none")
      {
        if (pset == nullptr && psol == nullptr)
        {
          throw std::runtime_error("unsupported: the splitting has no "
                                   "preconditioner");
        }
        check(ARKodeSetPreconditioner(mem_, pset, psol),
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
    check(ARKStepReInit(mem_, fe_, fi_, t, y), "ARKStepReInit");
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
    if (imex_)
    {
      ARKodeGetNumRhsEvals(mem_, 0, &s.nfe_e);
      s.split = split_.name;
    }
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
    ARKodeGetNumExpSteps(mem_, &s.nst_expstab);
    ARKodeGetNumAccSteps(mem_, &s.nst_acc);
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

  std::string name() const override { return imex_ ? "arkode-imex" : "arkode"; }

  std::string method() const override
  {
    std::string m = !table_.empty() ? table_
                    : imex_         ? default_ark_name(order_)
                                    : default_dirk_name(order_);
    if (optimal_) { m += "+optimal"; }
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
  double h0_    = 0.0;
  Problem* p_   = nullptr;
  bool roots_   = false;
  int order_    = 4;
  int stages_   = 0;
  bool linear_  = false;
  bool optimal_ = false;
  std::string table_; // --ark-itable (and --ark-etable)
  bool imex_ = false;
  RhsSplit split_;
  FIRKRhsFn fe_ = nullptr, fi_ = nullptr;
};

} // namespace

std::unique_ptr<Integrator> make_arkode_integrator()
{
  return std::unique_ptr<Integrator>(new ArkodeIntegrator(false));
}

std::unique_ptr<Integrator> make_arkode_imex_integrator()
{
  return std::unique_ptr<Integrator>(new ArkodeIntegrator(true));
}

} // namespace benchmark
