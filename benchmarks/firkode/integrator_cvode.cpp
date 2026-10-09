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
 * CVODE (BDF) integrator wrapper of the benchmark harness. The
 * problem callbacks are passed unchanged: the CVODE function-pointer
 * types are structurally identical to the FIRKODE ones.
 * -----------------------------------------------------------------*/

#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>

#include "integrator.hpp"

namespace benchmark {

namespace {

class CvodeIntegrator : public Integrator
{
public:
  ~CvodeIntegrator() override { teardown(); }

  void setup(Problem& p, const Options& o, N_Vector y, SUNContext ctx) override
  {
    const ProblemSpec& sp = p.spec();
    p_                    = &p;
    h0_                   = (o.mode == Mode::Evolve) ? o.h0 : 0.0;
    if (sp.has_mass)
    {
      throw std::runtime_error(
        "unsupported: CVODE has no mass-matrix interface");
    }
    if (o.mode == Mode::FixedStep)
    {
      throw std::runtime_error("unsupported: CVODE has no fixed-step mode");
    }

    mem_ = CVodeCreate(CV_BDF, ctx);
    check_ptr(mem_, "CVodeCreate");
    check(CVodeInit(mem_, p.rhs(), sp.t0, y), "CVodeInit");
    check(CVodeSetUserData(mem_, p.user_data()), "CVodeSetUserData");

    if (!o.atol_vec.empty())
    {
      N_Vector avec = N_VClone(y);
      for (sunindextype i = 0; i < sp.N; i++)
      {
        NV_Ith_S(avec,
                 i) = (sunrealtype)o.atol_vec[(size_t)i % o.atol_vec.size()];
      }
      check(CVodeSVtolerances(mem_, (sunrealtype)o.rtol, avec),
            "CVodeSVtolerances");
      N_VDestroy(avec);
    }
    else
    {
      check(CVodeSStolerances(mem_, (sunrealtype)o.rtol, (sunrealtype)o.atol),
            "CVodeSStolerances");
    }

    if (o.order > 0)
    {
      check(CVodeSetMaxOrd(mem_, o.order), "CVodeSetMaxOrd");
      maxord_ = o.order;
    }
    check(CVodeSetMaxNumSteps(mem_, o.maxsteps), "CVodeSetMaxNumSteps");
    if (o.maxnef > 0)
    {
      check(CVodeSetMaxErrTestFails(mem_, o.maxnef), "CVodeSetMaxErrTestFails");
    }
    if (o.stald)
    {
      check(CVodeSetStabLimDet(mem_, SUNTRUE), "CVodeSetStabLimDet");
      stald_ = true;
    }
    if (o.tstop)
    {
      check(CVodeSetStopTime(mem_, (sunrealtype)o.tf), "CVodeSetStopTime");
    }
    const double h0 = (o.mode == Mode::OneStep) ? o.h_one : o.h0;
    if (h0 > 0.0)
    {
      check(CVodeSetInitStep(mem_, (sunrealtype)h0), "CVodeSetInitStep");
    }
    if (o.msbp > 0)
    {
      check(CVodeSetLSetupFrequency(mem_, o.msbp), "CVodeSetLSetupFrequency");
    }
    if (o.dgmax > 0.0)
    {
      check(CVodeSetDeltaGammaMaxLSetup(mem_, (sunrealtype)o.dgmax),
            "CVodeSetDeltaGammaMaxLSetup");
    }
    if (o.maxcor > 0)
    {
      check(CVodeSetMaxNonlinIters(mem_, o.maxcor), "CVodeSetMaxNonlinIters");
    }
    if (o.nlscoef > 0.0)
    {
      check(CVodeSetNonlinConvCoef(mem_, (sunrealtype)o.nlscoef),
            "CVodeSetNonlinConvCoef");
    }

    const LinSolKind kind = parse_ls(o.ls, sp.ls_default);
    ls_                   = make_linsys(kind, sp, y, o, ctx);
    check(CVodeSetLinearSolver(mem_, ls_.LS, ls_.A), "CVodeSetLinearSolver");
    const bool direct = is_direct(kind);
    if (kind == LinSolKind::KLU && p.jac() == nullptr)
    {
      throw std::runtime_error("unsupported: KLU needs a user Jacobian");
    }
    if (direct)
    {
      if (p.jac() != nullptr)
      {
        check(CVodeSetJacFn(mem_, p.jac()), "CVodeSetJacFn");
      }
    }
    else
    {
      if (p.jtimes() != nullptr)
      {
        check(CVodeSetJacTimes(mem_, nullptr, p.jtimes()), "CVodeSetJacTimes");
      }
      if (sp.has_prec && ls_.prec != "none")
      {
        check(CVodeSetPreconditioner(mem_, p.psetup(), p.psolve()),
              "CVodeSetPreconditioner");
      }
    }
    if (o.epslin > 0.0)
    {
      check(CVodeSetEpsLin(mem_, (sunrealtype)o.epslin), "CVodeSetEpsLin");
    }
    if (o.msbj > 0)
    {
      check(CVodeSetJacEvalFrequency(mem_, o.msbj), "CVodeSetJacEvalFrequency");
    }
    if (o.dgmax_jbad > 0.0)
    {
      check(CVodeSetDeltaGammaMaxBadJac(mem_, (sunrealtype)o.dgmax_jbad),
            "CVodeSetDeltaGammaMaxBadJac");
    }

    if (o.mode == Mode::Roots && sp.nroots > 0)
    {
      check(CVodeRootInit(mem_, sp.nroots, p.roots()), "CVodeRootInit");
      roots_ = true;
    }
  }

  void reinit(sunrealtype t, N_Vector y) override
  {
    check(CVodeReInit(mem_, t, y), "CVodeReInit");
    if (h0_ > 0.0)
    {
      check(CVodeSetInitStep(mem_, (sunrealtype)h0_), "CVodeSetInitStep");
    }
  }

  void set_stop_time(sunrealtype tstop) override
  {
    check(CVodeSetStopTime(mem_, tstop), "CVodeSetStopTime");
  }

  int evolve(sunrealtype tout, N_Vector y, sunrealtype& tret, bool one_step) override
  {
    return CVode(mem_, tout, y, &tret, one_step ? CV_ONE_STEP : CV_NORMAL);
  }

  bool is_success(int flag) const override
  {
    return flag == CV_SUCCESS || flag == CV_TSTOP_RETURN;
  }

  bool is_root(int flag) const override { return flag == CV_ROOT_RETURN; }

  void root_info(std::vector<int>& dirs) override
  {
    dirs.assign((size_t)p_->spec().nroots, 0);
    check(CVodeGetRootInfo(mem_, dirs.data()), "CVodeGetRootInfo");
  }

  int get_dky(sunrealtype t, int k, N_Vector dky) override
  {
    return CVodeGetDky(mem_, t, k, dky);
  }

  int est_local_errors(N_Vector ele) override
  {
    return CVodeGetEstLocalErrors(mem_, ele);
  }

  int err_weights(N_Vector ewt) override
  {
    return CVodeGetErrWeights(mem_, ewt);
  }

  sunrealtype last_step() override
  {
    sunrealtype h = 0;
    CVodeGetLastStep(mem_, &h);
    return h;
  }

  int last_order() override
  {
    int q = -1;
    CVodeGetLastOrder(mem_, &q);
    return q;
  }

  void harvest(RunStats& s) override
  {
    CVodeGetNumSteps(mem_, &s.nst);
    CVodeGetNumRhsEvals(mem_, &s.nfe);
    CVodeGetNumLinSolvSetups(mem_, &s.nsetups);
    CVodeGetNumErrTestFails(mem_, &s.netf);
    CVodeGetNumStabLimOrderReds(mem_, &s.nstab);
    CVodeGetNumNonlinSolvIters(mem_, &s.nni);
    CVodeGetNumNonlinSolvConvFails(mem_, &s.nncf);
    CVodeGetNumStepSolveFails(mem_, &s.nsf);
    CVodeGetNumJacEvals(mem_, &s.nje);
    CVodeGetNumLinRhsEvals(mem_, &s.nfe_ls);
    CVodeGetNumPrecEvals(mem_, &s.npe);
    CVodeGetNumPrecSolves(mem_, &s.nps);
    CVodeGetNumLinIters(mem_, &s.nli);
    CVodeGetNumLinConvFails(mem_, &s.nlcf);
    CVodeGetNumJtimesEvals(mem_, &s.njtv);
    if (roots_) { CVodeGetNumGEvals(mem_, &s.ngeval); }
  }

  std::string name() const override { return "cvode"; }

  std::string method() const override { return stald_ ? "BDF+STALD" : "BDF"; }

  int order() const override { return maxord_; }

  int stages() const override { return 0; }

  std::string ls_name() const override { return to_string(ls_.kind); }

  std::string prec_name() const override { return ls_.prec; }

  void teardown() override
  {
    roots_ = false;
    if (mem_ != nullptr) { CVodeFree(&mem_); }
    ls_.destroy();
  }

private:
  void* mem_ = nullptr;
  LinSys ls_;
  Problem* p_ = nullptr;
  bool roots_ = false;
  bool stald_ = false;
  double h0_  = 0.0;
  int maxord_ = 5;
};

} // namespace

std::unique_ptr<Integrator> make_cvode_integrator()
{
  return std::unique_ptr<Integrator>(new CvodeIntegrator());
}

} // namespace benchmark
