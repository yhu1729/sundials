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
 * IDA (BDF) integrator wrapper of the benchmark harness. The ODE
 * system M y' = f(t,y) is posed as the residual F(t,y,y') = M y' - f
 * with Jacobian dF/dy + cj dF/dy' = cj M - J_f. The problem callbacks
 * are wrapped in adapters that translate the IDA signatures.
 * -----------------------------------------------------------------*/

#include <ida/ida.h>
#include <ida/ida_ls.h>

#include "integrator.hpp"

namespace benchmark {

namespace {

constexpr sunrealtype ONE = SUN_RCONST(1.0);

struct IdaAdapter
{
  Problem* p        = nullptr;
  void* udata       = nullptr;
  SUNMatrix M       = nullptr; // constant mass matrix (band) or NULL
  SUNMatrix Mcj     = nullptr; // scratch of the same shape as M
  N_Vector tmp      = nullptr;
  long nmass_mult   = 0;
  long nmass_solves = 0;
};

int ida_res(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
{
  IdaAdapter* a = static_cast<IdaAdapter*>(user_data);
  int retval    = a->p->rhs()(t, yy, rr, a->udata);
  if (retval != 0) { return retval; }
  if (a->M != nullptr)
  {
    if (SUNMatMatvec(a->M, yp, a->tmp) != 0) { return -1; }
    a->nmass_mult++;
    N_VLinearSum(ONE, a->tmp, -ONE, rr, rr);
  }
  else { N_VLinearSum(ONE, yp, -ONE, rr, rr); }
  return 0;
}

// Negate every stored entry of a dense, band or sparse matrix
void negate_matrix(SUNMatrix J)
{
  sunrealtype* data = nullptr;
  sunindextype n    = 0;
  switch (SUNMatGetID(J))
  {
  case SUNMATRIX_DENSE:
    data = SUNDenseMatrix_Data(J);
    n    = SUNDenseMatrix_LData(J);
    break;
  case SUNMATRIX_BAND:
    data = SUNBandMatrix_Data(J);
    n    = SUNBandMatrix_LData(J);
    break;
  default: // sparse: the entries in use
    data = SUNSparseMatrix_Data(J);
    n    = SUNSparseMatrix_IndexPointers(J)[SUNSparseMatrix_NP(J)];
    break;
  }
  for (sunindextype i = 0; i < n; i++) { data[i] = -data[i]; }
}

// Add cj to the diagonal of a CSR matrix; the problem's sparsity pattern
// must store every diagonal entry (-1 otherwise)
int add_to_csr_diagonal(SUNMatrix J, sunrealtype cj)
{
  const sunindextype N     = SUNSparseMatrix_NP(J);
  const sunindextype* ptrs = SUNSparseMatrix_IndexPointers(J);
  const sunindextype* cols = SUNSparseMatrix_IndexValues(J);
  sunrealtype* data        = SUNSparseMatrix_Data(J);
  for (sunindextype i = 0; i < N; i++)
  {
    sunindextype k = ptrs[i];
    while (k < ptrs[i + 1] && cols[k] != i) { k++; }
    if (k == ptrs[i + 1]) { return -1; }
    data[k] += cj;
  }
  return 0;
}

int ida_jac(sunrealtype t, sunrealtype cj, N_Vector yy, N_Vector, N_Vector,
            SUNMatrix J, void* user_data, N_Vector t1, N_Vector t2, N_Vector t3)
{
  IdaAdapter* a = static_cast<IdaAdapter*>(user_data);
  int retval    = a->p->jac()(t, yy, nullptr, J, a->udata, t1, t2, t3);
  if (retval != 0) { return retval; }
  negate_matrix(J);
  if (a->M == nullptr)
  {
    const sunindextype N = N_VGetLength(yy);
    switch (SUNMatGetID(J))
    {
    case SUNMATRIX_DENSE:
      for (sunindextype i = 0; i < N; i++) { SM_ELEMENT_D(J, i, i) += cj; }
      break;
    case SUNMATRIX_BAND:
      for (sunindextype i = 0; i < N; i++) { SM_ELEMENT_B(J, i, i) += cj; }
      break;
    default: return add_to_csr_diagonal(J, cj);
    }
  }
  else
  {
    if (SUNMatCopy(a->M, a->Mcj) != 0) { return -1; }
    sunrealtype* d       = SUNBandMatrix_Data(a->Mcj);
    const sunindextype n = SUNBandMatrix_LData(a->Mcj);
    for (sunindextype i = 0; i < n; i++) { d[i] *= cj; }
    if (SUNMatScaleAdd(ONE, J, a->Mcj) != 0) { return -1; }
  }
  return 0;
}

int ida_psetup(sunrealtype t, N_Vector yy, N_Vector, N_Vector, sunrealtype cj,
               void* user_data)
{
  IdaAdapter* a       = static_cast<IdaAdapter*>(user_data);
  sunbooleantype jcur = SUNFALSE;
  return a->p->psetup()(t, yy, nullptr, SUNFALSE, &jcur, ONE / cj, a->udata);
}

int ida_psolve(sunrealtype t, N_Vector yy, N_Vector, N_Vector, N_Vector r,
               N_Vector z, sunrealtype cj, sunrealtype delta, void* user_data)
{
  IdaAdapter* a = static_cast<IdaAdapter*>(user_data);
  int retval    = a->p->psolve()(t, yy, nullptr, r, z, ONE / cj, delta,
                              SUN_PREC_LEFT, a->udata);
  if (retval != 0) { return retval; }
  N_VScale(ONE / cj, z, z);
  return 0;
}

int ida_jtimes(sunrealtype t, N_Vector yy, N_Vector, N_Vector, N_Vector v,
               N_Vector Jv, sunrealtype cj, void* user_data, N_Vector tmp1,
               N_Vector tmp2)
{
  IdaAdapter* a = static_cast<IdaAdapter*>(user_data);
  int retval    = a->p->jtimes()(v, Jv, t, yy, nullptr, a->udata, tmp1);
  if (retval != 0) { return retval; }
  if (a->M != nullptr)
  {
    if (SUNMatMatvec(a->M, v, tmp2) != 0) { return -1; }
    a->nmass_mult++;
    N_VLinearSum(cj, tmp2, -ONE, Jv, Jv);
  }
  else { N_VLinearSum(cj, v, -ONE, Jv, Jv); }
  return 0;
}

int ida_roots(sunrealtype t, N_Vector y, N_Vector, sunrealtype* gout,
              void* user_data)
{
  IdaAdapter* a = static_cast<IdaAdapter*>(user_data);
  return a->p->roots()(t, y, gout, a->udata);
}

class IdaIntegrator : public Integrator
{
public:
  ~IdaIntegrator() override { teardown(); }

  void setup(Problem& p, const Options& o, N_Vector y, SUNContext ctx) override
  {
    const ProblemSpec& sp = p.spec();
    p_                    = &p;
    if (o.mode == Mode::FixedStep)
    {
      throw std::runtime_error("unsupported: IDA has no fixed-step mode");
    }
    ad_.p     = &p;
    ad_.udata = p.user_data();
    ad_.tmp   = N_VClone(y);

    ctx_ = ctx;
    h0_  = (o.mode == Mode::Evolve) ? o.h0 : 0.0;
    yp_  = N_VClone(y);
    if (sp.has_mass)
    {
      ad_.M   = SUNBandMatrix(sp.N, sp.mass_mu, sp.mass_ml, ctx);
      ad_.Mcj = SUNBandMatrix(sp.N, sp.mass_mu, sp.mass_ml, ctx);
      check_ptr(ad_.M, "SUNBandMatrix(M)");
      N_Vector t1 = N_VClone(y), t2 = N_VClone(y), t3 = N_VClone(y);
      check(p.mass()(sp.t0, ad_.M, ad_.udata, t1, t2, t3), "mass");
      N_VDestroy(t1);
      N_VDestroy(t2);
      N_VDestroy(t3);
    }
    consistent_yp(sp.t0, y);

    mem_ = IDACreate(ctx);
    check_ptr(mem_, "IDACreate");
    check(IDAInit(mem_, ida_res, sp.t0, y, yp_), "IDAInit");
    check(IDASetUserData(mem_, &ad_), "IDASetUserData");

    if (!o.atol_vec.empty())
    {
      N_Vector avec = N_VClone(y);
      for (sunindextype i = 0; i < sp.N; i++)
      {
        NV_Ith_S(avec,
                 i) = (sunrealtype)o.atol_vec[(size_t)i % o.atol_vec.size()];
      }
      check(IDASVtolerances(mem_, (sunrealtype)o.rtol, avec),
            "IDASVtolerances");
      N_VDestroy(avec);
    }
    else
    {
      check(IDASStolerances(mem_, (sunrealtype)o.rtol, (sunrealtype)o.atol),
            "IDASStolerances");
    }

    if (o.order > 0)
    {
      check(IDASetMaxOrd(mem_, o.order), "IDASetMaxOrd");
      maxord_ = o.order;
    }
    check(IDASetMaxNumSteps(mem_, o.maxsteps), "IDASetMaxNumSteps");
    if (o.maxnef > 0)
    {
      check(IDASetMaxErrTestFails(mem_, o.maxnef), "IDASetMaxErrTestFails");
    }
    if (o.tstop)
    {
      check(IDASetStopTime(mem_, (sunrealtype)o.tf), "IDASetStopTime");
    }
    const double h0 = (o.mode == Mode::OneStep) ? o.h_one : o.h0;
    if (h0 > 0.0)
    {
      check(IDASetInitStep(mem_, (sunrealtype)h0), "IDASetInitStep");
    }
    if (o.maxcor > 0)
    {
      check(IDASetMaxNonlinIters(mem_, o.maxcor), "IDASetMaxNonlinIters");
    }
    if (o.nlscoef > 0.0)
    {
      check(IDASetNonlinConvCoef(mem_, (sunrealtype)o.nlscoef),
            "IDASetNonlinConvCoef");
    }

    const LinSolKind kind = parse_ls(o.ls, sp.ls_default);
    if (kind == LinSolKind::KLU && p.jac() == nullptr)
    {
      throw std::runtime_error("unsupported: KLU needs a user Jacobian");
    }
    if (kind == LinSolKind::KLU && sp.has_mass)
    {
      throw std::runtime_error(
        "unsupported: the IDA adapter has no sparse mass matrix");
    }
    ls_ = make_linsys(kind, sp, y, o, ctx, sp.has_mass ? sp.mass_mu : 0,
                      sp.has_mass ? sp.mass_ml : 0);
    check(IDASetLinearSolver(mem_, ls_.LS, ls_.A), "IDASetLinearSolver");
    const bool direct = is_direct(kind);
    if (direct)
    {
      if (p.jac() != nullptr)
      {
        check(IDASetJacFn(mem_, ida_jac), "IDASetJacFn");
      }
    }
    else
    {
      if (p.jtimes() != nullptr)
      {
        check(IDASetJacTimes(mem_, nullptr, ida_jtimes), "IDASetJacTimes");
      }
      if (sp.has_prec && ls_.prec != "none")
      {
        check(IDASetPreconditioner(mem_, ida_psetup, ida_psolve),
              "IDASetPreconditioner");
      }
    }
    if (o.epslin > 0.0)
    {
      check(IDASetEpsLin(mem_, (sunrealtype)o.epslin), "IDASetEpsLin");
    }

    if (o.mode == Mode::Roots && sp.nroots > 0)
    {
      check(IDARootInit(mem_, sp.nroots, ida_roots), "IDARootInit");
      roots_ = true;
    }
  }

  // yp = M^{-1} f(t, y), the consistent derivative for IDAInit and IDAReInit
  void consistent_yp(sunrealtype t, N_Vector y)
  {
    check(p_->rhs()(t, y, yp_, ad_.udata), "rhs(yp)");
    if (ad_.M != nullptr)
    {
      const ProblemSpec& sp = p_->spec();
      SUNMatrix Mlu         = SUNBandMatrix(sp.N, sp.mass_mu, sp.mass_ml, ctx_);
      check_ptr(Mlu, "SUNBandMatrix(Mlu)");
      check(SUNMatCopy(ad_.M, Mlu), "SUNMatCopy");
      SUNLinearSolver mls = SUNLinSol_Band(y, Mlu, ctx_);
      check_ptr(mls, "SUNLinSol_Band(M)");
      check(SUNLinSolInitialize(mls), "SUNLinSolInitialize");
      check(SUNLinSolSetup(mls, Mlu), "SUNLinSolSetup(M)");
      N_Vector f0 = N_VClone(y);
      N_VScale(ONE, yp_, f0);
      check(SUNLinSolSolve(mls, Mlu, yp_, f0, SUN_RCONST(0.0)),
            "SUNLinSolSolve(M)");
      ad_.nmass_solves++;
      N_VDestroy(f0);
      SUNLinSolFree(mls);
      SUNMatDestroy(Mlu);
    }
  }

  void reinit(sunrealtype t, N_Vector y) override
  {
    consistent_yp(t, y);
    check(IDAReInit(mem_, t, y, yp_), "IDAReInit");
    if (h0_ > 0.0)
    {
      check(IDASetInitStep(mem_, (sunrealtype)h0_), "IDASetInitStep");
    }
  }

  void set_stop_time(sunrealtype tstop) override
  {
    check(IDASetStopTime(mem_, tstop), "IDASetStopTime");
  }

  int evolve(sunrealtype tout, N_Vector y, sunrealtype& tret, bool one_step) override
  {
    return IDASolve(mem_, tout, &tret, y, yp_,
                    one_step ? IDA_ONE_STEP : IDA_NORMAL);
  }

  bool is_success(int flag) const override
  {
    return flag == IDA_SUCCESS || flag == IDA_TSTOP_RETURN;
  }

  bool is_root(int flag) const override { return flag == IDA_ROOT_RETURN; }

  void root_info(std::vector<int>& dirs) override
  {
    dirs.assign((size_t)p_->spec().nroots, 0);
    check(IDAGetRootInfo(mem_, dirs.data()), "IDAGetRootInfo");
  }

  int get_dky(sunrealtype t, int k, N_Vector dky) override
  {
    return IDAGetDky(mem_, t, k, dky);
  }

  int est_local_errors(N_Vector ele) override
  {
    return IDAGetEstLocalErrors(mem_, ele);
  }

  int err_weights(N_Vector ewt) override { return IDAGetErrWeights(mem_, ewt); }

  sunrealtype last_step() override
  {
    sunrealtype h = 0;
    IDAGetLastStep(mem_, &h);
    return h;
  }

  int last_order() override
  {
    int q = -1;
    IDAGetLastOrder(mem_, &q);
    return q;
  }

  void harvest(RunStats& s) override
  {
    IDAGetNumSteps(mem_, &s.nst);
    IDAGetNumResEvals(mem_, &s.nfe);
    IDAGetNumLinSolvSetups(mem_, &s.nsetups);
    IDAGetNumErrTestFails(mem_, &s.netf);
    IDAGetNumNonlinSolvIters(mem_, &s.nni);
    IDAGetNumNonlinSolvConvFails(mem_, &s.nncf);
    IDAGetNumStepSolveFails(mem_, &s.nsf);
    IDAGetNumJacEvals(mem_, &s.nje);
    IDAGetNumLinResEvals(mem_, &s.nfe_ls);
    IDAGetNumPrecEvals(mem_, &s.npe);
    IDAGetNumPrecSolves(mem_, &s.nps);
    IDAGetNumLinIters(mem_, &s.nli);
    IDAGetNumLinConvFails(mem_, &s.nlcf);
    IDAGetNumJtimesEvals(mem_, &s.njtv);
    if (p_->spec().has_mass)
    {
      s.nmass_mult   = ad_.nmass_mult;
      s.nmass_solves = ad_.nmass_solves;
      s.nmass_setups = 1;
    }
    if (roots_) { IDAGetNumGEvals(mem_, &s.ngeval); }
  }

  std::string name() const override { return "ida"; }

  std::string method() const override { return "IDA-BDF"; }

  int order() const override { return maxord_; }

  int stages() const override { return 0; }

  std::string ls_name() const override { return to_string(ls_.kind); }

  std::string prec_name() const override { return ls_.prec; }

  void teardown() override
  {
    roots_ = false;
    if (mem_ != nullptr) { IDAFree(&mem_); }
    ls_.destroy();
    if (yp_ != nullptr) { N_VDestroy(yp_); }
    if (ad_.tmp != nullptr) { N_VDestroy(ad_.tmp); }
    if (ad_.M != nullptr) { SUNMatDestroy(ad_.M); }
    if (ad_.Mcj != nullptr) { SUNMatDestroy(ad_.Mcj); }
    yp_ = nullptr;
    ad_ = IdaAdapter();
  }

private:
  void* mem_ = nullptr;
  LinSys ls_;
  N_Vector yp_ = nullptr;
  IdaAdapter ad_;
  Problem* p_     = nullptr;
  bool roots_     = false;
  double h0_      = 0.0;
  int maxord_     = 5;
  SUNContext ctx_ = nullptr;
};

} // namespace

std::unique_ptr<Integrator> make_ida_integrator()
{
  return std::unique_ptr<Integrator>(new IdaIntegrator());
}

} // namespace benchmark
