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
 * Problem interface of the FIRKODE benchmark harness.
 *
 * Every problem is the ODE system M y' = f(t,y) with a constant mass
 * matrix M (the identity unless has_mass is set). The callbacks are
 * typed with the FIRKODE function-pointer types; CVODE and ARKODE use
 * structurally identical types so the same functions are passed to
 * them unchanged, and the IDA integrator wraps them in adapters.
 *
 * Callback contract: callbacks never dereference the fy argument and
 * preconditioner solves ignore the lr argument (block-diagonal
 * preconditioners are side-agnostic). Both are needed by the IDA
 * adapters, which pass fy = NULL and lr = SUN_PREC_LEFT.
 * -----------------------------------------------------------------*/

#ifndef FIRKODE_BENCHMARK_PROBLEM_HPP_
#define FIRKODE_BENCHMARK_PROBLEM_HPP_

#include <memory>
#include <string>
#include <vector>

#include <firkode/firkode.h>
#include <firkode/firkode_ls.h>
#include <sundials/sundials_core.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>

#include "benchmark_cli.hpp"

namespace benchmark {

enum class LinSolKind
{
  Dense,
  Band,
  KLU, // sparse direct (SUNMATRIX_SPARSE in CSR form); needs a user Jacobian
  SPGMR,
  PCG
};

inline std::string to_string(LinSolKind k)
{
  switch (k)
  {
  case LinSolKind::Dense: return "dense";
  case LinSolKind::Band: return "band";
  case LinSolKind::KLU: return "klu";
  case LinSolKind::SPGMR: return "spgmr";
  case LinSolKind::PCG: return "pcg";
  }
  return "?";
}

struct ProblemSpec
{
  std::string name;
  sunindextype N = 0;
  sunrealtype t0 = 0, tf = 0;
  sunrealtype rtol_default = 0, atol_default = 0;
  LinSolKind ls_default = LinSolKind::Dense;
  sunindextype mu = 0, ml = 0;           // Jacobian band widths
  sunindextype mass_mu = 0, mass_ml = 0; // mass matrix band widths
  sunindextype nnz = 0;                  // Jacobian nonzeros (sparse solvers)
  bool has_jac     = false;
  bool has_jtimes  = false;
  bool has_prec    = false;
  bool has_mass    = false;
  bool has_exact   = false; // exact solution of the ODE system
  bool linear      = false; // f linear in y with constant Jacobian
  bool symmetric   = false; // Jacobian symmetric (PCG admissible)
  int nroots       = 0;
};

class Problem
{
public:
  virtual ~Problem() = default;

  const ProblemSpec& spec() const { return spec_; }

  // Names of the generic "--key value" parameters this problem reads
  virtual std::vector<std::string> param_keys() const { return {}; }

  // Apply problem parameters from the command line (grid sizes etc.)
  virtual void configure(const Options&) {}

  // Fill y with the initial condition
  virtual void init(N_Vector y0) const = 0;

  // Pointer passed as user_data to the callbacks
  virtual void* user_data() = 0;

  virtual FIRKRhsFn rhs() const = 0;

  virtual FIRKLsJacFn jac() const { return nullptr; }

  virtual FIRKLsJacTimesVecFn jtimes() const { return nullptr; }

  virtual FIRKLsPrecSetupFn psetup() const { return nullptr; }

  virtual FIRKLsPrecSolveFn psolve() const { return nullptr; }

  virtual FIRKLsMassFn mass() const { return nullptr; }

  virtual FIRKRootFn roots() const { return nullptr; }

  // Exact solution of the ODE system (only if spec().has_exact)
  virtual int exact(sunrealtype, N_Vector) const { return -1; }

  // Applied to the state at every restart boundary (--restart-every), e.g. a
  // split-operator update, before the integrator is re-initialized
  virtual void restart_update(sunrealtype, N_Vector) const {}

  // Let the problem adjust its preconditioner to the chosen stage count
  // and options (matrix-free heat2d needs the psolve variant)
  virtual void prepare(const Options&) {}

protected:
  ProblemSpec spec_;
};

// Write a small dense Jacobian (row-major a[N*N]) into a dense, band or
// sparse SUNMatrix; the sparse pattern is the full N x N block, so the
// pattern never changes (KLU keeps its symbolic factorization).
inline int write_full_jacobian(SUNMatrix J, const sunrealtype* a, sunindextype N)
{
  switch (SUNMatGetID(J))
  {
  case SUNMATRIX_DENSE:
    for (sunindextype i = 0; i < N; i++)
    {
      for (sunindextype j = 0; j < N; j++)
      {
        SM_ELEMENT_D(J, i, j) = a[i * N + j];
      }
    }
    return 0;
  case SUNMATRIX_BAND:
    for (sunindextype i = 0; i < N; i++)
    {
      for (sunindextype j = 0; j < N; j++)
      {
        if (j - i <= SM_UBAND_B(J) && i - j <= SM_LBAND_B(J))
        {
          SM_ELEMENT_B(J, i, j) = a[i * N + j];
        }
      }
    }
    return 0;
  case SUNMATRIX_SPARSE:
  {
    if (SUNSparseMatrix_NNZ(J) < N * N) { return -1; }
    sunindextype* ptrs = SUNSparseMatrix_IndexPointers(J);
    sunindextype* vals = SUNSparseMatrix_IndexValues(J);
    sunrealtype* data  = SUNSparseMatrix_Data(J);
    const bool csr     = SUNSparseMatrix_SparseType(J) == SUN_CSR_MAT;
    sunindextype k     = 0;
    for (sunindextype r = 0; r < N; r++)
    {
      ptrs[r] = k;
      for (sunindextype c = 0; c < N; c++, k++)
      {
        vals[k] = c;
        data[k] = csr ? a[r * N + c] : a[c * N + r];
      }
    }
    ptrs[N] = k;
    return 0;
  }
  default: return -1;
  }
}

std::unique_ptr<Problem> make_problem(const std::string& name);

} // namespace benchmark

#endif // FIRKODE_BENCHMARK_PROBLEM_HPP_
