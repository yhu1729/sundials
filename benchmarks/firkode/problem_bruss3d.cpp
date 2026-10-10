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
 * 3-D Brusselator on the unit cube (family B12, the factorization-
 * dominated regime):
 *
 *   u_t = a + u^2 v - (b + 1) u + alpha Lap u,   a = 1, b = 3.4,
 *   v_t = b u - u^2 v + alpha Lap v,
 *
 * vertex-centred n x n x n grid, 7-point Laplacian with homogeneous
 * Neumann conditions by reflection, interleaved state (u, v) at every
 * node, N = 2 n^3. The Jacobian is written in sparse (CSR) form for KLU
 * (sorted column indices, every diagonal entry present, pattern fixed),
 * and into a dense or band matrix (mu = ml = 2 n^2 + 1) for small grids.
 * Initial state: the steady state (a, b/a) with a smooth cosine
 * perturbation. Parameters: --n (default 16), --alpha (default 1, which
 * gives h rho(J) of about 7e2 at n = 16 and rtol 1e-6, growing like n^2),
 * --pert (amplitude, default 0.1).
 * -----------------------------------------------------------------*/

#include <cmath>
#include <stdexcept>
#include <vector>

#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>

#include "problem.hpp"

namespace benchmark {

namespace {

constexpr sunrealtype ZERO    = SUN_RCONST(0.0);
constexpr sunrealtype ONE     = SUN_RCONST(1.0);
constexpr sunrealtype TWO     = SUN_RCONST(2.0);
constexpr sunrealtype PI      = SUN_RCONST(3.14159265358979323846);
constexpr sunrealtype BRUSS_A = SUN_RCONST(1.0);
constexpr sunrealtype BRUSS_B = SUN_RCONST(3.4);

struct BrussData
{
  long n            = 16;
  sunrealtype alpha = SUN_RCONST(1.0);
  sunrealtype pert  = SUN_RCONST(0.1);
  sunrealtype dx    = ONE / SUN_RCONST(15.0);
};

inline long node(const BrussData& d, long i, long j, long k)
{
  return i + d.n * (j + d.n * k);
}

// neighbours of a node along one axis with reflection: the pair of
// (index, weight) contributions to the Laplacian of node i
inline void axis_neighbours(long i, long n, long& lo, long& hi)
{
  lo = (i > 0) ? i - 1 : i + 1;
  hi = (i < n - 1) ? i + 1 : i - 1;
}

int bruss_rhs(sunrealtype, N_Vector y, N_Vector ydot, void* user_data)
{
  const BrussData* d   = static_cast<BrussData*>(user_data);
  const sunrealtype* Y = N_VGetArrayPointer(y);
  sunrealtype* F       = N_VGetArrayPointer(ydot);
  const long n         = d->n;
  const sunrealtype c  = d->alpha / (d->dx * d->dx);
  for (long k = 0; k < n; k++)
  {
    for (long j = 0; j < n; j++)
    {
      for (long i = 0; i < n; i++)
      {
        const long m = node(*d, i, j, k);
        long ilo, ihi, jlo, jhi, klo, khi;
        axis_neighbours(i, n, ilo, ihi);
        axis_neighbours(j, n, jlo, jhi);
        axis_neighbours(k, n, klo, khi);
        const sunrealtype u = Y[2 * m], v = Y[2 * m + 1];
        sunrealtype lu = -SUN_RCONST(6.0) * u, lv = -SUN_RCONST(6.0) * v;
        for (long nb :
             {node(*d, ilo, j, k), node(*d, ihi, j, k), node(*d, i, jlo, k),
              node(*d, i, jhi, k), node(*d, i, j, klo), node(*d, i, j, khi)})
        {
          lu += Y[2 * nb];
          lv += Y[2 * nb + 1];
        }
        F[2 * m]     = BRUSS_A + u * u * v - (BRUSS_B + ONE) * u + c * lu;
        F[2 * m + 1] = BRUSS_B * u - u * u * v + c * lv;
      }
    }
  }
  return 0;
}

// the column pattern of one node: sorted distinct neighbour nodes
void neighbour_nodes(const BrussData& d, long i, long j, long k,
                     std::vector<long>& cols, std::vector<sunrealtype>& wts)
{
  const long n = d.n;
  long ilo, ihi, jlo, jhi, klo, khi;
  axis_neighbours(i, n, ilo, ihi);
  axis_neighbours(j, n, jlo, jhi);
  axis_neighbours(k, n, klo, khi);
  const long nbs[6] = {node(d, i, j, klo), node(d, i, jlo, k),
                       node(d, ilo, j, k), node(d, ihi, j, k),
                       node(d, i, jhi, k), node(d, i, j, khi)};
  cols.clear();
  wts.clear();
  for (long nb : nbs)
  {
    // reflected neighbours coincide with an interior one: merge weights
    bool merged = false;
    for (size_t q = 0; q < cols.size(); q++)
    {
      if (cols[q] == nb)
      {
        wts[q] += ONE;
        merged = true;
      }
    }
    if (!merged)
    {
      cols.push_back(nb);
      wts.push_back(ONE);
    }
  }
  // sort by node index (insertion sort of at most 6 entries)
  for (size_t a = 1; a < cols.size(); a++)
  {
    for (size_t b = a; b > 0 && cols[b - 1] > cols[b]; b--)
    {
      std::swap(cols[b - 1], cols[b]);
      std::swap(wts[b - 1], wts[b]);
    }
  }
}

int bruss_jac(sunrealtype, N_Vector y, N_Vector, SUNMatrix J, void* user_data,
              N_Vector, N_Vector, N_Vector)
{
  const BrussData* d   = static_cast<BrussData*>(user_data);
  const sunrealtype* Y = N_VGetArrayPointer(y);
  const long n         = d->n;
  const sunrealtype c  = d->alpha / (d->dx * d->dx);
  const bool sparse    = SUNMatGetID(J) == SUNMATRIX_SPARSE;
  sunindextype *ptrs = nullptr, *vals = nullptr;
  sunrealtype* data = nullptr;
  bool csr          = true;
  if (sparse)
  {
    ptrs = SUNSparseMatrix_IndexPointers(J);
    vals = SUNSparseMatrix_IndexValues(J);
    data = SUNSparseMatrix_Data(J);
    csr  = SUNSparseMatrix_SparseType(J) == SUN_CSR_MAT;
    if (!csr) { return -1; } // the harness builds CSR matrices
  }
  else { SUNMatZero(J); }
  auto set = [&](sunindextype r, sunindextype col, sunrealtype v)
  {
    if (SUNMatGetID(J) == SUNMATRIX_DENSE) { SM_ELEMENT_D(J, r, col) = v; }
    else { SM_ELEMENT_B(J, r, col) = v; }
  };
  std::vector<long> cols;
  std::vector<sunrealtype> wts;
  sunindextype nz = 0;
  for (long k = 0; k < n; k++)
  {
    for (long j = 0; j < n; j++)
    {
      for (long i = 0; i < n; i++)
      {
        const long m        = node(*d, i, j, k);
        const sunrealtype u = Y[2 * m], v = Y[2 * m + 1];
        neighbour_nodes(*d, i, j, k, cols, wts);
        // reaction block
        const sunrealtype juu = TWO * u * v - (BRUSS_B + ONE) -
                                SUN_RCONST(6.0) * c;
        const sunrealtype juv = u * u;
        const sunrealtype jvu = BRUSS_B - TWO * u * v;
        const sunrealtype jvv = -u * u - SUN_RCONST(6.0) * c;
        for (int comp = 0; comp < 2; comp++)
        {
          const sunindextype row = 2 * m + comp;
          if (sparse) { ptrs[row] = nz; }
          // entries in increasing column order: neighbours below m, the
          // 2 x 2 block, neighbours above m
          size_t q = 0;
          for (; q < cols.size() && cols[q] < m; q++)
          {
            const sunindextype col = 2 * cols[q] + comp;
            if (sparse) { vals[nz] = col, data[nz++] = c * wts[q]; }
            else { set(row, col, c * wts[q]); }
          }
          const sunrealtype b0 = (comp == 0) ? juu : jvu;
          const sunrealtype b1 = (comp == 0) ? juv : jvv;
          if (sparse)
          {
            vals[nz] = 2 * m, data[nz++] = b0;
            vals[nz] = 2 * m + 1, data[nz++] = b1;
          }
          else
          {
            set(row, 2 * m, b0);
            set(row, 2 * m + 1, b1);
          }
          for (; q < cols.size(); q++)
          {
            const sunindextype col = 2 * cols[q] + comp;
            if (sparse) { vals[nz] = col, data[nz++] = c * wts[q]; }
            else { set(row, col, c * wts[q]); }
          }
        }
      }
    }
  }
  if (sparse) { ptrs[2 * n * n * n] = nz; }
  return 0;
}

class Bruss3D : public Problem
{
public:
  Bruss3D()
  {
    spec_.name         = "bruss3d";
    spec_.t0           = ZERO;
    spec_.tf           = SUN_RCONST(10.0);
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-8);
    spec_.ls_default   = LinSolKind::KLU;
    spec_.has_jac      = true;
  }

  std::vector<std::string> param_keys() const override
  {
    return {"alpha", "pert"};
  }

  void configure(const Options& o) override
  {
    d_.n = (o.n > 1) ? o.n : 16;
    if (d_.n < 2) { throw std::runtime_error("bruss3d: --n must be >= 2"); }
    d_.alpha     = (sunrealtype)param(o, "alpha", 1.0);
    d_.pert      = (sunrealtype)param(o, "pert", 0.1);
    d_.dx        = ONE / (sunrealtype)(d_.n - 1);
    const long n = d_.n;
    spec_.N      = 2 * n * n * n;
    spec_.mu = spec_.ml = (sunindextype)(2 * n * n + 1);
    // nonzeros: 2 per row for the reaction block plus the distinct
    // neighbours of the node
    std::vector<long> cols;
    std::vector<sunrealtype> wts;
    sunindextype nnz = 0;
    for (long k = 0; k < n; k++)
    {
      for (long j = 0; j < n; j++)
      {
        for (long i = 0; i < n; i++)
        {
          neighbour_nodes(d_, i, j, k, cols, wts);
          nnz += 2 * (2 + (sunindextype)cols.size());
        }
      }
    }
    spec_.nnz = nnz;
  }

  void init(N_Vector y0) const override
  {
    sunrealtype* Y = N_VGetArrayPointer(y0);
    const long n   = d_.n;
    for (long k = 0; k < n; k++)
    {
      for (long j = 0; j < n; j++)
      {
        for (long i = 0; i < n; i++)
        {
          const long m         = node(d_, i, j, k);
          const sunrealtype x  = d_.dx * (sunrealtype)i,
                            yy = d_.dx * (sunrealtype)j,
                            z  = d_.dx * (sunrealtype)k;
          const sunrealtype s  = std::cos(PI * x) * std::cos(PI * yy) *
                                std::cos(PI * z);
          Y[2 * m]     = BRUSS_A + d_.pert * s;
          Y[2 * m + 1] = BRUSS_B / BRUSS_A - d_.pert * s;
        }
      }
    }
  }

  void* user_data() override { return &d_; }

  FIRKRhsFn rhs() const override { return bruss_rhs; }

  FIRKLsJacFn jac() const override { return bruss_jac; }

private:
  BrussData d_;
};

} // namespace

std::unique_ptr<Problem> make_bruss3d()
{
  return std::unique_ptr<Problem>(new Bruss3D());
}

} // namespace benchmark
