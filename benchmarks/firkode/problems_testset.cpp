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
 * Published stiff test problems (families B10 and B11), transcribed
 * from their primary sources (scripts/fetch_testset.py downloads and
 * checks them; the reference solutions are in ../references):
 *
 *   pollu    Pollution problem, 20 species, 25 reactions, t in [0, 60]
 *            (Test Set for IVP Solvers, pollu.f: feval, jeval, solut).
 *            --nb copies the system nb times into a block-diagonal
 *            system (band mu = ml = 19) with perturbed initial states;
 *            --tau > 0 applies the relaxation y <- y + (1 - exp(-dt/tau))
 *            (y_bg - y) towards the initial state at every restart of a
 *            --restart-every run (a Lie-splitting workload).
 *   ringmod  Ring modulator, ODE form (Cs = 2e-12), 15 equations,
 *            t in [0, 1e-3] (ringmod.f: feval with its overflow guard,
 *            jeval, solut).
 *   beam     Elastic beam, 40 nodes, 80 equations, t in [0, 5]
 *            (beam.f of the Test Set, identical to E. Hairer's
 *            equation.f; no analytic Jacobian: every integrator uses
 *            difference quotients).
 *   plate    Rectangular plate with a moving load, 8 x 5 grid, 80
 *            equations, t in [0, 7] (E. Hairer's plate/equation.f with
 *            the driver's constants OMEGA = 1000, STIFFN = 100,
 *            WEIGHT = 200; analytic Jacobian JPLATE).
 *   detest_b5
 *            Problem B5 of the stiff DETEST set (Enright, Hull and
 *            Lindberg, BIT 15, 1975), linear with eigenvalues
 *            -10 +- 100 i, -4, -1, -0.5, -0.1, 6 equations,
 *            t in [0, 20], y(0) = 1 (FCN, PDERV and IVALU of STDTST,
 *            ACM TOMS Algorithm 648). The error is measured against the
 *            closed-form solution: the end-point values of STDTST's
 *            EVALU are a numerical approximation of it (y1 and y2 are
 *            given as about 4e-14, the exact values are about 1e-87).
 * -----------------------------------------------------------------*/

#include <cmath>
#include <vector>

#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>

#include "problem.hpp"

namespace benchmark {

namespace {

constexpr sunrealtype ZERO = SUN_RCONST(0.0);
constexpr sunrealtype ONE  = SUN_RCONST(1.0);
constexpr sunrealtype TWO  = SUN_RCONST(2.0);
constexpr sunrealtype PI   = SUN_RCONST(3.14159265358979323846);

// fractional part of j * golden ratio: a fixed, well-spread sequence
sunrealtype spread(long j)
{
  const double x = (double)(j + 1) * 0.6180339887498949;
  return (sunrealtype)(x - std::floor(x));
}

// write a dense block (row-major a[n*n]) at rows/columns off..off+n-1 of
// a dense or band matrix
void set_block(SUNMatrix J, sunindextype off, const sunrealtype* a, sunindextype n)
{
  for (sunindextype i = 0; i < n; i++)
  {
    for (sunindextype j = 0; j < n; j++)
    {
      const sunrealtype v = a[i * n + j];
      if (SUNMatGetID(J) == SUNMATRIX_DENSE)
      {
        SM_ELEMENT_D(J, off + i, off + j) = v;
      }
      else { SM_ELEMENT_B(J, off + i, off + j) = v; }
    }
  }
}

/* ------------------------------------------------------------------
 * POLLU
 * ------------------------------------------------------------------*/

constexpr int POLLU_N = 20;

constexpr sunrealtype POLLU_K[25] = {SUN_RCONST(.35e0),   SUN_RCONST(.266e2),
                                     SUN_RCONST(.123e5),  SUN_RCONST(.86e-3),
                                     SUN_RCONST(.82e-3),  SUN_RCONST(.15e5),
                                     SUN_RCONST(.13e-3),  SUN_RCONST(.24e5),
                                     SUN_RCONST(.165e5),  SUN_RCONST(.9e4),
                                     SUN_RCONST(.22e-1),  SUN_RCONST(.12e5),
                                     SUN_RCONST(.188e1),  SUN_RCONST(.163e5),
                                     SUN_RCONST(.48e7),   SUN_RCONST(.35e-3),
                                     SUN_RCONST(.175e-1), SUN_RCONST(.1e9),
                                     SUN_RCONST(.444e12), SUN_RCONST(.124e4),
                                     SUN_RCONST(.21e1),   SUN_RCONST(.578e1),
                                     SUN_RCONST(.474e-1), SUN_RCONST(.178e4),
                                     SUN_RCONST(.312e1)};

struct PolluData
{
  long nb         = 1;
  sunrealtype tau = ZERO; // relaxation time of the restart update (0: none)
  sunrealtype dt  = ZERO; // restart interval length
  std::vector<sunrealtype> ybg; // background (initial) state, N values
};

// f for one 20-species block; y and f point at the block
void pollu_block_rhs(const sunrealtype* y, sunrealtype* f)
{
  const sunrealtype* k = POLLU_K;
  sunrealtype r[25];
  // 1-based reaction and species indices of the source, shifted here
  r[0]  = k[0] * y[0];
  r[1]  = k[1] * y[1] * y[3];
  r[2]  = k[2] * y[4] * y[1];
  r[3]  = k[3] * y[6];
  r[4]  = k[4] * y[6];
  r[5]  = k[5] * y[6] * y[5];
  r[6]  = k[6] * y[8];
  r[7]  = k[7] * y[8] * y[5];
  r[8]  = k[8] * y[10] * y[1];
  r[9]  = k[9] * y[10] * y[0];
  r[10] = k[10] * y[12];
  r[11] = k[11] * y[9] * y[1];
  r[12] = k[12] * y[13];
  r[13] = k[13] * y[0] * y[5];
  r[14] = k[14] * y[2];
  r[15] = k[15] * y[3];
  r[16] = k[16] * y[3];
  r[17] = k[17] * y[15];
  r[18] = k[18] * y[15];
  r[19] = k[19] * y[16] * y[5];
  r[20] = k[20] * y[18];
  r[21] = k[21] * y[18];
  r[22] = k[22] * y[0] * y[3];
  r[23] = k[23] * y[18] * y[0];
  r[24] = k[24] * y[19];

  f[0] = -r[0] - r[9] - r[13] - r[22] - r[23] + r[1] + r[2] + r[8] + r[10] +
         r[11] + r[21] + r[24];
  f[1]  = -r[1] - r[2] - r[8] - r[11] + r[0] + r[20];
  f[2]  = -r[14] + r[0] + r[16] + r[18] + r[21];
  f[3]  = -r[1] - r[15] - r[16] - r[22] + r[14];
  f[4]  = -r[2] + r[3] + r[3] + r[5] + r[6] + r[12] + r[19];
  f[5]  = -r[5] - r[7] - r[13] - r[19] + r[2] + r[17] + r[17];
  f[6]  = -r[3] - r[4] - r[5] + r[12];
  f[7]  = r[3] + r[4] + r[5] + r[6];
  f[8]  = -r[6] - r[7];
  f[9]  = -r[11] + r[6] + r[8];
  f[10] = -r[8] - r[9] + r[7] + r[10];
  f[11] = r[8];
  f[12] = -r[10] + r[9];
  f[13] = -r[12] + r[11];
  f[14] = r[13];
  f[15] = -r[17] - r[18] + r[15];
  f[16] = -r[19];
  f[17] = r[19];
  f[18] = -r[20] - r[21] - r[23] + r[22] + r[24];
  f[19] = -r[24] + r[23];
}

// Jacobian of one block, row-major a[20*20]; entries as in jeval
void pollu_block_jac(const sunrealtype* y, sunrealtype* a)
{
  const sunrealtype* k = POLLU_K;
  for (int i = 0; i < POLLU_N * POLLU_N; i++) { a[i] = ZERO; }
  auto J = [&](int i, int j) -> sunrealtype&
  { return a[(i - 1) * POLLU_N + (j - 1)]; };
  auto Y = [&](int i) { return y[i - 1]; };
  // k1..k25 are k[0]..k[24]
  J(1, 1)  = -k[0] - k[9] * Y(11) - k[13] * Y(6) - k[22] * Y(4) - k[23] * Y(19);
  J(1, 11) = -k[9] * Y(1) + k[8] * Y(2);
  J(1, 6)  = -k[13] * Y(1);
  J(1, 4)  = -k[22] * Y(1) + k[1] * Y(2);
  J(1, 19) = -k[23] * Y(1) + k[21];
  J(1, 2)  = k[1] * Y(4) + k[8] * Y(11) + k[2] * Y(5) + k[11] * Y(10);
  J(1, 13) = k[10];
  J(1, 20) = k[24];
  J(1, 5)  = k[2] * Y(2);
  J(1, 10) = k[11] * Y(2);

  J(2, 4)  = -k[1] * Y(2);
  J(2, 5)  = -k[2] * Y(2);
  J(2, 11) = -k[8] * Y(2);
  J(2, 10) = -k[11] * Y(2);
  J(2, 19) = k[20];
  J(2, 1)  = k[0];
  J(2, 2)  = -k[1] * Y(4) - k[2] * Y(5) - k[8] * Y(11) - k[11] * Y(10);

  J(3, 1)  = k[0];
  J(3, 4)  = k[16];
  J(3, 16) = k[18];
  J(3, 19) = k[21];
  J(3, 3)  = -k[14];

  J(4, 4) = -k[1] * Y(2) - k[15] - k[16] - k[22] * Y(1);
  J(4, 2) = -k[1] * Y(4);
  J(4, 1) = -k[22] * Y(4);
  J(4, 3) = k[14];

  J(5, 5)  = -k[2] * Y(2);
  J(5, 2)  = -k[2] * Y(5);
  J(5, 7)  = TWO * k[3] + k[5] * Y(6);
  J(5, 6)  = k[5] * Y(7) + k[19] * Y(17);
  J(5, 9)  = k[6];
  J(5, 14) = k[12];
  J(5, 17) = k[19] * Y(6);

  J(6, 6)  = -k[5] * Y(7) - k[7] * Y(9) - k[13] * Y(1) - k[19] * Y(17);
  J(6, 7)  = -k[5] * Y(6);
  J(6, 9)  = -k[7] * Y(6);
  J(6, 1)  = -k[13] * Y(6);
  J(6, 17) = -k[19] * Y(6);
  J(6, 2)  = k[2] * Y(5);
  J(6, 5)  = k[2] * Y(2);
  J(6, 16) = TWO * k[17];

  J(7, 7)  = -k[3] - k[4] - k[5] * Y(6);
  J(7, 6)  = -k[5] * Y(7);
  J(7, 14) = k[12];

  J(8, 7) = k[3] + k[4] + k[5] * Y(6);
  J(8, 6) = k[5] * Y(7);
  J(8, 9) = k[6];

  J(9, 9) = -k[6] - k[7] * Y(6);
  J(9, 6) = -k[7] * Y(9);

  J(10, 10) = -k[11] * Y(2);
  J(10, 2)  = -k[11] * Y(10) + k[8] * Y(11);
  J(10, 9)  = k[6];
  J(10, 11) = k[8] * Y(2);

  J(11, 11) = -k[8] * Y(2) - k[9] * Y(1);
  J(11, 2)  = -k[8] * Y(11);
  J(11, 1)  = -k[9] * Y(11);
  J(11, 9)  = k[7] * Y(6);
  J(11, 6)  = k[7] * Y(9);
  J(11, 13) = k[10];

  J(12, 11) = k[8] * Y(2);
  J(12, 2)  = k[8] * Y(11);

  J(13, 13) = -k[10];
  J(13, 11) = k[9] * Y(1);
  J(13, 1)  = k[9] * Y(11);

  J(14, 14) = -k[12];
  J(14, 10) = k[11] * Y(2);
  J(14, 2)  = k[11] * Y(10);

  J(15, 1) = k[13] * Y(6);
  J(15, 6) = k[13] * Y(1);

  J(16, 16) = -k[17] - k[18];
  J(16, 4)  = k[15];

  J(17, 17) = -k[19] * Y(6);
  J(17, 6)  = -k[19] * Y(17);

  J(18, 17) = k[19] * Y(6);
  J(18, 6)  = k[19] * Y(17);

  J(19, 19) = -k[20] - k[21] - k[23] * Y(1);
  J(19, 1)  = -k[23] * Y(19) + k[22] * Y(4);
  J(19, 4)  = k[22] * Y(1);
  J(19, 20) = k[24];

  J(20, 20) = -k[24];
  J(20, 1)  = k[23] * Y(19);
  J(20, 19) = k[23] * Y(1);
}

int pollu_rhs(sunrealtype, N_Vector y, N_Vector ydot, void* user_data)
{
  PolluData* d         = static_cast<PolluData*>(user_data);
  const sunrealtype* Y = N_VGetArrayPointer(y);
  sunrealtype* F       = N_VGetArrayPointer(ydot);
  for (long b = 0; b < d->nb; b++)
  {
    pollu_block_rhs(Y + POLLU_N * b, F + POLLU_N * b);
  }
  return 0;
}

int pollu_jac(sunrealtype, N_Vector y, N_Vector, SUNMatrix J, void* user_data,
              N_Vector, N_Vector, N_Vector)
{
  PolluData* d         = static_cast<PolluData*>(user_data);
  const sunrealtype* Y = N_VGetArrayPointer(y);
  sunrealtype a[POLLU_N * POLLU_N];
  if (SUNMatGetID(J) == SUNMATRIX_SPARSE)
  {
    // block-diagonal pattern of full 20 x 20 blocks (CSR or CSC alike)
    sunindextype* ptrs = SUNSparseMatrix_IndexPointers(J);
    sunindextype* vals = SUNSparseMatrix_IndexValues(J);
    sunrealtype* data  = SUNSparseMatrix_Data(J);
    const bool csr     = SUNSparseMatrix_SparseType(J) == SUN_CSR_MAT;
    sunindextype nz    = 0;
    for (long b = 0; b < d->nb; b++)
    {
      pollu_block_jac(Y + POLLU_N * b, a);
      for (int r = 0; r < POLLU_N; r++)
      {
        ptrs[POLLU_N * b + r] = nz;
        for (int c = 0; c < POLLU_N; c++, nz++)
        {
          vals[nz] = POLLU_N * b + c;
          data[nz] = csr ? a[r * POLLU_N + c] : a[c * POLLU_N + r];
        }
      }
    }
    ptrs[POLLU_N * d->nb] = nz;
    return 0;
  }
  SUNMatZero(J);
  for (long b = 0; b < d->nb; b++)
  {
    pollu_block_jac(Y + POLLU_N * b, a);
    set_block(J, POLLU_N * b, a, POLLU_N);
  }
  return 0;
}

class Pollu : public Problem
{
public:
  Pollu()
  {
    spec_.name         = "pollu";
    spec_.t0           = ZERO;
    spec_.tf           = SUN_RCONST(60.0);
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-6);
    spec_.has_jac      = true;
  }

  std::vector<std::string> param_keys() const override { return {"nb", "tau"}; }

  void configure(const Options& o) override
  {
    d_.nb = param(o, "nb", 1L);
    if (d_.nb < 1) { throw std::runtime_error("pollu: --nb must be >= 1"); }
    d_.tau          = (sunrealtype)param(o, "tau", 0.0);
    const double tf = (o.tf > 0.0) ? o.tf : (double)spec_.tf;
    d_.dt           = (sunrealtype)((tf - (double)spec_.t0) / o.restart_every);
    const sunindextype N = (sunindextype)(POLLU_N * d_.nb);
    d_.ybg.assign((size_t)N, ZERO);
    static const sunrealtype base[POLLU_N] = {ZERO,
                                              SUN_RCONST(0.2),
                                              ZERO,
                                              SUN_RCONST(0.04),
                                              ZERO,
                                              ZERO,
                                              SUN_RCONST(0.1),
                                              SUN_RCONST(0.3),
                                              SUN_RCONST(0.01),
                                              ZERO,
                                              ZERO,
                                              ZERO,
                                              ZERO,
                                              ZERO,
                                              ZERO,
                                              ZERO,
                                              SUN_RCONST(0.007),
                                              ZERO,
                                              ZERO,
                                              ZERO};
    for (long b = 0; b < d_.nb; b++)
    {
      for (int i = 0; i < POLLU_N; i++)
      {
        // copies beyond the first start from perturbed concentrations
        const sunrealtype p               = (b == 0)
                                              ? ONE
                                              : ONE + SUN_RCONST(0.1) *
                                          (TWO * spread(POLLU_N * b + i) - ONE);
        d_.ybg[(size_t)(POLLU_N * b + i)] = base[i] * p;
      }
    }
    spec_.N          = N;
    spec_.mu         = POLLU_N - 1;
    spec_.ml         = POLLU_N - 1;
    spec_.nnz        = (sunindextype)(POLLU_N * POLLU_N) * (sunindextype)d_.nb;
    spec_.ls_default = (d_.nb == 1) ? LinSolKind::Dense : LinSolKind::Band;
  }

  void init(N_Vector y0) const override
  {
    sunrealtype* Y = N_VGetArrayPointer(y0);
    for (size_t i = 0; i < d_.ybg.size(); i++) { Y[i] = d_.ybg[i]; }
  }

  void* user_data() override { return &d_; }

  FIRKRhsFn rhs() const override { return pollu_rhs; }

  FIRKLsJacFn jac() const override { return pollu_jac; }

  // Lie splitting with an exact relaxation towards the background state
  void restart_update(sunrealtype, N_Vector y) const override
  {
    if (d_.tau <= ZERO) { return; }
    const sunrealtype w = ONE - std::exp(-d_.dt / d_.tau);
    sunrealtype* Y      = N_VGetArrayPointer(y);
    for (size_t i = 0; i < d_.ybg.size(); i++)
    {
      Y[i] += w * (d_.ybg[i] - Y[i]);
    }
  }

private:
  PolluData d_;
};

/* ------------------------------------------------------------------
 * RINGMOD (ODE case, Cs = 2e-12)
 * ------------------------------------------------------------------*/

namespace rm {
constexpr sunrealtype C     = SUN_RCONST(1.6e-8);
constexpr sunrealtype CS    = SUN_RCONST(2e-12);
constexpr sunrealtype CP    = SUN_RCONST(1e-8);
constexpr sunrealtype R     = SUN_RCONST(25e3);
constexpr sunrealtype RP    = SUN_RCONST(50.0);
constexpr sunrealtype LH    = SUN_RCONST(4.45);
constexpr sunrealtype LS1   = SUN_RCONST(2e-3);
constexpr sunrealtype LS2   = SUN_RCONST(5e-4);
constexpr sunrealtype LS3   = SUN_RCONST(5e-4);
constexpr sunrealtype RG1   = SUN_RCONST(36.3);
constexpr sunrealtype RG2   = SUN_RCONST(17.3);
constexpr sunrealtype RG3   = SUN_RCONST(17.3);
constexpr sunrealtype RI    = SUN_RCONST(5e1);
constexpr sunrealtype RC    = SUN_RCONST(6e2);
constexpr sunrealtype GAMMA = SUN_RCONST(40.67286402e-9);
constexpr sunrealtype DELTA = SUN_RCONST(17.7493332);
} // namespace rm

int ringmod_rhs(sunrealtype t, N_Vector yv, N_Vector ydot, void*)
{
  using namespace rm;
  const sunrealtype* Y   = N_VGetArrayPointer(yv);
  sunrealtype* f         = N_VGetArrayPointer(ydot);
  auto y                 = [&](int i) { return Y[i - 1]; };
  const sunrealtype uin1 = SUN_RCONST(0.5) * std::sin(SUN_RCONST(2e3) * PI * t);
  const sunrealtype uin2 = TWO * std::sin(SUN_RCONST(2e4) * PI * t);
  const sunrealtype ud1  = +y(3) - y(5) - y(7) - uin2;
  const sunrealtype ud2  = -y(4) + y(6) - y(7) - uin2;
  const sunrealtype ud3  = +y(4) + y(5) + y(7) + uin2;
  const sunrealtype ud4  = -y(3) - y(6) + y(7) + uin2;
  // overflow guard of the source: a recoverable failure (IEEE double)
  if (DELTA * std::max(std::max(ud1, ud2), std::max(ud3, ud4)) > SUN_RCONST(300.0))
  {
    return 1;
  }
  const sunrealtype qud1 = GAMMA * (std::exp(DELTA * ud1) - ONE);
  const sunrealtype qud2 = GAMMA * (std::exp(DELTA * ud2) - ONE);
  const sunrealtype qud3 = GAMMA * (std::exp(DELTA * ud3) - ONE);
  const sunrealtype qud4 = GAMMA * (std::exp(DELTA * ud4) - ONE);

  f[0] = (y(8) - SUN_RCONST(0.5) * y(10) + SUN_RCONST(0.5) * y(11) + y(14) -
          y(1) / R) /
         C;
  f[1] = (y(9) - SUN_RCONST(0.5) * y(12) + SUN_RCONST(0.5) * y(13) + y(15) -
          y(2) / R) /
         C;
  f[2]  = (y(10) - qud1 + qud4) / CS;
  f[3]  = (-y(11) + qud2 - qud3) / CS;
  f[4]  = (y(12) + qud1 - qud3) / CS;
  f[5]  = (-y(13) - qud2 + qud4) / CS;
  f[6]  = (-y(7) / RP + qud1 + qud2 - qud3 - qud4) / CP;
  f[7]  = -y(1) / LH;
  f[8]  = -y(2) / LH;
  f[9]  = (SUN_RCONST(0.5) * y(1) - y(3) - RG2 * y(10)) / LS2;
  f[10] = (-SUN_RCONST(0.5) * y(1) + y(4) - RG3 * y(11)) / LS3;
  f[11] = (SUN_RCONST(0.5) * y(2) - y(5) - RG2 * y(12)) / LS2;
  f[12] = (-SUN_RCONST(0.5) * y(2) + y(6) - RG3 * y(13)) / LS3;
  f[13] = (-y(1) + uin1 - (RI + RG1) * y(14)) / LS1;
  f[14] = (-y(2) - (RC + RG1) * y(15)) / LS1;
  return 0;
}

int ringmod_jac(sunrealtype t, N_Vector yv, N_Vector, SUNMatrix Jm, void*,
                N_Vector, N_Vector, N_Vector)
{
  using namespace rm;
  const sunrealtype* Y = N_VGetArrayPointer(yv);
  auto y               = [&](int i) { return Y[i - 1]; };
  sunrealtype a[15 * 15];
  for (int i = 0; i < 15 * 15; i++) { a[i] = ZERO; }
  auto J = [&](int i, int j) -> sunrealtype&
  { return a[(i - 1) * 15 + (j - 1)]; };
  const sunrealtype uin2  = TWO * std::sin(SUN_RCONST(2e4) * PI * t);
  const sunrealtype ud1   = +y(3) - y(5) - y(7) - uin2;
  const sunrealtype ud2   = -y(4) + y(6) - y(7) - uin2;
  const sunrealtype ud3   = +y(4) + y(5) + y(7) + uin2;
  const sunrealtype ud4   = -y(3) - y(6) + y(7) + uin2;
  const sunrealtype qpud1 = GAMMA * DELTA * std::exp(DELTA * ud1);
  const sunrealtype qpud2 = GAMMA * DELTA * std::exp(DELTA * ud2);
  const sunrealtype qpud3 = GAMMA * DELTA * std::exp(DELTA * ud3);
  const sunrealtype qpud4 = GAMMA * DELTA * std::exp(DELTA * ud4);

  J(1, 1)   = -ONE / (C * R);
  J(1, 8)   = ONE / C;
  J(1, 10)  = -SUN_RCONST(0.5) / C;
  J(1, 11)  = -J(1, 10);
  J(1, 14)  = J(1, 8);
  J(2, 2)   = J(1, 1);
  J(2, 9)   = J(1, 8);
  J(2, 12)  = J(1, 10);
  J(2, 13)  = J(1, 11);
  J(2, 15)  = J(1, 14);
  J(3, 3)   = (-qpud1 - qpud4) / CS;
  J(3, 5)   = qpud1 / CS;
  J(3, 6)   = -qpud4 / CS;
  J(3, 7)   = (qpud1 + qpud4) / CS;
  J(3, 10)  = ONE / CS;
  J(4, 4)   = (-qpud2 - qpud3) / CS;
  J(4, 5)   = -qpud3 / CS;
  J(4, 6)   = qpud2 / CS;
  J(4, 7)   = (-qpud2 - qpud3) / CS;
  J(4, 11)  = -ONE / CS;
  J(5, 3)   = qpud1 / CS;
  J(5, 4)   = -qpud3 / CS;
  J(5, 5)   = (-qpud1 - qpud3) / CS;
  J(5, 7)   = (-qpud1 - qpud3) / CS;
  J(5, 12)  = ONE / CS;
  J(6, 3)   = -qpud4 / CS;
  J(6, 4)   = qpud2 / CS;
  J(6, 6)   = (-qpud2 - qpud4) / CS;
  J(6, 7)   = (qpud2 + qpud4) / CS;
  J(6, 13)  = -ONE / CS;
  J(7, 3)   = (qpud1 + qpud4) / CP;
  J(7, 4)   = (-qpud2 - qpud3) / CP;
  J(7, 5)   = (-qpud1 - qpud3) / CP;
  J(7, 6)   = (qpud2 + qpud4) / CP;
  J(7, 7)   = (-qpud1 - qpud2 - qpud3 - qpud4 - ONE / RP) / CP;
  J(8, 1)   = -ONE / LH;
  J(9, 2)   = J(8, 1);
  J(10, 1)  = SUN_RCONST(0.5) / LS2;
  J(10, 3)  = -ONE / LS2;
  J(10, 10) = -RG2 / LS2;
  J(11, 1)  = -SUN_RCONST(0.5) / LS3;
  J(11, 4)  = ONE / LS3;
  J(11, 11) = -RG3 / LS3;
  J(12, 2)  = J(10, 1);
  J(12, 5)  = J(10, 3);
  J(12, 12) = J(10, 10);
  J(13, 2)  = J(11, 1);
  J(13, 6)  = J(11, 4);
  J(13, 13) = J(11, 11);
  J(14, 1)  = -ONE / LS1;
  J(14, 14) = -(RI + RG1) / LS1;
  J(15, 2)  = J(14, 1);
  J(15, 15) = -(RC + RG1) / LS1;
  return write_full_jacobian(Jm, a, 15);
}

class Ringmod : public Problem
{
public:
  Ringmod()
  {
    spec_.name         = "ringmod";
    spec_.N            = 15;
    spec_.t0           = ZERO;
    spec_.tf           = SUN_RCONST(1.0e-3);
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-6);
    spec_.ls_default   = LinSolKind::Dense;
    spec_.has_jac      = true;
    spec_.nnz          = 15 * 15;
  }

  void init(N_Vector y0) const override { N_VConst(ZERO, y0); }

  void* user_data() override { return nullptr; }

  FIRKRhsFn rhs() const override { return ringmod_rhs; }

  FIRKLsJacFn jac() const override { return ringmod_jac; }
};

/* ------------------------------------------------------------------
 * BEAM (FTIGE of E. Hairer, as in beam.f of the Test Set)
 * ------------------------------------------------------------------*/

constexpr int BEAM_N = 40;

int beam_rhs(sunrealtype t, N_Vector yv, N_Vector ydot, void*)
{
  const int n              = BEAM_N;
  const sunrealtype nsq    = (sunrealtype)(n * n);
  const sunrealtype nquatr = nsq * nsq;
  const sunrealtype* TH    = N_VGetArrayPointer(yv);
  sunrealtype* F           = N_VGetArrayPointer(ydot);
  auto th = [&](int i) { return TH[i - 1]; }; // 1-based as in the source
  sunrealtype u[BEAM_N + 1], v[BEAM_N + 1], w[BEAM_N + 1];
  sunrealtype alpha[BEAM_N + 1], beta[BEAM_N + 1], sth[BEAM_N + 1],
    cth[BEAM_N + 1];
  for (int i = 2; i <= n; i++)
  {
    const sunrealtype thdiff = th(i) - th(i - 1);
    sth[i]                   = std::sin(thdiff);
    cth[i]                   = std::cos(thdiff);
  }
  if (t > SUN_RCONST(3.14159265358979324))
  {
    v[1] = (-SUN_RCONST(3.0) * th(1) + th(2)) * nquatr;
    for (int i = 2; i <= n - 1; i++)
    {
      v[i] = (th(i - 1) - TWO * th(i) + th(i + 1)) * nquatr;
    }
    v[n] = (th(n - 1) - th(n)) * nquatr;
  }
  else
  {
    const sunrealtype fabs = SUN_RCONST(1.5) * std::sin(t) * std::sin(t);
    const sunrealtype fx = -fabs, fy = fabs;
    v[1] = (-SUN_RCONST(3.0) * th(1) + th(2)) * nquatr +
           nsq * (fy * std::cos(th(1)) - fx * std::sin(th(1)));
    for (int i = 2; i <= n - 1; i++)
    {
      v[i] = (th(i - 1) - TWO * th(i) + th(i + 1)) * nquatr +
             nsq * (fy * std::cos(th(i)) - fx * std::sin(th(i)));
    }
    v[n] = (th(n - 1) - th(n)) * nquatr +
           nsq * (fy * std::cos(th(n)) - fx * std::sin(th(n)));
  }
  // W = D V
  w[1] = sth[2] * v[2];
  for (int i = 2; i <= n - 1; i++)
  {
    w[i] = -sth[i] * v[i - 1] + sth[i + 1] * v[i + 1];
  }
  w[n] = -sth[n] * v[n - 1];
  for (int i = 1; i <= n; i++) { w[i] += th(n + i) * th(n + i); }
  // solve C W = W
  alpha[1] = ONE;
  for (int i = 2; i <= n; i++)
  {
    alpha[i]    = TWO;
    beta[i - 1] = -cth[i];
  }
  alpha[n] = SUN_RCONST(3.0);
  for (int i = n - 1; i >= 1; i--)
  {
    const sunrealtype q = beta[i] / alpha[i + 1];
    w[i] -= w[i + 1] * q;
    alpha[i] -= beta[i] * q;
  }
  w[1] /= alpha[1];
  for (int i = 2; i <= n; i++)
  {
    w[i] = (w[i] - beta[i - 1] * w[i - 1]) / alpha[i];
  }
  // U = C V + D W
  u[1] = v[1] - cth[2] * v[2] + sth[2] * w[2];
  for (int i = 2; i <= n - 1; i++)
  {
    u[i] = TWO * v[i] - cth[i] * v[i - 1] - cth[i + 1] * v[i + 1] -
           sth[i] * w[i - 1] + sth[i + 1] * w[i + 1];
  }
  u[n] = SUN_RCONST(3.0) * v[n] - cth[n] * v[n - 1] - sth[n] * w[n - 1];
  for (int i = 1; i <= n; i++)
  {
    F[i - 1]     = th(n + i);
    F[n + i - 1] = u[i];
  }
  return 0;
}

class Beam : public Problem
{
public:
  Beam()
  {
    spec_.name         = "beam";
    spec_.N            = 2 * BEAM_N;
    spec_.t0           = ZERO;
    spec_.tf           = SUN_RCONST(5.0);
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-6);
    spec_.ls_default   = LinSolKind::Dense;
  }

  void init(N_Vector y0) const override { N_VConst(ZERO, y0); }

  void* user_data() override { return nullptr; }

  FIRKRhsFn rhs() const override { return beam_rhs; }
};

/* ------------------------------------------------------------------
 * PLATE (FPLATE / JPLATE of E. Hairer with the driver's constants)
 * ------------------------------------------------------------------*/

namespace pl {
constexpr int NX             = 8;
constexpr int NY             = 5;
constexpr int NXM1           = NX - 1;
constexpr int NYM1           = NY - 1;
constexpr int NDEMI          = NX * NY;
constexpr int NACHS1         = 2;
constexpr int NACHS2         = 4;
constexpr sunrealtype OMEGA  = SUN_RCONST(1000.0);
constexpr sunrealtype STIFFN = SUN_RCONST(100.0);
constexpr sunrealtype WEIGHT = SUN_RCONST(200.0);
constexpr sunrealtype DELX   = TWO / (sunrealtype)(NX + 1);
constexpr sunrealtype USH4   = ONE / (DELX * DELX * DELX * DELX);
constexpr sunrealtype FAC    = STIFFN * USH4;
} // namespace pl

int plate_rhs(sunrealtype x, N_Vector yv, N_Vector ydot, void*)
{
  using namespace pl;
  const sunrealtype* Y = N_VGetArrayPointer(yv);
  sunrealtype* F       = N_VGetArrayPointer(ydot);
  auto y               = [&](int k) { return Y[k - 1]; };
  for (int i = 1; i <= NX; i++)
  {
    for (int j = 1; j <= NY; j++)
    {
      const int k    = i + NX * (j - 1);
      F[k - 1]       = y(k + NDEMI);
      sunrealtype uc = SUN_RCONST(16.0) * y(k);
      if (i > 1) { uc += y(k) - SUN_RCONST(8.0) * y(k - 1); }
      if (i < NX) { uc += y(k) - SUN_RCONST(8.0) * y(k + 1); }
      if (j > 1) { uc += y(k) - SUN_RCONST(8.0) * y(k - NX); }
      if (j < NY) { uc += y(k) - SUN_RCONST(8.0) * y(k + NX); }
      if (i > 1 && j > 1) { uc += TWO * y(k - NX - 1); }
      if (i < NX && j > 1) { uc += TWO * y(k - NX + 1); }
      if (i > 1 && j < NY) { uc += TWO * y(k + NX - 1); }
      if (i < NX && j < NY) { uc += TWO * y(k + NX + 1); }
      if (i > 2) { uc += y(k - 2); }
      if (i < NXM1) { uc += y(k + 2); }
      if (j > 2) { uc += y(k - 2 * NX); }
      if (j < NYM1) { uc += y(k + 2 * NX); }
      sunrealtype force = ZERO;
      if (j == NACHS1 || j == NACHS2)
      {
        const sunrealtype xi = (sunrealtype)i * DELX;
        const sunrealtype d2 = x - xi - TWO, d5 = x - xi - SUN_RCONST(5.0);
        force = std::exp(-SUN_RCONST(5.0) * d2 * d2) +
                std::exp(-SUN_RCONST(5.0) * d5 * d5);
      }
      F[k + NDEMI - 1] = -OMEGA * y(k + NDEMI) - FAC * uc + force * WEIGHT;
    }
  }
  return 0;
}

int plate_jac(sunrealtype, N_Vector, N_Vector, SUNMatrix Jm, void*, N_Vector,
              N_Vector, N_Vector)
{
  using namespace pl;
  constexpr int N = 2 * NDEMI;
  static sunrealtype a[N * N];
  for (int i = 0; i < N * N; i++) { a[i] = ZERO; }
  auto J = [&](int i, int j) -> sunrealtype& { return a[(i - 1) * N + (j - 1)]; };
  for (int i = 1; i <= NX; i++)
  {
    for (int j = 1; j <= NY; j++)
    {
      const int k     = i + NX * (j - 1);
      J(k, k + NDEMI) = ONE;
      J(k + NDEMI, k) = -FAC * SUN_RCONST(16.0);
      if (i > 1)
      {
        J(k + NDEMI, k) -= FAC;
        J(k + NDEMI, k - 1) = FAC * SUN_RCONST(8.0);
      }
      if (i < NX)
      {
        J(k + NDEMI, k) -= FAC;
        J(k + NDEMI, k + 1) = FAC * SUN_RCONST(8.0);
      }
      if (j > 1)
      {
        J(k + NDEMI, k) -= FAC;
        J(k + NDEMI, k - NX) = FAC * SUN_RCONST(8.0);
      }
      if (j < NY)
      {
        J(k + NDEMI, k) -= FAC;
        J(k + NDEMI, k + NX) = FAC * SUN_RCONST(8.0);
      }
      if (i > 1 && j > 1) { J(k + NDEMI, k - NX - 1) = -FAC * TWO; }
      if (i < NX && j > 1) { J(k + NDEMI, k - NX + 1) = -FAC * TWO; }
      if (i > 1 && j < NY) { J(k + NDEMI, k + NX - 1) = -FAC * TWO; }
      if (i < NX && j < NY) { J(k + NDEMI, k + NX + 1) = -FAC * TWO; }
      if (i > 2) { J(k + NDEMI, k - 2) = -FAC; }
      if (i < NXM1) { J(k + NDEMI, k + 2) = -FAC; }
      if (j > 2) { J(k + NDEMI, k - 2 * NX) = -FAC; }
      if (j < NYM1) { J(k + NDEMI, k + 2 * NX) = -FAC; }
      J(k + NDEMI, k + NDEMI) = -OMEGA;
    }
  }
  return write_full_jacobian(Jm, a, N);
}

class Plate : public Problem
{
public:
  Plate()
  {
    spec_.name         = "plate";
    spec_.N            = 2 * pl::NDEMI;
    spec_.t0           = ZERO;
    spec_.tf           = SUN_RCONST(7.0);
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-9); // the driver uses atol = 1e-3 rtol
    spec_.ls_default   = LinSolKind::Dense;
    spec_.has_jac      = true;
    spec_.linear       = true;
    spec_.nnz          = 4 * pl::NDEMI * pl::NDEMI;
  }

  void init(N_Vector y0) const override { N_VConst(ZERO, y0); }

  void* user_data() override { return nullptr; }

  FIRKRhsFn rhs() const override { return plate_rhs; }

  FIRKLsJacFn jac() const override { return plate_jac; }
};

/* ------------------------------------------------------------------
 * DETEST B5 (STDTST: FCN and PDERV, labels 200 and 300, with
 * BPARM(4) = 100; IVALU, label 220)
 * ------------------------------------------------------------------*/

constexpr sunrealtype DETEST_B5_ALPHA = SUN_RCONST(100.0);

int detest_b5_rhs(sunrealtype, N_Vector yv, N_Vector ydot, void*)
{
  const sunrealtype* y = N_VGetArrayPointer(yv);
  sunrealtype* f       = N_VGetArrayPointer(ydot);
  f[0]                 = -SUN_RCONST(10.0) * y[0] + DETEST_B5_ALPHA * y[1];
  f[1]                 = -DETEST_B5_ALPHA * y[0] - SUN_RCONST(10.0) * y[1];
  f[2]                 = -SUN_RCONST(4.0) * y[2];
  f[3]                 = -y[3];
  f[4]                 = -SUN_RCONST(0.5) * y[4];
  f[5]                 = -SUN_RCONST(0.1) * y[5];
  return 0;
}

int detest_b5_jac(sunrealtype, N_Vector, N_Vector, SUNMatrix Jm, void*,
                  N_Vector, N_Vector, N_Vector)
{
  constexpr sunindextype N = 6;
  sunrealtype a[N * N]     = {}; // row-major
  a[0 * N + 0]             = -SUN_RCONST(10.0);
  a[0 * N + 1]             = DETEST_B5_ALPHA;
  a[1 * N + 0]             = -DETEST_B5_ALPHA;
  a[1 * N + 1]             = -SUN_RCONST(10.0);
  a[2 * N + 2]             = -SUN_RCONST(4.0);
  a[3 * N + 3]             = -ONE;
  a[4 * N + 4]             = -SUN_RCONST(0.5);
  a[5 * N + 5]             = -SUN_RCONST(0.1);
  return write_full_jacobian(Jm, a, N);
}

class DetestB5 : public Problem
{
public:
  DetestB5()
  {
    spec_.name         = "detest_b5";
    spec_.N            = 6;
    spec_.t0           = ZERO;
    spec_.tf           = SUN_RCONST(20.0);
    spec_.rtol_default = SUN_RCONST(1.0e-6);
    spec_.atol_default = SUN_RCONST(1.0e-6);
    spec_.ls_default   = LinSolKind::Dense;
    spec_.mu           = 1;
    spec_.ml           = 1;
    spec_.has_jac      = true;
    spec_.has_exact    = true;
    spec_.linear       = true;
    spec_.nnz          = 6 * 6;
  }

  void init(N_Vector y0) const override { N_VConst(ONE, y0); }

  void* user_data() override { return nullptr; }

  FIRKRhsFn rhs() const override { return detest_b5_rhs; }

  FIRKLsJacFn jac() const override { return detest_b5_jac; }

  // (y1, y2) = exp(-10 t) R(alpha t) (1, 1) with the rotation
  // R(s) = [[cos s, sin s], [-sin s, cos s]]; the other components decay
  int exact(sunrealtype t, N_Vector yv) const override
  {
    sunrealtype* y      = N_VGetArrayPointer(yv);
    const sunrealtype e = std::exp(-SUN_RCONST(10.0) * t);
    const sunrealtype c = std::cos(DETEST_B5_ALPHA * t);
    const sunrealtype s = std::sin(DETEST_B5_ALPHA * t);
    y[0]                = e * (c + s);
    y[1]                = e * (c - s);
    y[2]                = std::exp(-SUN_RCONST(4.0) * t);
    y[3]                = std::exp(-t);
    y[4]                = std::exp(-SUN_RCONST(0.5) * t);
    y[5]                = std::exp(-SUN_RCONST(0.1) * t);
    return 0;
  }
};

} // namespace

std::unique_ptr<Problem> make_pollu()
{
  return std::unique_ptr<Problem>(new Pollu());
}

std::unique_ptr<Problem> make_ringmod()
{
  return std::unique_ptr<Problem>(new Ringmod());
}

std::unique_ptr<Problem> make_beam()
{
  return std::unique_ptr<Problem>(new Beam());
}

std::unique_ptr<Problem> make_plate()
{
  return std::unique_ptr<Problem>(new Plate());
}

std::unique_ptr<Problem> make_detest_b5()
{
  return std::unique_ptr<Problem>(new DetestB5());
}

} // namespace benchmark
