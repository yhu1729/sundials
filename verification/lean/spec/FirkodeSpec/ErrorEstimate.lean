/- -----------------------------------------------------------------
   Programmer(s): Yifan Hu @ UMBC
   -----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025-2026, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------- -/

import FirkodeSpec.Collocation

/-!
# The error estimate for a general number of stages

`firk_derive_ld` stores the weights `e = −γ0 l(0)ᵀA⁻¹`, where `l_j(0)` is the
Lagrange basis at `τ = 0`, and `firkErrorEstimate` forms
`γ h f(t_n, y_n) + Σ_i e_i M Z_i` with `γ = h γ0` (here with `M = I`). This
file shows, for any distinct nonzero nodes:

* `errWeights_eq`: `e_j = −γ0 l_j(0)/c_j`, which needs no matrix inverse;
* `estimate_eq`: with the stage increments `Z = h (A ⊗ I) F`, the estimate is
  `γ0 h (f(t_n, y_n) − Σ_j l_j(0) F_j)`, i.e. `γ0 h` times the difference
  between `f` at `t_n` and the interpolant of the stage derivatives there;
* `bracket_lt` and `bracket_X_pow`: that difference vanishes for polynomials
  of degree `< s` and equals `(−1)^s Π c_j ≠ 0` for `τ^s`, so the estimate
  has order `p = s`.
-/

namespace Firkode.Colloc

open Polynomial Finset

variable {F : Type*} [Field F] {s : ℕ} {c : Fin s → F}

/-- `Σ_j l_j(0) p(c_j) = p(0)` for `deg p < s`: interpolation evaluated at 0. -/
theorem interp_at_zero (hc : Function.Injective c) {p : F[X]} (hp : p.degree < s) :
    ∑ j, (ell c j).eval 0 * p.eval (c j) = p.eval 0 := by
  conv_rhs => rw [interp_exact hc hp]
  simp [eval_finsetSum, mul_comm]

/-- `Σ_j l_j(0) c_j^s = −M(0)` for `s ≥ 1`: `τ^s − M(τ)` has degree `< s`. -/
theorem interp_at_zero_X_pow (hc : Function.Injective c) (hs : 0 < s) :
    ∑ j, (ell c j).eval 0 * c j ^ s = -(nodePoly c).eval 0 := by
  have hdeg : ((X : F[X]) ^ s - nodePoly c).degree < s := by
    have := degree_sub_lt (p := (X : F[X]) ^ s) (q := nodePoly c)
      (by rw [degree_X_pow, degree_eq_natDegree (nodePoly_monic (c := c)).ne_zero,
        natDegree_nodePoly])
      (monic_X_pow s).ne_zero
      (by rw [leadingCoeff_X_pow, (nodePoly_monic (c := c)).leadingCoeff])
    rwa [degree_X_pow] at this
  have := interp_at_zero hc hdeg
  simp only [eval_sub, eval_pow, eval_X, nodePoly_eval_node, sub_zero] at this
  rw [this, zero_pow hs.ne', zero_sub]

/-- The bracket of the error estimate, `g(0) − Σ_j l_j(0) g(c_j)`. -/
noncomputable def bracket (c : Fin s → F) (g : F[X]) : F :=
  g.eval 0 - ∑ j, (ell c j).eval 0 * g.eval (c j)

/-- The bracket vanishes for polynomials of degree `< s`. -/
theorem bracket_lt (hc : Function.Injective c) {p : F[X]} (hp : p.degree < s) :
    bracket c p = 0 := by
  rw [bracket, interp_at_zero hc hp, sub_self]

/-- The bracket of `τ^s` is `M(0) = (−1)^s Π_j c_j`, nonzero for nonzero nodes:
    the error estimate has order exactly `s`. -/
theorem bracket_X_pow (hc : Function.Injective c) (hs : 0 < s) :
    bracket c (X ^ s) = (-1) ^ s * ∏ j, c j := by
  rw [bracket]
  simp only [eval_pow, eval_X]
  rw [interp_at_zero_X_pow hc hs, nodePoly, Lagrange.eval_nodal, zero_pow hs.ne', zero_sub,
    neg_neg]
  simp only [zero_sub]
  rw [Finset.prod_neg, card_univ, Fintype.card_fin]

/-- `l_s(0) = −M(0)/M′(1)` when the last node is `c_s = 1`. -/
theorem ell_last_eval_zero (last : Fin s) (hlast : c last = 1) :
    (ell c last).eval 0 = -(nodePoly c).eval 0 / (derivative (nodePoly c)).eval 1 := by
  classical
  rw [← hlast, nodePoly, Lagrange.eval_nodal_derivative_eval_node_eq (mem_univ last),
    Lagrange.eval_nodal, Lagrange.eval_nodal, ← mul_prod_erase univ _ (mem_univ last), hlast,
    ell, Lagrange.basis, eval_prod]
  simp only [Lagrange.basisDivisor, eval_mul, eval_C, eval_sub, eval_X, zero_sub, hlast]
  rw [prod_mul_distrib, prod_inv_distrib]
  ring

variable [CharZero F]

/-- `Σ_j (l_j(0)/c_j) a_jk = l_k(0)`. -/
theorem ell0_div_vecMul (hc : Function.Injective c) (h0 : ∀ i, c i ≠ 0) (k : Fin s) :
    ∑ j, (ell c j).eval 0 / c j * collocA c j k = (ell c k).eval 0 := by
  set H := divX (antideriv (ell c k))
  have hsplit : antideriv (ell c k) = X * H := by
    have := X_mul_divX_add (antideriv (ell c k))
    rw [coeff_antideriv_zero, map_zero, add_zero] at this
    exact this.symm
  have hA : ∀ j, collocA c j k = c j * H.eval (c j) := by
    intro j; simp [collocA, hsplit]
  have hHdeg : H.degree < s := by
    have h1 : H.natDegree ≤ s - 1 := by
      rw [natDegree_divX_eq_natDegree_tsub_one]
      have := natDegree_antideriv_le (ell c k)
      rw [natDegree_ell hc] at this
      omega
    have hs : 0 < s := Fin.pos k
    calc H.degree ≤ H.natDegree := degree_le_natDegree
      _ ≤ (s - 1 : ℕ) := by exact_mod_cast h1
      _ < s := by exact_mod_cast Nat.sub_lt hs one_pos
  have hH0 : H.eval 0 = (ell c k).eval 0 := by
    rw [← coeff_zero_eq_eval_zero, coeff_divX, coeff_antideriv_succ, coeff_zero_eq_eval_zero]
    simp
  rw [← hH0, ← interp_at_zero hc hHdeg]
  refine sum_congr rfl fun j _ => ?_
  rw [hA]
  field_simp [h0 j]

/-- `e = −γ0 l(0)ᵀA⁻¹ = −γ0 (l_j(0)/c_j)_j`. -/
theorem errWeights_eq (hc : Function.Injective c) (h0 : ∀ i, c i ≠ 0) (γ0 : F) :
    -γ0 • Matrix.vecMul (fun j => (ell c j).eval 0) (collocA c)⁻¹ =
      fun j => -γ0 * ((ell c j).eval 0 / c j) := by
  have hA : IsUnit (collocA c).det := isUnit_iff_ne_zero.mpr (collocA_det_ne_zero hc h0)
  have hw : Matrix.vecMul (fun j => (ell c j).eval 0 / c j) (collocA c) =
      fun k => (ell c k).eval 0 := by
    ext k; exact ell0_div_vecMul hc h0 k
  rw [← hw, Matrix.vecMul_vecMul, Matrix.mul_nonsing_inv _ hA, Matrix.vecMul_one]
  ext j; simp

/-- The unfiltered error estimate. With the stage increments
    `Z_i = h Σ_j a_ij F_j` and the weights `e_j = −γ0 l_j(0)/c_j`,
    `γ0 h f₀ + Σ_i e_i Z_i = γ0 h (f₀ − Σ_j l_j(0) F_j)`. -/
theorem estimate_eq (hc : Function.Injective c) (h0 : ∀ i, c i ≠ 0) (γ0 h : F)
    {V : Type*} [AddCommGroup V] [Module F V] (f₀ : V) (Fs : Fin s → V) :
    (γ0 * h) • f₀ + ∑ i, (-γ0 * ((ell c i).eval 0 / c i)) • (h • ∑ j, collocA c i j • Fs j) =
      (γ0 * h) • (f₀ - ∑ j, (ell c j).eval 0 • Fs j) := by
  have key : ∀ j, ∑ i, -γ0 * ((ell c i).eval 0 / c i) * (h * collocA c i j) =
      -(γ0 * h) * (ell c j).eval 0 := by
    intro j
    rw [← ell0_div_vecMul hc h0 j, mul_sum]
    exact sum_congr rfl fun i _ => by ring
  have hsum : ∑ i, (-γ0 * ((ell c i).eval 0 / c i)) • (h • ∑ j, collocA c i j • Fs j) =
      -(γ0 * h) • ∑ j, (ell c j).eval 0 • Fs j := by
    calc ∑ i, (-γ0 * ((ell c i).eval 0 / c i)) • (h • ∑ j, collocA c i j • Fs j)
        = ∑ i, ∑ j, (-γ0 * ((ell c i).eval 0 / c i) * (h * collocA c i j)) • Fs j := by
          simp only [smul_sum, smul_smul]
      _ = ∑ j, ∑ i, (-γ0 * ((ell c i).eval 0 / c i) * (h * collocA c i j)) • Fs j := sum_comm
      _ = ∑ j, (-(γ0 * h) * (ell c j).eval 0) • Fs j := by
          simp only [← sum_smul, key]
      _ = -(γ0 * h) • ∑ j, (ell c j).eval 0 • Fs j := by
          rw [smul_sum]; simp only [smul_smul]
  rw [hsum, smul_sub, neg_smul, sub_eq_add_neg]

end Firkode.Colloc
