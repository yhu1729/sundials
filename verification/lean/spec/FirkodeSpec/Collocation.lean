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

import Mathlib.LinearAlgebra.Lagrange
import Mathlib.LinearAlgebra.Vandermonde
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.LinearAlgebra.Matrix.Nondegenerate
import Mathlib.Algebra.Polynomial.Div
import Mathlib.Data.Nat.Factorial.BigOperators
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination

/-!
# Collocation methods for a general number of stages

For distinct nodes `c_1, …, c_s` the collocation method has coefficients
`a_ij = ∫_0^{c_i} l_j`, with `l_j` the Lagrange basis on the nodes. Here the
integral is the algebraic antiderivative `antideriv` (zero constant term),
which avoids measure theory; `Integral.lean` identifies it with the interval
integral over the reals.

Main results, for any field of characteristic zero:
* `collocA_exact`: `Σ_j a_ij p(c_j) = ∫_0^{c_i} p` for `deg p < s`; in
  particular the simplifying assumption C(s) (`collocA_C`).
* `collocA_mod`: `Σ_j a_ij p(c_j) = ∫_0^{c_i} (p mod M)` for every `p`, with `M`
  the node polynomial; hence B(2s−1) when `M` is orthogonal on `[0, 1]` to the
  polynomials of degree `≤ s − 2` (`collocA_B`), and the quadrature error in
  degree `2s − 1` (`collocA_B2s`, `collocA_B_remainder`).
* `collocA_unique`: C(s) determines the coefficients.
* `det_collocA`: `det A = Π c_i / s!`, so `A` is invertible for nonzero nodes.
-/

namespace Firkode.Colloc

open Polynomial Finset

variable {F : Type*} [Field F]

/-! ## The antiderivative -/

/-- The antiderivative with zero constant term. -/
noncomputable def antideriv : F[X] →ₗ[F] F[X] :=
  Polynomial.lsum fun n => ((n : F) + 1)⁻¹ • (monomial (n + 1) : F →ₗ[F] F[X])

theorem antideriv_monomial (n : ℕ) (a : F) :
    antideriv (monomial n a) = monomial (n + 1) (a / (n + 1)) := by
  rw [antideriv, lsum_apply, sum_monomial_index _ _ (by simp)]
  simp [smul_monomial, div_eq_inv_mul]

theorem antideriv_X_pow (n : ℕ) :
    antideriv (X ^ n : F[X]) = C ((n : F) + 1)⁻¹ * X ^ (n + 1) := by
  rw [← monomial_one_right_eq_X_pow, antideriv_monomial, C_mul_X_pow_eq_monomial]
  simp [div_eq_inv_mul]

theorem eval_antideriv_X_pow (n : ℕ) (x : F) :
    (antideriv (X ^ n : F[X])).eval x = x ^ (n + 1) / (n + 1) := by
  simp [antideriv_X_pow, div_eq_inv_mul]

theorem coeff_antideriv_succ (p : F[X]) (n : ℕ) :
    (antideriv p).coeff (n + 1) = p.coeff n / (n + 1) := by
  induction p using Polynomial.induction_on' with
  | add p q hp hq => simp [hp, hq, add_div]
  | monomial m a =>
    rw [antideriv_monomial, coeff_monomial, coeff_monomial]
    by_cases h : m = n
    · subst h; simp
    · simp [h]

theorem coeff_antideriv_zero (p : F[X]) : (antideriv p).coeff 0 = 0 := by
  induction p using Polynomial.induction_on' with
  | add p q hp hq => simp [hp, hq]
  | monomial m a => simp [antideriv_monomial, coeff_monomial]

theorem eval_antideriv_zero (p : F[X]) : (antideriv p).eval 0 = 0 := by
  rw [← coeff_zero_eq_eval_zero, coeff_antideriv_zero]

theorem derivative_antideriv [CharZero F] (p : F[X]) : derivative (antideriv p) = p := by
  induction p using Polynomial.induction_on' with
  | add p q hp hq => simp [hp, hq]
  | monomial n a =>
    rw [antideriv_monomial, derivative_monomial]
    have : ((n : F) + 1) ≠ 0 := by exact_mod_cast Nat.succ_ne_zero n
    simp [div_mul_cancel₀ _ this]

theorem natDegree_antideriv_le (p : F[X]) : (antideriv p).natDegree ≤ p.natDegree + 1 := by
  rw [natDegree_le_iff_coeff_eq_zero]
  intro N hN
  cases N with
  | zero => simp [coeff_antideriv_zero]
  | succ n =>
    rw [coeff_antideriv_succ, coeff_eq_zero_of_natDegree_lt (by omega), zero_div]

/-- `∫_0^x p` as a sum over the coefficients, as `firk_collocation_coeffs`
    computes it. -/
theorem eval_antideriv_eq_sum (p : F[X]) (x : F) {N : ℕ} (hN : p.natDegree < N) :
    (antideriv p).eval x = ∑ m ∈ range N, p.coeff m * x ^ (m + 1) / (m + 1) := by
  conv_lhs => rw [p.as_sum_range' N hN]
  simp [map_sum, eval_finsetSum, antideriv_monomial, div_eq_mul_inv, mul_comm, mul_left_comm]

/-! ## Lagrange basis and collocation coefficients -/

variable {s : ℕ} (c : Fin s → F)

/-- The Lagrange basis polynomial `l_j` on the nodes `c`. -/
noncomputable def ell (j : Fin s) : F[X] := Lagrange.basis univ c j

/-- The node polynomial `M(τ) = Π_i (τ − c_i)`. -/
noncomputable def nodePoly : F[X] := Lagrange.nodal univ c

/-- The collocation coefficients `a_ij = ∫_0^{c_i} l_j`. -/
noncomputable def collocA : Matrix (Fin s) (Fin s) F :=
  Matrix.of fun i j => (antideriv (ell c j)).eval (c i)

variable {c}

theorem ell_eval_node (hc : Function.Injective c) (i j : Fin s) :
    (ell c j).eval (c i) = if i = j then 1 else 0 := by
  split_ifs with h
  · subst h; exact Lagrange.eval_basis_self hc.injOn (mem_univ _)
  · exact Lagrange.eval_basis_of_ne (Ne.symm h) (mem_univ _)

theorem natDegree_ell (hc : Function.Injective c) (j : Fin s) : (ell c j).natDegree = s - 1 := by
  simpa [ell] using Lagrange.natDegree_basis hc.injOn (mem_univ j)

/-- Interpolation is exact below degree `s`: `p = Σ_j p(c_j) l_j`. -/
theorem interp_exact (hc : Function.Injective c) {p : F[X]} (hp : p.degree < s) :
    p = ∑ j, C (p.eval (c j)) * ell c j := by
  have := Lagrange.eq_interpolate (s := univ) (v := c) hc.injOn (by simpa using hp)
  simpa [Lagrange.interpolate_apply, ell] using this

/-- Collocation integrates exactly below degree `s`:
    `Σ_j a_ij p(c_j) = ∫_0^{c_i} p`. -/
theorem collocA_exact (hc : Function.Injective c) {p : F[X]} (hp : p.degree < s) (i : Fin s) :
    ∑ j, collocA c i j * p.eval (c j) = (antideriv p).eval (c i) := by
  conv_rhs => rw [interp_exact hc hp]
  rw [map_sum, eval_finsetSum]
  refine sum_congr rfl fun j _ => ?_
  rw [← smul_eq_C_mul, map_smul, eval_smul, smul_eq_mul, collocA, Matrix.of_apply, mul_comm]

/-- The simplifying assumption C(s): `Σ_j a_ij c_j^k = c_i^(k+1)/(k+1)` for `k < s`. -/
theorem collocA_C (hc : Function.Injective c) (i : Fin s) {k : ℕ} (hk : k < s) :
    ∑ j, collocA c i j * c j ^ k = c i ^ (k + 1) / (k + 1) := by
  have := collocA_exact hc (p := X ^ k) (by rw [degree_X_pow]; exact_mod_cast hk) i
  simpa [eval_antideriv_X_pow] using this

theorem nodePoly_monic : (nodePoly c).Monic := Lagrange.nodal_monic

theorem natDegree_nodePoly : (nodePoly c).natDegree = s := by
  simp [nodePoly, Lagrange.natDegree_nodal]

theorem nodePoly_eval_node (i : Fin s) : (nodePoly c).eval (c i) = 0 :=
  Lagrange.eval_nodal_at_node (mem_univ i)

/-- `Σ_j a_ij p(c_j) = ∫_0^{c_i} (p mod M)` for every polynomial `p`. -/
theorem collocA_mod (hc : Function.Injective c) (p : F[X]) (i : Fin s) :
    ∑ j, collocA c i j * p.eval (c j) = (antideriv (p %ₘ nodePoly c)).eval (c i) := by
  have hdeg : (p %ₘ nodePoly c).degree < s := by
    have := degree_modByMonic_lt p (nodePoly_monic (c := c))
    rwa [degree_eq_natDegree (nodePoly_monic (c := c)).ne_zero, natDegree_nodePoly] at this
  rw [← collocA_exact hc hdeg i]
  refine sum_congr rfl fun j _ => ?_
  conv_lhs => rw [← modByMonic_add_div p (nodePoly c)]
  simp [nodePoly_eval_node]

/-- B(2s−1) in polynomial form: if `c_s = 1` and the node polynomial is
    orthogonal on `[0, 1]` to the polynomials of degree `≤ s − 2`, the weights
    `b_j = a_sj` integrate every polynomial of degree `≤ 2s − 2` exactly. -/
theorem collocA_B_poly (hc : Function.Injective c) (last : Fin s) (hlast : c last = 1)
    (horth : ∀ j < s - 1, (antideriv (X ^ j * nodePoly c)).eval 1 = 0)
    {p : F[X]} (hp : p.natDegree < 2 * s - 1) :
    ∑ j, collocA c last j * p.eval (c j) = (antideriv p).eval 1 := by
  rw [collocA_mod hc p last, hlast]
  by_cases hps : p.natDegree < s
  · rw [(modByMonic_eq_self_iff (nodePoly_monic (c := c))).mpr
      (degree_lt_degree (by rwa [natDegree_nodePoly]))]
  have hqdeg := natDegree_divByMonic p (nodePoly_monic (c := c))
  rw [natDegree_nodePoly] at hqdeg
  set q := p /ₘ nodePoly c
  have hsplit : p %ₘ nodePoly c = p - nodePoly c * q := by
    rw [eq_sub_iff_add_eq, modByMonic_add_div]
  have hq : q.natDegree < s - 1 := by omega
  have hzero : (antideriv (nodePoly c * q)).eval 1 = 0 := by
    conv_lhs => rw [q.as_sum_range' (s - 1) hq]
    simp only [mul_sum, map_sum, eval_finsetSum]
    refine sum_eq_zero fun m hm => ?_
    rw [← C_mul_X_pow_eq_monomial, mul_left_comm, mul_comm (nodePoly c) (X ^ m), ← smul_eq_C_mul,
      map_smul, eval_smul, horth m (mem_range.mp hm), smul_zero]
  rw [hsplit, map_sub, eval_sub, hzero, sub_zero]

/-- B(2s−1): `Σ_j b_j c_j^k = 1/(k+1)` for `k ≤ 2s − 2`, under the hypotheses
    of `collocA_B_poly`. -/
theorem collocA_B (hc : Function.Injective c) (last : Fin s) (hlast : c last = 1)
    (horth : ∀ j < s - 1, (antideriv (X ^ j * nodePoly c)).eval 1 = 0)
    {k : ℕ} (hk : k < 2 * s - 1) :
    ∑ j, collocA c last j * c j ^ k = 1 / (k + 1) := by
  have := collocA_B_poly hc last hlast horth (p := X ^ k) (by rwa [natDegree_X_pow])
  simpa [eval_antideriv_X_pow] using this

/-- The B(2s) residual: under the hypotheses of `collocA_B_poly`,
    `Σ_j b_j c_j^(2s−1) = 1/(2s) − ∫_0^1 τ^(s−1) M(τ) dτ`. -/
theorem collocA_B2s (hc : Function.Injective c) (last : Fin s) (hlast : c last = 1)
    (horth : ∀ j < s - 1, (antideriv (X ^ j * nodePoly c)).eval 1 = 0) :
    ∑ j, collocA c last j * c j ^ (2 * s - 1) =
      1 / (2 * s : F) - (antideriv (X ^ (s - 1) * nodePoly c)).eval 1 := by
  have hs : 1 ≤ s := Nat.one_le_iff_ne_zero.mpr fun h => by subst h; exact last.elim0
  set p : F[X] := X ^ (2 * s - 1) - X ^ (s - 1) * nodePoly c
  have hdeg2 : (X ^ (s - 1) * nodePoly c : F[X]).natDegree = 2 * s - 1 := by
    rw [(monic_X_pow _).natDegree_mul (nodePoly_monic (c := c)), natDegree_X_pow,
      natDegree_nodePoly]
    omega
  have hp : p.natDegree < 2 * s - 1 := by
    have hlt : p.degree < ((X : F[X]) ^ (2 * s - 1)).degree := by
      apply degree_sub_lt
      · rw [degree_eq_natDegree ((monic_X_pow _).mul (nodePoly_monic (c := c))).ne_zero, hdeg2,
          degree_X_pow]
      · exact (monic_X_pow _).ne_zero
      · rw [leadingCoeff_X_pow, ((monic_X_pow _).mul (nodePoly_monic (c := c))).leadingCoeff]
    rw [degree_X_pow] at hlt
    by_cases hp0 : p = 0
    · rw [hp0, natDegree_zero]; omega
    · exact_mod_cast (degree_eq_natDegree hp0) ▸ hlt
  have h := collocA_B_poly hc last hlast horth hp
  have hnode : ∀ j, p.eval (c j) = c j ^ (2 * s - 1) := by
    intro j; simp [p, nodePoly_eval_node]
  simp only [hnode] at h
  rw [h, map_sub, eval_sub, eval_antideriv_X_pow]
  have h2 : ((2 * s - 1 : ℕ) : F) + 1 = 2 * s := by
    rw [Nat.cast_sub (by omega)]; push_cast; ring
  rw [h2, one_pow]

/-- Quadrature with remainder: under the hypotheses of `collocA_B_poly`,
    `Σ_j b_j p(c_j) = ∫_0^1 p − p_{2s−1} ∫_0^1 τ^(s−1) M(τ) dτ` for every `p`
    of degree `≤ 2s − 1`, with `p_{2s−1}` its coefficient of `τ^(2s−1)`. -/
theorem collocA_B_remainder (hc : Function.Injective c) (last : Fin s) (hlast : c last = 1)
    (horth : ∀ j < s - 1, (antideriv (X ^ j * nodePoly c)).eval 1 = 0)
    {p : F[X]} (hp : p.natDegree ≤ 2 * s - 1) :
    ∑ j, collocA c last j * p.eval (c j) =
      (antideriv p).eval 1 - p.coeff (2 * s - 1) * (antideriv (X ^ (s - 1) * nodePoly c)).eval 1 := by
  have hs : 1 ≤ s := Nat.one_le_iff_ne_zero.mpr fun h => by subst h; exact last.elim0
  set a := p.coeff (2 * s - 1)
  set q := p - C a * X ^ (2 * s - 1)
  have hq : q.natDegree < 2 * s - 1 := by
    have hle : q.natDegree ≤ 2 * s - 2 := by
      rw [natDegree_le_iff_coeff_eq_zero]
      intro N hN
      have hN' : 2 * s - 2 < N := by exact_mod_cast hN
      simp only [q, coeff_sub, coeff_C_mul_X_pow]
      rcases eq_or_lt_of_le (show 2 * s - 1 ≤ N by omega) with h | h
      · rw [← h, if_pos rfl, sub_self]
      · rw [coeff_eq_zero_of_natDegree_lt (by omega), if_neg (by omega), sub_zero]
    omega
  have hsplit : p = q + C a * X ^ (2 * s - 1) := by simp [q]
  have hB := collocA_B_poly hc last hlast horth hq
  have hB2 := collocA_B2s hc last hlast horth
  have h2 : ((2 * s - 1 : ℕ) : F) + 1 = 2 * s := by
    rw [Nat.cast_sub (by omega)]; push_cast; ring
  have hint : (antideriv p).eval 1 = (antideriv q).eval 1 + a / (2 * s) := by
    conv_lhs => rw [hsplit]
    rw [map_add, eval_add, ← smul_eq_C_mul, map_smul, eval_smul, eval_antideriv_X_pow, h2,
      one_pow, smul_eq_mul, mul_one_div]
  have hsum : ∑ j, collocA c last j * p.eval (c j) =
      ∑ j, collocA c last j * q.eval (c j) + a * ∑ j, collocA c last j * c j ^ (2 * s - 1) := by
    conv_lhs => rw [hsplit]
    simp only [eval_add, eval_mul, eval_C, eval_pow, eval_X, mul_add, Finset.sum_add_distrib,
      Finset.mul_sum]
    congr 1
    exact Finset.sum_congr rfl fun j _ => by ring
  rw [hsum, hB, hB2, hint]
  ring

/-- `∫_0^x p′ = p(x) − p(0)`. -/
theorem antideriv_derivative [CharZero F] (p : F[X]) : antideriv (derivative p) = p - C (p.coeff 0) := by
  ext n
  cases n with
  | zero => simp [coeff_antideriv_zero]
  | succ n =>
    rw [coeff_antideriv_succ, coeff_derivative, coeff_sub, coeff_C_succ, sub_zero]
    have : ((n : F) + 1) ≠ 0 := by exact_mod_cast Nat.succ_ne_zero n
    rw [mul_div_assoc, div_self this, mul_one]

/-- C(s) determines the coefficients. -/
theorem collocA_unique (hc : Function.Injective c) {A : Matrix (Fin s) (Fin s) F}
    (hA : ∀ i, ∀ k < s, ∑ j, A i j * c j ^ k = c i ^ (k + 1) / (k + 1)) : A = collocA c := by
  have hV : IsUnit (Matrix.vandermonde c).det :=
    isUnit_iff_ne_zero.mpr (Matrix.det_vandermonde_ne_zero_iff.mpr hc)
  have h : A * Matrix.vandermonde c = collocA c * Matrix.vandermonde c := by
    ext i k
    simp only [Matrix.mul_apply, Matrix.vandermonde_apply]
    rw [hA i k k.isLt, collocA_C hc i k.isLt]
  calc A = A * Matrix.vandermonde c * (Matrix.vandermonde c)⁻¹ := by
        rw [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hV, Matrix.mul_one]
    _ = collocA c := by
        rw [h, Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hV, Matrix.mul_one]

/-! ## Changing the field -/

section Map

variable {G : Type*} [Field G] (f : F →+* G)

theorem nodePoly_map (c : Fin s → F) : nodePoly (f ∘ c) = (nodePoly c).map f := by
  simp [nodePoly, Lagrange.nodal, Polynomial.map_prod]

theorem antideriv_map (p : F[X]) : antideriv (p.map f) = (antideriv p).map f := by
  ext n
  cases n with
  | zero => simp [coeff_antideriv_zero]
  | succ n => simp [coeff_antideriv_succ, map_div₀]

/-- The collocation coefficients on the image of the nodes are the images of
    the coefficients. -/
theorem collocA_map (hc : Function.Injective c) : collocA (f ∘ c) = (collocA c).map f := by
  symm
  refine collocA_unique (f.injective.comp hc) fun i k hk => ?_
  have := congrArg f (collocA_C hc i hk)
  simpa [map_sum, map_mul, map_pow, map_div₀] using this

end Map

/-- `A V = diag(c) V diag(1/(k+1))` with `V` the Vandermonde matrix of the nodes. -/
theorem collocA_mul_vandermonde (hc : Function.Injective c) :
    collocA c * Matrix.vandermonde c =
      Matrix.diagonal c * Matrix.vandermonde c *
        Matrix.diagonal (fun k : Fin s => ((k : F) + 1)⁻¹) := by
  ext i k
  rw [Matrix.mul_apply]
  simp only [Matrix.vandermonde_apply, Matrix.mul_diagonal, Matrix.diagonal_mul]
  rw [collocA_C hc i k.isLt, pow_succ]
  ring

/-- `det A = Π_i c_i / s!`. -/
theorem det_collocA [CharZero F] (hc : Function.Injective c) :
    (collocA c).det = (∏ i, c i) / (s.factorial : F) := by
  have hV : (Matrix.vandermonde c).det ≠ 0 := Matrix.det_vandermonde_ne_zero_iff.mpr hc
  have h := congrArg Matrix.det (collocA_mul_vandermonde hc)
  rw [Matrix.det_mul, Matrix.det_mul, Matrix.det_mul, Matrix.det_diagonal, Matrix.det_diagonal]
    at h
  have hfact : (∏ k : Fin s, (((k : ℕ) : F) + 1)⁻¹) = ((s.factorial : ℕ) : F)⁻¹ := by
    rw [Fin.prod_univ_eq_prod_range (fun k => ((k : F) + 1)⁻¹) s, prod_inv_distrib,
      ← prod_range_add_one_eq_factorial]
    push_cast
    rfl
  rw [hfact] at h
  field_simp at h ⊢
  exact h

theorem collocA_det_ne_zero [CharZero F] (hc : Function.Injective c) (h0 : ∀ i, c i ≠ 0) :
    (collocA c).det ≠ 0 := by
  rw [det_collocA hc]
  exact div_ne_zero (prod_ne_zero_iff.mpr fun i _ => h0 i) (by exact_mod_cast s.factorial_ne_zero)

/-- For a stiffly accurate table (`b = A[last]`), `d = bᵀA⁻¹ = e_last`. -/
theorem row_mul_inv {A : Matrix (Fin s) (Fin s) F} (hA : IsUnit A.det) (r : Fin s) :
    Matrix.vecMul (A.row r) A⁻¹ = Pi.single r 1 := by
  calc Matrix.vecMul (A.row r) A⁻¹ = Matrix.vecMul (Matrix.vecMul (Pi.single r 1) A) A⁻¹ := by
        rw [Matrix.single_one_vecMul]
    _ = Matrix.vecMul (Pi.single r 1) (A * A⁻¹) := Matrix.vecMul_vecMul _ _ _
    _ = Pi.single r 1 := by rw [Matrix.mul_nonsing_inv _ hA, Matrix.vecMul_one]

/-! ## The simplifying assumptions D, and which conditions fail -/

/-- D(s−1): if the node polynomial is orthogonal to the polynomials of degree
    `≤ s − 2` (B(2s−1)), then `Σ_i b_i c_i^(k−1) a_ij = b_j (1 − c_j^k)/k` for
    `1 ≤ k ≤ s − 1`. The difference `r_j` of the two sides satisfies
    `Σ_j r_j c_j^m = 0` for `m < s` by C(s) and B(2s−1). -/
theorem collocA_D [CharZero F] (hc : Function.Injective c) (last : Fin s) (hlast : c last = 1)
    (horth : ∀ j < s - 1, (antideriv (X ^ j * nodePoly c)).eval 1 = 0) {k : ℕ} (hk1 : 1 ≤ k)
    (hk : k ≤ s - 1) (j : Fin s) :
    ∑ i, collocA c last i * c i ^ (k - 1) * collocA c i j = collocA c last j * (1 - c j ^ k) / k := by
  set r : Fin s → F := fun j =>
    ∑ i, collocA c last i * c i ^ (k - 1) * collocA c i j - collocA c last j * (1 - c j ^ k) / k
  suffices hr : r = 0 by
    have := congrFun hr j
    simp only [r, Pi.zero_apply, sub_eq_zero] at this
    exact this
  refine Matrix.eq_zero_of_vecMul_eq_zero (M := Matrix.vandermonde c)
    (Matrix.det_vandermonde_ne_zero_iff.mpr hc) ?_
  funext m
  have hkF : (k : F) ≠ 0 := by exact_mod_cast (show k ≠ 0 by omega)
  have hm1 : ((m : ℕ) : F) + 1 ≠ 0 := by exact_mod_cast Nat.succ_ne_zero m
  have hkm : ((k : ℕ) : F) + m + 1 ≠ 0 := by exact_mod_cast Nat.succ_ne_zero (k + m)
  have h1 : ∑ j, (∑ i, collocA c last i * c i ^ (k - 1) * collocA c i j) * c j ^ (m : ℕ) =
      1 / (m + 1) * (1 / ((k - 1 + (m + 1) : ℕ) + 1)) := by
    simp only [Finset.sum_mul]
    rw [Finset.sum_comm]
    have : ∀ i, ∑ j, collocA c last i * c i ^ (k - 1) * collocA c i j * c j ^ (m : ℕ) =
        collocA c last i * c i ^ (k - 1) * (c i ^ ((m : ℕ) + 1) / ((m : ℕ) + 1)) := by
      intro i
      rw [← collocA_C hc i m.isLt, Finset.mul_sum]
      exact Finset.sum_congr rfl fun j _ => by ring
    simp only [this]
    have hB := collocA_B hc last hlast horth (k := k - 1 + (m + 1)) (by omega)
    rw [← hB, Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [pow_add (c i) (k - 1) ((m : ℕ) + 1)]
    ring
  have h2 : ∑ j, collocA c last j * (1 - c j ^ k) / k * c j ^ (m : ℕ) =
      1 / k * (1 / ((m : ℕ) + 1) - 1 / ((k + m : ℕ) + 1)) := by
    have hB1 := collocA_B hc last hlast horth (k := m) (by omega)
    have hB2 := collocA_B hc last hlast horth (k := k + m) (by omega)
    rw [← hB1, ← hB2, ← Finset.sum_sub_distrib, Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [pow_add (c i) k m]
    ring
  have hcast : ((k - 1 + ((m : ℕ) + 1) : ℕ) : F) + 1 = (k : F) + m + 1 := by
    push_cast [Nat.cast_sub hk1]; ring
  rw [hcast] at h1
  simp only [Matrix.vecMul, dotProduct, Matrix.vandermonde_apply, Pi.zero_apply, r, sub_mul,
    Finset.sum_sub_distrib, h1, h2]
  push_cast
  field_simp
  ring

/-- C(s+1) fails for nonzero nodes. Otherwise `∫_0^{c_i} M = 0` at every node,
    so `∫_0^τ M` and `τ M(τ)/(s+1)` agree at `0, c_1, …, c_s` and are equal;
    their coefficients of `τ` then give `M(0) = M(0)/(s+1)`, i.e. `M(0) = 0`. -/
theorem collocA_not_C_succ [CharZero F] (hc : Function.Injective c) (h0 : ∀ i, c i ≠ 0)
    (hs : 1 ≤ s) : ¬∀ i, ∑ j, collocA c i j * c j ^ s = c i ^ (s + 1) / (s + 1) := by
  intro H
  have hmon := (nodePoly_monic (c := c))
  have hMdeg : (nodePoly c).natDegree = s := natDegree_nodePoly
  have hMs : (nodePoly c).coeff s = 1 := by
    have := hmon.coeff_natDegree; rwa [hMdeg] at this
  have hdeg : ((X : F[X]) ^ s - nodePoly c).degree < s := by
    have := degree_sub_lt (p := (X : F[X]) ^ s) (q := nodePoly c)
      (by rw [degree_X_pow, degree_eq_natDegree hmon.ne_zero, hMdeg]) (monic_X_pow s).ne_zero
      (by rw [leadingCoeff_X_pow, hmon.leadingCoeff])
    rwa [degree_X_pow] at this
  have hM : ∀ i, (antideriv (nodePoly c)).eval (c i) = 0 := by
    intro i
    have h1 := collocA_exact hc hdeg i
    simp only [eval_sub, eval_pow, eval_X, nodePoly_eval_node, sub_zero] at h1
    rw [H i, map_sub, eval_sub, eval_antideriv_X_pow] at h1
    linear_combination h1
  set G := antideriv (nodePoly c) - C ((s : F) + 1)⁻¹ * (X * nodePoly c)
  have hG : G = 0 := by
    let f : Option (Fin s) → F := fun o => o.elim 0 c
    have hf : Function.Injective f := by
      intro a b hab
      cases a <;> cases b <;> simp only [f, Option.elim] at hab
      · rfl
      · exact absurd hab.symm (h0 _)
      · exact absurd hab (h0 _)
      · rw [hc hab]
    refine eq_zero_of_natDegree_lt_card_of_eval_eq_zero G hf (fun o => ?_) ?_
    · cases o with
      | none => simp [f, G, eval_antideriv_zero]
      | some i => simp [f, G, hM i, nodePoly_eval_node]
    · rw [Fintype.card_option, Fintype.card_fin, Nat.lt_succ_iff, natDegree_le_iff_coeff_eq_zero]
      intro N hN
      obtain ⟨n, rfl⟩ : ∃ n, N = n + 1 := ⟨N - 1, by omega⟩
      have hn : s ≤ n := by omega
      simp only [G, coeff_sub, coeff_C_mul, coeff_X_mul, coeff_antideriv_succ]
      rcases eq_or_lt_of_le hn with rfl | hlt
      · rw [hMs]; ring
      · rw [coeff_eq_zero_of_natDegree_lt (by omega)]; ring
  have h1 := congrArg (fun p => p.coeff 1) hG
  simp only [G, coeff_sub, coeff_C_mul, coeff_X_mul, coeff_zero, coeff_antideriv_succ] at h1
  have hM0 : (nodePoly c).coeff 0 ≠ 0 := by
    rw [coeff_zero_eq_eval_zero, nodePoly, Lagrange.eval_nodal]
    exact prod_ne_zero_iff.mpr fun i _ => by rw [zero_sub, neg_ne_zero]; exact h0 i
  have hs1 : (s : F) + 1 ≠ 0 := by exact_mod_cast Nat.succ_ne_zero s
  have hs0 : (s : F) ≠ 0 := by exact_mod_cast (show s ≠ 0 by omega)
  apply hM0
  field_simp at h1
  have : (nodePoly c).coeff 0 * s = 0 := by linear_combination h1
  exact (mul_eq_zero.mp this).resolve_right hs0

/-- D(s) fails when B(2s) does: with C(s) and B(s), D(s) would force
    `Σ_j b_j c_j^(2s−1) = 1/(2s)`. -/
theorem collocA_not_D [CharZero F] (hc : Function.Injective c) (last : Fin s) (hlast : c last = 1)
    (horth : ∀ j < s - 1, (antideriv (X ^ j * nodePoly c)).eval 1 = 0)
    (hres : (antideriv (X ^ (s - 1) * nodePoly c)).eval 1 ≠ 0) :
    ¬∀ j, ∑ i, collocA c last i * c i ^ (s - 1) * collocA c i j =
      collocA c last j * (1 - c j ^ s) / s := by
  intro H
  have hs : 1 ≤ s := Nat.one_le_iff_ne_zero.mpr fun h => by subst h; exact last.elim0
  have hs0 : (s : F) ≠ 0 := by exact_mod_cast (show s ≠ 0 by omega)
  have hcast : ((s - 1 : ℕ) : F) + 1 = s := by rw [Nat.cast_sub hs]; push_cast; ring
  have hL : ∑ j, (∑ i, collocA c last i * c i ^ (s - 1) * collocA c i j) * c j ^ (s - 1) =
      (∑ i, collocA c last i * c i ^ (2 * s - 1)) / s := by
    simp only [Finset.sum_mul]
    rw [Finset.sum_comm, div_eq_mul_inv, Finset.sum_mul]
    refine Finset.sum_congr rfl fun i _ => ?_
    have hC := collocA_C hc i (k := s - 1) (by omega)
    rw [hcast, show s - 1 + 1 = s by omega] at hC
    rw [show ∑ j, collocA c last i * c i ^ (s - 1) * collocA c i j * c j ^ (s - 1) =
      collocA c last i * c i ^ (s - 1) * ∑ j, collocA c i j * c j ^ (s - 1) by
        rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun j _ => by ring, hC,
      show 2 * s - 1 = (s - 1) + s by omega, pow_add]
    ring
  have hR : ∑ j, collocA c last j * (1 - c j ^ s) / s * c j ^ (s - 1) =
      (∑ j, collocA c last j * c j ^ (s - 1) - ∑ j, collocA c last j * c j ^ (2 * s - 1)) / s := by
    rw [← Finset.sum_sub_distrib, div_eq_mul_inv, Finset.sum_mul]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [show 2 * s - 1 = s + (s - 1) by omega, pow_add]
    ring
  have hB := collocA_B hc last hlast horth (k := s - 1) (by omega)
  rw [hcast] at hB
  have hB2 := collocA_B2s hc last hlast horth
  have hsum := congrArg (fun r => ∑ j, r j * c j ^ (s - 1)) (funext H)
  simp only at hsum
  rw [hL, hR, hB, hB2] at hsum
  field_simp at hsum
  have h4 : (4 * s : F) * (antideriv (X ^ (s - 1) * nodePoly c)).eval 1 = 0 := by
    linear_combination -hsum
  exact hres ((mul_eq_zero.mp h4).resolve_left (by
    have : (4 : F) ≠ 0 := by norm_num
    exact mul_ne_zero this hs0))

end Firkode.Colloc
