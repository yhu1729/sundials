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

import FirkodeSpec.Integral
import FirkodeSpec.Nodes
import Mathlib.Analysis.Polynomial.Basic

/-!
# A-stability and L-stability

Let a collocation method on real distinct nodes with `c_s = 1` have a node
polynomial `M` that is orthogonal on `[0, 1]` to the polynomials of degree
`≤ s − 2` (B(2s−1)), with `J = ∫_0^1 τ^(s−1) M(τ) dτ < 0`. Apply it to
`y′ = zy` with complex `z` and let `u` be the collocation polynomial. Then
`w = u ū` is `|u|²` on the reals and `w′(c_j) = 2 Re z |Y_j|²`. The quadrature
error formula `collocA_B_remainder` applied to `w′`, which has degree
`2s − 1` and leading coefficient `2s |u_s|²`, gives the energy identity
(`energy_identity`)

  `|Y_s|² − |y0|² = 2 Re z Σ_j b_j |Y_j|² + 2s |u_s|² J`.

The weights are `b_j ≥ 0` (`collocA_last_nonneg`), so for `Re z ≤ 0` both terms
are `≤ 0`. This gives:

* `isUnit_det_of_re_nonpos`: `I − zA` is invertible, because a nonzero
  solution of `(I − zA)Y = 0` would force `u_s = 0` and then `u′ = zu`;
* `norm_stabR_le_one`: `|R(z)| ≤ 1`, i.e. A-stability;
* `stabR_tendsto_cobounded`: `R(z) → 0` as `|z| → ∞`, since `deg N < deg D`
  in Nørsett's formula; with A-stability this is L-stability.

The Radau IIA methods satisfy the hypotheses for s ≤ 9 (`radau_A_stable`,
`radau_L_stable`): orthogonality is `NodePoly.nodePoly_orthogonal`, and
`J = −(s!(s−1)!)²/((2s)!(2s−1)!)` is `NodePoly.nodePoly_B2s_residual`.
No complex analysis is used, and neither is the E-polynomial nor the
Routh–Hurwitz criterion.
-/

namespace Firkode.Colloc

open Polynomial Finset ComplexConjugate
open scoped Matrix

variable {s : ℕ} {c : Fin s → ℝ}

/-- The nodes as complex numbers. -/
noncomputable abbrev nodesC (c : Fin s → ℝ) : Fin s → ℂ := Complex.ofRealHom ∘ c

theorem nodesC_injective (hc : Function.Injective c) : Function.Injective (nodesC c) :=
  Complex.ofReal_injective.comp hc

theorem collocA_nodesC (hc : Function.Injective c) (i j : Fin s) :
    collocA (nodesC c) i j = ((collocA c i j : ℝ) : ℂ) := by
  rw [nodesC, collocA_map Complex.ofRealHom hc]
  rfl

theorem antideriv_nodesC (k : ℕ) :
    (antideriv (X ^ k * nodePoly (nodesC c))).eval 1 =
      (((antideriv (X ^ k * nodePoly c)).eval 1 : ℝ) : ℂ) := by
  rw [nodesC, nodePoly_map, ← Polynomial.map_X Complex.ofRealHom, ← Polynomial.map_pow,
    ← Polynomial.map_mul, antideriv_map, eval_one_map]
  rfl

/-- `ū(t) = conj u(t)` at real points, for `ū = u.map conj`. -/
theorem eval_map_conj (p : ℂ[X]) (t : ℝ) :
    (p.map (starRingEnd ℂ)).eval (t : ℂ) = conj (p.eval (t : ℂ)) := by
  rw [eval_map]
  conv_lhs => rw [← Complex.conj_ofReal t]
  exact eval₂_at_apply _ _

/-- The energy identity for `y′ = zy`: under the stage equations
    `Y_i = y0 + z Σ_j a_ij Y_j`, with `u` the collocation polynomial,
    `|Y_s|² − |y0|² = 2 Re z Σ_j b_j |Y_j|² + 2s |u_s|² ∫_0^1 τ^(s−1) M`. -/
theorem energy_identity (hc : Function.Injective c) (last : Fin s) (hlast : c last = 1)
    (horth : ∀ j < s - 1, (antideriv (X ^ j * nodePoly c)).eval 1 = 0) (y0 z : ℂ)
    (Y : Fin s → ℂ) (hY : ∀ i, Y i = y0 + z * ∑ j, collocA (nodesC c) i j * Y j) :
    Complex.normSq (Y last) - Complex.normSq y0 =
      2 * z.re * ∑ j, collocA c last j * Complex.normSq (Y j) +
        2 * s * Complex.normSq ((stagePoly (nodesC c) y0 z Y).coeff s) *
          (antideriv (X ^ (s - 1) * nodePoly c)).eval 1 := by
  have hs : 1 ≤ s := Nat.one_le_iff_ne_zero.mpr fun h => by subst h; exact last.elim0
  have hcC := nodesC_injective hc
  have hlastC : nodesC c last = 1 := by simp [nodesC, hlast]
  have horthC : ∀ j < s - 1, (antideriv (X ^ j * nodePoly (nodesC c))).eval 1 = 0 :=
    fun j hj => by rw [antideriv_nodesC, horth j hj, Complex.ofReal_zero]
  set u := stagePoly (nodesC c) y0 z Y
  set ub := u.map (starRingEnd ℂ)
  have hu : u.natDegree ≤ s := natDegree_stagePoly_le hcC y0 z Y
  have hub : ub.natDegree ≤ s := natDegree_map_le.trans hu
  have hwdeg : (derivative (u * ub)).natDegree ≤ 2 * s - 1 := by
    have := natDegree_mul_le (p := u) (q := ub)
    exact (natDegree_derivative_le _).trans (by omega)
  have hquad := collocA_B_remainder hcC last hlastC horthC hwdeg
  -- the integrand at the nodes
  have hnode : ∀ j, (derivative (u * ub)).eval (nodesC c j) =
      ((2 * z.re * Complex.normSq (Y j) : ℝ) : ℂ) := by
    intro j
    have h1 : u.eval (nodesC c j) = Y j := stagePoly_eval_node hY j
    have h2 : (derivative u).eval (nodesC c j) = z * Y j :=
      derivative_stagePoly_eval_node hcC y0 z Y j
    have h3 : ub.eval (nodesC c j) = conj (Y j) := by rw [← h1]; exact eval_map_conj u (c j)
    have h4 : (derivative ub).eval (nodesC c j) = conj (z * Y j) := by
      rw [← h2, derivative_map]; exact eval_map_conj _ (c j)
    rw [derivative_mul, eval_add, eval_mul, eval_mul, h1, h2, h3, h4, map_mul]
    have : z * Y j * conj (Y j) + Y j * (conj z * conj (Y j)) = (z + conj z) * (Y j * conj (Y j)) := by
      ring
    rw [this, Complex.add_conj, Complex.mul_conj]
    push_cast
    ring
  -- the integral of `w′` is `|u(1)|² − |u(0)|²`
  have hint : (antideriv (derivative (u * ub))).eval 1 =
      ((Complex.normSq (Y last) - Complex.normSq y0 : ℝ) : ℂ) := by
    have hu1 : u.eval 1 = Y last := by rw [← hlastC]; exact stagePoly_eval_node hY last
    have hu0 : u.eval 0 = y0 := stagePoly_eval_zero y0 z Y
    have hub1 : ub.eval 1 = conj (Y last) := by
      rw [← hu1]; exact_mod_cast eval_map_conj u 1
    have hub0 : ub.eval 0 = conj y0 := by
      rw [← hu0]; exact_mod_cast eval_map_conj u 0
    rw [antideriv_derivative, eval_sub, eval_C, coeff_zero_eq_eval_zero, eval_mul, eval_mul, hu1,
      hu0, hub1, hub0, Complex.mul_conj, Complex.mul_conj]
    push_cast
    ring
  -- the leading coefficient of `w′`
  have hcoeff : (derivative (u * ub)).coeff (2 * s - 1) =
      ((2 * s * Complex.normSq (u.coeff s) : ℝ) : ℂ) := by
    rw [coeff_derivative, show 2 * s - 1 + 1 = s + s by omega,
      coeff_mul_add_eq_of_natDegree_le hu hub, coeff_map, Complex.mul_conj,
      show ((2 * s - 1 : ℕ) : ℂ) + 1 = 2 * s by rw [Nat.cast_sub (by omega)]; push_cast; ring]
    push_cast
    ring
  have hsum : ∑ j, collocA (nodesC c) last j * (derivative (u * ub)).eval (nodesC c j) =
      ((2 * z.re * ∑ j, collocA c last j * Complex.normSq (Y j) : ℝ) : ℂ) := by
    simp only [collocA_nodesC hc, hnode]
    push_cast
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun j _ => by ring
  rw [hsum, hint, hcoeff, antideriv_nodesC] at hquad
  have key : ((Complex.normSq (Y last) - Complex.normSq y0 : ℝ) : ℂ) =
      ((2 * z.re * ∑ j, collocA c last j * Complex.normSq (Y j) +
        2 * s * Complex.normSq (u.coeff s) * (antideriv (X ^ (s - 1) * nodePoly c)).eval 1 : ℝ) :
          ℂ) := by
    push_cast at hquad ⊢
    linear_combination -hquad
  exact_mod_cast key

/-- For `Re z ≤ 0` the stage system is uniquely solvable: `I − zA` is
    invertible. -/
theorem isUnit_det_of_re_nonpos (hc : Function.Injective c) (last : Fin s) (hlast : c last = 1)
    (horth : ∀ j < s - 1, (antideriv (X ^ j * nodePoly c)).eval 1 = 0)
    (hJ : (antideriv (X ^ (s - 1) * nodePoly c)).eval 1 < 0) {z : ℂ} (hz : z.re ≤ 0) :
    IsUnit (1 - z • collocA (nodesC c)).det := by
  have hs : 1 ≤ s := Nat.one_le_iff_ne_zero.mpr fun h => by subst h; exact last.elim0
  rcases eq_or_ne z 0 with rfl | hz0
  · simp
  rw [isUnit_iff_ne_zero]
  intro hdet
  obtain ⟨Y, hY0, hYk⟩ := Matrix.exists_mulVec_eq_zero_iff.mpr hdet
  have hY : ∀ i, Y i = 0 + z * ∑ j, collocA (nodesC c) i j * Y j := by
    intro i
    have := congrFun hYk i
    rw [Matrix.sub_mulVec, Matrix.one_mulVec, Matrix.smul_mulVec] at this
    simp only [Pi.sub_apply, Pi.smul_apply, smul_eq_mul, Pi.zero_apply] at this
    rw [show ∑ j, collocA (nodesC c) i j * Y j = (collocA (nodesC c) *ᵥ Y) i from rfl]
    linear_combination this
  have hE := energy_identity hc last hlast horth 0 z Y hY
  set u := stagePoly (nodesC c) 0 z Y
  have hb : 0 ≤ ∑ j, collocA c last j * Complex.normSq (Y j) := Finset.sum_nonneg fun j _ =>
    mul_nonneg (collocA_last_nonneg hc last hlast horth j) (Complex.normSq_nonneg _)
  have h1 : 2 * z.re * ∑ j, collocA c last j * Complex.normSq (Y j) ≤ 0 :=
    mul_nonpos_of_nonpos_of_nonneg (by linarith) hb
  rw [Complex.normSq_zero, sub_zero] at hE
  -- the leading coefficient vanishes
  have hus0 : u.coeff s = 0 := by
    by_contra hne
    have hpos : 0 < 2 * (s : ℝ) * Complex.normSq (u.coeff s) := by
      have := Complex.normSq_pos.mpr hne
      positivity
    have := mul_neg_of_pos_of_neg hpos hJ
    linarith [Complex.normSq_nonneg (Y last)]
  have hcC := nodesC_injective hc
  have hudeg : u.natDegree < s := by
    have hle := natDegree_stagePoly_le hcC 0 z Y
    by_contra h
    have heq : u.natDegree = s := le_antisymm hle (not_lt.mp h)
    have hu0 : u ≠ 0 := by intro h0; rw [h0, natDegree_zero] at heq; omega
    exact hu0 (leadingCoeff_eq_zero.mp (by rw [leadingCoeff, heq, hus0]))
  -- so `u′ − z u` has degree `< s` and vanishes at the nodes
  have hq : derivative u - C z * u = 0 := by
    refine Polynomial.eq_zero_of_natDegree_lt_card_of_eval_eq_zero _ hcC (fun i => ?_) ?_
    · rw [eval_sub, eval_mul, eval_C, derivative_stagePoly_eval_node hcC,
        stagePoly_eval_node hY, sub_self]
    · rw [Fintype.card_fin]
      calc (derivative u - C z * u).natDegree
          ≤ max (derivative u).natDegree (C z * u).natDegree := natDegree_sub_le _ _
        _ < s := max_lt ((natDegree_derivative_le u).trans_lt (by omega))
            ((natDegree_C_mul_le _ _).trans_lt hudeg)
  have hu0 : u = 0 := eq_zero_of_derivative_eq hz0 (sub_eq_zero.mp hq)
  apply hY0
  funext i
  rw [← stagePoly_eval_node hY i, Pi.zero_apply]
  change u.eval _ = 0
  rw [hu0, eval_zero]

/-- A-stability: `|R(z)| ≤ 1` for `Re z ≤ 0`. -/
theorem norm_stabR_le_one (hc : Function.Injective c) (last : Fin s) (hlast : c last = 1)
    (horth : ∀ j < s - 1, (antideriv (X ^ j * nodePoly c)).eval 1 = 0)
    (hJ : (antideriv (X ^ (s - 1) * nodePoly c)).eval 1 < 0) {z : ℂ} (hz : z.re ≤ 0) :
    ‖stabR (collocA (nodesC c)) ((collocA (nodesC c)).row last) z‖ ≤ 1 := by
  have hunit := isUnit_det_of_re_nonpos hc last hlast horth hJ hz
  set A := collocA (nodesC c)
  set Y : Fin s → ℂ := (1 - z • A)⁻¹ *ᵥ fun _ => 1
  have hY : ∀ i, Y i = 1 + z * ∑ j, A i j * Y j := stage_eq_of_isUnit hunit 1
  have hR : stabR A (A.row last) z = Y last := by rw [hY last]; rfl
  have hE := energy_identity hc last hlast horth 1 z Y hY
  have hb : 0 ≤ ∑ j, collocA c last j * Complex.normSq (Y j) := Finset.sum_nonneg fun j _ =>
    mul_nonneg (collocA_last_nonneg hc last hlast horth j) (Complex.normSq_nonneg _)
  have h1 : 2 * z.re * ∑ j, collocA c last j * Complex.normSq (Y j) ≤ 0 :=
    mul_nonpos_of_nonpos_of_nonneg (by linarith) hb
  have h2 : 2 * (s : ℝ) * Complex.normSq ((stagePoly (nodesC c) 1 z Y).coeff s) *
      (antideriv (X ^ (s - 1) * nodePoly c)).eval 1 ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos (mul_nonneg (by positivity) (Complex.normSq_nonneg _)) hJ.le
  rw [Complex.normSq_one] at hE
  rw [hR, ← sq_le_one_iff₀ (norm_nonneg _), ← Complex.normSq_eq_norm_sq]
  linarith

/-- `R(z) → 0` as `|z| → ∞`: with A-stability, L-stability. -/
theorem stabR_tendsto_cobounded (hc : Function.Injective c) (h0 : ∀ i, c i ≠ 0) (last : Fin s)
    (hlast : c last = 1) :
    Filter.Tendsto (stabR (collocA (nodesC c)) ((collocA (nodesC c)).row last))
      (Bornology.cobounded ℂ) (nhds 0) := by
  have hcC := nodesC_injective hc
  have h0C : ∀ i, nodesC c i ≠ 0 := fun i => by simpa [nodesC] using h0 i
  have hlastC : nodesC c last = 1 := by simp [nodesC, hlast]
  have hdeg := degree_norsettNPoly_lt (c := nodesC c) last hlastC h0C
  have hD : norsettDPoly (nodesC c) ≠ 0 := by
    intro h; rw [h, degree_zero] at hdeg; exact not_lt_bot hdeg
  have hQ : (collocA (nodesC c)).charpolyRev ≠ 0 := by
    intro h
    have := Matrix.eval_charpolyRev (M := collocA (nodesC c))
    rw [h, eval_zero] at this
    exact zero_ne_one this
  have hlim := (Polynomial.isLittleO_cobounded_of_degree_lt hdeg).tendsto_div_nhds_zero
  refine hlim.congr' ?_
  have hb1 := Bornology.isBounded_def.mp (Polynomial.finite_setOf_isRoot hD).isBounded
  have hb2 := Bornology.isBounded_def.mp (Polynomial.finite_setOf_isRoot hQ).isBounded
  filter_upwards [hb1, hb2] with z hzD hzQ
  simp only [Set.mem_compl_iff, Set.mem_setOf_eq, IsRoot.def] at hzD hzQ
  have hz : IsUnit (1 - z • collocA (nodesC c)).det := by
    rw [← Colloc.eval_charpolyRev]; exact isUnit_iff_ne_zero.mpr hzQ
  have h := norsett hcC last hlastC z hz
  rw [eval_norsettNPoly, eval_norsettDPoly, ← h, mul_div_cancel_right₀]
  rwa [← eval_norsettDPoly]

/-! ## The Radau IIA methods -/

theorem radau_orth (hs : s ∈ NodePoly.stages) :
    ∀ j < s - 1, (antideriv (X ^ j * nodePoly (Nodes.radauNodes s))).eval 1 = 0 := fun j hj => by
  rw [Nodes.radauNodes_nodePoly hs, ← Radau.ratPoly_X, ← Radau.ratPoly_pow, ← Radau.ratPoly_mul,
    ← Radau.integral01_eq, NodePoly.nodePoly_orthogonal s hs j hj, Rat.cast_zero]

theorem radau_J_neg (hs : s ∈ NodePoly.stages) :
    (antideriv (X ^ (s - 1) * nodePoly (Nodes.radauNodes s))).eval 1 < 0 := by
  rw [Nodes.radauNodes_nodePoly hs, ← Radau.ratPoly_X, ← Radau.ratPoly_pow, ← Radau.ratPoly_mul,
    ← Radau.integral01_eq, NodePoly.nodePoly_B2s_residual s hs]
  push_cast [Radau.fact_eq]
  exact div_neg_of_neg_of_pos (neg_neg_of_pos (by positivity)) (by positivity)

/-- The Radau IIA methods are A-stable (s ≤ 9): for `Re z ≤ 0`, `I − zA` is
    invertible and `|R(z)| ≤ 1`. -/
theorem radau_A_stable (hs : s ∈ NodePoly.stages) (last : Fin s) (hlast : (last : ℕ) = s - 1)
    {z : ℂ} (hz : z.re ≤ 0) :
    IsUnit (1 - z • collocA (nodesC (Nodes.radauNodes s))).det ∧
      ‖stabR (collocA (nodesC (Nodes.radauNodes s)))
        ((collocA (nodesC (Nodes.radauNodes s))).row last) z‖ ≤ 1 :=
  ⟨isUnit_det_of_re_nonpos (Nodes.radauNodes_injective hs) last
      (Nodes.radauNodes_last hs last hlast) (radau_orth hs) (radau_J_neg hs) hz,
    norm_stabR_le_one (Nodes.radauNodes_injective hs) last (Nodes.radauNodes_last hs last hlast)
      (radau_orth hs) (radau_J_neg hs) hz⟩

/-- The Radau IIA methods are L-stable (s ≤ 9): A-stable, and `R(z) → 0` as
    `|z| → ∞`. -/
theorem radau_L_stable (hs : s ∈ NodePoly.stages) (last : Fin s) (hlast : (last : ℕ) = s - 1) :
    Filter.Tendsto (stabR (collocA (nodesC (Nodes.radauNodes s)))
      ((collocA (nodesC (Nodes.radauNodes s))).row last)) (Bornology.cobounded ℂ) (nhds 0) :=
  stabR_tendsto_cobounded (Nodes.radauNodes_injective hs)
    (fun i => (Nodes.radauNodes_pos hs i).ne') last (Nodes.radauNodes_last hs last hlast)

end Firkode.Colloc
