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

import FirkodeSpec.Nodes
import Mathlib.LinearAlgebra.Matrix.SchurComplement
import Mathlib.LinearAlgebra.Matrix.ToLinearEquiv
import Mathlib.RingTheory.Coprime.Basic

/-!
# The characteristic polynomial of the Radau IIA matrix (s ≤ 9)

For a stiffly accurate table the stability function is
`R(z) = det(I − z(A − 𝟙bᵀ)) / det(I − zA)` (`det_stab`, the matrix
determinant lemma), and Nørsett's formula gives `R = N/D`. Clearing
denominators, `P · D = Q · N` as polynomials, with `Q(z) = det(I − zA)` and
`P(z) = det(I − z(A − 𝟙bᵀ))` (`charpolyRev_norsett`). For the Radau nodes
`N` and `D` are multiples of the Padé polynomials, which are coprime
(`NodePoly.pade_coprime`), and `deg Q ≤ s`, `Q(0) = 1`; so

* `det(I − zA)` is the Padé denominator `Q_s(z)` (`radau_charpolyRev`,
  `radau_det`) and `det(I − z(A − 𝟙bᵀ))` the numerator
  (`radau_charpolyRev_num`);
* the real eigenvalues of `A` are the real roots of `det(μI − A) = μ^s Q_s(1/μ)`
  (`radau_isEigen_iff`).

With the root counts of `RootCount` this proves the comment of
`firk_gamma0`: for odd `s` the matrix `A⁻¹` has exactly one real eigenvalue
`U1`, and `firk_gamma0[s]` is `1/U1` to within half a unit in its last printed
digit (`radau_gamma0_odd`, `radau_gamma0_odd_inv`); for even `s` neither `A`
nor `A⁻¹` has a real eigenvalue (`radau_even_no_real_eigenvalue`).
-/

namespace Firkode.Charpoly

open Polynomial Finset Matrix Radau Nodes Colloc

variable {s : ℕ}

/-! ## The stability function as a quotient of determinants -/

/-- The matrix `𝟙 bᵀ`. -/
noncomputable abbrev onesRow {F : Type*} [Field F] (b : Fin s → F) : Matrix (Fin s) (Fin s) F :=
  replicateCol (Fin 1) (fun _ => (1 : F)) * replicateRow (Fin 1) b

/-- The matrix determinant lemma for the stability function:
    `det(I − z(A − 𝟙bᵀ)) = det(I − zA) R(z)` wherever `I − zA` is
    invertible. -/
theorem det_stab {F : Type*} [Field F] (A : Matrix (Fin s) (Fin s) F) (b : Fin s → F) (z : F)
    (hz : IsUnit (1 - z • A).det) :
    (1 - z • (A - onesRow b)).det = (1 - z • A).det * stabR A b z := by
  have hsplit : 1 - z • (A - onesRow b) =
      (1 - z • A) + replicateCol (Fin 1) (fun _ => z) * replicateRow (Fin 1) b := by
    ext i j
    simp only [onesRow, Matrix.sub_apply, Matrix.add_apply, Matrix.smul_apply, Matrix.mul_apply,
      replicateCol_apply, replicateRow_apply, smul_eq_mul, Finset.univ_unique, Finset.sum_singleton]
    ring
  rw [hsplit, Matrix.det_add_mul _ _ hz, Matrix.det_fin_one]
  congr 1
  simp only [Matrix.add_apply, Matrix.one_apply_eq, Matrix.mul_apply, replicateRow_apply,
    replicateCol_apply, stabR, dotProduct, mulVec, Finset.mul_sum, Finset.sum_mul]
  rw [Finset.sum_comm]
  exact congrArg _ (Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by ring)

/-- Nørsett's formula with the denominators cleared: `P · D = Q · N` for
    collocation on distinct nodes with `c_s = 1`. -/
theorem charpolyRev_norsett {c : Fin s → ℝ} (hc : Function.Injective c) (last : Fin s)
    (hlast : c last = 1) :
    (collocA c - onesRow ((collocA c).row last)).charpolyRev * norsettDPoly c =
      (collocA c).charpolyRev * norsettNPoly c := by
  have hQ : (collocA c).charpolyRev ≠ 0 := by
    intro h
    have := Matrix.eval_charpolyRev (M := collocA c)
    rw [h, eval_zero] at this
    exact zero_ne_one this
  refine Polynomial.eq_of_infinite_eval_eq _ _
    ((Polynomial.finite_setOf_isRoot hQ).infinite_compl.mono fun z hz => ?_)
  simp only [Set.mem_compl_iff, Set.mem_setOf_eq, IsRoot.def] at hz
  have hunit : IsUnit (1 - z • collocA c).det := by
    rw [← Colloc.eval_charpolyRev]; exact isUnit_iff_ne_zero.mpr hz
  simp only [Set.mem_setOf_eq, eval_mul, Colloc.eval_charpolyRev, eval_norsettDPoly,
    eval_norsettNPoly]
  rw [det_stab _ _ z hunit, mul_assoc, norsett hc last hlast z hunit]

/-! ## The Radau tables -/

theorem radau_norsettDPoly (hs : s ∈ NodePoly.stages) :
    norsettDPoly (radauNodes s) = C (fact s : ℝ) * ratPoly (padeDen s) :=
  Polynomial.funext fun z => by
    rw [eval_norsettDPoly, radau_norsettD hs (radauNodes_nodePoly hs), eval_mul, eval_C]

theorem radau_norsettNPoly (hs : s ∈ NodePoly.stages) :
    norsettNPoly (radauNodes s) = C (fact s : ℝ) * ratPoly (padeNum s) :=
  Polynomial.funext fun z => by
    rw [eval_norsettNPoly, radau_norsettN hs (radauNodes_nodePoly hs), eval_mul, eval_C]

/-- The Padé numerator and denominator are coprime (`NodePoly.pade_coprime`). -/
theorem pade_isCoprime (hs : s ∈ NodePoly.stages) :
    IsCoprime (ratPoly (padeNum s)) (ratPoly (padeDen s)) := by
  have h := congrArg ratPoly (NodePoly.pade_coprime s hs)
  simp only [ratPoly, toPoly_trim, toPoly_add, toPoly_mul, Polynomial.map_add,
    Polynomial.map_mul] at h
  exact ⟨ratPoly (NodePoly.padeBezout s).1, ratPoly (NodePoly.padeBezout s).2, by
    simpa [ratPoly] using h⟩

theorem natDegree_padeDen (hs : s ∈ NodePoly.stages) : (ratPoly (padeDen s)).natDegree = s := by
  have hlen : (padeDen s).length = s + 1 := by simp [padeDen]
  have hle : (ratPoly (padeDen s)).natDegree < s + 1 := by
    rw [← hlen]
    exact natDegree_map_le.trans_lt (natDegree_toPoly_lt (List.ne_nil_of_length_eq_add_one hlen))
  refine le_antisymm (by omega) (le_natDegree_of_ne_zero ?_)
  rw [RuleOfSigns.coeff_ratPoly]
  exact_mod_cast (NodePoly.padeDen_ends s hs).2

theorem eval_zero_padeDen (hs : s ∈ NodePoly.stages) : (ratPoly (padeDen s)).eval 0 = 1 := by
  rw [← coeff_zero_eq_eval_zero, RuleOfSigns.coeff_ratPoly, (NodePoly.padeDen_ends s hs).1,
    Rat.cast_one]

theorem natDegree_charpolyRev_le (A : Matrix (Fin s) (Fin s) ℝ) : A.charpolyRev.natDegree ≤ s := by
  rw [← Matrix.reverse_charpoly]
  exact (Polynomial.reverse_natDegree_le _).trans
    (by rw [Matrix.charpoly_natDegree_eq_dim, Fintype.card_fin])

/-- `P · Q_s = Q · P_s` for the Radau tables. -/
theorem radau_charpolyRev_mul (hs : s ∈ NodePoly.stages) (last : Fin s)
    (hlast : (last : ℕ) = s - 1) :
    (collocA (radauNodes s) - onesRow ((collocA (radauNodes s)).row last)).charpolyRev *
        ratPoly (padeDen s) =
      (collocA (radauNodes s)).charpolyRev * ratPoly (padeNum s) := by
  have h := charpolyRev_norsett (radauNodes_injective hs) last (radauNodes_last hs last hlast)
  rw [radau_norsettDPoly hs, radau_norsettNPoly hs, mul_left_comm, mul_left_comm _ (C _)] at h
  exact mul_left_cancel₀ (C_ne_zero.mpr (fact_ne_zero_real s)) h

/-- For s ≤ 9, `det(I − zA)` is the denominator `Q_s` of the `(s−1, s)` Padé
    approximant of `e^z`. -/
theorem radau_charpolyRev (hs : s ∈ NodePoly.stages) :
    (collocA (radauNodes s)).charpolyRev = ratPoly (padeDen s) := by
  have hs0 : 0 < s := by simp [NodePoly.stages] at hs; omega
  have hQ0 : (collocA (radauNodes s)).charpolyRev.eval 0 = 1 := Matrix.eval_charpolyRev
  have hQ : (collocA (radauNodes s)).charpolyRev ≠ 0 := by
    intro h; rw [h, eval_zero] at hQ0; exact zero_ne_one hQ0
  have hid := radau_charpolyRev_mul hs ⟨s - 1, by omega⟩ rfl
  have hdvd : ratPoly (padeDen s) ∣ (collocA (radauNodes s)).charpolyRev * ratPoly (padeNum s) :=
    ⟨_, by rw [← hid, mul_comm]⟩
  obtain ⟨w, hw⟩ := (pade_isCoprime hs).symm.dvd_of_dvd_mul_right hdvd
  have hD0 : ratPoly (padeDen s) ≠ 0 := fun h => by
    have := eval_zero_padeDen hs; rw [h, eval_zero] at this; exact zero_ne_one this
  have hw0 : w ≠ 0 := by rintro rfl; rw [mul_zero] at hw; exact hQ hw
  have hdeg := natDegree_charpolyRev_le (collocA (radauNodes s))
  rw [hw, natDegree_mul hD0 hw0, natDegree_padeDen hs] at hdeg
  have hwc : w = C (w.coeff 0) := eq_C_of_natDegree_eq_zero (by omega)
  have h1 : w.coeff 0 = 1 := by
    rw [hw, eval_mul, eval_zero_padeDen hs, one_mul, hwc, eval_C] at hQ0
    exact hQ0
  rw [hw, hwc, h1, C_1, mul_one]

/-- `det(I − z(A − 𝟙bᵀ))` is the Padé numerator `P_s`, so `R = P_s/Q_s`. -/
theorem radau_charpolyRev_num (hs : s ∈ NodePoly.stages) (last : Fin s)
    (hlast : (last : ℕ) = s - 1) :
    (collocA (radauNodes s) - onesRow ((collocA (radauNodes s)).row last)).charpolyRev =
      ratPoly (padeNum s) := by
  have h := radau_charpolyRev_mul hs last hlast
  rw [radau_charpolyRev hs, mul_comm (ratPoly (padeDen s))] at h
  have hD0 : ratPoly (padeDen s) ≠ 0 := fun h0 => by
    have := eval_zero_padeDen hs; rw [h0, eval_zero] at this; exact zero_ne_one this
  exact mul_right_cancel₀ hD0 h

theorem radau_det (hs : s ∈ NodePoly.stages) (z : ℝ) :
    (1 - z • collocA (radauNodes s)).det = (ratPoly (padeDen s)).eval z := by
  rw [← Colloc.eval_charpolyRev, radau_charpolyRev hs]

/-! ## Eigenvalues -/

/-- `μ` is a real eigenvalue of `M`. -/
def IsEigen (M : Matrix (Fin s) (Fin s) ℝ) (μ : ℝ) : Prop :=
  ∃ v : Fin s → ℝ, v ≠ 0 ∧ M *ᵥ v = μ • v

theorem isEigen_iff_det (M : Matrix (Fin s) (Fin s) ℝ) (μ : ℝ) :
    IsEigen M μ ↔ (μ • (1 : Matrix (Fin s) (Fin s) ℝ) - M).det = 0 := by
  rw [← Matrix.exists_mulVec_eq_zero_iff]
  simp only [IsEigen, sub_mulVec, smul_mulVec, one_mulVec, sub_eq_zero]
  exact exists_congr fun v => and_congr_right fun _ => eq_comm

/-- The eigenvalues of `M⁻¹` are the inverses of those of `M`. -/
theorem isEigen_inv_iff {M : Matrix (Fin s) (Fin s) ℝ} (hM : IsUnit M.det) (l : ℝ) :
    IsEigen M⁻¹ l ↔ l ≠ 0 ∧ IsEigen M l⁻¹ := by
  constructor
  · rintro ⟨v, hv, h⟩
    have hback : M *ᵥ (M⁻¹ *ᵥ v) = v := by rw [mulVec_mulVec, mul_nonsing_inv _ hM, one_mulVec]
    have hl : l ≠ 0 := by
      rintro rfl
      rw [h, zero_smul, mulVec_zero] at hback
      exact hv hback.symm
    refine ⟨hl, v, hv, ?_⟩
    rw [h, mulVec_smul] at hback
    calc M *ᵥ v = l⁻¹ • l • M *ᵥ v := by rw [smul_smul, inv_mul_cancel₀ hl, one_smul]
      _ = l⁻¹ • v := by rw [hback]
  · rintro ⟨hl, v, hv, h⟩
    refine ⟨v, hv, ?_⟩
    have : M⁻¹ *ᵥ (M *ᵥ v) = v := by rw [mulVec_mulVec, nonsing_inv_mul _ hM, one_mulVec]
    rw [h, mulVec_smul] at this
    calc M⁻¹ *ᵥ v = l • l⁻¹ • M⁻¹ *ᵥ v := by rw [smul_smul, mul_inv_cancel₀ hl, one_smul]
      _ = l • v := by rw [this]

theorem get_reverse (l : Poly ℚ) {k : ℕ} (hk : k < l.length) :
    get l.reverse k = get l (l.length - 1 - k) := by
  simp only [get, List.getD_eq_getElem?_getD, List.getElem?_reverse hk]

/-- `x^n p(1/x)` for the reversed coefficient list. -/
theorem eval_ratPoly_reverse (l : Poly ℚ) {x : ℝ} (hx : x ≠ 0) :
    (ratPoly l.reverse).eval x = x ^ (l.length - 1) * (ratPoly l).eval x⁻¹ := by
  rw [ratPoly_eval, ratPoly_eval, List.length_reverse, Finset.mul_sum, ← Finset.sum_range_reflect]
  refine Finset.sum_congr rfl fun j hj => ?_
  have hj' := Finset.mem_range.mp hj
  rw [get_reverse l (by omega), show l.length - 1 - (l.length - 1 - j) = j by omega,
    pow_sub₀ x hx (by omega : j ≤ l.length - 1), inv_pow]
  ring

theorem charPoly_eval_zero_ne (hs : s ∈ NodePoly.stages) :
    (ratPoly (Gamma0.charPoly s)).eval 0 ≠ 0 := by
  have hlen : (padeDen s).length = s + 1 := by simp [padeDen]
  rw [← coeff_zero_eq_eval_zero, RuleOfSigns.coeff_ratPoly, Gamma0.charPoly,
    get_reverse _ (by omega), hlen, show s + 1 - 1 - 0 = s by omega]
  exact_mod_cast (NodePoly.padeDen_ends s hs).2

/-- For s ≤ 9 the real eigenvalues of `A` are the real roots of
    `det(μI − A) = μ^s Q_s(1/μ)`, the polynomial of the `firk_gamma0`
    certificates. -/
theorem radau_isEigen_iff (hs : s ∈ NodePoly.stages) (μ : ℝ) :
    IsEigen (collocA (radauNodes s)) μ ↔ (ratPoly (Gamma0.charPoly s)).eval μ = 0 := by
  have hlen : (padeDen s).length = s + 1 := by simp [padeDen]
  rw [isEigen_iff_det]
  by_cases hμ : μ = 0
  · subst hμ
    have hA := collocA_det_ne_zero (radauNodes_injective hs) fun i => (radauNodes_pos hs i).ne'
    have h1 : ((0 : ℝ) • (1 : Matrix (Fin s) (Fin s) ℝ) - collocA (radauNodes s)).det ≠ 0 := by
      rw [zero_smul, zero_sub, det_neg]
      exact mul_ne_zero (pow_ne_zero _ (by norm_num)) hA
    exact ⟨fun h => absurd h h1, fun h => absurd h (charPoly_eval_zero_ne hs)⟩
  · have hsplit : μ • (1 : Matrix (Fin s) (Fin s) ℝ) - collocA (radauNodes s) =
        μ • (1 - μ⁻¹ • collocA (radauNodes s)) := by
      rw [smul_sub, smul_smul, mul_inv_cancel₀ hμ, one_smul]
    rw [hsplit, det_smul, Fintype.card_fin, radau_det hs, Gamma0.charPoly,
      eval_ratPoly_reverse _ hμ, hlen, Nat.add_sub_cancel]

/-! ## The `firk_gamma0` literals -/

theorem mem_stages_of_odd {s : ℕ} (hs : s ∈ [3, 5, 7, 9]) : s ∈ NodePoly.stages := by
  simp only [List.mem_cons, List.not_mem_nil, or_false] at hs
  rcases hs with rfl | rfl | rfl | rfl <;> decide

theorem mem_stages_of_even {s : ℕ} (hs : s ∈ [2, 4, 6, 8]) : s ∈ NodePoly.stages := by
  simp only [List.mem_cons, List.not_mem_nil, or_false] at hs
  rcases hs with rfl | rfl | rfl | rfl <;> decide

/-- For s = 3, 5, 7, 9 the Radau IIA matrix `A` has exactly one real
    eigenvalue `μ`, and `firk_gamma0[s]` is within half a unit in its last
    printed digit of it. -/
theorem radau_gamma0_odd :
    ∀ s ∈ [3, 5, 7, 9], ∃ μ : ℝ, IsEigen (collocA (radauNodes s)) μ ∧
      |μ - (Gamma0.literal s : ℝ)| < (Gamma0.ulp s : ℝ) / 2 ∧
        ∀ ν : ℝ, IsEigen (collocA (radauNodes s)) ν → ν = μ := by
  intro s hs
  have hst := mem_stages_of_odd hs
  obtain ⟨μ, h1, h2, h3, h4⟩ := RootCount.gamma0_odd_unique s hs
  refine ⟨μ, (radau_isEigen_iff hst μ).mpr h3, ?_,
    fun ν hν => h4 ν ((radau_isEigen_iff hst ν).mp hν)⟩
  push_cast at h1 h2
  rw [abs_sub_lt_iff]
  constructor <;> linarith

/-- The comment of `firk_gamma0`: for odd s, `A⁻¹` has a unique real eigenvalue
    `U1`, and `firk_gamma0[s]` is `1/U1` to within half a unit in its last
    printed digit. -/
theorem radau_gamma0_odd_inv :
    ∀ s ∈ [3, 5, 7, 9], ∃ U1 : ℝ, IsEigen (collocA (radauNodes s))⁻¹ U1 ∧
      (∀ l : ℝ, IsEigen (collocA (radauNodes s))⁻¹ l → l = U1) ∧
        |1 / U1 - (Gamma0.literal s : ℝ)| < (Gamma0.ulp s : ℝ) / 2 := by
  intro s hs
  have hst := mem_stages_of_odd hs
  have hA : IsUnit (collocA (radauNodes s)).det := isUnit_iff_ne_zero.mpr
    (collocA_det_ne_zero (radauNodes_injective hst) fun i => (radauNodes_pos hst i).ne')
  obtain ⟨μ, hμ, hlit, huniq⟩ := radau_gamma0_odd s hs
  have hμ0 : μ ≠ 0 := by
    rintro rfl
    exact charPoly_eval_zero_ne hst ((radau_isEigen_iff hst 0).mp hμ)
  refine ⟨μ⁻¹, (isEigen_inv_iff hA _).mpr ⟨inv_ne_zero hμ0, by rwa [inv_inv]⟩, fun l hl => ?_, ?_⟩
  · obtain ⟨hl0, hl⟩ := (isEigen_inv_iff hA l).mp hl
    rw [← inv_inv l, huniq _ hl]
  · rwa [one_div, inv_inv]

/-- For even s, neither `A` nor `A⁻¹` has a real eigenvalue, as the comment of
    `firk_gamma0` says. -/
theorem radau_even_no_real_eigenvalue :
    ∀ s ∈ [2, 4, 6, 8], ∀ μ : ℝ, ¬IsEigen (collocA (radauNodes s)) μ ∧
      ¬IsEigen (collocA (radauNodes s))⁻¹ μ := by
  intro s hs μ
  have hst := mem_stages_of_even hs
  have hA : IsUnit (collocA (radauNodes s)).det := isUnit_iff_ne_zero.mpr
    (collocA_det_ne_zero (radauNodes_injective hst) fun i => (radauNodes_pos hst i).ne')
  have hno : ∀ ν, ¬IsEigen (collocA (radauNodes s)) ν := fun ν h =>
    RootCount.gamma0_even_no_root s hs ν ((radau_isEigen_iff hst ν).mp h)
  exact ⟨hno μ, fun h => hno _ ((isEigen_inv_iff hA μ).mp h).2⟩

end Firkode.Charpoly
