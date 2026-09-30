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

import FirkodeCert.NodePoly
import FirkodeSpec.Bridge
import FirkodeSpec.Collocation
import FirkodeSpec.ErrorEstimate
import FirkodeSpec.Stability
import Mathlib.Data.Rat.BigOperators
import Mathlib.RingTheory.Polynomial.ShiftedLegendre

/-!
# The Radau IIA tables of FIRKODE over the reals

* For s ≤ 3 the closed-form tables of `firk_radau_closed_form`, read as real
  numbers, are the collocation tables on their nodes (`closedForm_eq_collocA`),
  so every theorem of `Collocation` applies to them.
* For s ≤ 9 the tables that `firk_radau_build` computes have nodes that are
  the roots of the Radau node polynomial `M_s` of `FirkodeCert.NodePoly`.
  For real nodes with this node polynomial (`Nodes` proves that they exist),
  the table satisfies B(2s−1) (`radau_B`) and the B(2s) residual is
  `(s!(s−1)!)²/((2s)!(2s−1)!)` (`radau_B2s_residual`).
-/

namespace Firkode

open Polynomial Finset Colloc

/-! ## Loops of `FIRKodeTable_CheckOrder` -/

theorem lastTrue_go_spec (cond : Nat → Bool) :
    ∀ fuel k, 1 ≤ k → ∀ j, k ≤ j → j ≤ Mirror.lastTrue.go cond fuel k → cond j = true := by
  intro fuel
  induction fuel with
  | zero => intro k hk j h1 h2; simp [Mirror.lastTrue.go] at h2; omega
  | succ fuel ih =>
    intro k hk j h1 h2
    simp only [Mirror.lastTrue.go] at h2
    split at h2
    · rcases Nat.eq_or_lt_of_le h1 with rfl | h
      · assumption
      · exact ih (k + 1) (by omega) j h h2
    · omega

/-- If `lastTrue cond kmax = m` then `cond 1, …, cond m` all hold. -/
theorem lastTrue_spec (cond : Nat → Bool) (kmax : Nat) :
    ∀ k, 1 ≤ k → k ≤ Mirror.lastTrue cond kmax → cond k = true :=
  fun k h1 h2 => lastTrue_go_spec cond kmax 1 le_rfl k h1 h2

theorem toReal_sum_map_range (n : Nat) (f : Nat → Q6) :
    (sum ((List.range n).map f)).toReal = ∑ j ∈ Finset.range n, (f j).toReal := by
  rw [sum_toReal, List.map_map]
  induction n with
  | zero => simp
  | succ n ih => simp [List.range_succ, Finset.sum_range_succ, ← ih]

/-! ## Closed-form tables, s ≤ 3 -/

namespace Radau

open Mirror Generated Tables

/-- The nodes of the closed-form table for `s` stages, as real numbers. -/
noncomputable def closedC (s : Nat) : Fin s → ℝ := fun i => (get (closedFormC s) i).toReal

/-- The coefficients of the closed-form table for `s` stages, as real numbers. -/
noncomputable def closedA (s : Nat) : Matrix (Fin s) (Fin s) ℝ :=
  Matrix.of fun i j => (get2 (closedFormA s) i j).toReal

theorem closedForm_length : ∀ s ∈ closedStages, (closedFormC s).length = s := by decide +kernel

theorem closedForm_distinct :
    ∀ s ∈ closedStages, ∀ i < s, ∀ j < s, i ≠ j → get (closedFormC s) i ≠ get (closedFormC s) j := by
  decide +kernel

theorem closedC_injective {s : Nat} (hs : s ∈ closedStages) : Function.Injective (closedC s) := by
  intro i j h
  by_contra hij
  exact closedForm_distinct s hs i i.isLt j j.isLt (fun h' => hij (Fin.ext h'))
    (Q6.toReal_injective h)

/-- The C(k) conditions that `FIRKodeTable_CheckOrder` verifies for the
    closed-form tables, as equalities of real numbers. -/
theorem closedForm_C_real {s : Nat} (hs : s ∈ closedStages) (i : Fin s) {k : Nat} (hk : k < s) :
    ∑ j, closedA s i j * closedC s j ^ k = closedC s i ^ (k + 1) / (k + 1) := by
  have hq := closedForm_checkOrder s hs
  have heta : lastTrue (condC (closedForm s)) maxCheckOrder = s := by
    simp only [checkOrder, Prod.mk.injEq] at hq; exact hq.2.1
  have hC := lastTrue_spec (condC (closedForm s)) maxCheckOrder (k + 1) (by omega)
    (by rw [heta]; omega)
  have hlen := closedForm_length s hs
  simp only [condC, List.all_eq_true, List.mem_range, beq_iff_eq] at hC
  have hi := congrArg Q6.toReal (hC i (by simp [closedForm, radauTable, hlen]))
  simp only [closedForm, radauTable, hlen, toReal_sum_map_range, Q6.toReal_mul, Q6.toReal_div,
    npow_toReal, Q6.toReal_natCast, Nat.add_sub_cancel] at hi
  simp only [closedA, closedC, Matrix.of_apply]
  rw [Fin.sum_univ_eq_sum_range (fun j => (get2 (closedFormA s) i j).toReal *
    (get (closedFormC s) j).toReal ^ k) s]
  push_cast at hi
  exact hi

/-- For s ≤ 3 the closed-form table is the collocation table on its nodes. -/
theorem closedForm_eq_collocA {s : Nat} (hs : s ∈ closedStages) :
    closedA s = collocA (closedC s) :=
  collocA_unique (closedC_injective hs) fun i _ hk => closedForm_C_real hs i hk

end Radau

/-! ## Radau node polynomial, s ≤ 9 -/

namespace Radau

open NodePoly

theorem fact_eq (n : Nat) : fact n = n.factorial := by
  induction n with
  | zero => rfl
  | succ n ih => simp [fact, ih, Nat.factorial_succ]

/-- A rational coefficient list as a real polynomial. -/
noncomputable def ratPoly (p : Poly ℚ) : ℝ[X] := (toPoly p).map (Rat.castHom ℝ)

theorem ratPoly_mul (p q : Poly ℚ) : ratPoly (Poly.mul p q) = ratPoly p * ratPoly q := by
  simp [ratPoly, Polynomial.map_mul]

theorem ratPoly_pow (p : Poly ℚ) (n : Nat) : ratPoly (Poly.pow p n) = ratPoly p ^ n := by
  simp [ratPoly, Polynomial.map_pow]

theorem ratPoly_X : ratPoly (Poly.X : Poly ℚ) = X := by simp [ratPoly]

theorem natDegree_toPoly_lt {R : Type} [CommRing R] {p : Poly R} (hp : p ≠ []) :
    (toPoly p).natDegree < p.length := by
  have : (toPoly p).natDegree ≤ p.length - 1 := by
    rw [natDegree_le_iff_coeff_eq_zero]
    intro N hN
    have : p.length ≤ N := by omega
    rw [← toPoly_coeff]
    simp [Poly.coeff, List.getD_eq_getElem?_getD, List.getElem?_eq_none this]
  have : 0 < p.length := List.length_pos_iff.mpr hp
  omega

/-- The certificate integral `integral01` is `∫_0^1` of the polynomial. -/
theorem integral01_eq (p : Poly ℚ) : ((integral01 p : ℚ) : ℝ) = (antideriv (ratPoly p)).eval 1 := by
  rcases eq_or_ne p [] with rfl | hp
  · simp [integral01, sum, ratPoly]
  have hdeg : (ratPoly p).natDegree < p.length :=
    (natDegree_map_le).trans_lt (natDegree_toPoly_lt hp)
  rw [eval_antideriv_eq_sum _ _ hdeg, integral01, sum_eq, Rat.cast_list_sum, List.map_map,
    ← sum_eq (R := ℝ), sum_map_range]
  refine Finset.sum_congr rfl fun k _ => ?_
  simp [ratPoly, coeff_map, ← toPoly_coeff, Poly.coeff, get]

variable {s : Nat} {c : Fin s → ℝ}

/-- If the nodes are the roots of the Radau node polynomial `M_s` and end at
    1, the weights satisfy B(2s−1): `Σ_j b_j c_j^k = 1/(k+1)` for `k ≤ 2s − 2`. -/
theorem radau_B (hs : s ∈ stages) (hc : Function.Injective c)
    (hM : nodePoly c = ratPoly (NodePoly.nodePoly s)) (last : Fin s) (hlast : c last = 1)
    {k : Nat} (hk : k < 2 * s - 1) :
    ∑ j, collocA c last j * c j ^ k = 1 / (k + 1) := by
  refine collocA_B hc last hlast (fun j hj => ?_) hk
  rw [hM, ← ratPoly_X, ← ratPoly_pow, ← ratPoly_mul, ← integral01_eq,
    nodePoly_orthogonal s hs j hj, Rat.cast_zero]

/-- The B(2s) residual of the Radau IIA table,
    `Σ_j b_j c_j^(2s−1) − 1/(2s) = (s!(s−1)!)²/((2s)!(2s−1)!)`: about 9.4e−11
    for s = 9, the value that the tolerance of `FIRKodeTable_CheckOrder` must
    resolve. -/
theorem radau_B2s_residual (hs : s ∈ stages) (hc : Function.Injective c)
    (hM : nodePoly c = ratPoly (NodePoly.nodePoly s)) (last : Fin s) (hlast : c last = 1) :
    ∑ j, collocA c last j * c j ^ (2 * s - 1) - 1 / (2 * s) =
      ((s.factorial * (s - 1).factorial) ^ 2 : ℝ) /
        ((2 * s).factorial * (2 * s - 1).factorial) := by
  have horth : ∀ j < s - 1, (antideriv (X ^ j * nodePoly c)).eval 1 = 0 := fun j hj => by
    rw [hM, ← ratPoly_X, ← ratPoly_pow, ← ratPoly_mul, ← integral01_eq,
      nodePoly_orthogonal s hs j hj, Rat.cast_zero]
  have hres : (antideriv (X ^ (s - 1) * nodePoly c)).eval 1 =
      -((s.factorial * (s - 1).factorial) ^ 2 : ℝ) /
        ((2 * s).factorial * (2 * s - 1).factorial) := by
    rw [hM, ← ratPoly_X, ← ratPoly_pow, ← ratPoly_mul, ← integral01_eq,
      nodePoly_B2s_residual s hs]
    simp [fact_eq]
  rw [collocA_B2s hc last hlast horth, hres]
  ring

theorem choose_eq (n k : Nat) : NodePoly.choose n k = n.choose k := by
  induction n generalizing k with
  | zero => cases k <;> simp [NodePoly.choose]
  | succ n ih => cases k <;> simp [NodePoly.choose, Nat.choose_succ_succ, ih]

/-- `firk_radau_poly` at `x = 2c − 1` is `(−1)^s (P̃_s + P̃_{s−1})`, with
    Mathlib's shifted Legendre polynomials `P̃_n(c) = P_n(1 − 2c)`: the
    recurrence in the C routine gives `P_s(2c − 1) − P_{s−1}(2c − 1)`, the
    definition of the Radau IIA nodes in `Mathematics.rst`. -/
theorem radauPoly_eq_shiftedLegendre (hs : s ∈ stages) :
    ratPoly (Poly.comp (Mirror.radauPoly s) [-1, 2]) =
      C ((-1 : ℝ) ^ s) * (shiftedLegendre s + shiftedLegendre (s - 1)).map (Int.castRingHom ℝ) := by
  rw [ratPoly, ← toPoly_trim, radauPoly_shiftedLegendre s hs]
  ext k
  rw [coeff_map, ← toPoly_coeff, coeff_C_mul, coeff_map, coeff_add, coeff_shiftedLegendre,
    coeff_shiftedLegendre]
  simp only [Poly.coeff, List.getD_eq_getElem?_getD, List.getElem?_map]
  by_cases hk : k < s + 1
  · rw [List.getElem?_range hk]
    simp [shLegendreCoeff, choose_eq]
  · rw [List.getElem?_eq_none (by simp; omega)]
    simp [Nat.choose_eq_zero_of_lt (by omega : s < k), Nat.choose_eq_zero_of_lt (by omega : s - 1 < k)]

theorem ratPoly_eval_ratCast (p : Poly ℚ) (q : ℚ) :
    (ratPoly p).eval (q : ℝ) = ((Poly.eval p q : ℚ) : ℝ) := by
  have hq : (q : ℝ) = Rat.castHom ℝ q := rfl
  rw [ratPoly, eval_map, hq, eval₂_hom, ← toPoly_eval]
  rfl

theorem ratPoly_deriv (p : Poly ℚ) : derivative (ratPoly p) = ratPoly (Poly.deriv p) := by
  rw [ratPoly, ratPoly, toPoly_deriv, derivative_map]

/-- `e_s = (−1)^s γ0/s` for the Radau nodes: with `e_j = −γ0 l_j(0)/c_j`
    (`errWeights_eq`), this is the check `e_s = (−1)^s γ0/s` of
    `firk_test_tables.c`. -/
theorem radau_e_last (hs : s ∈ stages) (hM : nodePoly c = ratPoly (NodePoly.nodePoly s))
    (last : Fin s) (hlast : c last = 1) (γ0 : ℝ) :
    -γ0 * ((ell c last).eval 0 / c last) = (-1) ^ s * γ0 / s := by
  have h := congrArg (fun q : ℚ => (q : ℝ)) (nodePoly_interior_product s hs)
  simp only [Rat.cast_div, Rat.cast_mul, Rat.cast_pow, Rat.cast_neg, Rat.cast_one,
    Rat.cast_natCast] at h
  set a := ((Poly.eval (NodePoly.nodePoly s) 0 : ℚ) : ℝ)
  set b := ((Poly.eval (Poly.deriv (NodePoly.nodePoly s)) 1 : ℚ) : ℝ)
  have e0 : (ratPoly (NodePoly.nodePoly s)).eval 0 = a := by
    simpa using ratPoly_eval_ratCast (NodePoly.nodePoly s) 0
  have e1 : (ratPoly (Poly.deriv (NodePoly.nodePoly s))).eval 1 = b := by
    simpa using ratPoly_eval_ratCast (Poly.deriv (NodePoly.nodePoly s)) 1
  have hab : a / b = (-1) ^ s / s := by
    have h2 : ((-1 : ℝ) ^ s) ^ 2 = 1 := by rw [← pow_mul, mul_comm, pow_mul]; simp
    calc a / b = ((-1) ^ s) ^ 2 * (a / b) := by rw [h2, one_mul]
      _ = (-1) ^ s * ((-1) ^ s * a / b) := by ring
      _ = (-1) ^ s * (1 / s) := by rw [h]
      _ = (-1) ^ s / s := by ring
  rw [ell_last_eval_zero last hlast, hlast, hM, ratPoly_deriv, e0, e1, div_one, neg_div,
    mul_neg, neg_mul, neg_neg, hab]
  ring

theorem ratPoly_derivN (p : Poly ℚ) (k : Nat) :
    ratPoly (derivN p k) = derivative^[k] (ratPoly p) := by
  induction k with
  | zero => rfl
  | succ k ih =>
    rw [derivN, ← ratPoly_deriv, ih, Function.iterate_succ_apply']

theorem eval_toPoly_eq_sum {R : Type} [CommRing R] (L : Poly R) (x : R) :
    (toPoly L).eval x = ∑ j ∈ Finset.range L.length, get L j * x ^ j := by
  induction L with
  | nil => simp
  | cons a L ih =>
    rw [toPoly_cons, eval_add, eval_C, eval_mul, eval_X, ih, List.length_cons,
      Finset.sum_range_succ', Finset.mul_sum]
    simp only [get, List.getD_cons_succ, List.getD_cons_zero, pow_zero, mul_one, pow_succ]
    rw [add_comm]
    exact congrArg (· + a) (Finset.sum_congr rfl fun _ _ => by ring)

theorem get_map_cast (L : List ℚ) (j : Nat) :
    get (L.map (Rat.castHom ℝ)) j = ((get L j : ℚ) : ℝ) := by
  simp only [get, List.getD_eq_getElem?_getD, List.getElem?_map]
  cases L[j]? <;> simp

theorem ratPoly_eval (L : Poly ℚ) (z : ℝ) :
    (ratPoly L).eval z = ∑ j ∈ Finset.range L.length, ((get L j : ℚ) : ℝ) * z ^ j := by
  rw [ratPoly, ← toPoly_map, eval_toPoly_eq_sum, List.length_map]
  exact Finset.sum_congr rfl fun j _ => by rw [get_map_cast]

theorem ratPoly_smul (a : ℚ) (L : Poly ℚ) : ratPoly (Poly.smul a L) = C (a : ℝ) * ratPoly L := by
  simp [ratPoly, Polynomial.map_mul]

theorem radau_norsett_term (hM : nodePoly c = ratPoly (NodePoly.nodePoly s)) (x : ℚ) (j : Nat) :
    ((Poly.eval (derivN (NodePoly.nodePoly s) (s - j)) x : ℚ) : ℝ) =
      (derivative^[s - j] (nodePoly c)).eval (x : ℝ) := by
  rw [hM, ← ratPoly_derivN, ratPoly_eval_ratCast]

/-- For the Radau nodes, the denominator of Nørsett's formula is
    `s! · padeDen s`. -/
theorem radau_norsettD (hs : s ∈ stages) (hM : nodePoly c = ratPoly (NodePoly.nodePoly s))
    (z : ℝ) : norsettD c z = (fact s : ℝ) * (ratPoly (padeDen s)).eval z := by
  have h1 : (ratPoly (norsettDen (NodePoly.nodePoly s) s)).eval z = norsettD c z := by
    rw [ratPoly_eval, norsettD, norsettDen, List.length_map, List.length_range]
    refine Finset.sum_congr rfl fun j hj => ?_
    rw [get_map_range _ _ (Finset.mem_range.mp hj), radau_norsett_term hM 0 j, Rat.cast_zero]
  rw [← h1, (nodePoly_norsett s hs).2, ratPoly_smul, eval_mul, eval_C, Rat.cast_natCast]

/-- For the Radau nodes, the numerator of Nørsett's formula is
    `s! · padeNum s`. -/
theorem radau_norsettN (hs : s ∈ stages) (hM : nodePoly c = ratPoly (NodePoly.nodePoly s))
    (z : ℝ) : norsettN c z = (fact s : ℝ) * (ratPoly (padeNum s)).eval z := by
  have h1 : (ratPoly (norsettNum (NodePoly.nodePoly s) s)).eval z = norsettN c z := by
    rw [ratPoly_eval, norsettN, norsettNum, List.length_map, List.length_range]
    refine Finset.sum_congr rfl fun j hj => ?_
    rw [get_map_range _ _ (Finset.mem_range.mp hj), radau_norsett_term hM 1 j, Rat.cast_one]
  rw [← h1, ratPoly, ← toPoly_trim, ← ratPoly, (nodePoly_norsett s hs).1, ratPoly_smul, eval_mul,
    eval_C, Rat.cast_natCast]

theorem fact_ne_zero_real (s : Nat) : (fact s : ℝ) ≠ 0 := by
  rw [fact_eq]; exact_mod_cast s.factorial_ne_zero

/-- The stability function of the Radau IIA table is the `(s−1, s)` Padé
    approximant of `e^z`: `R(z) Q(z) = P(z)` with `P = padeNum s` and
    `Q = padeDen s`, wherever `I − zA` is invertible. -/
theorem radau_pade (hs : s ∈ stages) (hc : Function.Injective c)
    (hM : nodePoly c = ratPoly (NodePoly.nodePoly s)) (last : Fin s) (hlast : c last = 1) (z : ℝ)
    (hz : IsUnit (1 - z • collocA c).det) :
    stabR (collocA c) ((collocA c).row last) z * (ratPoly (padeDen s)).eval z =
      (ratPoly (padeNum s)).eval z := by
  have hnorsett := norsett hc last hlast z hz
  rw [radau_norsettD hs hM, radau_norsettN hs hM] at hnorsett
  rw [mul_left_comm] at hnorsett
  exact mul_left_cancel₀ (fact_ne_zero_real s) hnorsett

end Radau

end Firkode
