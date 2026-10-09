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
import Mathlib.Tactic.LinearCombination
import Mathlib.Analysis.Polynomial.Basic
import Mathlib.LinearAlgebra.Matrix.Charpoly.Coeff

/-!
# The stability function of a collocation method

The stability function of a Butcher table is `R(z) = 1 + z bᵀ(I − zA)⁻¹𝟙`.
For collocation on nodes with `c_s = 1` and node polynomial `M`, Nørsett's
theorem (`norsett`) states `R(z) D(z) = N(z)` with

  `N(z) = Σ_j M^(s−j)(1) z^j`,  `D(z) = Σ_j M^(s−j)(0) z^j`.

The proof follows the collocation polynomial `u` of `y′ = zy`, `y(0) = 1`:
`u′ − zu` vanishes at the nodes, so it is a multiple of `M`, and the only
polynomial with `w′ = zw` for `z ≠ 0` is `w = 0`, which forces
`u = α Σ_j z^(s−j) M^(j)`.

Since `c_s = 1` the top coefficient `M(1)` of `N` vanishes while that of `D`,
`M(0)`, does not; so `R(z) → 0` as `z → −∞` (`stabR_tendsto_atBot`), the
stiff decay of L-stable methods.
-/

namespace Firkode.Colloc

open Polynomial Finset Matrix

variable {F : Type*} [Field F] {s : ℕ} {c : Fin s → F}

/-- The stability function `R(z) = 1 + z bᵀ(I − zA)⁻¹𝟙`. -/
noncomputable def stabR (A : Matrix (Fin s) (Fin s) F) (b : Fin s → F) (z : F) : F :=
  1 + z * (b ⬝ᵥ ((1 - z • A)⁻¹ *ᵥ fun _ => 1))

variable (c) in
/-- Numerator of Nørsett's formula, `N(z) = Σ_j M^(s−j)(1) z^j`. -/
noncomputable def norsettN (z : F) : F :=
  ∑ j ∈ range (s + 1), (derivative^[s - j] (nodePoly c)).eval 1 * z ^ j

variable (c) in
/-- Denominator of Nørsett's formula, `D(z) = Σ_j M^(s−j)(0) z^j`. -/
noncomputable def norsettD (z : F) : F :=
  ∑ j ∈ range (s + 1), (derivative^[s - j] (nodePoly c)).eval 0 * z ^ j

/-- `w′ = z w` with `z ≠ 0` has only the solution `w = 0` among polynomials. -/
theorem eq_zero_of_derivative_eq {z : F} (hz : z ≠ 0) {w : F[X]}
    (h : derivative w = C z * w) : w = 0 := by
  by_contra hw
  have h1 : (derivative w).coeff w.natDegree = 0 := by
    rw [coeff_derivative, coeff_eq_zero_of_natDegree_lt (Nat.lt_succ_self _), zero_mul]
  rw [h, coeff_C_mul] at h1
  exact mul_ne_zero hz (leadingCoeff_ne_zero.mpr hw) h1

variable (c) in
/-- `V(τ) = Σ_j z^(s−j) M^(j)(τ)`. -/
noncomputable def norsettV (z : F) : F[X] :=
  ∑ j ∈ range (s + 1), C (z ^ (s - j)) * derivative^[j] (nodePoly c)

/-- `V′ − zV = −z^(s+1) M`, by telescoping. -/
theorem derivative_norsettV (z : F) :
    derivative (norsettV c z) - C z * norsettV c z = -(C (z ^ (s + 1)) * nodePoly c) := by
  set f : ℕ → F[X] := fun j => C (z ^ (s + 1 - j)) * derivative^[j] (nodePoly c)
  have htel : ∑ j ∈ range (s + 1), (f (j + 1) - f j) = f (s + 1) - f 0 := sum_range_sub f (s + 1)
  have hlast : f (s + 1) = 0 := by
    simp only [f]
    rw [iterate_derivative_eq_zero (by rw [natDegree_nodePoly]; omega), mul_zero]
  have hterm : ∀ j ∈ range (s + 1), f (j + 1) - f j =
      C (z ^ (s - j)) * derivative^[j + 1] (nodePoly c) -
        C z * (C (z ^ (s - j)) * derivative^[j] (nodePoly c)) := by
    intro j hj
    have hj' : j ≤ s := Nat.lt_succ_iff.mp (mem_range.mp hj)
    simp only [f, show s + 1 - (j + 1) = s - j by omega, show s + 1 - j = s - j + 1 by omega,
      pow_succ, C_mul]
    ring
  rw [sum_congr rfl hterm, hlast, zero_sub] at htel
  simp only [f, Nat.sub_zero, Function.iterate_zero, id] at htel
  rw [norsettV, derivative_sum, mul_sum, ← sum_sub_distrib, ← htel]
  refine sum_congr rfl fun j _ => ?_
  rw [derivative_C_mul, Function.iterate_succ_apply']

theorem norsettV_eval_zero (z : F) : (norsettV c z).eval 0 = norsettD c z := by
  rw [norsettV, eval_finsetSum, norsettD, ← sum_range_reflect]
  refine sum_congr rfl fun j hj => ?_
  have hj' : j < s + 1 := mem_range.mp hj
  simp only [eval_mul, eval_C, show s + 1 - 1 - j = s - j by omega,
    show s - (s - j) = j by omega]
  ring

theorem norsettV_eval_one (z : F) : (norsettV c z).eval 1 = norsettN c z := by
  rw [norsettV, eval_finsetSum, norsettN, ← sum_range_reflect]
  refine sum_congr rfl fun j hj => ?_
  have hj' : j < s + 1 := mem_range.mp hj
  simp only [eval_mul, eval_C, show s + 1 - 1 - j = s - j by omega,
    show s - (s - j) = j by omega]
  ring

theorem natDegree_antideriv_ell_le (hc : Function.Injective c) (j : Fin s) :
    (antideriv (ell c j)).natDegree ≤ s := by
  have := natDegree_antideriv_le (ell c j)
  rw [natDegree_ell hc] at this
  have : 0 < s := Fin.pos j
  omega

/-! ## The collocation polynomial of `y′ = zy` -/

variable (c) in
/-- The collocation polynomial of `y′ = zy`, `y(0) = y0`, with stage values
    `Y`: `u(τ) = y0 + z Σ_j Y_j ∫_0^τ l_j`. -/
noncomputable def stagePoly (y0 z : F) (Y : Fin s → F) : F[X] :=
  C y0 + C z * ∑ j, C (Y j) * antideriv (ell c j)

theorem stagePoly_eval_zero (y0 z : F) (Y : Fin s → F) : (stagePoly c y0 z Y).eval 0 = y0 := by
  simp [stagePoly, eval_finsetSum, eval_antideriv_zero]

/-- Under the stage equations `Y_i = y0 + z Σ_j a_ij Y_j`, `u(c_i) = Y_i`. -/
theorem stagePoly_eval_node {y0 z : F} {Y : Fin s → F}
    (hY : ∀ i, Y i = y0 + z * ∑ j, collocA c i j * Y j) (i : Fin s) :
    (stagePoly c y0 z Y).eval (c i) = Y i := by
  rw [hY i]
  simp [stagePoly, eval_finsetSum, collocA, mul_comm]

theorem derivative_stagePoly [CharZero F] (y0 z : F) (Y : Fin s → F) :
    derivative (stagePoly c y0 z Y) = C z * ∑ j, C (Y j) * ell c j := by
  simp [stagePoly, derivative_antideriv]

/-- `u′(c_i) = z Y_i`: `u` satisfies the collocation equations. -/
theorem derivative_stagePoly_eval_node [CharZero F] (hc : Function.Injective c) (y0 z : F)
    (Y : Fin s → F) (i : Fin s) : (derivative (stagePoly c y0 z Y)).eval (c i) = z * Y i := by
  simp [derivative_stagePoly, eval_finsetSum, ell_eval_node hc]

theorem natDegree_stagePoly_le (hc : Function.Injective c) (y0 z : F) (Y : Fin s → F) :
    (stagePoly c y0 z Y).natDegree ≤ s := by
  refine (natDegree_add_le _ _).trans (max_le (by simp) ?_)
  refine (natDegree_C_mul_le _ _).trans (natDegree_sum_le_of_forall_le _ _ fun j _ => ?_)
  exact (natDegree_C_mul_le _ _).trans (natDegree_antideriv_ell_le hc j)

/-- `Y = (I − zA)⁻¹ y0 𝟙` solves the stage equations when `I − zA` is
    invertible. -/
theorem stage_eq_of_isUnit {A : Matrix (Fin s) (Fin s) F} {z : F} (hz : IsUnit (1 - z • A).det)
    (y0 : F) (i : Fin s) :
    ((1 - z • A)⁻¹ *ᵥ fun _ => y0) i = y0 + z * ∑ j, A i j * ((1 - z • A)⁻¹ *ᵥ fun _ => y0) j := by
  set Y : Fin s → F := (1 - z • A)⁻¹ *ᵥ fun _ => y0
  have h1 : (1 - z • A) *ᵥ Y = fun _ => y0 := by
    simp only [Y, mulVec_mulVec, mul_nonsing_inv _ hz, one_mulVec]
  rw [sub_mulVec, one_mulVec, smul_mulVec] at h1
  have := congrFun h1 i
  simp only [Pi.sub_apply, Pi.smul_apply, smul_eq_mul] at this
  rw [show ∑ j, A i j * Y j = (A *ᵥ Y) i from rfl]
  linear_combination this

/-- Nørsett's theorem: for collocation on distinct nodes with `c_s = 1` and
    weights `b = A[s]`, `R(z) D(z) = N(z)` wherever `I − zA` is invertible. -/
theorem norsett [CharZero F] (hc : Function.Injective c) (last : Fin s) (hlast : c last = 1)
    (z : F) (hz : IsUnit (1 - z • collocA c).det) :
    stabR (collocA c) ((collocA c).row last) z * norsettD c z = norsettN c z := by
  classical
  have hMdeg : (nodePoly c).natDegree = s := natDegree_nodePoly
  -- z = 0: both sides are the constant `M^(s)`
  rcases eq_or_ne z 0 with rfl | hz0
  · have hconst : (derivative^[s] (nodePoly c)).natDegree = 0 := by
      have := natDegree_iterate_derivative (nodePoly c) s
      omega
    have hC := eq_C_of_natDegree_eq_zero hconst
    simp only [stabR, norsettD, norsettN, zero_mul, add_zero, one_mul, sum_range_succ', pow_succ,
      mul_zero, sum_const_zero, zero_add, pow_zero, mul_one, Nat.sub_zero]
    rw [hC, eval_C, eval_C]
  set A := collocA c
  set Y : Fin s → F := (1 - z • A)⁻¹ *ᵥ fun _ => 1
  -- the stage equations `Y = 𝟙 + z A Y`
  have hY : ∀ i, Y i = 1 + z * ∑ j, A i j * Y j := stage_eq_of_isUnit hz 1
  -- the collocation polynomial `u`
  set u : F[X] := stagePoly c 1 z Y
  have hu0 : u.eval 0 = 1 := stagePoly_eval_zero 1 z Y
  have hunode : ∀ i, u.eval (c i) = Y i := stagePoly_eval_node hY
  have hu1 : u.eval 1 = stabR A (A.row last) z := by
    have := hunode last
    rw [hlast] at this
    rw [this, hY last]
    rfl
  have hdu : derivative u = C z * ∑ j, C (Y j) * ell c j := derivative_stagePoly 1 z Y
  have hudeg : u.natDegree ≤ s := natDegree_stagePoly_le hc 1 z Y
  -- `q = u′ − zu` vanishes at the nodes, so `q = K M`
  set q := derivative u - C z * u
  have hqnode : ∀ i, q.eval (c i) = 0 := by
    intro i
    simp [q, hdu, eval_finsetSum, ell_eval_node hc, hunode]
  have hqdeg : q.natDegree ≤ s :=
    (natDegree_sub_le _ _).trans (max_le ((natDegree_derivative_le u).trans (by omega))
      ((natDegree_C_mul_le _ _).trans hudeg))
  set K := q.coeff s
  have hqM : q = C K * nodePoly c := by
    refine Polynomial.eq_of_degree_sub_lt_of_eval_index_eq (s := univ) (v := c) hc.injOn ?_ ?_
    · have hlc : (nodePoly c).coeff s = 1 := by
        have := (nodePoly_monic (c := c)).coeff_natDegree
        rwa [hMdeg] at this
      rw [card_univ, Fintype.card_fin, degree_lt_iff_coeff_zero]
      intro m hm
      rw [coeff_sub, coeff_C_mul]
      rcases eq_or_lt_of_le hm with rfl | hm
      · rw [hlc, mul_one, sub_self]
      · rw [coeff_eq_zero_of_natDegree_lt (by omega), coeff_eq_zero_of_natDegree_lt (by omega),
          mul_zero, sub_zero]
    · intro i _
      simp [hqnode, nodePoly_eval_node]
  -- `w = u + (K/z^(s+1)) V` solves `w′ = zw`, hence vanishes
  set α := K / z ^ (s + 1)
  have hw : derivative (u + C α * norsettV c z) = C z * (u + C α * norsettV c z) := by
    have hV := derivative_norsettV (c := c) z
    have hz1 : z ^ (s + 1) ≠ 0 := pow_ne_zero _ hz0
    have hq : derivative u - C z * u = C K * nodePoly c := hqM
    rw [derivative_add, derivative_C_mul]
    have hαz : C α * C (z ^ (s + 1)) = C K := by
      rw [← C_mul]; congr 1
      show K / z ^ (s + 1) * z ^ (s + 1) = K
      exact div_mul_cancel₀ K hz1
    linear_combination hq + C α * hV - nodePoly c * hαz
  have hw0 := eq_zero_of_derivative_eq hz0 hw
  have hu : u = -(C α * norsettV c z) := eq_neg_of_add_eq_zero_left hw0
  have h0 : 1 = -(α * norsettD c z) := by
    rw [← hu0, hu, eval_neg, eval_mul, eval_C, norsettV_eval_zero]
  have h1 : stabR A (A.row last) z = -(α * norsettN c z) := by
    rw [← hu1, hu, eval_neg, eval_mul, eval_C, norsettV_eval_one]
  rw [h1]
  linear_combination (-(norsettN c z)) * h0

/-! ## Stiff decay -/

variable (c) in
/-- `N` as a polynomial in `z`. -/
noncomputable def norsettNPoly : F[X] :=
  ∑ j ∈ range (s + 1), C ((derivative^[s - j] (nodePoly c)).eval 1) * X ^ j

variable (c) in
/-- `D` as a polynomial in `z`. -/
noncomputable def norsettDPoly : F[X] :=
  ∑ j ∈ range (s + 1), C ((derivative^[s - j] (nodePoly c)).eval 0) * X ^ j

theorem eval_norsettNPoly (z : F) : (norsettNPoly c).eval z = norsettN c z := by
  simp [norsettNPoly, norsettN, eval_finsetSum]

theorem eval_norsettDPoly (z : F) : (norsettDPoly c).eval z = norsettD c z := by
  simp [norsettDPoly, norsettD, eval_finsetSum]

theorem coeff_norsettPoly (x : F) (m : ℕ) :
    (∑ j ∈ range (s + 1), C ((derivative^[s - j] (nodePoly c)).eval x) * X ^ j).coeff m =
      if m < s + 1 then (derivative^[s - m] (nodePoly c)).eval x else 0 := by
  simp only [finsetSum_coeff, coeff_C_mul_X_pow]
  split_ifs with hm
  · rw [sum_eq_single m (fun j _ hj => if_neg (Ne.symm hj)) (fun h => absurd (mem_range.mpr hm) h),
      if_pos rfl]
  · exact sum_eq_zero fun j hj => if_neg fun h => hm (by rw [h]; exact mem_range.mp hj)

/-- With `c_s = 1` and nonzero nodes, `deg N < deg D = s`. -/
theorem degree_norsettNPoly_lt (last : Fin s) (hlast : c last = 1) (h0 : ∀ i, c i ≠ 0) :
    (norsettNPoly c).degree < (norsettDPoly c).degree := by
  have hD : (norsettDPoly c).coeff s ≠ 0 := by
    rw [norsettDPoly, coeff_norsettPoly, if_pos (Nat.lt_succ_self s), Nat.sub_self,
      Function.iterate_zero, id, nodePoly, Lagrange.eval_nodal]
    exact prod_ne_zero_iff.mpr fun i _ => by rw [zero_sub, neg_ne_zero]; exact h0 i
  have hDdeg : (norsettDPoly c).degree = s := by
    refine le_antisymm ?_ (le_degree_of_ne_zero hD)
    rw [degree_le_iff_coeff_zero]
    intro m hm
    rw [norsettDPoly, coeff_norsettPoly, if_neg (by exact_mod_cast (not_lt.mpr (Nat.succ_le_of_lt
      (by exact_mod_cast hm))))]
  rw [hDdeg, degree_lt_iff_coeff_zero]
  intro m hm
  rw [norsettNPoly, coeff_norsettPoly]
  split_ifs with h
  · have hms : m = s := by omega
    subst hms
    rw [Nat.sub_self, Function.iterate_zero, id, ← hlast, nodePoly_eval_node]
  · rfl

/-- `det(I − zA)` as a polynomial in `z`. -/
theorem eval_charpolyRev (A : Matrix (Fin s) (Fin s) F) (z : F) :
    A.charpolyRev.eval z = (1 - z • A).det := by
  rw [Matrix.charpolyRev, ← coe_evalRingHom, RingHom.map_det]
  congr 1
  ext i j
  simp only [RingHom.mapMatrix_apply, Matrix.map_apply, Matrix.sub_apply, Matrix.one_apply,
    Matrix.smul_apply, smul_eq_mul, apply_ite, coe_evalRingHom, eval_sub, eval_mul, eval_X, eval_C,
    eval_one, eval_zero]

/-- Stiff decay: `R(z) → 0` as `z → −∞`, for collocation on distinct nonzero
    nodes with `c_s = 1`. -/
theorem stabR_tendsto_atBot {c : Fin s → ℝ} (hc : Function.Injective c) (h0 : ∀ i, c i ≠ 0)
    (last : Fin s) (hlast : c last = 1) :
    Filter.Tendsto (stabR (collocA c) ((collocA c).row last)) Filter.atBot (nhds 0) := by
  have hdeg := degree_norsettNPoly_lt last hlast h0
  have hDne : norsettDPoly c ≠ 0 := by
    intro h; rw [h, degree_zero] at hdeg; exact not_lt_bot hdeg
  have hdetne : (collocA c).charpolyRev ≠ 0 := by
    intro h
    have := Matrix.eval_charpolyRev (M := collocA c)
    rw [h, eval_zero] at this
    exact zero_ne_one this
  -- `N/D → 0` at `−∞`, through `z ↦ −z`
  have hlim : Filter.Tendsto (fun z => (norsettNPoly c).eval z / (norsettDPoly c).eval z)
      Filter.atBot (nhds 0) := by
    have hdeg' : ((norsettNPoly c).comp (-X)).degree < ((norsettDPoly c).comp (-X)).degree := by
      rwa [degree_comp_neg_X, degree_comp_neg_X]
    have := (div_tendsto_atTop_zero_of_degree_lt _ _ hdeg').comp Filter.tendsto_neg_atBot_atTop
    refine this.congr fun z => ?_
    simp [eval_comp]
  refine hlim.congr' ?_
  filter_upwards [eventually_atBot_not_isRoot _ hDne, eventually_atBot_not_isRoot _ hdetne]
    with z hD hdet
  have hz : IsUnit (1 - z • collocA c).det := by
    rw [← eval_charpolyRev]; exact isUnit_iff_ne_zero.mpr hdet
  have h := norsett hc last hlast z hz
  rw [eval_norsettNPoly, eval_norsettDPoly, ← h, mul_div_cancel_right₀]
  rwa [← eval_norsettDPoly]

end Firkode.Colloc
