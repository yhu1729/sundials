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

import FirkodeCert.Mirror
import FirkodeSpec.ListPoly
import FirkodeSpec.DenseOutput
import FirkodeSpec.ErrorEstimate

/-!
# The transcribed C routines compute the collocation quantities

`FirkodeCert.Mirror` transcribes the C table routines into exact arithmetic.
Here, for any list of distinct nodes (so for every stage count, including the
tables FIRKODE computes for s ≥ 4), the transcriptions are shown to compute
the objects of `Collocation`, `ErrorEstimate` and `DenseOutput`:

* `lagrangeCoeffs_eq` (`firk_lagrange_coeffs`): the Lagrange basis `l_j`;
* `collocationA_eq` (`firk_collocation_coeffs`): the coefficients
  `a_ij = ∫_0^{c_i} l_j`;
* `derive_P_eq` (`firk_derive_ld`): `P[m][j]` is coefficient `m + 1` of
  `L_j(θ) = θ l_j(θ)/c_j`;
* `denseWeights_eq` (`FIRKodeGetDky` for `k = 0`, `firkPredict`): the weights
  `L_j(θ) − d_j`;
* `dkyWeights_eq` (`FIRKodeGetDky` for `k ≥ 1`): the derivatives
  `L_j^(k)(θ)`.
-/

namespace Firkode

open Polynomial Finset Colloc Mirror

variable {F : Type} [Field F]

/-- The nodes of a list, indexed by `Fin`. -/
def nodesOf (cl : List F) : Fin cl.length → F := fun i => get cl i

/-- The Lagrange basis as a product over `Fin s` with the factor `j` replaced
    by 1. -/
theorem ell_eq_prod {s : Nat} (c : Fin s → F) (j : Fin s) :
    ell c j = ∏ k, if k = j then 1 else Lagrange.basisDivisor (c j) (c k) := by
  classical
  rw [ell, Lagrange.basis, ← mul_prod_erase univ _ (mem_univ j), if_pos rfl, one_mul]
  exact prod_congr rfl fun k hk => by rw [if_neg (ne_of_mem_erase hk)]

theorem toPoly_lagrange_foldl (cl : List F) (j : Nat) (n : Nat) :
    toPoly ((List.range n).foldl
        (fun (r : Poly F) (k : Nat) =>
          if k = j then r
          else Poly.smul (1 / (get cl j - get cl k)) (Poly.sub (Poly.mulX r) (Poly.smul (get cl k) r)))
        [1]) =
      ∏ k ∈ range n, if k = j then 1 else Lagrange.basisDivisor (get cl j) (get cl k) := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [List.range_succ, List.foldl_append, List.foldl_cons, List.foldl_nil, prod_range_succ, ← ih]
    split_ifs with h
    · simp
    · simp only [toPoly_smul, toPoly_sub, toPoly_mulX, Lagrange.basisDivisor, one_div]
      ring

/-- `firk_lagrange_coeffs` computes the Lagrange basis polynomial `l_j`. -/
theorem lagrangeCoeffs_eq (cl : List F) (j : Fin cl.length) :
    toPoly (lagrangeCoeffs cl j) = ell (nodesOf cl) j := by
  rw [lagrangeCoeffs, toPoly_lagrange_foldl, ell_eq_prod,
    ← Fin.prod_univ_eq_prod_range (fun k => if k = (j : Nat) then (1 : F[X])
      else Lagrange.basisDivisor (get cl j) (get cl k)) cl.length]
  refine prod_congr rfl fun k _ => ?_
  simp [nodesOf, Fin.ext_iff]

/-- `firk_collocation_coeffs` computes `a_ij = ∫_0^{c_i} l_j`. -/
theorem collocationA_eq (cl : List F) (hc : Function.Injective (nodesOf cl))
    (i j : Fin cl.length) :
    get2 (collocationA cl) i j = collocA (nodesOf cl) i j := by
  rw [collocationA, get2_map_range _ _ _ i.isLt j.isLt, sum_map_range, collocA, Matrix.of_apply,
    eval_antideriv_eq_sum _ _ (N := cl.length) (by
      rw [natDegree_ell hc]; exact Nat.sub_lt (Fin.pos j) one_pos)]
  refine sum_congr rfl fun m _ => ?_
  rw [← lagrangeCoeffs_eq, ← toPoly_coeff, npow_eq]
  simp only [nodesOf, Nat.cast_add, Nat.cast_one]

/-- `firk_derive_ld` computes `P[m][j]`, coefficient `m + 1` of `L_j`. -/
theorem derive_P_eq (T : Table F) (γ0 : F) (Ainv : List (List F)) {m : Nat}
    (hm : m < T.c.length) (j : Fin T.c.length) :
    get2 (derive T γ0 Ainv).P m j = (denseL (nodesOf T.c) j).coeff (m + 1) := by
  rw [derive, get2_map_range _ _ _ hm j.isLt, denseL_coeff_succ, ← lagrangeCoeffs_eq,
    ← toPoly_coeff]
  rfl

/-- `firk_derive_ld` on the computed table of `firk_radau_build` gives
    `d = bᵀA⁻¹ = e_s`, when `Ainv` holds the exact inverse. -/
theorem derive_d_eq [CharZero F] (cl : List F) (hc : Function.Injective (nodesOf cl))
    (h0 : ∀ i, nodesOf cl i ≠ 0) (γ0 : F) (Ainv : List (List F))
    (hinv : ∀ i j : Fin cl.length, get2 Ainv i j = (collocA (nodesOf cl))⁻¹ i j)
    (last : Fin cl.length) (hlast : (last : Nat) = cl.length - 1) (j : Fin cl.length) :
    get (derive (radauTable cl (collocationA cl)) γ0 Ainv).d j = if last = j then 1 else 0 := by
  have hA : IsUnit (collocA (nodesOf cl)).det :=
    isUnit_iff_ne_zero.mpr (collocA_det_ne_zero hc h0)
  have h2 := congrFun (row_mul_inv hA last) j
  simp only [Matrix.vecMul, dotProduct, Matrix.row] at h2
  simp only [derive, radauTable]
  rw [get_map_range _ _ j.isLt, sum_map_range,
    ← Fin.sum_univ_eq_sum_range
      (fun i => get ((collocationA cl).getD (cl.length - 1) []) i * get2 Ainv i j) cl.length]
  have hterm : ∀ i : Fin cl.length,
      get ((collocationA cl).getD (cl.length - 1) []) i * get2 Ainv i j =
        collocA (nodesOf cl) last i * (collocA (nodesOf cl))⁻¹ i j := by
    intro i
    rw [hinv, ← collocationA_eq cl hc, hlast]
    rfl
  rw [sum_congr rfl (fun i _ => hterm i), h2, Pi.single_apply]
  simp [eq_comm]

/-- `firk_derive_ld` on the computed table of `firk_radau_build` gives the
    error-estimate weights `e_j = −γ0 l_j(0)/c_j`, when `Ainv` holds the exact
    inverse. -/
theorem derive_e_eq [CharZero F] (cl : List F) (hc : Function.Injective (nodesOf cl))
    (h0 : ∀ i, nodesOf cl i ≠ 0) (γ0 : F) (Ainv : List (List F))
    (hinv : ∀ i j : Fin cl.length, get2 Ainv i j = (collocA (nodesOf cl))⁻¹ i j)
    (j : Fin cl.length) :
    get (derive (radauTable cl (collocationA cl)) γ0 Ainv).e j =
      -γ0 * ((ell (nodesOf cl) j).eval 0 / nodesOf cl j) := by
  have h := congrFun (errWeights_eq hc h0 γ0) j
  simp only [Pi.smul_apply, smul_eq_mul, Matrix.vecMul, dotProduct] at h
  rw [← h]
  simp only [derive, radauTable]
  set l0 := (List.range cl.length).map fun i => Poly.coeff (lagrangeCoeffs cl i) 0
  rw [get_map_range _ _ j.isLt, sum_map_range,
    ← Fin.sum_univ_eq_sum_range (fun i => get l0 i * get2 Ainv i j) cl.length, mul_comm]
  congr 1
  refine sum_congr rfl fun i _ => ?_
  rw [get_map_range _ _ i.isLt, hinv, toPoly_coeff, lagrangeCoeffs_eq, coeff_zero_eq_eval_zero]

theorem foldl_add_map (l : List Nat) (a : F) (f : Nat → F) :
    l.foldl (fun acc m => acc + f m) a = a + (l.map f).sum := by
  induction l generalizing a with
  | nil => simp
  | cons x l ih => simp [ih, add_assoc]

variable {s : Nat} {c : Fin s → F}

/-- The weights of `FIRKodeGetDky` (`k = 0`) and `firkPredict`, `L_j(θ) − d_j`,
    for any `P` holding the coefficients of the `L_j`. -/
theorem denseWeights_eq (hc : Function.Injective c) (D : Derived F)
    (hP : ∀ m < s, ∀ j : Fin s, get2 D.P m j = (denseL c j).coeff (m + 1)) (θ : F) (j : Fin s) :
    get (denseWeights D s θ) j = (denseL c j).eval θ - get D.d j := by
  rw [denseWeights, get_map_range _ _ j.isLt, foldl_add_map, ← sum_eq, sum_map_range,
    eval_eq_sum_range' (n := s + 1) (Nat.lt_succ_of_le (natDegree_denseL_le hc j)),
    sum_range_succ', denseL_coeff_zero]
  have : ∀ m ∈ range s, npow θ (m + 1) * get2 D.P m j = (denseL c j).coeff (m + 1) * θ ^ (m + 1) :=
    fun m hm => by rw [hP m (mem_range.mp hm), npow_eq, mul_comm]
  rw [sum_congr rfl this]
  ring

theorem foldl_descFactorial (m k : Nat) :
    (List.range k).foldl (fun f t => f * ((m + 1 - t : Nat) : F)) 1 = ((m + 1).descFactorial k : F) := by
  induction k with
  | zero => simp
  | succ k ih =>
    rw [List.range_succ, List.foldl_append, List.foldl_cons, List.foldl_nil, ih,
      Nat.descFactorial_succ]
    push_cast
    ring

theorem sum_range'_eq_sum_Ico (a n : Nat) (f : Nat → F) :
    ((List.range' a n).map f).sum = ∑ m ∈ Ico a (a + n), f m := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [List.range'_concat, List.map_append, List.sum_append, ih, ← Nat.add_assoc,
      sum_Ico_succ_top (by omega)]
    simp

/-- The weights of `FIRKodeGetDky` for `1 ≤ k ≤ s`: the `k`-th derivatives
    `L_j^(k)(θ)`, for any `P` holding the coefficients of the `L_j`. -/
theorem dkyWeights_eq [CharZero F] (hc : Function.Injective c) (D : Derived F)
    (hP : ∀ m < s, ∀ j : Fin s, get2 D.P m j = (denseL c j).coeff (m + 1)) {k : Nat}
    (hk : 1 ≤ k) (hks : k ≤ s) (θ : F) (j : Fin s) :
    get (dkyWeights D s k θ) j = (derivative^[k] (denseL c j)).eval θ := by
  rw [dkyWeights, get_map_range _ _ j.isLt, foldl_add_map, zero_add, sum_range'_eq_sum_Ico,
    iterate_derivative_denseL_eval hc j hk hks, show k - 1 + (s - (k - 1)) = s by omega]
  refine sum_congr rfl fun m hm => ?_
  rw [foldl_descFactorial, npow_eq, hP m (mem_Ico.mp hm).2 j, denseL_coeff_succ]

end Firkode
