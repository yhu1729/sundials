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
# Dense output for a general number of stages

FIRKODE's dense output is `u(t_n + θh) = y_n + Σ_j L_j(θ) Z_j` with
`L_j(θ) = θ l_j(θ)/c_j`, the Lagrange basis on the nodes `{0, c_1, …, c_s}`;
`firk_derive_ld` stores `P[m][j] = l_j`'s coefficient `m` over `c_j`, which is
coefficient `m + 1` of `L_j`. For distinct nonzero nodes:

* `denseL_eval_zero`, `denseL_eval_node`: `L_j(0) = 0` and `L_j(c_i) = δ_ij`, so
  the dense output passes through `y_n` and the stage values (`dense_interp`);
* `denseL_eval_one`: `L_j(1) = δ_js` when `c_s = 1`, so it ends at
  `y_{n+1} = y_n + Z_s`;
* `antideriv_eq_sum_denseL` and `denseL_deriv_collocA`: `∫_0^θ l_k =
  Σ_j a_jk L_j(θ)`, hence `L′(c) A = I`: the dense output satisfies the
  collocation equations `h u′(t_n + c_i h) = Σ_j L_j′(c_i) Z_j = h F_i`;
* `iterate_derivative_denseL_eval`: the `k`-th derivative formula of
  `FIRKodeGetDky`.
-/

namespace Firkode.Colloc

open Polynomial Finset

variable {F : Type*} [Field F] {s : ℕ} {c : Fin s → F}

variable (c) in
/-- `L_j(θ) = θ l_j(θ)/c_j`. -/
noncomputable def denseL (j : Fin s) : F[X] := C (c j)⁻¹ * X * ell c j

theorem denseL_eval_zero (j : Fin s) : (denseL c j).eval 0 = 0 := by
  simp [denseL]

theorem denseL_eval_node (hc : Function.Injective c) (h0 : ∀ i, c i ≠ 0) (i j : Fin s) :
    (denseL c j).eval (c i) = if i = j then 1 else 0 := by
  simp only [denseL, eval_mul, eval_C, eval_X, ell_eval_node hc]
  split_ifs with h
  · subst h; field_simp [h0 i]
  · simp

/-- The dense-output coefficients of `firk_derive_ld`:
    coefficient `m + 1` of `L_j` is coefficient `m` of `l_j` over `c_j`. -/
theorem denseL_coeff_succ (j : Fin s) (m : ℕ) :
    (denseL c j).coeff (m + 1) = (ell c j).coeff m / c j := by
  simp [denseL, mul_assoc, coeff_C_mul, coeff_X_mul, div_eq_inv_mul]

theorem denseL_coeff_zero (j : Fin s) : (denseL c j).coeff 0 = 0 := by
  simp [denseL, mul_assoc]

theorem natDegree_denseL_le (hc : Function.Injective c) (j : Fin s) :
    (denseL c j).natDegree ≤ s := by
  have hs : 0 < s := Fin.pos j
  calc (denseL c j).natDegree ≤ (C (c j)⁻¹ * X).natDegree + (ell c j).natDegree :=
        natDegree_mul_le
    _ ≤ 1 + (s - 1) := by
        gcongr
        · exact (natDegree_C_mul_le _ _).trans natDegree_X_le
        · exact (natDegree_ell hc j).le
    _ = s := by omega

/-- The dense output `u(θ) = y_n + Σ_j L_j(θ) Z_j` takes the value `y_n` at
    `θ = 0` and the stage value `y_n + Z_i` at `θ = c_i`. -/
theorem dense_interp (hc : Function.Injective c) (h0 : ∀ i, c i ≠ 0)
    {V : Type*} [AddCommGroup V] [Module F V] (yn : V) (Z : Fin s → V) :
    yn + ∑ j, (denseL c j).eval 0 • Z j = yn ∧
      ∀ i, yn + ∑ j, (denseL c j).eval (c i) • Z j = yn + Z i := by
  refine ⟨by simp [denseL_eval_zero], fun i => ?_⟩
  simp [denseL_eval_node hc h0]

/-- `L_j(1) = δ_js` when the last node is `c_s = 1`: the dense output ends at
    `y_n + Z_s = y_{n+1}`. -/
theorem denseL_eval_one (hc : Function.Injective c) (h0 : ∀ i, c i ≠ 0) (last : Fin s)
    (hlast : c last = 1) (j : Fin s) : (denseL c j).eval 1 = if last = j then 1 else 0 := by
  have := denseL_eval_node hc h0 last j
  rwa [hlast] at this

/-- `∫_0^θ l_k = Σ_j a_jk L_j(θ)`: both sides have degree `≤ s` and agree at
    `0` and at the nodes. -/
theorem antideriv_eq_sum_denseL (hc : Function.Injective c) (h0 : ∀ i, c i ≠ 0) (k : Fin s) :
    antideriv (ell c k) = ∑ j, C (collocA c j k) * denseL c j := by
  classical
  have hs : 0 < s := Fin.pos k
  set pts : Finset F := insert 0 (univ.image c)
  have h0mem : (0 : F) ∉ univ.image c := by
    simp only [mem_image, mem_univ, true_and, not_exists]
    exact fun i h => h0 i h
  have hcard : pts.card = s + 1 := by
    rw [card_insert_of_notMem h0mem, card_image_of_injective _ hc, card_univ, Fintype.card_fin]
  refine Polynomial.eq_of_degree_sub_lt_of_eval_finset_eq pts ?_ ?_
  · rw [hcard]
    refine (degree_sub_le _ _).trans_lt (max_lt ?_ ?_)
    · calc (antideriv (ell c k)).degree ≤ (antideriv (ell c k)).natDegree := degree_le_natDegree
        _ ≤ ((s - 1 + 1 : ℕ) : WithBot ℕ) := by
            exact_mod_cast (natDegree_antideriv_le _).trans (by rw [natDegree_ell hc])
        _ < (s + 1 : ℕ) := by exact_mod_cast (by omega : s - 1 + 1 < s + 1)
    · refine (degree_sum_le _ _).trans_lt ?_
      rw [Finset.sup_lt_iff (WithBot.bot_lt_coe _)]
      intro j _
      calc (C (collocA c j k) * denseL c j).degree
          ≤ (C (collocA c j k) * denseL c j).natDegree := degree_le_natDegree
        _ ≤ (s : WithBot ℕ) := by
            exact_mod_cast (natDegree_C_mul_le _ _).trans (natDegree_denseL_le hc j)
        _ < (s + 1 : ℕ) := by exact_mod_cast Nat.lt_succ_self s
  · intro x hx
    rcases mem_insert.mp hx with rfl | hx
    · simp [eval_antideriv_zero, eval_finsetSum, denseL_eval_zero]
    · obtain ⟨i, -, rfl⟩ := mem_image.mp hx
      simp [eval_finsetSum, denseL_eval_node hc h0, collocA]

/-- `L′(c) A = I`: `Σ_j L_j′(c_i) a_jk = δ_ik`. -/
theorem denseL_deriv_collocA [CharZero F] (hc : Function.Injective c) (h0 : ∀ i, c i ≠ 0) (i k : Fin s) :
    ∑ j, (derivative (denseL c j)).eval (c i) * collocA c j k = if i = k then 1 else 0 := by
  have h := congrArg (fun p => (derivative p).eval (c i)) (antideriv_eq_sum_denseL hc h0 k)
  simp only [derivative_antideriv, derivative_sum, derivative_C_mul, eval_finsetSum, eval_mul,
    eval_C] at h
  rw [← ell_eval_node hc i k, h]
  exact sum_congr rfl fun j _ => mul_comm _ _

/-- The `k`-th derivative formula of `FIRKodeGetDky`, for `1 ≤ k ≤ s`:
    `L_j^(k)(θ) = Σ_{m=k−1}^{s−1} (m+1)!/(m+1−k)! θ^(m+1−k) P[m][j]`, with
    `P[m][j] = (coefficient m of l_j)/c_j`. -/
theorem iterate_derivative_denseL_eval (hc : Function.Injective c) (j : Fin s) {k : ℕ}
    (hk : 1 ≤ k) (hks : k ≤ s) (θ : F) :
    (derivative^[k] (denseL c j)).eval θ =
      ∑ m ∈ Ico (k - 1) s,
        ((m + 1).descFactorial k : F) * θ ^ (m + 1 - k) * ((ell c j).coeff m / c j) := by
  have hdeg : (derivative^[k] (denseL c j)).natDegree < s - k + 1 := by
    have := natDegree_iterate_derivative (denseL c j) k
    have := natDegree_denseL_le hc j
    omega
  rw [eval_eq_sum_range' hdeg, sum_Ico_eq_sum_range, show s - (k - 1) = s - k + 1 by omega]
  refine sum_congr rfl fun x _ => ?_
  have e1 : k - 1 + x + 1 = x + k := by omega
  have e2 : k - 1 + x + 1 - k = x := by omega
  rw [coeff_iterate_derivative, nsmul_eq_mul, e2, ← denseL_coeff_succ, e1]
  ring

end Firkode.Colloc
