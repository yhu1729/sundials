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

import FirkodeSpec.AStability
import FirkodeSpec.Mirror

/-!
# `FIRKodeTable_CheckOrder` on the Radau IIA tables (s ≤ 9)

`Mirror.checkOrder` transcribes `FIRKodeTable_CheckOrder` with its tolerance
replaced by exact equality. On the table that `firk_radau_build` computes
from real nodes (`radauTable cl (collocationA cl)`) it finds exactly the
simplifying assumptions of a collocation method with B(2s−1) and a nonzero
B(2s) error (`checkOrder_colloc`):

* B(k) holds for `k ≤ 2s − 1` (`collocA_B`) and fails for `k = 2s`;
* C(k) holds for `k ≤ s` (`collocA_C`) and fails for `k = s + 1`
  (`collocA_not_C_succ`);
* D(k) holds for `k ≤ s − 1` (`collocA_D`) and fails for `k = s`
  (`collocA_not_D`).

So for the Radau IIA tables of every stage count FIRKODE supports the routine
returns `(qB, η, ζ, q) = (2s − 1, s, s − 1, 2s − 1)` in exact arithmetic
(`radau_checkOrder`). That these conditions give order `2s − 1` is Butcher's
theorem (cited).
-/

namespace Firkode.CheckOrder

open Polynomial Finset Mirror Generated Colloc

/-- `lastTrue cond kmax = m` when `cond 1, …, cond m` hold and `cond (m + 1)`
    fails, for `m < kmax`. -/
theorem lastTrue_go_eq (cond : ℕ → Bool) (m : ℕ) (hfalse : cond (m + 1) = false) :
    ∀ fuel k, k ≤ m + 1 → m + 1 < k + fuel → (∀ j, k ≤ j → j ≤ m → cond j = true) →
      lastTrue.go cond fuel k = m := by
  intro fuel
  induction fuel with
  | zero => intro k hk hf; omega
  | succ fuel ih =>
    intro k hk hf htrue
    simp only [lastTrue.go]
    rcases eq_or_lt_of_le hk with rfl | hlt
    · rw [hfalse]; simp
    · rw [htrue k le_rfl (by omega), if_pos rfl]
      exact ih (k + 1) (by omega) (by omega) fun j h1 h2 => htrue j (by omega) h2

theorem lastTrue_eq (cond : ℕ → Bool) (kmax m : ℕ) (hm : m + 1 ≤ kmax)
    (htrue : ∀ j, 1 ≤ j → j ≤ m → cond j = true) (hfalse : cond (m + 1) = false) :
    lastTrue cond kmax = m :=
  lastTrue_go_eq cond m hfalse kmax 1 (by omega) (by omega) htrue

theorem sum_range_eq_fin (n : ℕ) (f : ℕ → ℝ) :
    sum ((List.range n).map f) = ∑ i : Fin n, f i := by
  rw [sum_map_range, Fin.sum_univ_eq_sum_range]

section Table

variable (cl : List ℝ) (last : Fin cl.length) (hlastv : (last : ℕ) = cl.length - 1)
include hlastv

theorem get_b (hc : Function.Injective (nodesOf cl)) (i : Fin cl.length) :
    get (radauTable cl (collocationA cl)).b i = collocA (nodesOf cl) last i := by
  simp only [radauTable]
  rw [← collocationA_eq cl hc last i, hlastv]
  rfl

end Table

/-- The conditions of `FIRKodeTable_CheckOrder` in terms of the collocation
    coefficients. -/
theorem condB_iff (cl : List ℝ) (hc : Function.Injective (nodesOf cl)) (last : Fin cl.length)
    (hlastv : (last : ℕ) = cl.length - 1) (k : ℕ) :
    condB (radauTable cl (collocationA cl)) k = true ↔
      ∑ i, collocA (nodesOf cl) last i * nodesOf cl i ^ (k - 1) = 1 / (k : ℝ) := by
  simp only [condB, beq_iff_eq]
  rw [show (radauTable cl (collocationA cl)).c = cl from rfl, sum_range_eq_fin]
  simp only [get_b cl last hlastv hc, npow_eq]
  rfl

theorem condC_iff (cl : List ℝ) (hc : Function.Injective (nodesOf cl)) (k : ℕ) :
    condC (radauTable cl (collocationA cl)) k = true ↔
      ∀ i, ∑ j, collocA (nodesOf cl) i j * nodesOf cl j ^ (k - 1) = nodesOf cl i ^ k / (k : ℝ) := by
  have hrow : ∀ i : Fin cl.length,
      sum ((List.range cl.length).map fun j => get2 (collocationA cl) i j * npow (get cl j) (k - 1)) =
        ∑ j, collocA (nodesOf cl) i j * nodesOf cl j ^ (k - 1) := fun i => by
    rw [sum_range_eq_fin]
    exact Finset.sum_congr rfl fun j _ => by rw [collocationA_eq cl hc i j, npow_eq]; rfl
  simp only [condC, List.all_eq_true, List.mem_range, beq_iff_eq]
  change (∀ i, i < cl.length → sum ((List.range cl.length).map fun j =>
      get2 (collocationA cl) i j * npow (get cl j) (k - 1)) = npow (get cl i) k / (k : ℝ)) ↔ _
  constructor
  · intro h i
    have := h i i.isLt
    rw [hrow i, npow_eq] at this
    exact this
  · intro h i hi
    have := hrow ⟨i, hi⟩
    simp only at this
    rw [this, npow_eq]
    exact h ⟨i, hi⟩

theorem condD_iff (cl : List ℝ) (hc : Function.Injective (nodesOf cl)) (last : Fin cl.length)
    (hlastv : (last : ℕ) = cl.length - 1) (k : ℕ) :
    condD (radauTable cl (collocationA cl)) k = true ↔
      ∀ j, ∑ i, collocA (nodesOf cl) last i * nodesOf cl i ^ (k - 1) * collocA (nodesOf cl) i j =
        collocA (nodesOf cl) last j * (1 - nodesOf cl j ^ k) / (k : ℝ) := by
  have hcol : ∀ j : Fin cl.length,
      sum ((List.range cl.length).map fun i => get (radauTable cl (collocationA cl)).b i *
          npow (get cl i) (k - 1) * get2 (collocationA cl) i j) =
        ∑ i, collocA (nodesOf cl) last i * nodesOf cl i ^ (k - 1) * collocA (nodesOf cl) i j :=
    fun j => by
      rw [sum_range_eq_fin]
      exact Finset.sum_congr rfl fun i _ => by
        rw [get_b cl last hlastv hc i, collocationA_eq cl hc i j, npow_eq]; rfl
  simp only [condD, List.all_eq_true, List.mem_range, beq_iff_eq]
  change (∀ j, j < cl.length → sum ((List.range cl.length).map fun i =>
      get (radauTable cl (collocationA cl)).b i * npow (get cl i) (k - 1) *
        get2 (collocationA cl) i j) =
      get (radauTable cl (collocationA cl)).b j * (1 - npow (get cl j) k) / (k : ℝ)) ↔ _
  constructor
  · intro h j
    have := h j j.isLt
    rw [hcol j, get_b cl last hlastv hc j, npow_eq] at this
    exact this
  · intro h j hj
    have := hcol ⟨j, hj⟩
    simp only at this
    rw [this, get_b cl last hlastv hc ⟨j, hj⟩, npow_eq]
    exact h ⟨j, hj⟩

/-- `FIRKodeTable_CheckOrder`, in exact arithmetic, on the table built from
    distinct nonzero nodes with `c_s = 1`, B(2s−1) and a nonzero B(2s) error,
    returns `(2s − 1, s, s − 1, 2s − 1)`, provided `2s ≤ FIRK_MAX_CHECK_ORDER`
    so that the loops reach the failing conditions. -/
theorem checkOrder_colloc (cl : List ℝ) (hc : Function.Injective (nodesOf cl))
    (h0 : ∀ i, nodesOf cl i ≠ 0) (last : Fin cl.length) (hlastv : (last : ℕ) = cl.length - 1)
    (hlast : nodesOf cl last = 1)
    (horth : ∀ j < cl.length - 1, (antideriv (X ^ j * nodePoly (nodesOf cl))).eval 1 = 0)
    (hres : (antideriv (X ^ (cl.length - 1) * nodePoly (nodesOf cl))).eval 1 ≠ 0)
    (hlen : 2 * cl.length ≤ maxCheckOrder) :
    checkOrder (radauTable cl (collocationA cl)) maxCheckOrder =
      (2 * cl.length - 1, cl.length, cl.length - 1, 2 * cl.length - 1) := by
  have hL : 1 ≤ cl.length := by have := last.isLt; omega
  have hqB : lastTrue (condB (radauTable cl (collocationA cl))) maxCheckOrder =
      2 * cl.length - 1 := by
    refine lastTrue_eq _ _ _ (by omega) (fun k hk1 hk => ?_) ?_
    · rw [condB_iff cl hc last hlastv, collocA_B hc last hlast horth (k := k - 1) (by omega)]
      push_cast [Nat.cast_sub hk1]
      ring
    · rw [Bool.eq_false_iff]
      intro h
      rw [condB_iff cl hc last hlastv, show 2 * cl.length - 1 + 1 - 1 = 2 * cl.length - 1 by omega,
        collocA_B2s hc last hlast horth, show 2 * cl.length - 1 + 1 = 2 * cl.length by omega] at h
      apply hres
      push_cast at h
      linear_combination -h
  have hη : lastTrue (condC (radauTable cl (collocationA cl))) maxCheckOrder = cl.length := by
    refine lastTrue_eq _ _ _ (by omega) (fun k hk1 hk => ?_) ?_
    · rw [condC_iff cl hc]
      intro i
      have := collocA_C hc i (k := k - 1) (by omega)
      rw [show k - 1 + 1 = k by omega] at this
      rw [this]
      push_cast [Nat.cast_sub hk1]
      ring
    · rw [Bool.eq_false_iff]
      intro h
      rw [condC_iff cl hc, show cl.length + 1 - 1 = cl.length by omega] at h
      apply collocA_not_C_succ hc h0 hL
      intro i
      have := h i
      push_cast at this
      exact this
  have hζ : lastTrue (condD (radauTable cl (collocationA cl))) maxCheckOrder = cl.length - 1 := by
    refine lastTrue_eq _ _ _ (by omega) (fun k hk1 hk => ?_) ?_
    · rw [condD_iff cl hc last hlastv]
      exact collocA_D hc last hlast horth hk1 hk
    · rw [Bool.eq_false_iff]
      intro h
      rw [condD_iff cl hc last hlastv, show cl.length - 1 + 1 = cl.length by omega] at h
      exact collocA_not_D hc last hlast horth hres h
  simp only [checkOrder, hqB, hη, hζ, Prod.mk.injEq, true_and]
  omega

theorem nodePoly_comp_equiv {n m : ℕ} (e : Fin n ≃ Fin m) (f : Fin m → ℝ) :
    nodePoly (f ∘ e) = nodePoly f := by
  simp only [nodePoly, Lagrange.nodal]
  exact Fintype.prod_equiv e _ _ fun i => rfl

/-- For the Radau IIA tables of every stage count FIRKODE supports,
    `FIRKodeTable_CheckOrder` returns `(2s − 1, s, s − 1, 2s − 1)` in exact
    arithmetic: B(2s−1), C(s) and D(s−1) hold, and B(2s), C(s+1) and D(s)
    fail. -/
theorem radau_checkOrder {s : ℕ} (hs : s ∈ NodePoly.stages) :
    checkOrder (radauTable (List.ofFn (Nodes.radauNodes s))
        (collocationA (List.ofFn (Nodes.radauNodes s)))) maxCheckOrder =
      (2 * s - 1, s, s - 1, 2 * s - 1) := by
  have hs9 : 1 ≤ s ∧ s ≤ 9 := by simp [NodePoly.stages, maxStages] at hs; omega
  set f := Nodes.radauNodes s
  have hlen : (List.ofFn f).length = s := List.length_ofFn
  have hnodes : nodesOf (List.ofFn f) = f ∘ finCongr hlen := by
    funext i
    have hi : (i : ℕ) < s := by have := i.isLt; omega
    simp only [nodesOf, get, List.getD_eq_getElem?_getD, List.getElem?_ofFn, Function.comp,
      finCongr_apply, dif_pos hi, Option.getD_some]
    rfl
  have hc : Function.Injective (nodesOf (List.ofFn f)) := by
    rw [hnodes]; exact (Nodes.radauNodes_injective hs).comp (finCongr hlen).injective
  have hM : nodePoly (nodesOf (List.ofFn f)) = nodePoly f := by
    rw [hnodes]; exact nodePoly_comp_equiv (finCongr hlen) f
  have key := checkOrder_colloc (List.ofFn f) hc
    (fun i => by rw [hnodes]; exact (Nodes.radauNodes_pos hs _).ne')
    ⟨s - 1, by omega⟩ (by simp [hlen])
    (by rw [hnodes]; exact Nodes.radauNodes_last hs _ (by simp))
    (fun j hj => by rw [hM]; exact radau_orth hs j (by omega))
    (by rw [hM, hlen]; exact (radau_J_neg hs).ne)
    (by rw [hlen]; simp [maxCheckOrder]; omega)
  rw [hlen] at key
  exact key

end Firkode.CheckOrder
