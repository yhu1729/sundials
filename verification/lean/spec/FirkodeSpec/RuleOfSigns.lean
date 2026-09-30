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

import FirkodeCert.Gamma0
import FirkodeSpec.Radau
import Mathlib.Algebra.Polynomial.RuleOfSigns

/-!
# Descartes' rule of signs on intervals

`FirkodeCert.Gamma0` counts sign variations of coefficient lists with
`Gamma0.signVariations`. Here that count is shown to be Mathlib's
`Polynomial.signVariations` of the real polynomial (`signVariations_ratPoly`),
so Descartes' rule (`Polynomial.roots_countP_pos_le_signVariations`) bounds
the positive roots. Transforming a polynomial first moves an interval to
`(0, ∞)`:

* `mobius p a b` has a positive root `x` for each root `(a + b x)/(1 + x)`
  of `p` in `(a, b)` (`eval_mobius`);
* `reflect p` is `p(−x)` (`eval_reflect`), for the roots below 0;
* `shift p a` is `p(x + a)` (`eval_shift`), for the roots above `a`.

No sign variation means no root in the interval, and one variation means at
most one root (`no_root_Ioo`, `root_Ioo_unique`, `no_root_neg`,
`no_root_Ioi`).
-/

namespace Firkode.RuleOfSigns

open Polynomial Radau

/-! ## Counting sign variations -/

section Lists

variable {α : Type*} [DecidableEq α]

/-- The number of adjacent unequal entries. -/
def adjNe (L : List α) : ℕ := ((L.zip L.tail).filter fun ab => decide (ab.1 ≠ ab.2)).length

theorem adjNe_cons_cons (a b : α) (L : List α) :
    adjNe (a :: b :: L) = (if a ≠ b then 1 else 0) + adjNe (b :: L) := by
  by_cases h : a = b <;> simp [adjNe, h, Nat.add_comm]

theorem length_destutter' (a : α) (L : List α) :
    (L.destutter' (· ≠ ·) a).length = adjNe (a :: L) + 1 := by
  induction L generalizing a with
  | nil => simp [adjNe]
  | cons b L ih =>
    rw [List.destutter'_cons, adjNe_cons_cons]
    by_cases h : a = b
    · subst h; simp [ih]
    · simp [h, ih]; omega

/-- `destutter` keeps one entry per run, so its length is one more than the
    number of changes between adjacent entries. -/
theorem length_destutter_sub_one (L : List α) :
    (L.destutter (· ≠ ·)).length - 1 = adjNe L := by
  cases L with
  | nil => simp [adjNe]
  | cons a L => rw [List.destutter_cons', length_destutter', Nat.add_sub_cancel]

end Lists

theorem mul_neg_iff_sign_ne {a b : ℚ} (ha : a ≠ 0) (hb : b ≠ 0) :
    a * b < 0 ↔ SignType.sign a ≠ SignType.sign b := by
  rcases ha.lt_or_gt with ha | ha <;> rcases hb.lt_or_gt with hb | hb <;>
    simp [sign_neg, sign_pos, ha, hb, mul_neg_of_neg_of_pos, mul_neg_of_pos_of_neg,
      mul_pos_of_neg_of_neg, le_of_lt, not_lt_of_ge]

theorem sign_ratCast (q : ℚ) : SignType.sign (q : ℝ) = SignType.sign q :=
  (Rat.cast_strictMono (K := ℝ)).sign_comp (F := ℚ →+* ℝ) (f := Rat.castHom ℝ) q

theorem coeff_ratPoly (l : Poly ℚ) (k : ℕ) : (ratPoly l).coeff k = ((get l k : ℚ) : ℝ) := by
  rw [ratPoly, coeff_map, ← toPoly_coeff]
  rfl

theorem map_get_range (l : List ℚ) : (List.range l.length).map (get l) = l := by
  refine List.ext_getElem (by simp) fun i h1 h2 => ?_
  simp [get, List.getD_eq_getElem?_getD, List.getElem?_eq_getElem h2]

/-- The nonzero coefficients of `ratPoly l`, highest power first, are the
    nonzero entries of `l.reverse`. -/
theorem coeffList_filter (l : Poly ℚ) :
    (ratPoly l).coeffList.filter (fun x => decide (x ≠ 0)) =
      (l.reverse.map (Rat.cast : ℚ → ℝ)).filter (fun x => decide (x ≠ 0)) := by
  have hrev : l.reverse.map (Rat.cast : ℚ → ℝ) =
      (List.range l.length).reverse.map (ratPoly l).coeff := by
    conv_lhs => rw [← map_get_range l]
    rw [List.map_reverse, List.map_map, List.map_reverse]
    congr 2
    funext k
    simp [coeff_ratPoly]
  rw [hrev]
  by_cases hP : ratPoly l = 0
  · rw [hP, coeffList_zero, List.filter_nil]
    symm
    rw [List.filter_eq_nil_iff]
    intro x hx
    simp only [List.mem_map] at hx
    obtain ⟨k, _, rfl⟩ := hx
    simp
  · have hlen : (ratPoly l).natDegree + 1 ≤ l.length := by
      rcases eq_or_ne l [] with rfl | hne
      · simp [ratPoly] at hP
      · have := (natDegree_map_le (f := Rat.castHom ℝ) (p := toPoly l)).trans_lt
          (natDegree_toPoly_lt hne)
        rw [← ratPoly] at this
        omega
    obtain ⟨m, hm⟩ : ∃ m, l.length = (ratPoly l).natDegree + 1 + m :=
      ⟨_, (Nat.add_sub_cancel' hlen).symm⟩
    rw [coeffList, withBotSucc_degree_eq_natDegree_add_one hP, hm,
      @List.range_add ((ratPoly l).natDegree + 1) m, List.reverse_append, List.map_append,
      List.filter_append]
    have htail : (((List.range m).map ((ratPoly l).natDegree + 1 + ·)).reverse.map
        (ratPoly l).coeff).filter (fun x => decide (x ≠ 0)) = [] := by
      rw [List.filter_eq_nil_iff]
      intro x hx
      simp only [List.mem_map, List.mem_reverse, List.mem_range] at hx
      obtain ⟨k, ⟨i, _, rfl⟩, rfl⟩ := hx
      simp [coeff_eq_zero_of_natDegree_lt (by omega :
        (ratPoly l).natDegree < (ratPoly l).natDegree + 1 + i)]
    rw [htail, List.nil_append]

/-- The sign-variation count of `FirkodeCert.Gamma0` is Mathlib's. -/
theorem signVariations_ratPoly (l : Poly ℚ) :
    Gamma0.signVariations l = (ratPoly l).signVariations := by
  set nz := l.reverse.filter fun q => decide (q ≠ 0)
  have hnz : ∀ q ∈ nz, q ≠ 0 := fun q hq => by simpa [nz] using (List.mem_filter.mp hq).2
  have hsigns : ((ratPoly l).coeffList.map SignType.sign).filter (fun x => decide (x ≠ 0)) =
      nz.map SignType.sign := by
    rw [List.filter_map]
    have : (((fun x : SignType => decide (x ≠ 0)) ∘ SignType.sign) : ℝ → Bool) =
        fun x => decide (x ≠ 0) := by
      funext x; simp
    rw [this, coeffList_filter, List.filter_map]
    have : (((fun x : ℝ => decide (x ≠ 0)) ∘ (Rat.cast : ℚ → ℝ)) : ℚ → Bool) =
        fun q => decide (q ≠ 0) := by
      funext q; simp
    rw [this, List.map_map]
    congr 1
    funext q
    exact sign_ratCast q
  rw [Polynomial.signVariations, hsigns, length_destutter_sub_one, adjNe, ← List.map_tail,
    List.zip_map, List.filter_map, List.length_map, Gamma0.signVariations]
  congr 1
  refine List.filter_congr fun ab hab => ?_
  obtain ⟨h1, h2⟩ := List.of_mem_zip hab
  have ha := hnz _ h1
  have hb := hnz _ (List.mem_of_mem_tail h2)
  simp [mul_neg_iff_sign_ne ha hb]

/-! ## Transformations moving an interval to `(0, ∞)` -/

theorem toPoly_foldl_add {R : Type} [CommRing R] (f : ℕ → Poly R) (n : ℕ) :
    toPoly ((List.range n).foldl (fun acc k => Poly.add acc (f k)) []) =
      ∑ k ∈ Finset.range n, toPoly (f k) := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [List.range_succ, List.foldl_append, List.foldl_cons, List.foldl_nil, toPoly_add, ih,
      Finset.sum_range_succ]

/-- `mobius p a b` at `x > −1` is `(1 + x)^n p((a + b x)/(1 + x))`, with
    `n + 1` the length of `p`. -/
theorem eval_mobius (p : Poly ℚ) (a b : ℚ) {x : ℝ} (hx : 1 + x ≠ 0) :
    (ratPoly (Gamma0.mobius p a b)).eval x =
      (1 + x) ^ (p.length - 1) * (ratPoly p).eval ((a + b * x) / (1 + x)) := by
  conv_rhs => rw [ratPoly_eval, Finset.mul_sum]
  conv_lhs => rw [ratPoly, Gamma0.mobius, toPoly_foldl_add, Polynomial.map_sum, eval_finsetSum]
  refine Finset.sum_congr rfl fun k hk => ?_
  have hk : k ≤ p.length - 1 := by have := Finset.mem_range.mp hk; omega
  simp only [toPoly_smul, toPoly_mul, toPoly_pow, toPoly_cons, toPoly_nil, Polynomial.map_mul,
    Polynomial.map_pow, Polynomial.map_add, Polynomial.map_C, Polynomial.map_X, eval_mul, eval_pow,
    eval_add, eval_C, eval_X, mul_zero, add_zero, Rat.cast_one, eq_ratCast]
  rw [div_pow, ← pow_mul_pow_sub (1 + x) hk]
  field_simp

/-- `shift p a` is `p(x + a)`. -/
theorem eval_shift (p : Poly ℚ) (a : ℚ) (x : ℝ) :
    (ratPoly (Gamma0.shift p a)).eval x = (ratPoly p).eval (x + a) := by
  simp [ratPoly, Gamma0.shift, Polynomial.map_comp, eval_comp, add_comm]

/-- `reflect p` is `p(−x)`. -/
theorem eval_reflect (p : Poly ℚ) (x : ℝ) :
    (ratPoly (NodePoly.reflect p)).eval x = (ratPoly p).eval (-x) := by
  rw [ratPoly_eval, ratPoly_eval, NodePoly.reflect, List.length_map, List.length_range]
  refine Finset.sum_congr rfl fun k hk => ?_
  rw [get_map_range _ _ (Finset.mem_range.mp hk), neg_pow]
  rcases Nat.mod_two_eq_zero_or_one k with h | h
  · rw [if_pos h, (Nat.even_iff.mpr h).neg_one_pow]; ring
  · rw [if_neg (by omega), (Nat.odd_iff.mpr h).neg_one_pow]; push_cast; ring

/-! ## Root counts -/

/-- Descartes: a nonzero polynomial without sign variations has no positive
    root. -/
theorem no_pos_root {P : ℝ[X]} (hP : P ≠ 0) (h : P.signVariations = 0) {x : ℝ} (hx : 0 < x) :
    P.eval x ≠ 0 := by
  intro hroot
  have hle := P.roots_countP_pos_le_signVariations
  have hpos : 0 < P.roots.countP (0 < ·) :=
    Multiset.countP_pos.mpr ⟨x, (mem_roots hP).mpr hroot, hx⟩
  omega

/-- Descartes: a nonzero polynomial with at most one sign variation has at
    most one positive root. -/
theorem pos_root_unique {P : ℝ[X]} (hP : P ≠ 0) (h : P.signVariations ≤ 1) {x y : ℝ}
    (hx : 0 < x) (hy : 0 < y) (hxr : P.eval x = 0) (hyr : P.eval y = 0) : x = y := by
  by_contra hxy
  have hle := P.roots_countP_pos_le_signVariations
  have hsub : ({x, y} : Multiset ℝ) ≤ P.roots.filter (0 < ·) := by
    rw [Multiset.le_iff_subset (by simp [hxy])]
    intro t ht
    simp only [Multiset.insert_eq_cons, Multiset.mem_cons, Multiset.mem_singleton] at ht
    rw [Multiset.mem_filter, mem_roots hP]
    rcases ht with rfl | rfl
    · exact ⟨hxr, hx⟩
    · exact ⟨hyr, hy⟩
  have h2 := Multiset.card_le_card hsub
  rw [← Multiset.countP_eq_card_filter] at h2
  simp at h2
  omega

theorem ratPoly_ne_zero_of_eval {p : Poly ℚ} {x : ℝ} (h : (ratPoly p).eval x ≠ 0) :
    ratPoly p ≠ 0 := fun h0 => h (by rw [h0, eval_zero])

theorem ratPoly_eval_ne_zero {p : Poly ℚ} {a : ℚ} (h : Poly.eval p a ≠ 0) :
    (ratPoly p).eval (a : ℝ) ≠ 0 := by
  rw [ratPoly_eval_ratCast]; exact_mod_cast h

/-- The point of `(0, ∞)` that `mobius p a b` associates with `t ∈ (a, b)`. -/
theorem mobius_point {a b t : ℝ} (hat : a < t) (htb : t < b) :
    0 < (t - a) / (b - t) ∧ 1 + (t - a) / (b - t) ≠ 0 ∧
      (a + b * ((t - a) / (b - t))) / (1 + (t - a) / (b - t)) = t := by
  have hbt : 0 < b - t := by linarith
  have hpos : 0 < (t - a) / (b - t) := div_pos (by linarith) hbt
  have h1 : 1 + (t - a) / (b - t) ≠ 0 := by positivity
  refine ⟨hpos, h1, ?_⟩
  rw [div_eq_iff h1]
  field_simp
  ring

theorem mobius_root {p : Poly ℚ} {a b : ℚ} {t : ℝ} (hat : (a : ℝ) < t) (htb : t < b)
    (ht : (ratPoly p).eval t = 0) :
    (ratPoly (Gamma0.mobius p a b)).eval ((t - a) / (b - t)) = 0 := by
  obtain ⟨_, h1, h2⟩ := mobius_point hat htb
  rw [eval_mobius p a b h1, h2, ht, mul_zero]

theorem mobius_ne_zero {p : Poly ℚ} {a b : ℚ} (ha : Poly.eval p a ≠ 0) :
    ratPoly (Gamma0.mobius p a b) ≠ 0 := by
  refine ratPoly_ne_zero_of_eval (x := 0) ?_
  rw [eval_mobius p a b (by norm_num)]
  simpa using ratPoly_eval_ne_zero ha

/-- No sign variation of `mobius p a b`: `p` has no root in `(a, b)`. -/
theorem no_root_Ioo {p : Poly ℚ} {a b : ℚ} (ha : Poly.eval p a ≠ 0)
    (h : Gamma0.signVariations (Gamma0.mobius p a b) = 0) {t : ℝ} (hat : (a : ℝ) < t)
    (htb : t < b) : (ratPoly p).eval t ≠ 0 := by
  intro ht
  rw [signVariations_ratPoly] at h
  exact no_pos_root (mobius_ne_zero ha) h (mobius_point hat htb).1 (mobius_root hat htb ht)

/-- At most one sign variation of `mobius p a b`: `p` has at most one root in
    `(a, b)`. -/
theorem root_Ioo_unique {p : Poly ℚ} {a b : ℚ} (ha : Poly.eval p a ≠ 0)
    (h : Gamma0.signVariations (Gamma0.mobius p a b) ≤ 1) {t u : ℝ} (hat : (a : ℝ) < t)
    (htb : t < b) (hau : (a : ℝ) < u) (hub : u < b) (ht : (ratPoly p).eval t = 0)
    (hu : (ratPoly p).eval u = 0) : t = u := by
  rw [signVariations_ratPoly] at h
  have hx := pos_root_unique (mobius_ne_zero ha) h (mobius_point hat htb).1
    (mobius_point hau hub).1 (mobius_root hat htb ht) (mobius_root hau hub hu)
  rw [← (mobius_point hat htb).2.2, ← (mobius_point hau hub).2.2, hx]

/-- No sign variation of `reflect p`: `p` has no negative root. -/
theorem no_root_neg {p : Poly ℚ} (h0 : Poly.eval p 0 ≠ 0)
    (h : Gamma0.signVariations (NodePoly.reflect p) = 0) {t : ℝ} (ht : t < 0) :
    (ratPoly p).eval t ≠ 0 := by
  have hne : ratPoly (NodePoly.reflect p) ≠ 0 := by
    refine ratPoly_ne_zero_of_eval (x := 0) ?_
    rw [eval_reflect, neg_zero]
    simpa using ratPoly_eval_ne_zero h0
  intro hroot
  rw [signVariations_ratPoly] at h
  exact no_pos_root hne h (neg_pos.mpr ht) (by rw [eval_reflect, neg_neg, hroot])

/-- No sign variation of `shift p a`: `p` has no root above `a`. -/
theorem no_root_Ioi {p : Poly ℚ} {a : ℚ} (ha : Poly.eval p a ≠ 0)
    (h : Gamma0.signVariations (Gamma0.shift p a) = 0) {t : ℝ} (ht : (a : ℝ) < t) :
    (ratPoly p).eval t ≠ 0 := by
  have hne : ratPoly (Gamma0.shift p a) ≠ 0 := by
    refine ratPoly_ne_zero_of_eval (x := 0) ?_
    rw [eval_shift, zero_add]
    exact ratPoly_eval_ne_zero ha
  intro hroot
  rw [signVariations_ratPoly] at h
  exact no_pos_root hne h (sub_pos.mpr ht) (by rw [eval_shift, sub_add_cancel, hroot])

end Firkode.RuleOfSigns
