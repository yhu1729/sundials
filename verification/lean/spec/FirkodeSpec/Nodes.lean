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

import FirkodeSpec.RootCount
import Mathlib.Data.Finset.Sort

/-!
# The Radau IIA nodes exist (s ≤ 9)

`firk_radau_nodes` scans `P_s − P_{s−1}` on a grid of `40 s` points of
`[−1, 1)`, and `NodePoly.nodes_grid` certifies that the scan finds `s − 1`
sign changes and no grid zeros. By the intermediate value theorem each of
these grid intervals contains a root, and together with the root `c = 1`
that gives `s` distinct roots of the monic node polynomial `M_s` of degree
`s`, hence all of them. So for every stage count FIRKODE supports:

* `radauNodes s` enumerates the roots of `M_s` in increasing order;
  `0 < c_1 < … < c_s = 1` (`radauNodes_strictMono`, `radauNodes_pos`,
  `radauNodes_last`), and `M_s` is their node polynomial
  (`radauNodes_nodePoly`);
* each grid interval where the scan of `firk_radau_nodes` sees a sign change
  contains exactly one node, and every interior node lies in one
  (`radauNodes_bracket`, `radauNodes_mem_bracket`), so the bisection of the C
  routine targets each node once;
* the results of `Radau` for the Radau nodes hold without hypotheses
  (`radauNodes_B`, `radauNodes_B2s_residual`, `radauNodes_e_last`,
  `radauNodes_pade`).
-/

namespace Firkode.Nodes

open Polynomial Finset Mirror Generated Radau

variable {s : ℕ}

/-! ## The two transcriptions of `firk_radau_poly` -/

/-- One step of the recurrence of `firk_radau_poly` on numbers. -/
def stepV (x : ℚ) (p : ℚ × ℚ) (k : ℕ) : ℚ × ℚ :=
  (p.2, (((2 * k + 1 : ℕ) : ℚ) * x * p.2 - (k : ℚ) * p.1) / ((k + 1 : ℕ) : ℚ))

/-- One step of the recurrence of `firk_radau_poly` on coefficient lists. -/
def stepP (p : Poly ℚ × Poly ℚ) (k : ℕ) : Poly ℚ × Poly ℚ :=
  (p.2, Poly.smul (1 / ((k + 1 : ℕ) : ℚ))
    (Poly.sub (Poly.smul ((2 * k + 1 : ℕ) : ℚ) (Poly.mulX p.2)) (Poly.smul (k : ℚ) p.1)))

/-- The recurrence of `firk_radau_poly` run on coefficient lists
    (`radauPoly`) evaluates to the recurrence run on numbers
    (`radauPolyAt`). -/
theorem radauPolyAt_eq_eval (s : ℕ) (x : ℚ) : radauPolyAt s x = Poly.eval (radauPoly s) x := by
  by_cases hs : s = 1
  · subst hs; simp [radauPolyAt, radauPoly, toPoly_eval]
  have hV : radauPolyAt s x = ((List.range' 1 (s - 1)).foldl (stepV x) (1, x)).2 -
      ((List.range' 1 (s - 1)).foldl (stepV x) (1, x)).1 := by
    simp only [radauPolyAt, if_neg hs]; rfl
  have hP : radauPoly s = Poly.sub ((List.range' 1 (s - 1)).foldl stepP (Poly.C 1, Poly.X)).2
      ((List.range' 1 (s - 1)).foldl stepP (Poly.C 1, Poly.X)).1 := by
    simp only [radauPoly, if_neg hs]; rfl
  let f : Poly ℚ × Poly ℚ → ℚ × ℚ := fun P => (Poly.eval P.1 x, Poly.eval P.2 x)
  have hstep : ∀ P k, stepV x (f P) k = f (stepP P k) := by
    intro P k
    simp only [f, stepV, stepP, toPoly_eval, toPoly_smul, toPoly_sub, toPoly_mulX, eval_mul,
      eval_C, eval_sub, eval_X, Prod.mk.injEq, true_and]
    push_cast
    ring
  have hinit : ((1 : ℚ), x) = f (Poly.C 1, Poly.X) := by simp [f, toPoly_eval]
  have hfold := List.foldl_hom f (l := List.range' 1 (s - 1)) (init := (Poly.C 1, Poly.X)) hstep
  rw [← hinit] at hfold
  rw [hV, hP, hfold]
  simp [f, toPoly_eval]

/-! ## The bracketing scan of `firk_radau_nodes` -/

/-- One step of the scan in `radauGridScan`. -/
def scanStep (acc : ℕ × ℕ) (f : ℚ × ℚ) : ℕ × ℕ :=
  if f.1 = 0 then (acc.1 + 1, acc.2 + 1) else if f.1 * f.2 < 0 then (acc.1 + 1, acc.2) else acc

theorem scan_foldl (L : List (ℚ × ℚ)) (acc : ℕ × ℕ) :
    L.foldl scanStep acc = (acc.1 + L.countP (fun f => decide (f.1 = 0 ∨ f.1 * f.2 < 0)),
      acc.2 + L.countP (fun f => decide (f.1 = 0))) := by
  induction L generalizing acc with
  | nil => simp
  | cons f L ih =>
    rw [List.foldl_cons, ih, List.countP_cons, List.countP_cons]
    unfold scanStep
    by_cases h1 : f.1 = 0
    · simp [h1]; omega
    · by_cases h2 : f.1 * f.2 < 0
      · simp [h1, h2]; omega
      · simp [h1, h2]

theorem zip_tail_map_range (g : ℕ → ℚ) (n : ℕ) :
    ((List.range n).map g).zip ((List.range n).map g).tail =
      (List.range (n - 1)).map fun j => (g j, g (j + 1)) := by
  refine List.ext_getElem (by simp) fun i h1 h2 => ?_
  simp

/-- Number of grid points of `firk_radau_nodes`. -/
abbrev gridN (s : ℕ) : ℕ := nodeGridPerStage * s

/-- The value of `firk_radau_poly` at grid point `j`. -/
def gridVal (s j : ℕ) : ℚ := radauPolyAt s (gridPoint (gridN s) j)

/-- The grid intervals `[x_j, x_{j+1}]` across which the scan of
    `firk_radau_nodes` sees a sign change (none for `s = 1`, where the C
    routine skips the scan). -/
def brackets (s : ℕ) : Finset ℕ :=
  if s = 1 then ∅ else (range (gridN s - 1)).filter fun j => gridVal s j * gridVal s (j + 1) < 0

theorem card_filter_range (p : ℕ → Prop) [DecidablePred p] (m : ℕ) :
    ((range m).filter p).card = (List.range m).countP (fun j => decide (p j)) := by
  rw [List.countP_eq_length_filter, ← List.toFinset_card_of_nodup ((List.nodup_range).filter _)]
  congr 1
  ext j
  simp

/-- The scan finds `s − 1` sign changes and no grid zeros. -/
theorem brackets_card (hs : s ∈ NodePoly.stages) : (brackets s).card = s - 1 := by
  by_cases hs1 : s = 1
  · simp [brackets, hs1]
  have h := (NodePoly.nodes_grid s hs).1
  have hv : gridValues nodeGridPerStage s = (List.range (gridN s)).map (gridVal s) := rfl
  simp only [radauGridScan, if_neg hs1, hv] at h
  change List.foldl scanStep (0, 0) (((List.range (gridN s)).map (gridVal s)).zip
    ((List.range (gridN s)).map (gridVal s)).tail) = _ at h
  rw [zip_tail_map_range, scan_foldl, List.countP_map, List.countP_map] at h
  simp only [Prod.mk.injEq, zero_add] at h
  obtain ⟨h1, h0⟩ := h
  rw [List.countP_eq_zero] at h0
  rw [brackets, if_neg hs1, card_filter_range, ← h1]
  refine List.countP_congr fun j hj => ?_
  have hj0 : gridVal s j ≠ 0 := by simpa using h0 j hj
  simp [hj0]

theorem gridPoint_lt {n j k : ℕ} (hn : 0 < n) (hjk : j < k) :
    (gridPoint n j : ℝ) < gridPoint n k := by
  have : (j : ℝ) < k := by exact_mod_cast hjk
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  simp only [gridPoint]
  push_cast
  gcongr

theorem mem_brackets {s j : ℕ} (hj : j ∈ brackets s) :
    s ≠ 1 ∧ j + 1 < gridN s ∧ gridVal s j * gridVal s (j + 1) < 0 := by
  by_cases hs1 : s = 1
  · simp [brackets, hs1] at hj
  simp only [brackets, if_neg hs1, mem_filter, mem_range] at hj
  exact ⟨hs1, by omega, hj.2⟩

theorem ratPoly_radauPoly_grid (s j : ℕ) :
    (ratPoly (radauPoly s)).eval (gridPoint (gridN s) j : ℝ) = gridVal s j := by
  rw [ratPoly_eval_ratCast, gridVal, radauPolyAt_eq_eval]

/-- Each bracket contains a root of `P_s − P_{s−1}`. -/
theorem exists_root_bracket {s j : ℕ} (hj : j ∈ brackets s) :
    ∃ x : ℝ, (gridPoint (gridN s) j : ℝ) < x ∧ x < gridPoint (gridN s) (j + 1) ∧
      (ratPoly (radauPoly s)).eval x = 0 := by
  obtain ⟨_, hjn, hsign⟩ := mem_brackets hj
  have hsign' : (ratPoly (radauPoly s)).eval (gridPoint (gridN s) j : ℝ) *
      (ratPoly (radauPoly s)).eval (gridPoint (gridN s) (j + 1) : ℝ) < 0 := by
    rw [ratPoly_radauPoly_grid, ratPoly_radauPoly_grid]; exact_mod_cast hsign
  obtain ⟨x, ⟨h1, h2⟩, h3⟩ := RootCount.exists_root_of_sign_change _
    (gridPoint_lt (by omega) (Nat.lt_succ_self j)).le hsign'
  exact ⟨x, h1, h2, h3⟩

theorem gridPoint_mono {n j k : ℕ} (hn : 0 < n) (hjk : j ≤ k) :
    (gridPoint n j : ℝ) ≤ gridPoint n k := by
  rcases eq_or_lt_of_le hjk with rfl | h
  · exact le_rfl
  · exact (gridPoint_lt hn h).le

theorem neg_one_le_gridPoint (n j : ℕ) : (-1 : ℝ) ≤ gridPoint n j := by
  have : (0 : ℝ) ≤ 2 * (j : ℝ) / n := by positivity
  simp only [gridPoint]; push_cast; linarith

theorem gridPoint_lt_one {n j : ℕ} (hjn : j < n) : (gridPoint n j : ℝ) < 1 := by
  have hn : (0 : ℝ) < n := by exact_mod_cast (show 0 < n by omega)
  have : 2 * (j : ℝ) / n < 2 := by
    rw [div_lt_iff₀ hn]; exact_mod_cast (show 2 * j < 2 * n by omega)
  simp only [gridPoint]; push_cast; linarith

/-- Distinct brackets are disjoint. -/
theorem bracket_unique {s j k : ℕ} (hj : j ∈ brackets s) {x : ℝ}
    (hxj : (gridPoint (gridN s) j : ℝ) < x ∧ x < gridPoint (gridN s) (j + 1))
    (hxk : (gridPoint (gridN s) k : ℝ) < x ∧ x < gridPoint (gridN s) (k + 1)) : j = k := by
  have hn : 0 < gridN s := by have := (mem_brackets hj).2.1; omega
  by_contra hne
  rcases lt_or_gt_of_ne hne with h | h
  · linarith [gridPoint_mono (n := gridN s) hn (show j + 1 ≤ k by omega)]
  · linarith [gridPoint_mono (n := gridN s) hn (show k + 1 ≤ j by omega)]

/-! ## The nodes -/

/-- A root `x ∈ (x_j, x_{j+1})` of `P_s − P_{s−1}` in bracket `j`. -/
noncomputable def bracketRoot (s j : ℕ) : ℝ :=
  if h : j ∈ brackets s then Classical.choose (exists_root_bracket h) else 0

theorem bracketRoot_spec {j : ℕ} (hj : j ∈ brackets s) :
    ((gridPoint (gridN s) j : ℝ) < bracketRoot s j ∧
      bracketRoot s j < gridPoint (gridN s) (j + 1)) ∧
      (ratPoly (radauPoly s)).eval (bracketRoot s j) = 0 := by
  rw [bracketRoot, dif_pos hj]
  obtain ⟨h1, h2, h3⟩ := Classical.choose_spec (exists_root_bracket hj)
  exact ⟨⟨h1, h2⟩, h3⟩

/-- The node `c = (x + 1)/2` in bracket `j`. -/
noncomputable def bracketNode (s j : ℕ) : ℝ := (bracketRoot s j + 1) / 2

theorem bracketNode_mem {j : ℕ} (hj : j ∈ brackets s) :
    0 < bracketNode s j ∧ bracketNode s j < 1 := by
  obtain ⟨⟨h1, h2⟩, _⟩ := bracketRoot_spec hj
  have := neg_one_le_gridPoint (gridN s) j
  have := gridPoint_lt_one (mem_brackets hj).2.1
  constructor <;> simp only [bracketNode] <;> linarith

theorem bracketNode_injOn : Set.InjOn (bracketNode s) (brackets s) := by
  intro j hj k hk hjk
  have h : bracketRoot s j = bracketRoot s k := by
    simp only [bracketNode] at hjk; linarith
  exact bracket_unique hj (bracketRoot_spec hj).1 (h ▸ (bracketRoot_spec hk).1)

/-- `M_s(c)` is a nonzero multiple of `(P_s − P_{s−1})(2c − 1)`. -/
theorem eval_nodePoly (s : ℕ) (c : ℝ) :
    (ratPoly (NodePoly.nodePoly s)).eval c =
      ((1 / (Poly.trim (Poly.comp (radauPoly s) [-1, 2])).getLastD 1 : ℚ) : ℝ) *
        (ratPoly (radauPoly s)).eval (2 * c - 1) := by
  simp only [NodePoly.nodePoly, ratPoly_smul, eval_mul, eval_C]
  congr 1
  rw [ratPoly, toPoly_trim, toPoly_comp, Polynomial.map_comp, eval_comp, ← ratPoly]
  congr 1
  simp [ratPoly]; ring

theorem bracketNode_root {j : ℕ} (hj : j ∈ brackets s) :
    (ratPoly (NodePoly.nodePoly s)).eval (bracketNode s j) = 0 := by
  rw [eval_nodePoly, bracketNode, show 2 * ((bracketRoot s j + 1) / 2) - 1 = bracketRoot s j by ring,
    (bracketRoot_spec hj).2, mul_zero]

theorem natDegree_nodePoly (hs : s ∈ NodePoly.stages) :
    (ratPoly (NodePoly.nodePoly s)).natDegree = s := by
  have hcoeff : (ratPoly (NodePoly.nodePoly s)).coeff s = 1 := by
    rw [RuleOfSigns.coeff_ratPoly]
    have := NodePoly.nodePoly_monic s hs
    exact_mod_cast this
  have hlen := (NodePoly.nodePoly_degree_endpoint s hs).1
  have hle : (ratPoly (NodePoly.nodePoly s)).natDegree < s + 1 := by
    rw [← hlen]
    exact natDegree_map_le.trans_lt (natDegree_toPoly_lt (List.ne_nil_of_length_eq_add_one hlen))
  exact le_antisymm (by omega) (le_natDegree_of_ne_zero (by rw [hcoeff]; exact one_ne_zero))

theorem monic_nodePoly (hs : s ∈ NodePoly.stages) : (ratPoly (NodePoly.nodePoly s)).Monic := by
  rw [Monic, leadingCoeff, natDegree_nodePoly hs, RuleOfSigns.coeff_ratPoly]
  exact_mod_cast NodePoly.nodePoly_monic s hs

/-- The real roots of `M_s`. -/
noncomputable def rootFinset (s : ℕ) : Finset ℝ := (ratPoly (NodePoly.nodePoly s)).roots.toFinset

/-- The roots of `M_s` are the bracket nodes and 1; there are `s` of them. -/
theorem rootFinset_eq (hs : s ∈ NodePoly.stages) :
    rootFinset s = (brackets s).image (bracketNode s) ∪ {1} ∧ (rootFinset s).card = s := by
  have hM := (monic_nodePoly hs).ne_zero
  have hs0 : 1 ≤ s := by simp [NodePoly.stages] at hs; omega
  have hsub : (brackets s).image (bracketNode s) ∪ {1} ⊆ rootFinset s := by
    intro x hx
    rw [rootFinset, Multiset.mem_toFinset, mem_roots hM, IsRoot.def]
    rcases mem_union.mp hx with hx | hx
    · obtain ⟨j, hj, rfl⟩ := mem_image.mp hx
      exact bracketNode_root hj
    · rw [mem_singleton.mp hx]
      have := (NodePoly.nodePoly_degree_endpoint s hs).2
      have h1 := ratPoly_eval_ratCast (NodePoly.nodePoly s) 1
      rw [this] at h1
      simpa using h1
  have hcard : ((brackets s).image (bracketNode s) ∪ {1}).card = s := by
    rw [card_union_of_disjoint, card_image_of_injOn bracketNode_injOn, brackets_card hs,
      card_singleton]
    · omega
    · rw [disjoint_singleton_right, mem_image]
      rintro ⟨j, hj, h⟩
      exact absurd h (bracketNode_mem hj).2.ne
  have hle : (rootFinset s).card ≤ s :=
    (Multiset.toFinset_card_le _).trans ((card_roots' _).trans (natDegree_nodePoly hs).le)
  have heq := eq_of_subset_of_card_le hsub (by omega)
  exact ⟨heq.symm, by rw [← heq, hcard]⟩

theorem mem_rootFinset (hs : s ∈ NodePoly.stages) {x : ℝ} (hx : x ∈ rootFinset s) :
    x = 1 ∨ ∃ j ∈ brackets s, x = bracketNode s j := by
  rw [(rootFinset_eq hs).1, mem_union, mem_singleton, mem_image] at hx
  rcases hx with ⟨j, hj, rfl⟩ | rfl
  · exact Or.inr ⟨j, hj, rfl⟩
  · exact Or.inl rfl

theorem rootFinset_bounds (hs : s ∈ NodePoly.stages) {x : ℝ} (hx : x ∈ rootFinset s) :
    0 < x ∧ x ≤ 1 := by
  rcases mem_rootFinset hs hx with rfl | ⟨j, hj, rfl⟩
  · norm_num
  · exact ⟨(bracketNode_mem hj).1, (bracketNode_mem hj).2.le⟩

/-- The Radau IIA nodes for `s` stages, in increasing order: the roots of
    `M_s` (for `s` in `stages`). -/
noncomputable def radauNodes (s : ℕ) : Fin s → ℝ :=
  if h : (rootFinset s).card = s then fun i => (rootFinset s).orderEmbOfFin h i else fun _ => 0

theorem radauNodes_eq (hs : s ∈ NodePoly.stages) (i : Fin s) :
    radauNodes s i = (rootFinset s).orderEmbOfFin (rootFinset_eq hs).2 i := by
  rw [radauNodes, dif_pos (rootFinset_eq hs).2]

theorem radauNodes_mem (hs : s ∈ NodePoly.stages) (i : Fin s) : radauNodes s i ∈ rootFinset s := by
  rw [radauNodes_eq hs]; exact orderEmbOfFin_mem _ _ i

theorem radauNodes_strictMono (hs : s ∈ NodePoly.stages) : StrictMono (radauNodes s) := by
  intro i j hij
  rw [radauNodes_eq hs, radauNodes_eq hs]
  exact ((rootFinset s).orderEmbOfFin (rootFinset_eq hs).2).strictMono hij

theorem radauNodes_injective (hs : s ∈ NodePoly.stages) : Function.Injective (radauNodes s) :=
  (radauNodes_strictMono hs).injective

theorem radauNodes_pos (hs : s ∈ NodePoly.stages) (i : Fin s) : 0 < radauNodes s i :=
  (rootFinset_bounds hs (radauNodes_mem hs i)).1

theorem radauNodes_le_one (hs : s ∈ NodePoly.stages) (i : Fin s) : radauNodes s i ≤ 1 :=
  (rootFinset_bounds hs (radauNodes_mem hs i)).2

/-- The last node is `c_s = 1`. -/
theorem radauNodes_last (hs : s ∈ NodePoly.stages) (i : Fin s) (hi : (i : ℕ) = s - 1) :
    radauNodes s i = 1 := by
  have hs0 : 0 < s := Fin.pos i
  have hi' : i = ⟨s - 1, Nat.sub_lt hs0 (Nat.succ_pos 0)⟩ := Fin.ext hi
  rw [radauNodes_eq hs, hi', orderEmbOfFin_last (rootFinset_eq hs).2 hs0]
  have h1 : (1 : ℝ) ∈ rootFinset s := by
    rw [(rootFinset_eq hs).1]; exact mem_union_right _ (mem_singleton_self 1)
  exact le_antisymm (max'_le _ _ _ fun x hx => (rootFinset_bounds hs hx).2) (le_max' _ _ h1)

theorem radauNodes_root (hs : s ∈ NodePoly.stages) (i : Fin s) :
    (ratPoly (NodePoly.nodePoly s)).eval (radauNodes s i) = 0 := by
  have := radauNodes_mem hs i
  rw [rootFinset, Multiset.mem_toFinset, mem_roots (monic_nodePoly hs).ne_zero] at this
  exact this

/-- `M_s` is the node polynomial of the Radau nodes. -/
theorem radauNodes_nodePoly (hs : s ∈ NodePoly.stages) :
    Colloc.nodePoly (radauNodes s) = ratPoly (NodePoly.nodePoly s) := by
  have hc := radauNodes_injective hs
  have h1 : (Colloc.nodePoly (radauNodes s)).Monic := Colloc.nodePoly_monic
  have h2 := monic_nodePoly hs
  have hd1 : (Colloc.nodePoly (radauNodes s)).natDegree = s := Colloc.natDegree_nodePoly
  refine Polynomial.eq_of_degree_sub_lt_of_eval_index_eq (s := univ) (v := radauNodes s)
    hc.injOn ?_ ?_
  · rw [card_univ, Fintype.card_fin]
    have hlt := degree_sub_lt (p := Colloc.nodePoly (radauNodes s))
      (q := ratPoly (NodePoly.nodePoly s))
      (by rw [degree_eq_natDegree h1.ne_zero, degree_eq_natDegree h2.ne_zero, hd1,
        natDegree_nodePoly hs])
      h1.ne_zero (by rw [h1.leadingCoeff, h2.leadingCoeff])
    rwa [degree_eq_natDegree h1.ne_zero, hd1] at hlt
  · intro i _
    rw [Colloc.nodePoly_eval_node, radauNodes_root hs i]

/-- Each bracket of the scan of `firk_radau_nodes` holds exactly one node
    (in the variable `x = 2c − 1` of the C routine). -/
theorem radauNodes_bracket (hs : s ∈ NodePoly.stages) {j : ℕ} (hj : j ∈ brackets s) :
    ∃! i : Fin s, (gridPoint (gridN s) j : ℝ) < 2 * radauNodes s i - 1 ∧
      2 * radauNodes s i - 1 < gridPoint (gridN s) (j + 1) := by
  have hin : ∀ i : Fin s, (gridPoint (gridN s) j : ℝ) < 2 * radauNodes s i - 1 →
      2 * radauNodes s i - 1 < gridPoint (gridN s) (j + 1) → radauNodes s i = bracketNode s j := by
    intro i h1 h2
    rcases mem_rootFinset hs (radauNodes_mem hs i) with h | ⟨k, hk, h⟩
    · have := gridPoint_lt_one (mem_brackets hj).2.1
      rw [h] at h2; linarith
    · rw [h] at h1 h2 ⊢
      have hjk := bracket_unique hj ⟨by simp only [bracketNode] at h1; linarith,
        by simp only [bracketNode] at h2; linarith⟩ (bracketRoot_spec hk).1
      rw [hjk]
  have hmem : bracketNode s j ∈ rootFinset s := by
    rw [(rootFinset_eq hs).1]; exact mem_union_left _ (mem_image_of_mem _ hj)
  have hrange := range_orderEmbOfFin (rootFinset s) (rootFinset_eq hs).2
  obtain ⟨i, hi⟩ : ∃ i, (rootFinset s).orderEmbOfFin (rootFinset_eq hs).2 i = bracketNode s j := by
    have : bracketNode s j ∈ Set.range ((rootFinset s).orderEmbOfFin (rootFinset_eq hs).2) := by
      rw [hrange]; exact hmem
    exact this
  rw [← radauNodes_eq hs] at hi
  obtain ⟨⟨h1, h2⟩, _⟩ := bracketRoot_spec hj
  refine ⟨i, ⟨?_, ?_⟩, fun i' ⟨h1', h2'⟩ => radauNodes_injective hs ?_⟩
  · rw [hi, bracketNode]; linarith
  · rw [hi, bracketNode]; linarith
  · rw [hin i' h1' h2', hi]

/-- Every interior node lies in a bracket of the scan. -/
theorem radauNodes_mem_bracket (hs : s ∈ NodePoly.stages) (i : Fin s) (hi : radauNodes s i < 1) :
    ∃ j ∈ brackets s, (gridPoint (gridN s) j : ℝ) < 2 * radauNodes s i - 1 ∧
      2 * radauNodes s i - 1 < gridPoint (gridN s) (j + 1) := by
  rcases mem_rootFinset hs (radauNodes_mem hs i) with h | ⟨j, hj, h⟩
  · exact absurd h hi.ne
  · obtain ⟨⟨h1, h2⟩, _⟩ := bracketRoot_spec hj
    refine ⟨j, hj, ?_, ?_⟩ <;> rw [h, bracketNode] <;> linarith

/-! ## The results of `Radau` for the Radau nodes -/

/-- B(2s−1): `Σ_j b_j c_j^k = 1/(k+1)` for `k ≤ 2s − 2`. -/
theorem radauNodes_B (hs : s ∈ NodePoly.stages) (last : Fin s) (hlast : (last : ℕ) = s - 1)
    {k : ℕ} (hk : k < 2 * s - 1) :
    ∑ j, Colloc.collocA (radauNodes s) last j * radauNodes s j ^ k = 1 / (k + 1) :=
  radau_B hs (radauNodes_injective hs) (radauNodes_nodePoly hs) last
    (radauNodes_last hs last hlast) hk

/-- The B(2s) residual `(s!(s−1)!)²/((2s)!(2s−1)!)`. -/
theorem radauNodes_B2s_residual (hs : s ∈ NodePoly.stages) (last : Fin s)
    (hlast : (last : ℕ) = s - 1) :
    ∑ j, Colloc.collocA (radauNodes s) last j * radauNodes s j ^ (2 * s - 1) - 1 / (2 * s) =
      ((s.factorial * (s - 1).factorial) ^ 2 : ℝ) /
        ((2 * s).factorial * (2 * s - 1).factorial) :=
  radau_B2s_residual hs (radauNodes_injective hs) (radauNodes_nodePoly hs) last
    (radauNodes_last hs last hlast)

/-- `e_s = (−1)^s γ0/s`. -/
theorem radauNodes_e_last (hs : s ∈ NodePoly.stages) (last : Fin s)
    (hlast : (last : ℕ) = s - 1) (γ0 : ℝ) :
    -γ0 * ((Colloc.ell (radauNodes s) last).eval 0 / radauNodes s last) = (-1) ^ s * γ0 / s :=
  radau_e_last hs (radauNodes_nodePoly hs) last (radauNodes_last hs last hlast) γ0

/-- The stability function is the `(s−1, s)` Padé approximant of `e^z`. -/
theorem radauNodes_pade (hs : s ∈ NodePoly.stages) (last : Fin s) (hlast : (last : ℕ) = s - 1)
    (z : ℝ) (hz : IsUnit (1 - z • Colloc.collocA (radauNodes s)).det) :
    Colloc.stabR (Colloc.collocA (radauNodes s)) ((Colloc.collocA (radauNodes s)).row last) z *
        (ratPoly (padeDen s)).eval z = (ratPoly (padeNum s)).eval z :=
  radau_pade hs (radauNodes_injective hs) (radauNodes_nodePoly hs) last
    (radauNodes_last hs last hlast) z hz

end Firkode.Nodes
