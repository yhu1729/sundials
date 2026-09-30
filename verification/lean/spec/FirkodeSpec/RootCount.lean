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
import FirkodeSpec.RuleOfSigns
import Mathlib.Analysis.Calculus.Deriv.Polynomial
import Mathlib.Topology.Order.IntermediateValue

/-!
# Real roots behind the `firk_gamma0` certificates

`FirkodeCert.Gamma0` certifies, for odd `s`, that `p_s(μ) = μ^s Q_s(1/μ)`
changes sign across `[L − u/2, L + u/2]`, where `L` is the literal
`firk_gamma0[s]` and `u` the unit in its last printed digit. By the
intermediate value theorem the real polynomial has a root there
(`gamma0_odd_root_exists`): the literal is within `u/2` of a real root.
With Descartes' rule of signs on the intervals of the certificate
(`RuleOfSigns`), this root is the only real root (`uniqueRootIn_spec`,
`gamma0_odd_unique`).
-/

namespace Firkode.RootCount

open Polynomial Gamma0 Radau

/-- A sign change of a real polynomial on `[lo, hi]` gives a root in
    `(lo, hi)`. -/
theorem exists_root_of_sign_change (p : ℝ[X]) {lo hi : ℝ} (hlh : lo ≤ hi)
    (h : p.eval lo * p.eval hi < 0) : ∃ x ∈ Set.Ioo lo hi, p.eval x = 0 := by
  have hcont : ContinuousOn (fun x => p.eval x) (Set.Icc lo hi) :=
    p.differentiable.continuous.continuousOn
  rcases mul_neg_iff.mp h with ⟨h1, h2⟩ | ⟨h1, h2⟩
  · obtain ⟨x, hx, hx0⟩ := intermediate_value_Ioo' hlh hcont ⟨h2, h1⟩
    exact ⟨x, hx, hx0⟩
  · obtain ⟨x, hx, hx0⟩ := intermediate_value_Ioo hlh hcont ⟨h1, h2⟩
    exact ⟨x, hx, hx0⟩

/-- For `s = 3, 5, 7, 9`, `det(μI − A) = μ^s Q_s(1/μ)` has a real root within
    half a unit in the last printed digit of `firk_gamma0[s]`. -/
theorem gamma0_odd_root_exists :
    ∀ s ∈ [3, 5, 7, 9], ∃ μ : ℝ,
      ((literal s - ulp s / 2 : ℚ) : ℝ) < μ ∧ μ < ((literal s + ulp s / 2 : ℚ) : ℝ) ∧
        (ratPoly (charPoly s)).eval μ = 0 := by
  intro s hs
  obtain ⟨hord, hsign⟩ := gamma0_odd_bracket s hs
  have hsign' : (ratPoly (charPoly s)).eval ((literal s - ulp s / 2 : ℚ) : ℝ) *
      (ratPoly (charPoly s)).eval ((literal s + ulp s / 2 : ℚ) : ℝ) < 0 := by
    rw [ratPoly_eval_ratCast, ratPoly_eval_ratCast, ← Rat.cast_mul, ← Rat.cast_zero, Rat.cast_lt]
    exact hsign
  obtain ⟨μ, ⟨h1, h2⟩, h3⟩ :=
    exists_root_of_sign_change _ (by exact_mod_cast hord.le) hsign'
  exact ⟨μ, h1, h2, h3⟩

/-- A passing `uniqueRootIn p lo hi` certificate: `p` has a root in
    `(lo, hi)`, and it is the only real root of `p`. -/
theorem uniqueRootIn_spec {p : Poly ℚ} {lo hi : ℚ} (h : uniqueRootIn p lo hi = true) :
    ∃ μ : ℝ, (lo : ℝ) < μ ∧ μ < hi ∧ (ratPoly p).eval μ = 0 ∧
      ∀ ν : ℝ, (ratPoly p).eval ν = 0 → ν = μ := by
  simp only [uniqueRootIn, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at h
  obtain ⟨⟨⟨⟨⟨⟨⟨⟨h0lo, hlohi, hhi2⟩, hp0, hplo, hphi, hp2⟩, hsign⟩, hneg⟩, hleft⟩, hmid⟩,
    hright⟩, htail⟩ := h
  have hsign' : (ratPoly p).eval (lo : ℝ) * (ratPoly p).eval (hi : ℝ) < 0 := by
    rw [ratPoly_eval_ratCast, ratPoly_eval_ratCast, ← Rat.cast_mul, ← Rat.cast_zero, Rat.cast_lt]
    exact hsign
  obtain ⟨μ, ⟨hμlo, hμhi⟩, hμ⟩ :=
    exists_root_of_sign_change _ (by exact_mod_cast hlohi.le) hsign'
  refine ⟨μ, hμlo, hμhi, hμ, fun ν hν => ?_⟩
  have h0lo' : (0 : ℝ) < lo := by exact_mod_cast h0lo
  have hhi2' : (hi : ℝ) < 2 := by exact_mod_cast hhi2
  have hne : ∀ q : ℚ, Poly.eval p q ≠ 0 → ν ≠ q := fun q hq hνq =>
    RuleOfSigns.ratPoly_eval_ne_zero hq (hνq ▸ hν)
  rcases lt_trichotomy ν 0 with hν0 | hν0 | hν0
  · exact absurd hν (RuleOfSigns.no_root_neg hp0 hneg hν0)
  · exact absurd (by simpa using hν0) (hne 0 hp0)
  rcases lt_trichotomy ν lo with hνlo | hνlo | hνlo
  · exact absurd hν (RuleOfSigns.no_root_Ioo hp0 hleft (by simpa using hν0) hνlo)
  · exact absurd hνlo (hne lo hplo)
  rcases lt_trichotomy ν hi with hνhi | hνhi | hνhi
  · exact RuleOfSigns.root_Ioo_unique hplo hmid.le hνlo hνhi hμlo hμhi hν hμ
  · exact absurd hνhi (hne hi hphi)
  rcases lt_trichotomy ν 2 with hν2 | hν2 | hν2
  · exact absurd hν (RuleOfSigns.no_root_Ioo hphi hright hνhi (by simpa using hν2))
  · exact absurd (by simpa using hν2) (hne 2 hp2)
  · exact absurd hν (RuleOfSigns.no_root_Ioi hp2 htail (by simpa using hν2))

/-- For `s = 3, 5, 7, 9`, `det(μI − A) = μ^s Q_s(1/μ)` has exactly one real
    root, and it lies within half a unit in the last printed digit of
    `firk_gamma0[s]`. -/
theorem gamma0_odd_unique :
    ∀ s ∈ [3, 5, 7, 9], ∃ μ : ℝ,
      ((literal s - ulp s / 2 : ℚ) : ℝ) < μ ∧ μ < ((literal s + ulp s / 2 : ℚ) : ℝ) ∧
        (ratPoly (charPoly s)).eval μ = 0 ∧
          ∀ ν : ℝ, (ratPoly (charPoly s)).eval ν = 0 → ν = μ :=
  fun s hs => uniqueRootIn_spec (gamma0_odd s hs)

/-- A passing `noRootIn p fuel a b` certificate with `p(a) ≠ 0`: `p` has no
    root in `(a, b)`. -/
theorem noRootIn_spec {p : Poly ℚ} :
    ∀ (fuel : ℕ) {a b : ℚ}, Poly.eval p a ≠ 0 → noRootIn p fuel a b = true →
      ∀ t : ℝ, (a : ℝ) < t → t < b → (ratPoly p).eval t ≠ 0 := by
  intro fuel
  induction fuel with
  | zero =>
    intro a b ha h t hat htb
    simp only [noRootIn, beq_iff_eq] at h
    exact RuleOfSigns.no_root_Ioo ha h hat htb
  | succ fuel ih =>
    intro a b ha h t hat htb
    simp only [noRootIn, Bool.or_eq_true, beq_iff_eq, Bool.and_eq_true, decide_eq_true_eq] at h
    rcases h with h | ⟨⟨hm, hl⟩, hr⟩
    · exact RuleOfSigns.no_root_Ioo ha h hat htb
    · rcases lt_trichotomy t (((a + b) / 2 : ℚ) : ℝ) with htm | htm | htm
      · exact ih ha hl t hat htm
      · rw [htm]; exact RuleOfSigns.ratPoly_eval_ne_zero hm
      · exact ih hm hr t htm htb

/-- A passing `noRealRoot p` certificate: `p` has no real root. -/
theorem noRealRoot_spec {p : Poly ℚ} (h : noRealRoot p = true) (t : ℝ) :
    (ratPoly p).eval t ≠ 0 := by
  simp only [noRealRoot, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at h
  obtain ⟨⟨⟨⟨hp0, hp2⟩, hneg⟩, hmid⟩, htail⟩ := h
  rcases lt_trichotomy t 0 with ht0 | ht0 | ht0
  · exact RuleOfSigns.no_root_neg hp0 hneg ht0
  · rw [ht0]; simpa using RuleOfSigns.ratPoly_eval_ne_zero hp0
  rcases lt_trichotomy t 2 with ht2 | ht2 | ht2
  · exact noRootIn_spec 6 hp0 hmid t (by simpa using ht0) (by simpa using ht2)
  · rw [ht2]; simpa using RuleOfSigns.ratPoly_eval_ne_zero hp2
  · exact RuleOfSigns.no_root_Ioi hp2 htail (by simpa using ht2)

/-- For `s = 2, 4, 6, 8`, `det(μI − A) = μ^s Q_s(1/μ)` has no real root. -/
theorem gamma0_even_no_root :
    ∀ s ∈ [2, 4, 6, 8], ∀ μ : ℝ, (ratPoly (charPoly s)).eval μ ≠ 0 :=
  fun s hs => noRealRoot_spec (gamma0_even_no_real s hs)

end Firkode.RootCount
