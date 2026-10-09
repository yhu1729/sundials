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

import FirkodeCert.Tables
import FirkodeSpec.ListPoly
import Mathlib.NumberTheory.Real.Irrational
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.IntervalCases
import Mathlib.Tactic.LinearCombination

/-!
# From ℚ(√6) to the reals

`Q6.toReal` sends `a + b√6` to the real number it denotes. It preserves the
field operations of `Q6` and is injective because √6 is irrational, so the
exact certificates of `FirkodeCert` transfer to statements about real numbers.
-/

namespace Firkode

namespace Q6

@[simp] theorem add_a (x y : Q6) : (x + y).a = x.a + y.a := rfl
@[simp] theorem add_b (x y : Q6) : (x + y).b = x.b + y.b := rfl
@[simp] theorem sub_a (x y : Q6) : (x - y).a = x.a - y.a := rfl
@[simp] theorem sub_b (x y : Q6) : (x - y).b = x.b - y.b := rfl
@[simp] theorem neg_a (x : Q6) : (-x).a = -x.a := rfl
@[simp] theorem neg_b (x : Q6) : (-x).b = -x.b := rfl
@[simp] theorem mul_a (x y : Q6) : (x * y).a = x.a * y.a + 6 * x.b * y.b := rfl
@[simp] theorem mul_b (x y : Q6) : (x * y).b = x.a * y.b + x.b * y.a := rfl
@[simp] theorem inv_a (x : Q6) : (inv x).a = x.a / x.norm := rfl
@[simp] theorem inv_b (x : Q6) : (inv x).b = -x.b / x.norm := rfl
@[simp] theorem div_eq (x y : Q6) : x / y = x * inv y := rfl
@[simp] theorem zero_a : (0 : Q6).a = 0 := rfl
@[simp] theorem zero_b : (0 : Q6).b = 0 := rfl
@[simp] theorem one_a : (1 : Q6).a = 1 := rfl
@[simp] theorem one_b : (1 : Q6).b = 0 := rfl
@[simp] theorem ofNat_a (n : Nat) [n.AtLeastTwo] : (no_index (OfNat.ofNat n) : Q6).a = n := rfl
@[simp] theorem ofNat_b (n : Nat) [n.AtLeastTwo] : (no_index (OfNat.ofNat n) : Q6).b = 0 := rfl
@[simp] theorem natCast_a (n : Nat) : (n : Q6).a = n := rfl
@[simp] theorem natCast_b (n : Nat) : (n : Q6).b = 0 := rfl

theorem ext' {x y : Q6} (ha : x.a = y.a) (hb : x.b = y.b) : x = y := by
  cases x; cases y; simp_all

/-- The real number `a + b√6`. -/
noncomputable def toReal (x : Q6) : ℝ := x.a + x.b * √6

theorem sqrt6_mul_self : √6 * √6 = (6 : ℝ) := Real.mul_self_sqrt (by norm_num)

theorem irrational_sqrt6 : Irrational (√6 : ℝ) := by
  have : ¬IsSquare 6 := by
    rintro ⟨r, hr⟩
    have : r ≤ 3 := by nlinarith
    interval_cases r <;> omega
  exact_mod_cast irrational_sqrt_natCast_iff.mpr this

@[simp] theorem toReal_add (x y : Q6) : (x + y).toReal = x.toReal + y.toReal := by
  simp [toReal]; ring

@[simp] theorem toReal_sub (x y : Q6) : (x - y).toReal = x.toReal - y.toReal := by
  simp [toReal]; ring

@[simp] theorem toReal_neg (x : Q6) : (-x).toReal = -x.toReal := by
  simp [toReal]; ring

@[simp] theorem toReal_mul (x y : Q6) : (x * y).toReal = x.toReal * y.toReal := by
  simp only [toReal, mul_a, mul_b]
  push_cast
  linear_combination (-((x.b : ℝ) * y.b)) * sqrt6_mul_self

@[simp] theorem toReal_ofNat (n : Nat) [n.AtLeastTwo] :
    (no_index (OfNat.ofNat n) : Q6).toReal = OfNat.ofNat n := by
  simp only [toReal, ofNat_a, ofNat_b, Rat.cast_natCast, Rat.cast_zero, zero_mul, add_zero]
  exact Nat.cast_ofNat

@[simp] theorem toReal_natCast (n : Nat) : (n : Q6).toReal = n := by
  simp [toReal]

@[simp] theorem toReal_zero : (0 : Q6).toReal = 0 := by simp [toReal]

@[simp] theorem toReal_one : (1 : Q6).toReal = 1 := by simp [toReal]

@[simp] theorem toReal_ofRat (q : Rat) : (ofRat q).toReal = q := by
  simp [toReal, ofRat]

@[simp] theorem toReal_r6 : r6.toReal = √6 := by
  simp [toReal, r6]

/-- `a + b√6 = 0` only for `a = b = 0`. -/
theorem toReal_eq_zero {x : Q6} : x.toReal = 0 ↔ x = 0 := by
  constructor
  · intro h
    by_cases hb : x.b = 0
    · have ha : x.a = 0 := by simpa [toReal, hb] using h
      exact ext' (by simp [ha]) (by simp [hb])
    · exfalso
      apply irrational_sqrt6.ne_rat (-x.a / x.b)
      have hb' : (x.b : ℝ) ≠ 0 := by exact_mod_cast hb
      push_cast
      field_simp
      linarith [show (x.a : ℝ) + x.b * √6 = 0 from h]
  · rintro rfl; simp

theorem toReal_injective : Function.Injective toReal := by
  intro x y h
  have : (x - y).toReal = 0 := by simp [h]
  have h0 := toReal_eq_zero.mp this
  exact ext' (by have := congrArg Q6.a h0; simp at this; linarith)
    (by have := congrArg Q6.b h0; simp at this; linarith)

/-- The conjugate `a − b√6`. -/
def conj (x : Q6) : Q6 := ⟨x.a, -x.b⟩

theorem norm_eq_toReal_mul_conj (x : Q6) : (x.norm : ℝ) = x.toReal * (conj x).toReal := by
  simp only [norm, toReal, conj]
  push_cast
  linear_combination ((x.b : ℝ) * x.b) * sqrt6_mul_self

@[simp] theorem toReal_inv (x : Q6) : (inv x).toReal = x.toReal⁻¹ := by
  by_cases hx : x = 0
  · subst hx; simp [toReal, norm]
  · have hc : conj x ≠ 0 := by
      intro h; apply hx
      have ha := congrArg Q6.a h; have hb := congrArg Q6.b h
      simp [conj] at ha hb
      exact ext' (by simp [ha]) (by simp [hb])
    have h1 : x.toReal ≠ 0 := fun h => hx (toReal_eq_zero.mp h)
    have h2 : (conj x).toReal ≠ 0 := fun h => hc (toReal_eq_zero.mp h)
    have hn : (x.norm : ℝ) ≠ 0 := by rw [norm_eq_toReal_mul_conj]; exact mul_ne_zero h1 h2
    have hconj : (inv x).toReal = (conj x).toReal / (x.norm : ℝ) := by
      simp only [toReal, inv_a, inv_b, conj]
      push_cast
      field_simp
    rw [hconj, norm_eq_toReal_mul_conj]
    field_simp

@[simp] theorem toReal_div (x y : Q6) : (x / y).toReal = x.toReal / y.toReal := by
  simp [div_eq, div_eq_mul_inv]

end Q6

/-! ## List helpers under `toReal` -/

theorem sum_toReal (l : List Q6) : (sum l).toReal = (l.map Q6.toReal).sum := by
  simp only [sum]
  suffices ∀ a : Q6, (l.foldl (· + ·) a).toReal = a.toReal + (l.map Q6.toReal).sum by
    simpa using this 0
  induction l with
  | nil => simp
  | cons x l ih => intro a; simp [ih, add_assoc]

theorem npow_toReal (x : Q6) (n : Nat) : (npow x n).toReal = x.toReal ^ n := by
  induction n with
  | zero => simp [npow]
  | succ n ih => simp [npow, ih, pow_succ]

theorem get_toReal (v : List Q6) (i : Nat) : (get v i).toReal = get (v.map Q6.toReal) i := by
  simp only [get, List.getD_eq_getElem?_getD, List.getElem?_map]
  cases v[i]? <;> simp

theorem get2_toReal (A : List (List Q6)) (i j : Nat) :
    (get2 A i j).toReal = get2 (A.map (·.map Q6.toReal)) i j := by
  simp only [get2, List.getD_eq_getElem?_getD, List.getElem?_map]
  cases A[i]? with
  | none => simp
  | some r => simpa [get] using get_toReal r j

end Firkode
