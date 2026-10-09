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

import FirkodeCert.Basic
import Mathlib.Algebra.Polynomial.Derivative
import Mathlib.Algebra.BigOperators.Fin

/-!
# Coefficient lists as polynomials

`toPoly` reads a `Firkode.Poly` coefficient list (ascending powers) as a
Mathlib polynomial, and the lemmas here show that the list operations of
`FirkodeCert.Basic` are the polynomial operations. The list helpers `sum`,
`prod`, `npow` and `get` are related to their Mathlib counterparts as well.
-/

namespace Firkode

open Polynomial

section Helpers
variable {R : Type} [CommRing R]

theorem foldl_add_eq (a : R) (l : List R) : l.foldl (· + ·) a = a + l.sum := by
  induction l generalizing a with
  | nil => simp
  | cons x l ih => simp [ih, add_assoc]

theorem sum_eq (l : List R) : sum l = l.sum := by
  simp [sum, foldl_add_eq]

theorem foldl_mul_eq (a : R) (l : List R) : l.foldl (· * ·) a = a * l.prod := by
  induction l generalizing a with
  | nil => simp
  | cons x l ih => simp [ih, mul_assoc]

theorem prod_eq (l : List R) : prod l = l.prod := by
  simp [prod, foldl_mul_eq]

theorem npow_eq (x : R) (n : Nat) : npow x n = x ^ n := by
  induction n with
  | zero => simp [npow]
  | succ n ih => simp [npow, ih, pow_succ]

theorem sum_map_range (n : Nat) (f : Nat → R) :
    sum ((List.range n).map f) = ∑ i ∈ Finset.range n, f i := by
  rw [sum_eq]
  induction n with
  | zero => simp
  | succ n ih => simp [List.range_succ, Finset.sum_range_succ, ih]

end Helpers

theorem get_map_range {α : Type} [OfNat α 0] (n : Nat) (f : Nat → α) {i : Nat} (hi : i < n) :
    get ((List.range n).map f) i = f i := by
  simp [get, List.getD_eq_getElem?_getD, hi]

theorem get2_map_range {α : Type} [OfNat α 0] (n m : Nat) (f : Nat → Nat → α) {i j : Nat}
    (hi : i < n) (hj : j < m) :
    get2 ((List.range n).map fun i => (List.range m).map (f i)) i j = f i j := by
  simp [get2, List.getD_eq_getElem?_getD, hi, hj]

section Poly

variable {R : Type} [CommRing R]

/-- The polynomial with coefficient list `p`, lowest power first. -/
noncomputable def toPoly : Poly R → R[X]
  | [] => 0
  | a :: p => C a + X * toPoly p

@[simp] theorem toPoly_nil : toPoly ([] : Poly R) = 0 := rfl

@[simp] theorem toPoly_cons (a : R) (p : Poly R) : toPoly (a :: p) = C a + X * toPoly p := rfl

@[simp] theorem toPoly_add (p q : Poly R) : toPoly (Poly.add p q) = toPoly p + toPoly q := by
  induction p generalizing q with
  | nil => simp [Poly.add]
  | cons a p ih =>
    cases q with
    | nil => simp [Poly.add]
    | cons b q => simp [Poly.add, ih]; ring

@[simp] theorem toPoly_smul (c : R) (p : Poly R) : toPoly (Poly.smul c p) = C c * toPoly p := by
  induction p with
  | nil => simp [Poly.smul]
  | cons a p ih => simp_all [Poly.smul]; ring

@[simp] theorem toPoly_neg (p : Poly R) : toPoly (Poly.neg p) = -toPoly p := by
  induction p with
  | nil => simp [Poly.neg]
  | cons a p ih => simp_all [Poly.neg]; ring

@[simp] theorem toPoly_sub (p q : Poly R) : toPoly (Poly.sub p q) = toPoly p - toPoly q := by
  simp [Poly.sub, sub_eq_add_neg]

@[simp] theorem toPoly_mulX (p : Poly R) : toPoly (Poly.mulX p) = X * toPoly p := by
  simp [Poly.mulX]

@[simp] theorem toPoly_mul (p q : Poly R) : toPoly (Poly.mul p q) = toPoly p * toPoly q := by
  induction p with
  | nil => simp [Poly.mul]
  | cons a p ih => simp [Poly.mul, ih]; ring

@[simp] theorem toPoly_C (a : R) : toPoly (Poly.C a) = C a := by simp [Poly.C]

@[simp] theorem toPoly_X : toPoly (Poly.X : Poly R) = X := by simp [Poly.X]

@[simp] theorem toPoly_pow (p : Poly R) (n : Nat) : toPoly (Poly.pow p n) = toPoly p ^ n := by
  induction n with
  | zero => simp [Poly.pow]
  | succ n ih => simp [Poly.pow, ih, pow_succ, mul_comm]

theorem toPoly_eval (p : Poly R) (x : R) : Poly.eval p x = (toPoly p).eval x := by
  induction p with
  | nil => simp [Poly.eval]
  | cons a p ih => simp [Poly.eval, List.foldr] at *; rw [← ih]

theorem toPoly_coeff (p : Poly R) (k : Nat) : Poly.coeff p k = (toPoly p).coeff k := by
  induction p generalizing k with
  | nil => simp [Poly.coeff]
  | cons a p ih =>
    cases k with
    | zero => simp [Poly.coeff]
    | succ k => simp [Poly.coeff, coeff_X_mul, ← ih]

@[simp] theorem toPoly_comp (p q : Poly R) : toPoly (Poly.comp p q) = (toPoly p).comp (toPoly q) := by
  induction p with
  | nil => simp [Poly.comp]
  | cons a p ih => simp [Poly.comp, List.foldr] at *; rw [← ih]

theorem toPoly_derivAux (k : Nat) (p : Poly R) :
    toPoly (Poly.derivAux k p) = C (k : R) * toPoly p + X * derivative (toPoly p) := by
  induction p generalizing k with
  | nil => simp [Poly.derivAux]
  | cons a p ih =>
    simp only [Poly.derivAux, toPoly_cons, ih, derivative_add, derivative_C, derivative_mul,
      derivative_X, zero_add, one_mul, C_mul, Nat.cast_add, Nat.cast_one, C_add, C_1]
    ring

@[simp] theorem toPoly_deriv (p : Poly R) : toPoly (Poly.deriv p) = derivative (toPoly p) := by
  cases p with
  | nil => simp [Poly.deriv]
  | cons a p => simp [Poly.deriv, toPoly_derivAux]

theorem toPoly_replicate_zero (n : Nat) : toPoly (List.replicate n (0 : R)) = 0 := by
  induction n with
  | zero => rfl
  | succ n ih => simp [List.replicate_succ, ih]

theorem toPoly_append_replicate_zero (l : Poly R) (n : Nat) :
    toPoly (l ++ List.replicate n 0) = toPoly l := by
  induction l with
  | nil => simp [toPoly_replicate_zero]
  | cons a l ih => simp [ih]

/-- Trailing zero coefficients do not change the polynomial. -/
theorem toPoly_trim [DecidableEq R] (p : Poly R) : toPoly (Poly.trim p) = toPoly p := by
  have h := List.takeWhile_append_dropWhile (p := fun a : R => decide (a = 0)) (l := p.reverse)
  have hz : ∀ x ∈ p.reverse.takeWhile (fun a : R => decide (a = 0)), x = 0 := fun x hx => by
    simpa using List.all_eq_true.mp List.all_takeWhile x hx
  obtain ⟨n, hn⟩ : ∃ n, (p.reverse.takeWhile fun a : R => decide (a = 0)) = List.replicate n 0 :=
    ⟨_, List.eq_replicate_iff.mpr ⟨rfl, hz⟩⟩
  conv_rhs => rw [← List.reverse_reverse p, ← h, hn, List.reverse_append, List.reverse_replicate]
  rw [toPoly_append_replicate_zero]
  rfl

theorem toPoly_map {S : Type} [CommRing S] (f : R →+* S) (p : Poly R) :
    toPoly (p.map f) = (toPoly p).map f := by
  induction p with
  | nil => simp
  | cons a p ih => simp [ih]

end Poly

end Firkode
