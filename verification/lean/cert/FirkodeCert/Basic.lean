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

/-!
# Exact arithmetic for the FIRKODE certificates

* `Q6` is the field ℚ(√6) with exact rational coordinates. The closed-form
  Radau IIA tables of `firk_radau_closed_form` (s ≤ 3) have entries in it.
* `Poly α` is a dense coefficient list in ascending powers.

Everything here is computable and reduces in the kernel, so the certificates
of this library are checked with `decide +kernel`.
-/

namespace Firkode

/-- `a + b·√6` with rational `a` and `b`. -/
structure Q6 where
  a : Rat
  b : Rat
  deriving DecidableEq, Repr

namespace Q6

instance {n : Nat} : OfNat Q6 n := ⟨⟨(n : Rat), 0⟩⟩
instance : NatCast Q6 := ⟨fun n => ⟨(n : Rat), 0⟩⟩
instance : Inhabited Q6 := ⟨0⟩
instance : Add Q6 := ⟨fun x y => ⟨x.a + y.a, x.b + y.b⟩⟩
instance : Sub Q6 := ⟨fun x y => ⟨x.a - y.a, x.b - y.b⟩⟩
instance : Neg Q6 := ⟨fun x => ⟨-x.a, -x.b⟩⟩
instance : Mul Q6 := ⟨fun x y => ⟨x.a * y.a + 6 * x.b * y.b, x.a * y.b + x.b * y.a⟩⟩

/-- Field norm `a² − 6b²`, nonzero for `x ≠ 0` because √6 is irrational. -/
def norm (x : Q6) : Rat := x.a * x.a - 6 * x.b * x.b

/-- Inverse through the conjugate; `inv 0 = 0`, as for `Rat`. -/
def inv (x : Q6) : Q6 := ⟨x.a / norm x, -x.b / norm x⟩

instance : Div Q6 := ⟨fun x y => x * inv y⟩

/-- Embedding of the rationals. -/
def ofRat (q : Rat) : Q6 := ⟨q, 0⟩

/-- √6, the value of `r6 = sqrtl(6.0L)` in `firk_radau_closed_form`. -/
def r6 : Q6 := ⟨0, 1⟩

/-- Exact test `0 < a + b√6`. -/
def isPos (x : Q6) : Bool :=
  if 0 ≤ x.a ∧ 0 ≤ x.b then decide (0 < x.a ∨ 0 < x.b)
  else if x.a ≤ 0 ∧ x.b ≤ 0 then false
  else if 0 < x.a then decide (6 * x.b * x.b < x.a * x.a)
  else decide (x.a * x.a < 6 * x.b * x.b)

instance : LT Q6 := ⟨fun x y => isPos (y - x) = true⟩
instance : LE Q6 := ⟨fun x y => x = y ∨ isPos (y - x) = true⟩

instance (x y : Q6) : Decidable (x < y) :=
  inferInstanceAs (Decidable (isPos (y - x) = true))
instance (x y : Q6) : Decidable (x ≤ y) :=
  inferInstanceAs (Decidable (x = y ∨ isPos (y - x) = true))

end Q6

/-- Factorial (core Lean has no `Nat.factorial`). -/
def fact : Nat → Nat
  | 0 => 1
  | n + 1 => (n + 1) * fact n

section Helpers
variable {α : Type} [Add α] [Mul α] [OfNat α 0] [OfNat α 1]

/-- Power by repeated multiplication, `x^(n+1) = x^n * x`, as the C loops
    accumulate it. -/
def npow (x : α) : Nat → α
  | 0 => 1
  | n + 1 => npow x n * x

/-- Left-to-right sum, as `sum += ...` in C. -/
def sum (l : List α) : α := l.foldl (· + ·) 0

/-- Entry `i` of a vector (0 when out of range). -/
def get (v : List α) (i : Nat) : α := v.getD i 0

/-- Entry `(i, j)` of a row-major matrix (0 when out of range). -/
def get2 (A : List (List α)) (i j : Nat) : α := (A.getD i []).getD j 0

/-- Product of the entries. -/
def prod (l : List α) : α := l.foldl (· * ·) 1

end Helpers

/-- Dense polynomial: coefficients in ascending powers. -/
abbrev Poly (α : Type) := List α

namespace Poly

variable {α : Type} [Add α] [Mul α] [Neg α] [OfNat α 0] [OfNat α 1]

def add : Poly α → Poly α → Poly α
  | [], q => q
  | p, [] => p
  | a :: p, b :: q => (a + b) :: add p q

def smul (c : α) (p : Poly α) : Poly α := p.map (c * ·)

def neg (p : Poly α) : Poly α := p.map (- ·)

def sub (p q : Poly α) : Poly α := add p (neg q)

/-- Multiplication by the indeterminate. -/
def mulX (p : Poly α) : Poly α := 0 :: p

def mul : Poly α → Poly α → Poly α
  | [], _ => []
  | a :: p, q => add (smul a q) (mulX (mul p q))

/-- Constant polynomial. -/
def C (a : α) : Poly α := [a]

/-- The indeterminate. -/
def X : Poly α := [0, 1]

def pow (p : Poly α) : Nat → Poly α
  | 0 => [1]
  | n + 1 => mul p (pow p n)

/-- Coefficient of the `k`-th power. -/
def coeff (p : Poly α) (k : Nat) : α := p.getD k 0

/-- Horner evaluation. -/
def eval (p : Poly α) (x : α) : α := p.foldr (fun a acc => a + x * acc) 0

/-- Composition `p ∘ q`. -/
def comp (p q : Poly α) : Poly α := p.foldr (fun a acc => add (C a) (mul q acc)) []

section Deriv
variable [NatCast α]

def derivAux : Nat → Poly α → Poly α
  | _, [] => []
  | k, a :: p => ((k : α) * a) :: derivAux (k + 1) p

/-- Formal derivative. -/
def deriv : Poly α → Poly α
  | [] => []
  | _ :: p => derivAux 1 p

end Deriv

/-- Drop trailing zero coefficients, so that equal polynomials have equal
    coefficient lists. -/
def trim [DecidableEq α] (p : Poly α) : Poly α :=
  (p.reverse.dropWhile (fun a => decide (a = 0))).reverse

end Poly

/-! ## Division with remainder and extended Euclid over ℚ

Used only to compute certificates: the identities they must satisfy are
checked by the kernel, so these routines need no correctness proof. -/

namespace Poly

/-- Division with remainder by a trimmed `q ≠ 0`: `(d, r)` with `p = d q + r`
    and `deg r < deg q`, in at most `fuel` steps. -/
def divModAux (q : Poly Rat) : Nat → Poly Rat → Poly Rat → Poly Rat × Poly Rat
  | 0, d, r => (d, trim r)
  | fuel + 1, d, r =>
    let r := trim r
    if r.length < q.length then (d, r)
    else
      let m : Poly Rat := List.replicate (r.length - q.length) 0 ++ [r.getLastD 0 / q.getLastD 1]
      divModAux q fuel (add d m) (sub r (mul m q))

/-- Division with remainder. -/
def divMod (p q : Poly Rat) : Poly Rat × Poly Rat := divModAux (trim q) p.length [] p

/-- Euclid's algorithm on `(r0, r1)` with the cofactors `(s0, s1)` and
    `(t0, t1)`. -/
def xgcdAux :
    Nat → Poly Rat → Poly Rat → Poly Rat → Poly Rat → Poly Rat → Poly Rat →
      Poly Rat × Poly Rat × Poly Rat
  | 0, r0, _, s0, _, t0, _ => (r0, s0, t0)
  | fuel + 1, r0, r1, s0, s1, t0, t1 =>
    if trim r1 = [] then (r0, s0, t0)
    else
      let dr := divMod r0 r1
      xgcdAux fuel r1 dr.2 s1 (sub s0 (mul dr.1 s1)) t1 (sub t0 (mul dr.1 t1))

/-- Extended Euclid: `(g, u, v)` with `u p + v q = g` a greatest common
    divisor. -/
def xgcd (p q : Poly Rat) : Poly Rat × Poly Rat × Poly Rat :=
  xgcdAux (p.length + q.length + 1) p q [1] [] [] [1]

end Poly

/-! ## Small dense matrices (row-major lists) -/

section Matrices
variable {α : Type} [Add α] [Sub α] [Mul α] [Div α] [Neg α] [OfNat α 0] [OfNat α 1]

/-- The rows `rows` and columns `cols` of `A`. -/
def submatrix (A : List (List α)) (rows cols : List Nat) : List (List α) :=
  rows.map fun i => cols.map fun j => get2 A i j

/-- Determinant of an `n × n` matrix by Laplace expansion along the first row. -/
def detN : Nat → List (List α) → α
  | 0, _ => 1
  | n + 1, A =>
    sum ((List.range (n + 1)).map fun j =>
      (if j % 2 = 0 then 1 else -1) * get (A.headD []) j * detN n (A.tail.map (·.eraseIdx j)))

def det (A : List (List α)) : α := detN A.length A

def ident (n : Nat) : List (List α) :=
  (List.range n).map fun i => (List.range n).map fun j => if i = j then 1 else 0

def matMul (A B : List (List α)) : List (List α) :=
  A.map fun row =>
    (List.range (B.headD []).length).map fun j =>
      sum ((List.range row.length).map fun k => get row k * get2 B k j)

/-- Inverse by the adjugate, `(A⁻¹)_ij = (−1)^(i+j) det(A without row j and
    column i) / det A`. -/
def inverse (A : List (List α)) : List (List α) :=
  let n := A.length
  let others (k : Nat) := (List.range n).filter (· ≠ k)
  (List.range n).map fun i => (List.range n).map fun j =>
    (if (i + j) % 2 = 0 then 1 else -1) * det (submatrix A (others j) (others i)) / det A

/-- The `k`-element subsets of `l`. -/
def subsets : List Nat → Nat → List (List Nat)
  | _, 0 => [[]]
  | [], _ + 1 => []
  | x :: xs, k + 1 => (subsets xs k).map (x :: ·) ++ subsets xs (k + 1)

/-- Coefficients of `det(I − zA)` in `z`: `(−1)^k` times the sum of the
    principal `k × k` minors of `A`. -/
def detIMinusZ (A : List (List α)) : Poly α :=
  (List.range (A.length + 1)).map fun k =>
    (if k % 2 = 0 then 1 else -1) *
      sum ((subsets (List.range A.length) k).map fun S => det (submatrix A S S))

/-- `|p(iy)|²` as a polynomial in `y`, for a polynomial `p` with real
    coefficients: the squares of the real and imaginary parts of `p(iy)`. -/
def absSqImag (p : Poly α) : Poly α :=
  let re : Poly α := (List.range p.length).map fun k =>
    if k % 4 = 0 then get p k else if k % 4 = 2 then -get p k else 0
  let im : Poly α := (List.range p.length).map fun k =>
    if k % 4 = 1 then get p k else if k % 4 = 3 then -get p k else 0
  Poly.add (Poly.mul re re) (Poly.mul im im)

end Matrices

/-! ## Padé approximants of the exponential -/

/-- Numerator of the `(s−1, s)` Padé approximant of `e^z`:
    `Σ_{m<s} (2s−1−m)! (s−1)! / ((2s−1)! m! (s−1−m)!) z^m`. -/
def padeNum (s : Nat) : Poly Rat :=
  (List.range s).map fun m =>
    ((fact (2 * s - 1 - m) * fact (s - 1) : Nat) : Rat) /
      ((fact (2 * s - 1) * fact m * fact (s - 1 - m) : Nat) : Rat)

/-- Denominator of the `(s−1, s)` Padé approximant of `e^z`:
    `Σ_{m≤s} (2s−1−m)! s! / ((2s−1)! m! (s−m)!) (−z)^m`. -/
def padeDen (s : Nat) : Poly Rat :=
  (List.range (s + 1)).map fun m =>
    (-1) ^ m * ((fact (2 * s - 1 - m) * fact s : Nat) : Rat) /
      ((fact (2 * s - 1) * fact m * fact (s - m) : Nat) : Rat)


end Firkode
