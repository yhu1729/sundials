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

import FirkodeCert.NodePoly
import FirkodeCert.Tables

/-!
# Certificates for the `firk_gamma0` literals

`firk_gamma0[s]` is the real shift `γ0` of the block preconditioner
`M − hγ0J` and of the filtered error estimate. It is defined as `1/U1`, with
`U1` the real eigenvalue of `A⁻¹`, for odd `s`, and as `1/Re λ`, with `λ` the
smallest-modulus eigenvalue of `A⁻¹`, for even `s`. The eigenvalues of `A`
are the roots of `det(μI − A) = μ^s Q_s(1/μ)`, where `Q_s(z) = det(I − zA)`
is the Padé denominator. That identity is proved from the tables for s ≤ 3
(`Tables.closedForm_pade`), and in `spec` for every s ≤ 9
(`Charpoly.radau_charpolyRev`). For even s ≥ 4 the certificates show that
there is no real eigenvalue; the values of γ0 are only checked numerically
(roadmap R6).
-/

namespace Firkode.Gamma0

open Generated NodePoly

/-- `det(μI − A) = μ^s Q_s(1/μ)`, the reversed Padé denominator. -/
def charPoly (s : Nat) : Poly Rat := (padeDen s).reverse

/-- The literal `firk_gamma0[s]`. -/
def literal (s : Nat) : Rat := get gamma0Literal s

/-- The unit in the last printed digit of `firk_gamma0[s]`. -/
def ulp (s : Nat) : Rat := 1 / 10 ^ get gamma0Digits s

/-- Number of sign changes in the sequence of nonzero coefficients, read from
    the highest power down as in Mathlib's `Polynomial.signVariations`. -/
def signVariations (p : Poly Rat) : Nat :=
  let nz := p.reverse.filter (· ≠ 0)
  ((nz.zip nz.tail).filter fun ab => ab.1 * ab.2 < 0).length

/-- `p(x + a)`: its positive roots are the roots of `p` above `a`. -/
def shift (p : Poly Rat) (a : Rat) : Poly Rat := Poly.comp p [a, 1]

/-- `(1 + x)^n p((a + b x)/(1 + x))` with `n = deg p`: its positive roots
    correspond to the roots of `p` in `(a, b)`. -/
def mobius (p : Poly Rat) (a b : Rat) : Poly Rat :=
  (List.range p.length).foldl
    (fun acc k => Poly.add acc
      (Poly.smul (get p k) (Poly.mul (Poly.pow [a, b] k) (Poly.pow [1, 1] (p.length - 1 - k)))))
    []

/-- Certificate that `p` has exactly one real root, and that it lies in
    `(lo, hi)`, for `0 < lo < hi < 2`. By Descartes' rule of signs (proved in
    Mathlib, applied to the transformed polynomials in `spec`) the numbers of
    roots in `(−∞, 0)`, `(0, lo)`, `(lo, hi)`, `(hi, 2)` and `(2, ∞)` are at
    most the sign variations `0, 0, 1, 0, 0`; `p` is nonzero at the break
    points, and its sign change on `[lo, hi]` gives the root
    (`RootCount.uniqueRootIn_spec`). -/
def uniqueRootIn (p : Poly Rat) (lo hi : Rat) : Bool :=
  decide (0 < lo ∧ lo < hi ∧ hi < 2) &&
    decide (Poly.eval p 0 ≠ 0 ∧ Poly.eval p lo ≠ 0 ∧ Poly.eval p hi ≠ 0 ∧ Poly.eval p 2 ≠ 0) &&
    decide (Poly.eval p lo * Poly.eval p hi < 0) &&
    signVariations (reflect p) == 0 &&
    signVariations (mobius p 0 lo) == 0 &&
    signVariations (mobius p lo hi) == 1 &&
    signVariations (mobius p hi 2) == 0 &&
    signVariations (shift p 2) == 0

/-- `p` has no root in `(a, b)`: the Möbius transform has no sign variation,
    or `p` is nonzero at the midpoint and both halves pass, with at most
    `fuel` bisections. -/
def noRootIn (p : Poly Rat) : Nat → Rat → Rat → Bool
  | 0, a, b => signVariations (mobius p a b) == 0
  | fuel + 1, a, b =>
    signVariations (mobius p a b) == 0 ||
      (decide (Poly.eval p ((a + b) / 2) ≠ 0) && noRootIn p fuel a ((a + b) / 2) &&
        noRootIn p fuel ((a + b) / 2) b)

/-- Certificate that `p` has no real root: `p` is nonzero at 0 and 2, has no
    root in `(0, 2)` by `noRootIn`, and `p(−x)` and `p(x + 2)` have no sign
    variation. -/
def noRealRoot (p : Poly Rat) : Bool :=
  decide (Poly.eval p 0 ≠ 0 ∧ Poly.eval p 2 ≠ 0) && signVariations (reflect p) == 0 &&
    noRootIn p 6 0 2 && signVariations (shift p 2) == 0

/-- s = 1: `A = (1)`, so `γ0 = 1`. -/
theorem gamma0_one : literal 1 = 1 ∧ Poly.eval (charPoly 1) (literal 1) = 0 := by
  decide +kernel

/-- s = 2: `det(μI − A) = μ² − t1 μ + t2` with `t1² < 4 t2`, a complex pair
    `μ` of equal modulus. The eigenvalues `λ = 1/μ` of `A⁻¹` have
    `Re λ = (t1/2)/t2`, so `γ0 = 1/Re λ = 2 t2/t1 = 1/2`. -/
theorem gamma0_two :
    let p := charPoly 2
    get p 1 ^ 2 < 4 * get p 0 * get p 2 ∧ literal 2 = 2 * get p 0 / -get p 1 := by
  decide +kernel

/-- Odd s: `det(μI − A)` has exactly one real root, the real eigenvalue of `A`
    (so `γ0 = 1/U1`), and `firk_gamma0[s]` is within half a unit in its last
    printed digit of it, i.e. correctly rounded. -/
theorem gamma0_odd :
    ∀ s ∈ [3, 5, 7, 9],
      uniqueRootIn (charPoly s) (literal s - ulp s / 2) (literal s + ulp s / 2) = true := by
  decide +kernel

/-- Even s: `det(μI − A)` has no real root, so neither `A` nor `A⁻¹` has a
    real eigenvalue, and `γ0 = 1/Re λ` is taken from a complex pair. -/
theorem gamma0_even_no_real : ∀ s ∈ [2, 4, 6, 8], noRealRoot (charPoly s) = true := by
  decide +kernel

/-- The sign change part of `gamma0_odd`, stated directly: `p_s` changes sign
    across `[L − u/2, L + u/2]`. -/
theorem gamma0_odd_bracket :
    ∀ s ∈ [3, 5, 7, 9], literal s - ulp s / 2 < literal s + ulp s / 2 ∧
      Poly.eval (charPoly s) (literal s - ulp s / 2) *
        Poly.eval (charPoly s) (literal s + ulp s / 2) < 0 := by
  decide +kernel

end Firkode.Gamma0
