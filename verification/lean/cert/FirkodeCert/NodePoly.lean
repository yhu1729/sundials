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

import FirkodeCert.Generated
import FirkodeCert.Mirror

/-!
# Certificates from the Radau node polynomial (s ≤ 9)

For s ≥ 4 FIRKODE computes the Radau IIA nodes numerically. The nodes are
irrational, but their symmetric functions are rational, so exact certificates
are still possible through the monic node polynomial `M_s(c) = Π_i (c − c_i)`,
obtained here from `firk_radau_poly` at `x = 2c − 1`. Each docstring says what
the certificate implies once combined with a theorem proved in `spec` or cited.
-/

namespace Firkode.NodePoly

open Mirror Generated

/-- The stage counts FIRKODE supports, `1, …, FIRK_MAX_STAGES`. -/
abbrev stages : List Nat := List.range' 1 maxStages

/-- The monic node polynomial `M_s(c) = Π_i (c − c_i)`: `firk_radau_poly` at
    `x = 2c − 1`, scaled to leading coefficient 1. -/
def nodePoly (s : Nat) : Poly Rat :=
  let p := Poly.trim (Poly.comp (radauPoly s) [-1, 2])
  Poly.smul (1 / p.getLastD 1) p

/-- `∫₀¹ p`. -/
def integral01 (p : Poly Rat) : Rat :=
  sum ((List.range p.length).map fun k => get p k / ((k + 1 : Nat) : Rat))

/-- The `k`-th derivative. -/
def derivN (p : Poly Rat) : Nat → Poly Rat
  | 0 => p
  | k + 1 => Poly.deriv (derivN p k)

/-- The polynomial coefficients from `radauPoly` agree with the evaluation
    `radauPolyAt` (both transcribe `firk_radau_poly`) at `s + 1` points, hence
    everywhere, as both have degree at most `s`. (`spec` proves the identity
    at every point directly: `Nodes.radauPolyAt_eq_eval`.) -/
theorem radauPoly_eval :
    ∀ s ∈ stages, ∀ j < s + 1, Poly.eval (radauPoly s) (j : Rat) = radauPolyAt s (j : Rat) := by
  decide +kernel

/-- `M_s` has degree `s` and vanishes at `c = 1`, the last node. -/
theorem nodePoly_degree_endpoint :
    ∀ s ∈ stages, (nodePoly s).length = s + 1 ∧ Poly.eval (nodePoly s) 1 = 0 := by
  decide +kernel

/-- `M_s` is monic: its coefficient of `c^s` is 1. -/
theorem nodePoly_monic : ∀ s ∈ stages, Poly.coeff (nodePoly s) s = 1 := by
  decide +kernel

/-- Binomial coefficient (core Lean has no `Nat.choose`). -/
def choose : Nat → Nat → Nat
  | _, 0 => 1
  | 0, _ + 1 => 0
  | n + 1, k + 1 => choose n k + choose n (k + 1)

/-- Coefficient `k` of Mathlib's shifted Legendre polynomial
    `P̃_n(x) = P_n(1 − 2x)`: `(−1)^k C(n, k) C(n + k, n)`. -/
def shLegendreCoeff (n k : Nat) : Int := (-1) ^ k * (choose n k : Int) * (choose (n + k) n : Int)

/-- `firk_radau_poly` at `x = 2c − 1` is `(−1)^s (P̃_s + P̃_{s−1})`: the
    recurrence of the C routine produces the Legendre polynomials. -/
theorem radauPoly_shiftedLegendre :
    ∀ s ∈ stages, Poly.trim (Poly.comp (radauPoly s) [-1, 2]) =
      (List.range (s + 1)).map fun k =>
        (((-1) ^ s * (shLegendreCoeff s k + shLegendreCoeff (s - 1) k) : Int) : Rat) := by
  decide +kernel

/-! ## Quadrature order -/

/-- `∫₀¹ c^j M_s(c) dc = 0` for `j ≤ s − 2`. Because `spec` proves that the
    collocation weights `b_i = ∫₀¹ l_i` satisfy
    `Σ_i b_i p(c_i) = ∫₀¹ (p mod M_s)`, this gives B(2s−1): the quadrature is
    exact up to degree `2s − 2`. -/
theorem nodePoly_orthogonal :
    ∀ s ∈ stages, ∀ j < s - 1,
      integral01 (Poly.mul (Poly.pow Poly.X j) (nodePoly s)) = 0 := by
  decide +kernel

/-- `∫₀¹ c^(s−1) M_s = −(s!(s−1)!)²/((2s)!(2s−1)!)`. By the same argument
    the B(2s) residual `Σ_i b_i c_i^(2s−1) − 1/(2s)` is the negative of this
    integral, about 9.4e−11 for s = 9: the smallest residual of a violated
    condition quoted in the tolerance comment of `FIRKodeTable_CheckOrder`. -/
theorem nodePoly_B2s_residual :
    ∀ s ∈ stages,
      integral01 (Poly.mul (Poly.pow Poly.X (s - 1)) (nodePoly s)) =
        -(((fact s * fact (s - 1)) ^ 2 : Nat) : Rat) /
          ((fact (2 * s) * fact (2 * s - 1) : Nat) : Rat) := by
  decide +kernel

/-! ## Stability function -/

/-- Numerator of Nørsett's formula for a collocation method with monic node
    polynomial `M` of degree `s`: `N(z) = Σ_j M^(s−j)(1) z^j`. -/
def norsettNum (M : Poly Rat) (s : Nat) : Poly Rat :=
  (List.range (s + 1)).map fun j => Poly.eval (derivN M (s - j)) 1

/-- Denominator of Nørsett's formula: `D(z) = Σ_j M^(s−j)(0) z^j`. -/
def norsettDen (M : Poly Rat) (s : Nat) : Poly Rat :=
  (List.range (s + 1)).map fun j => Poly.eval (derivN M (s - j)) 0

/-- Nørsett's formula on the Radau IIA node polynomial gives the `(s−1, s)`
    Padé approximant: `N = s! · padeNum` (degree `s − 1`) and
    `D = s! · padeDen`. With `spec`'s proof that the stability function of a
    collocation method is `N/D`, this covers every table FIRKODE builds. -/
theorem nodePoly_norsett :
    ∀ s ∈ stages,
      Poly.trim (norsettNum (nodePoly s) s) = Poly.smul (fact s : Rat) (padeNum s) ∧
        norsettDen (nodePoly s) s = Poly.smul (fact s : Rat) (padeDen s) := by
  decide +kernel

/-- `get (padeDen s) 0 = 1` and the coefficient of `z^s` is nonzero, so
    `padeDen s` has degree `s` and value 1 at 0. -/
theorem padeDen_ends : ∀ s ∈ stages, get (padeDen s) 0 = 1 ∧ get (padeDen s) s ≠ 0 := by
  decide +kernel

/-- Bézout cofactors `(u, v)` with `u · padeNum s + v · padeDen s = 1`. -/
def padeBezout (s : Nat) : Poly Rat × Poly Rat :=
  let g := Poly.xgcd (padeNum s) (padeDen s)
  let c := (Poly.trim g.1).headD 1
  (Poly.smul (1 / c) g.2.1, Poly.smul (1 / c) g.2.2)

/-- The numerator and the denominator of the `(s−1, s)` Padé approximant are
    coprime. With Nørsett's formula this identifies `det(I − zA)` with the
    Padé denominator (`spec`, `Charpoly.radau_charpolyRev`). -/
theorem pade_coprime :
    ∀ s ∈ stages, Poly.trim (Poly.add (Poly.mul (padeBezout s).1 (padeNum s))
      (Poly.mul (padeBezout s).2 (padeDen s))) = [1] := by
  decide +kernel

/-- `E(y) = |Q(iy)|² − |P(iy)|² = ((s−1)!/(2s−1)!)² y^(2s)` for the
    `(s−1, s)` Padé approximant `P/Q`, so `|R(iy)| ≤ 1`. -/
theorem pade_Epoly :
    ∀ s ∈ stages,
      Poly.trim (Poly.sub (absSqImag (padeDen s)) (absSqImag (padeNum s))) =
        List.replicate (2 * s) 0 ++ [((fact (s - 1) : Rat) / (fact (2 * s - 1) : Rat)) ^ 2] := by
  decide +kernel

/-- First column of the Routh array of `p` (ascending coefficients, highest
    coefficient first in the array). -/
def routhFirstColumn (p : Poly Rat) : List Rat :=
  let d := p.reverse
  let evens := (List.range ((d.length + 1) / 2)).map fun k => get d (2 * k)
  let odds := (List.range (d.length / 2)).map fun k => get d (2 * k + 1)
  go (d.length - 1) evens odds [get evens 0]
where
  go : Nat → List Rat → List Rat → List Rat → List Rat
    | 0, _, _, acc => acc.reverse
    | fuel + 1, r0, r1, acc =>
      let a := get r1 0
      let r2 := (List.range r0.length).map fun j => (a * get r0 (j + 1) - get r0 0 * get r1 (j + 1)) / a
      go fuel r1 r2 (a :: acc)

/-- `p(−z)`. -/
def reflect (p : Poly Rat) : Poly Rat :=
  (List.range p.length).map fun k => if k % 2 = 0 then get p k else -get p k

/-- The Routh array of `Q_s(−z)` has a positive first column, so by the
    Routh–Hurwitz criterion (not formalized) the poles of `R` lie in the open
    right half-plane. Together with `pade_Epoly` this is the classical
    A-stability argument. `spec` proves A-stability without it, by an energy
    identity (`Colloc.radau_A_stable`); this certificate is an independent
    check. -/
theorem pade_poles_right :
    ∀ s ∈ stages, ∀ x ∈ routhFirstColumn (reflect (padeDen s)), 0 < x := by
  decide +kernel

/-! ## Error-estimate weight `e_s` -/

/-- `Π_{k<s} c_k/(1 − c_k) = (−1)^s M_s(0)/M_s′(1) = 1/s` over the interior
    nodes. Hence `l_s(0) = (−1)^(s−1)/s`, and with `e_j = −γ0 l_j(0)/c_j`
    (proved in `spec`) `e_s = (−1)^s γ0/s`. -/
theorem nodePoly_interior_product :
    ∀ s ∈ stages,
      (-1) ^ s * Poly.eval (nodePoly s) 0 / Poly.eval (Poly.deriv (nodePoly s)) 1 =
        1 / (s : Rat) := by
  decide +kernel

/-! ## Root bracketing in `firk_radau_nodes` -/

/-- The bracketing scan of `firk_radau_nodes` finds `s − 1` sign changes and
    no grid zeros, so in exact arithmetic it never takes its failure exit and
    each interior node lies in its own grid interval, where bisection converges
    to it (`spec` derives the nodes from these brackets:
    `Nodes.radauNodes_bracket`). At every grid point
    `|P_s − P_{s−1}| ≥ 8 · 10⁻⁵`, far above the error of evaluating it in
    `long double`, so the signs the C code computes are the exact ones
    (rounding bound cited). Both parts are one theorem so that the kernel
    evaluates the grid values once. -/
theorem nodes_grid :
    ∀ s ∈ stages, radauGridScan nodeGridPerStage s = (s - 1, 0) ∧
      (s ≠ 1 → (gridValues nodeGridPerStage s).all fun f =>
        decide ((8 / 100000 : Rat) ≤ f ∨ f ≤ -(8 / 100000 : Rat))) := by
  decide +kernel

/-! ## Non-vacuity -/

/-- Legendre polynomial `P_s` (not a transcription). -/
def legendrePoly (s : Nat) : Poly Rat :=
  ((List.range' 1 (s - 1)).foldl
    (fun (p : Poly Rat × Poly Rat) (k : Nat) =>
      (p.2, Poly.smul (1 / ((k + 1 : Nat) : Rat))
        (Poly.sub (Poly.smul ((2 * k + 1 : Nat) : Rat) (Poly.mulX p.2)) (Poly.smul (k : Rat) p.1))))
    (Poly.C 1, if s = 0 then Poly.C 1 else Poly.X)).2

/-- The monic Gauss node polynomial, `P_s(2c − 1)` scaled. -/
def gaussNodePoly (s : Nat) : Poly Rat :=
  let p := Poly.trim (Poly.comp (legendrePoly s) [-1, 2])
  Poly.smul (1 / p.getLastD 1) p

/-- On the Gauss nodes Nørsett's formula does not give the `(s−1, s)` Padé
    approximant, so `nodePoly_norsett` is not vacuous. -/
theorem gauss_not_pade :
    ∀ s ∈ [2, 3, 4], norsettDen (gaussNodePoly s) s ≠ Poly.smul (fact s : Rat) (padeDen s) := by
  decide +kernel

end Firkode.NodePoly
