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
# Certificates for the closed-form Radau IIA tables (s ≤ 3)

The tables of `firk_radau_closed_form`, taken from `Generated`, checked in exact
arithmetic over ℚ(√6).
-/

namespace Firkode.Tables

open Q6 Mirror Generated

/-- Stage counts with a closed-form table. -/
abbrev closedStages : List Nat := [1, 2, 3]

/-- The closed-form table for `s` stages, as built by `firk_radau_build`. -/
def closedForm (s : Nat) : Table Q6 := radauTable (closedFormC s) (closedFormA s)

/-- The exact inverse of the closed-form coefficient matrix. -/
def closedFormAinv (s : Nat) : List (List Q6) := inverse (closedFormA s)

/-- The quantities of `firk_derive_ld` for `γ0 = 1`. They are linear in `γ0`,
    so the error-estimate weights here are `e/γ0`. -/
def closedFormDerived (s : Nat) : Derived Q6 := derive (closedForm s) 1 (closedFormAinv s)

/-! ## Nodes and coefficients -/

/-- The closed-form nodes are the Radau IIA nodes: `firk_radau_poly`
    vanishes at `x = 2c − 1` for every node `c`. -/
theorem closedForm_nodes_radau :
    ∀ s ∈ closedStages, ∀ c ∈ closedFormC s, radauPolyAt s (2 * c - 1) = 0 := by
  decide +kernel

/-- The nodes satisfy `0 < c_1 < … < c_s = 1`. -/
theorem closedForm_nodes_increasing :
    ∀ s ∈ closedStages, 0 < get (closedFormC s) 0 ∧
      (∀ i < s - 1, get (closedFormC s) i < get (closedFormC s) (i + 1)) ∧
      get (closedFormC s) (s - 1) = 1 := by
  decide +kernel

/-- The closed-form coefficients equal the collocation coefficients that
    `firk_collocation_coeffs` computes on the closed-form nodes, so both paths
    of `firk_radau_build` define the same tables for s ≤ 3. -/
theorem closedForm_eq_collocation :
    ∀ s ∈ closedStages, collocationA (closedFormC s) = closedFormA s := by
  decide +kernel

/-- `FIRKodeTable_CheckOrder` finds exactly B(2s−1), C(s) and D(s−1), and so
    returns the stored order `q = 2s − 1`. That these simplifying assumptions
    give order `2s − 1` is Butcher's theorem (cited, not formalized). -/
theorem closedForm_checkOrder :
    ∀ s ∈ closedStages,
      checkOrder (closedForm s) maxCheckOrder = (2 * s - 1, s, s - 1, 2 * s - 1) := by
  decide +kernel

/-- `det A = Π c_i / s! = (s−1)!/(2s−1)!`; in particular `A` is invertible. -/
theorem closedForm_det :
    ∀ s ∈ closedStages, det (closedFormA s) = prod (closedFormC s) / (fact s : Q6) ∧
      det (closedFormA s) = ofRat ((fact (s - 1) : Rat) / (fact (2 * s - 1) : Rat)) := by
  decide +kernel

theorem closedForm_inverse :
    ∀ s ∈ closedStages, matMul (closedFormA s) (closedFormAinv s) = ident s ∧
      matMul (closedFormAinv s) (closedFormA s) = ident s := by
  decide +kernel

/-- `d = bᵀA⁻¹ = e_s`, so the step update is `y_{n+1} = y_n + Z_s`, and the
    stability function satisfies `R(∞) = 1 − bᵀA⁻¹𝟙 = 0`. -/
theorem closedForm_d :
    ∀ s ∈ closedStages,
      (closedFormDerived s).d = (List.range s).map (fun j => if j = s - 1 then 1 else 0) ∧
        1 - sum (closedFormDerived s).d = 0 := by
  decide +kernel

/-! ## Error estimate

`firk_derive_ld` sets `e = −γ0 l(0)ᵀA⁻¹`. With the stage increments
`z = h (A ⊗ I) F`, the unfiltered estimate is
`γ0 h f(t_n, y_n) + Σ_i e_i z_i = γ0 h (f(t_n, y_n) − Σ_i l_i(0) F_i)`: `γ0 h`
times the difference between `f` at `t_n` and the interpolant of the stage
derivatives evaluated there. -/

/-- `l_i(0)`, the constant coefficient computed by `firk_lagrange_coeffs`. -/
def l0 (c : List Q6) (i : Nat) : Q6 := Poly.coeff (lagrangeCoeffs c i) 0

/-- `Σ_i l_i(0) c_i^k`: the interpolant of `τ^k` evaluated at `τ = 0`. -/
def moment (c : List Q6) (k : Nat) : Q6 :=
  sum ((List.range c.length).map fun i => l0 c i * npow (get c i) k)

/-- `eᵀA = −γ0 l(0)ᵀ`. -/
theorem closedForm_eA :
    ∀ s ∈ closedStages,
      (List.range s).map (fun k =>
          sum ((List.range s).map fun i => get (closedFormDerived s).e i * get2 (closedFormA s) i k))
        = (List.range s).map fun k => -l0 (closedFormC s) k := by
  decide +kernel

/-- The bracket `g(0) − Σ_i l_i(0) g(c_i)` vanishes for `g = τ^k`, `k < s`, and
    equals `(−1)^s Π c_i ≠ 0` for `g = τ^s`: the estimate has order `p = s`. -/
theorem closedForm_moments :
    ∀ s ∈ closedStages, (∀ k < s, moment (closedFormC s) k = if k = 0 then 1 else 0) ∧
      moment (closedFormC s) s = ofRat ((-1) ^ (s + 1)) * prod (closedFormC s) := by
  decide +kernel

/-- `e_j = −γ0 l_j(0)/c_j`, a closed form that needs no matrix inverse. -/
theorem closedForm_e :
    ∀ s ∈ closedStages, (closedFormDerived s).e =
      (List.range s).map fun j => -l0 (closedFormC s) j / get (closedFormC s) j := by
  decide +kernel

/-- `e_s = (−1)^s γ0/s`. -/
theorem closedForm_e_last :
    ∀ s ∈ closedStages,
      get (closedFormDerived s).e (s - 1) = ofRat ((-1) ^ s / (s : Rat)) := by
  decide +kernel

/-- For s = 3, `e/γ0` is RADAU5's `(DD1, DD2, DD3)`. -/
theorem radau3_e_radau5 :
    (closedFormDerived 3).e = [-(13 + 7 * r6) / 3, (-13 + 7 * r6) / 3, -1 / 3] := by
  decide +kernel

/-! ## Dense output -/

/-- `L_i(θ) = Σ_m θ^(m+1) P[m][i]`, the polynomial multiplying `Z_i` in the
    dense output `u(t_n + θh) = y_n + Σ_i L_i(θ) Z_i`. `L_i(0) = 0` holds by
    construction. -/
def denseL (D : Derived Q6) (s i : Nat) (θ : Q6) : Q6 :=
  sum ((List.range s).map fun m => npow θ (m + 1) * get2 D.P m i)

/-- `L_i(c_j) = δ_ij` and `L_i(1) = d_i`: the dense output passes through the
    stage values `Y_j = y_n + Z_j` and ends at `y_{n+1}`. -/
theorem closedForm_dense :
    ∀ s ∈ closedStages,
      (∀ i < s, ∀ j < s,
          denseL (closedFormDerived s) s i (get (closedFormC s) j) = if i = j then 1 else 0) ∧
        ∀ i < s, denseL (closedFormDerived s) s i 1 = get (closedFormDerived s).d i := by
  decide +kernel

/-- `L′(c) A = I`, with `L_j′(c_i)` given by the `k = 1` formula of
    `FIRKodeGetDky`: the dense output satisfies the collocation equations
    `h u′(t_n + c_i h) = Σ_j L_j′(c_i) Z_j = h F_i`. -/
theorem closedForm_dense_deriv :
    ∀ s ∈ closedStages,
      matMul ((List.range s).map fun i =>
          dkyWeights (closedFormDerived s) s 1 (get (closedFormC s) i)) (closedFormA s)
        = ident s := by
  decide +kernel

/-! ## Stability function

`R(z) = det(I − zA + z𝟙bᵀ) / det(I − zA)`. -/

/-- The numerator `det(I − z(A − 𝟙bᵀ))`. -/
def stabNum (T : Table Q6) : Poly Q6 :=
  detIMinusZ (T.A.map fun row => (List.range row.length).map fun j => get row j - get T.b j)

/-- `R` is the `(s−1, s)` Padé approximant of `e^z`. Its numerator has degree
    `s − 1` and its denominator degree `s`, so `R(∞) = 0`. -/
theorem closedForm_pade :
    ∀ s ∈ closedStages, Poly.trim (detIMinusZ (closedFormA s)) = (padeDen s).map ofRat ∧
      Poly.trim (stabNum (closedForm s)) = (padeNum s).map ofRat := by
  decide +kernel

/-- `E(y) = |Q(iy)|² − |P(iy)|² = (det A)² y^(2s) ≥ 0`, so `|R(iy)| ≤ 1` on the
    imaginary axis. With the poles of `R` in the right half-plane this gives
    A-stability (criterion not formalized; `spec` proves A-stability directly,
    `Colloc.radau_A_stable`). -/
theorem closedForm_Epoly :
    ∀ s ∈ closedStages,
      Poly.trim (Poly.sub (absSqImag (detIMinusZ (closedFormA s)))
          (absSqImag (stabNum (closedForm s))))
        = List.replicate (2 * s) 0 ++ [det (closedFormA s) * det (closedFormA s)] := by
  decide +kernel

end Firkode.Tables
