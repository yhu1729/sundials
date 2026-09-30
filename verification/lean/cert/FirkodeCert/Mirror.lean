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

/-!
# Lean transcriptions of the FIRKODE table routines

Each definition transcribes the C routine named in its docstring into exact
arithmetic over any type with the field operations: `Q6` and `Rat` here, the
reals in `spec`. The transcriptions follow the C loops (`foldl` for a `for`
loop) so that `spec` can prove them correct; floating-point rounding is not
modeled. `gen_firkode_tables.py --check` fails when one of these C routines
changes, which prompts a review of its transcription here.
-/

namespace Firkode.Mirror

variable {α : Type} [Add α] [Sub α] [Mul α] [Div α] [Neg α] [OfNat α 0] [OfNat α 1]
  [NatCast α]

/-! ## Nodes -/

/-- `firk_radau_poly`: `P_s(x) − P_{s−1}(x)`, with the Legendre polynomials
    from the recurrence `(k+1) P_{k+1} = (2k+1) x P_k − k P_{k−1}`. -/
def radauPolyAt (s : Nat) (x : α) : α :=
  if s = 1 then x - 1
  else
    let p := (List.range' 1 (s - 1)).foldl
      (fun (p : α × α) (k : Nat) =>
        (p.2, (((2 * k + 1 : Nat) : α) * x * p.2 - (k : α) * p.1) / ((k + 1 : Nat) : α)))
      (1, x)
    p.2 - p.1

/-- `firk_radau_poly` run on polynomials: the coefficients of
    `P_s(x) − P_{s−1}(x)`. -/
def radauPoly (s : Nat) : Poly α :=
  if s = 1 then Poly.sub Poly.X (Poly.C 1)
  else
    let p := (List.range' 1 (s - 1)).foldl
      (fun (p : Poly α × Poly α) (k : Nat) =>
        (p.2, Poly.smul (1 / ((k + 1 : Nat) : α))
          (Poly.sub (Poly.smul ((2 * k + 1 : Nat) : α) (Poly.mulX p.2))
            (Poly.smul (k : α) p.1))))
      (Poly.C 1, Poly.X)
    Poly.sub p.2 p.1

/-- Grid point `x_j = −1 + 2j/n` of `firk_radau_nodes`. -/
def gridPoint (n j : Nat) : Rat := -1 + 2 * (j : Rat) / (n : Rat)

/-- The values of `firk_radau_poly` at the grid points `x_0, …, x_{n−1}`,
    `n = perStage · s`, of `firk_radau_nodes`. -/
def gridValues (perStage s : Nat) : List Rat :=
  let n := perStage * s
  (List.range n).map fun j => radauPolyAt s (gridPoint n j)

/-- The bracketing scan of `firk_radau_nodes`. Returns `(roots, zeros)`:
    `roots` counts the grid intervals `[x_{j−1}, x_j]`, `j = 1, …, n − 1`,
    whose left end is a zero of `radauPolyAt s` or across which it changes
    sign, and `zeros` counts the zeros among those left ends. Like the C loop,
    it evaluates each grid value once and carries it to the next interval. The
    C routine fails unless `roots = s − 1`; it skips the scan for `s = 1`. -/
def radauGridScan (perStage s : Nat) : Nat × Nat :=
  if s = 1 then (0, 0)
  else
    let v := gridValues perStage s
    (v.zip v.tail).foldl
      (fun (acc : Nat × Nat) (f : Rat × Rat) =>
        if f.1 = 0 then (acc.1 + 1, acc.2 + 1)
        else if f.1 * f.2 < 0 then (acc.1 + 1, acc.2)
        else acc)
      (0, 0)

/-! ## Collocation coefficients -/

/-- `firk_lagrange_coeffs`: coefficients of
    `l_j(τ) = Π_{k≠j} (τ − c_k)/(c_j − c_k)`, built one factor at a time:
    multiply by `τ − c_k`, then scale by `1/(c_j − c_k)`. -/
def lagrangeCoeffs (c : List α) (j : Nat) : Poly α :=
  (List.range c.length).foldl
    (fun (r : Poly α) (k : Nat) =>
      if k = j then r
      else Poly.smul (1 / (get c j - get c k)) (Poly.sub (Poly.mulX r) (Poly.smul (get c k) r)))
    [1]

/-- `firk_collocation_coeffs`: `a_ij = Σ_m r_m c_i^(m+1)/(m+1)`, the integral
    of `l_j` over `[0, c_i]`, where `r = lagrangeCoeffs c j`. Row major. -/
def collocationA (c : List α) : List (List α) :=
  (List.range c.length).map fun i =>
    (List.range c.length).map fun j =>
      sum ((List.range c.length).map fun m =>
        Poly.coeff (lagrangeCoeffs c j) m * npow (get c i) (m + 1) / ((m + 1 : Nat) : α))

/-- A Butcher table; `A` is row major. -/
structure Table (α : Type) where
  c : List α
  A : List (List α)
  b : List α

/-- `firk_radau_build`: the table from the nodes and coefficients, with the
    weights `b_i = A[s−1][i]` (stiffly accurate). The C routine also stores the
    orders `q = 2s − 1` and `p = s`. -/
def radauTable (c : List α) (A : List (List α)) : Table α :=
  { c, A, b := A.getD (c.length - 1) [] }

/-- Quantities derived by `firk_derive_ld`. -/
structure Derived (α : Type) where
  /-- `d = bᵀA⁻¹` -/
  d : List α
  /-- error-estimate weights `e = −γ0 l(0)ᵀA⁻¹` -/
  e : List α
  /-- dense-output coefficients: `P[k][i]` is coefficient `k` of `l_i`
      divided by `c_i`, i.e. coefficient `k + 1` of `θ l_i(θ)/c_i` -/
  P : List (List α)

/-- `firk_derive_ld`, given `Ainv = A⁻¹` (which the C routine computes by LU
    factorization); `l_i(0)` is the constant coefficient of
    `lagrangeCoeffs c i`. -/
def derive (T : Table α) (gamma0 : α) (Ainv : List (List α)) : Derived α :=
  let s := T.c.length
  let l0 := (List.range s).map fun i => Poly.coeff (lagrangeCoeffs T.c i) 0
  { d := (List.range s).map fun j =>
      sum ((List.range s).map fun i => get T.b i * get2 Ainv i j)
    e := (List.range s).map fun j =>
      sum ((List.range s).map fun i => get l0 i * get2 Ainv i j) * -gamma0
    P := (List.range s).map fun k =>
      (List.range s).map fun i => Poly.coeff (lagrangeCoeffs T.c i) k / get T.c i }

/-! ## Order check -/

section Order
variable [DecidableEq α]

/-- `B(k)`: `Σ_i b_i c_i^(k−1) = 1/k`. -/
def condB (T : Table α) (k : Nat) : Bool :=
  sum ((List.range T.c.length).map fun i => get T.b i * npow (get T.c i) (k - 1))
    == 1 / (k : α)

/-- `C(k)`: `Σ_j a_ij c_j^(k−1) = c_i^k / k` for all `i`. -/
def condC (T : Table α) (k : Nat) : Bool :=
  (List.range T.c.length).all fun i =>
    sum ((List.range T.c.length).map fun j => get2 T.A i j * npow (get T.c j) (k - 1))
      == npow (get T.c i) k / (k : α)

/-- `D(k)`: `Σ_i b_i c_i^(k−1) a_ij = b_j (1 − c_j^k) / k` for all `j`. -/
def condD (T : Table α) (k : Nat) : Bool :=
  (List.range T.c.length).all fun j =>
    sum ((List.range T.c.length).map fun i =>
        get T.b i * npow (get T.c i) (k - 1) * get2 T.A i j)
      == get T.b j * (1 - npow (get T.c j) k) / (k : α)

/-- Largest `m ≤ kmax` such that `cond 1, …, cond m` all hold: a loop of
    `FIRKodeTable_CheckOrder`, with its tolerance replaced by exact equality. -/
def lastTrue (cond : Nat → Bool) (kmax : Nat) : Nat := go kmax 1
where
  go : Nat → Nat → Nat
    | 0, k => k - 1
    | fuel + 1, k => if cond k then go fuel (k + 1) else k - 1

/-- `FIRKodeTable_CheckOrder` in exact arithmetic, checking up to order
    `kmax` (`FIRK_MAX_CHECK_ORDER`): returns `(qB, η, ζ, q)` with
    `q = min(qB, 2η + 2, η + ζ + 1)`. -/
def checkOrder (T : Table α) (kmax : Nat) : Nat × Nat × Nat × Nat :=
  let qB := lastTrue (condB T) kmax
  let eta := lastTrue (condC T) kmax
  let zeta := lastTrue (condD T) kmax
  (qB, eta, zeta, min qB (min (2 * eta + 2) (eta + zeta + 1)))

end Order

/-! ## Dense output and prediction -/

/-- The weights of `Zprev_j` in `FIRKodeGetDky` for `k = 0`, and in
    `firkPredict`: `−d_j + Σ_{m<s} θ^(m+1) P[m][j]`. -/
def denseWeights (D : Derived α) (s : Nat) (θ : α) : List α :=
  (List.range s).map fun j =>
    (List.range s).foldl (fun acc m => acc + npow θ (m + 1) * get2 D.P m j) (-(get D.d j))

/-- `θ = 1 + (t − tcur)/hold` in `FIRKodeGetDky`. -/
def dkyTheta (t tcur hold : α) : α := 1 + (t - tcur) / hold

/-- `FIRKodeGetDky` for `k ≥ 1`: the weights
    `Σ_{m=k−1}^{s−1} (m+1)!/(m+1−k)! θ^(m+1−k) P[m][j]` of `Zprev_j`; the C
    routine then scales the result by `hold^(−k)`. -/
def dkyWeights (D : Derived α) (s k : Nat) (θ : α) : List α :=
  (List.range s).map fun j =>
    (List.range' (k - 1) (s - (k - 1))).foldl
      (fun acc m =>
        let fac := (List.range k).foldl (fun f t => f * ((m + 1 - t : Nat) : α)) 1
        acc + fac * npow θ (m + 1 - k) * get2 D.P m j)
      0

/-- `firkPredict`: the weights of `Zprev_j` in `Zpred_i`, evaluated at
    `θ_i = 1 + c_i h / hold`. -/
def predictWeights (T : Table α) (D : Derived α) (h hold : α) : List (List α) :=
  (List.range T.c.length).map fun i => denseWeights D T.c.length (1 + get T.c i * h / hold)

/-! ## User options -/

/-- `FIRKodeSetOrder`: the number of stages selected for order `ord`. Lean's
    `/` on `Int` rounds down and C's rounds toward zero; they agree here
    because the division only happens for `ord > 0`. -/
def setOrderStages (defaultStages : Nat) (ord : Int) : Int :=
  if ord ≤ 0 then defaultStages else ord / 2 + 1

end Firkode.Mirror
