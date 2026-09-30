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

import FirkodeSpec.Collocation
import Mathlib.Analysis.Calculus.Deriv.Polynomial
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

/-!
# The antiderivative is the interval integral

`Collocation` integrates polynomials with the algebraic antiderivative
`antideriv`. Over the reals it is the interval integral
(`integral_eval_eq_antideriv`, by the fundamental theorem of calculus), so
`∫_0^1 p² ≥ 0` (`antideriv_sq_eval_one_nonneg`). With the exactness of the
quadrature up to degree `2s − 2`, the weights of a collocation method whose
node polynomial is orthogonal to the polynomials of degree `≤ s − 2` are
`b_j = ∫_0^1 l_j² ≥ 0` (`collocA_last_nonneg`).
-/

namespace Firkode.Colloc

open Polynomial Finset

/-- `∫_a^b p = [antideriv p]_a^b`. -/
theorem integral_eval_eq_antideriv (p : ℝ[X]) (a b : ℝ) :
    ∫ x in a..b, p.eval x = (antideriv p).eval b - (antideriv p).eval a := by
  refine intervalIntegral.integral_eq_sub_of_hasDerivAt (f := fun x => (antideriv p).eval x)
    (fun x _ => ?_) (p.differentiable.continuous.intervalIntegrable a b)
  have := (antideriv p).hasDerivAt x
  rwa [derivative_antideriv] at this

/-- `∫_0^1 p² ≥ 0`. -/
theorem antideriv_sq_eval_one_nonneg (p : ℝ[X]) : 0 ≤ (antideriv (p ^ 2)).eval 1 := by
  have h := integral_eval_eq_antideriv (p ^ 2) 0 1
  rw [eval_antideriv_zero, sub_zero] at h
  rw [← h]
  exact intervalIntegral.integral_nonneg zero_le_one fun x _ => by
    rw [eval_pow]; exact sq_nonneg _

/-- The weights `b_j = a_sj` of a collocation method with `c_s = 1` whose node
    polynomial is orthogonal to the polynomials of degree `≤ s − 2` are
    nonnegative: `b_j = Σ_i b_i l_j(c_i)² = ∫_0^1 l_j²`. -/
theorem collocA_last_nonneg {s : ℕ} {c : Fin s → ℝ} (hc : Function.Injective c) (last : Fin s)
    (hlast : c last = 1) (horth : ∀ j < s - 1, (antideriv (X ^ j * nodePoly c)).eval 1 = 0)
    (j : Fin s) : 0 ≤ collocA c last j := by
  have hdeg : (ell c j ^ 2).natDegree < 2 * s - 1 := by
    rw [natDegree_pow, natDegree_ell hc]
    have := Fin.pos j
    omega
  have hsum : ∑ i, collocA c last i * (ell c j ^ 2).eval (c i) = collocA c last j := by
    rw [Finset.sum_eq_single j]
    · simp [ell_eval_node hc]
    · intro i _ hij
      simp [ell_eval_node hc, hij]
    · simp
  rw [← hsum, collocA_B_poly hc last hlast horth hdeg]
  exact antideriv_sq_eval_one_nonneg _

end Firkode.Colloc
