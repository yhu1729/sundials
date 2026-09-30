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
# User-option logic
-/

namespace Firkode.Api

open Mirror Generated

/-- `FIRKodeSetOrder` with `ord ≤ 0` selects `FIRK_DEFAULT_STAGES` stages. -/
theorem setOrder_default (ord : Int) (h : ord ≤ 0) :
    setOrderStages defaultStages ord = defaultStages := by
  simp [setOrderStages, h]

/-- For `ord > 0`, `FIRKodeSetOrder` selects the smallest stage count `s`
    whose order `2s − 1` is at least `ord`. -/
theorem setOrder_minimal (ord : Int) (h : 0 < ord) :
    ord ≤ 2 * setOrderStages defaultStages ord - 1 ∧
      2 * (setOrderStages defaultStages ord - 1) - 1 < ord := by
  simp only [setOrderStages, if_neg (show ¬ord ≤ 0 by omega)]
  omega

/-- For `ord > 0` the selected stage count is at most `FIRK_MAX_STAGES` exactly
    when `ord ≤ 2 · FIRK_MAX_STAGES − 1 = 17`. -/
theorem setOrder_supported (ord : Int) (h : 0 < ord) :
    setOrderStages defaultStages ord ≤ maxStages ↔ ord ≤ 2 * maxStages - 1 := by
  simp only [setOrderStages, if_neg (show ¬ord ≤ 0 by omega), maxStages]
  omega

/-- The stage count `ord / 2 + 1` of `FIRKodeSetOrder` is `⌈(ord + 1)/2⌉`,
    the formula in `User_callable.rst`, here `(ord + 2) / 2`. Unlike
    `ord + 2`, `ord / 2 + 1` cannot overflow `int`. -/
theorem setOrder_ceil (ord : Nat) : ord / 2 + 1 = (ord + 2) / 2 := by
  omega

end Firkode.Api
