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

import Lean
import FirkodeSpec.Mirror
import FirkodeSpec.Stability
import FirkodeSpec.Radau
import FirkodeSpec.RootCount
import FirkodeSpec.Nodes
import FirkodeSpec.Charpoly
import FirkodeSpec.Integral
import FirkodeSpec.AStability
import FirkodeSpec.CheckOrder

/-!
# Axiom audit

Checks that every theorem in the `Firkode` namespace, of `spec` and of the
`cert` modules it imports, depends only on the standard axioms `propext`,
`Classical.choice` and `Quot.sound`. In particular
no proof uses `sorry` (`sorryAx`) or `native_decide` (`Lean.ofReduceBool`,
`Lean.trustCompiler`), which would trust the compiler instead of the kernel.
-/

open Lean Elab Command

/-- Theorems Lean generates for structures and definitions (injectivity,
    `sizeOf`, equation lemmas); they are audited but not counted. -/
private def isGeneratedSpec (n : Name) : Bool :=
  match n with
  | .str _ s => s == "inj" || s == "injEq" || s == "sizeOf_spec" || s.startsWith "eq_"
  | _ => false

/-- Fails if a theorem in the `Firkode` namespace depends on an axiom other
    than `propext`, `Classical.choice` and `Quot.sound`; otherwise reports how
    many certificate theorems were audited. -/
elab "#audit_firkode_spec_axioms" : command => do
  let allowed : List Name := [``propext, ``Classical.choice, ``Quot.sound]
  let theorems := (← getEnv).constants.fold (init := []) fun acc n info =>
    match info with
    | .thmInfo _ => if (`Firkode).isPrefixOf n && !n.isInternal then n :: acc else acc
    | _ => acc
  for n in theorems do
    for ax in ← collectAxioms n do
      unless allowed.contains ax do
        throwError "{n} depends on the axiom {ax}"
  let count := (theorems.filter fun n => !isGeneratedSpec n).length
  logInfo m!"{count} theorems use only propext, Classical.choice and Quot.sound"

/-- info: 337 theorems use only propext, Classical.choice and Quot.sound -/
#guard_msgs in
#audit_firkode_spec_axioms
