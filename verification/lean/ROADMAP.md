# FIRKODE verification roadmap

These are components of FIRKODE that could be verified beyond what the `cert`
and `spec` packages cover. Difficulty: E easy, M medium, H hard. Items that
are done are listed under [Completed](#completed); an item done only in part
keeps its number and lists what is left.

## Mathematics

| # | Component | Approach and prerequisites | Diff. |
|---|---|---|---|
| R1 | B(2s−1) for every s (D(s−1) then follows, `Colloc.collocA_D`) | Orthogonality of the shifted Legendre polynomials on [0, 1], by integration by parts from Mathlib's Rodrigues formula `Polynomial.factorial_mul_shiftedLegendre_eq`. | M |
| R2 | The Radau nodes are real, simple, lie in (0, 1], and end with c_s = 1, for every s (done for s ≤ 9) | Rolle's theorem or real-rootedness of the node polynomial. | H |
| R3 | det(I − zA) equals the Padé (s−1, s) denominator for every s (done for s ≤ 9) | Companion-matrix structure of the collocation matrix, or the argument of `Charpoly` with the coprimality of the Padé polynomials proved for every s. | H |
| R4 | A-stability for every s (done for s ≤ 9) | The energy argument of `AStability` needs only B(2s−1), `c_s = 1` and `∫_0^1 τ^(s−1) M < 0`; for every s these follow from R1 and R2. | M/H |
| R5 | Butcher's theorem (B, C, D ⇒ order) and the rooted-tree order conditions | No Mathlib support. Prior art for reference only: lean-pool `RungeKuttaOrderConditions` (licence unknown) and uw-math-ai `OpenMath` (no licence file). | H |
| R6 | γ0 for even s = 4, 6, 8 (that there is no real eigenvalue is done) | Certified isolation of the complex eigenvalue pair of largest modulus (Kantorovich or Krawczyk tests, or the interval tactics of newer Mathlib). | M/H |
| R7 | The comment in `firk_gamma0` that γ0 "minimizes the stiff-limit deviation" | Formalize the stiff-limit preconditioned operator A/γ0. The claim holds eigenvalue by eigenvalue but numerically fails for the worst case over eigenvalues at s = 3 and 4, so the comment should be made precise. | M |
| R8 | The error estimate is O(h^{s+1}) for smooth right-hand sides | Taylor's theorem (`taylor_mean_remainder_lagrange`) on top of the polynomial exactness proved in `spec`. | M |
| R9 | Convergence: O(h^{2s−1}) for nonstiff problems; B-convergence for stiff problems | Gronwall and discrete Gronwall are in Mathlib; Runge–Kutta convergence theory is not. | H |
| R10 | Equivalence of the mass-matrix stage equations with the M⁻¹f form | Linear algebra. | E/M |
| R11 | The stacked Newton operator; the preconditioner is exact for s = 1; the effect of using γ from the last setup | Linear algebra. | M |

## Algorithmic logic

These items model C control flow over the rationals or reals.

| # | Component | Approach and prerequisites | Diff. |
|---|---|---|---|
| R12 | Nonlinear convergence test `firkNlsConvTest` | At iteration maxcor − 1 the convergence and predicted-failure tests are complementary, so iteration maxcor never runs. Also bounds on the returned step ratio. | E/M |
| R13 | Step-size selection `firkAdapt` and the rejection path | Bounds on the step ratio η, hmin/hmax handling, ImpGus with p = s, and fixed-step mode with a stop time (η from the tstop clip persists). | E/M |
| R14 | Rootfinding with the Illinois method | Bracket invariant, interval shrinkage, termination within the time tolerance. | M |
| R15 | Initial-step heuristic | hlb ≤ h0 ≤ hub. | E |
| R16 | The `FIRKodeTable_CheckOrder` tolerance separates satisfied from violated conditions (s ≤ 9 in double, s ≤ 3 in single precision); the exact-arithmetic result is done | The exact residual formulas of the node-polynomial certificates plus a rounding model. | M |

## Floating point and C code

| # | Component | Approach and prerequisites | Diff. |
|---|---|---|---|
| R17 | Rounding error of the `long double` table construction | Lean has no 80-bit extended-precision model. | H |
| R18 | Double-precision kernels such as dense-output evaluation | Core Lean `Float.Model` (Lean ≥ 4.33) with FloatSpec (a Flocq port); needs a toolchain upgrade. | H |
| R19 | The C code itself | No mature C-to-Lean path exists (CIRCE handles an integer subset only). Memory safety of the table routines could be checked with Frama-C or CBMC outside Lean. | — |

## Infrastructure

| # | Component | Approach and prerequisites | Diff. |
|---|---|---|---|
| R20 | CI job | A `check-lean.yml` workflow gated by `ci.yml` change detection on `verification/` and `src/firkode/firkode_tables.c`. Re-pin Mathlib to a release tag first. | E |
| R21 | `.lean` handler in `scripts/updateCopyright` | Same pattern as the `.py` handler. | E |
| R22 | Upgrade the toolchain from v4.32.0-rc1 to a stable release | Re-pin Mathlib and re-check the build cache. | E/M |

## Completed

| # | Component | Lean |
|---|---|---|
| R23 | Uniqueness of the real root behind each odd-s `firk_gamma0` literal | `RootCount.uniqueRootIn_spec`, `RootCount.gamma0_odd_unique`, from `RuleOfSigns` |
| R2, s ≤ 9 | The Radau nodes exist for every stage count FIRKODE supports, one in each bracket of the grid scan | `Nodes.radauNodes` with `radauNodes_nodePoly`, `radauNodes_strictMono`, `radauNodes_pos`, `radauNodes_last`, `radauNodes_bracket` |
| R3, s ≤ 9 | `det(I − zA)` is the Padé denominator, and the real eigenvalues of `A` are the real roots of `det(μI − A)`, for every stage count FIRKODE supports; hence the `firk_gamma0` comment for odd s, and no real eigenvalue for even s | `Charpoly.radau_charpolyRev`, `Charpoly.radau_isEigen_iff`, `Charpoly.radau_gamma0_odd_inv`, `Charpoly.radau_even_no_real_eigenvalue` |
| R4, s ≤ 9 | A-stability and L-stability of the Radau IIA methods, by an energy identity instead of the cited E-polynomial and Routh–Hurwitz criteria | `Colloc.energy_identity`, `Colloc.radau_A_stable`, `Colloc.radau_L_stable` |
| R1 in part | D(s−1) for every s from B(2s−1) and C(s); with the failures of B(2s), C(s+1) and D(s), `FIRKodeTable_CheckOrder` returns `(2s − 1, s, s − 1, 2s − 1)` in exact arithmetic for s ≤ 9 (also the exact part of R16) | `Colloc.collocA_D`, `Colloc.collocA_not_C_succ`, `Colloc.collocA_not_D`, `CheckOrder.radau_checkOrder` |
