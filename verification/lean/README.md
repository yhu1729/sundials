# FIRKODE Lean verification

Kernel-checked proofs of the mathematics that FIRKODE (`src/firkode`) relies
on: the Radau IIA tables and their order, the `firk_gamma0` shifts, the
error-estimate weights, the dense-output coefficients, and the stability
function.

The proofs use exact arithmetic. They do not model the C code or its
floating-point rounding, and they do not cover the behaviour of the nonlinear
solver, the linear solvers, the step controller or rootfinding.
[ROADMAP.md](ROADMAP.md) lists what could be verified later.

This directory is not part of the CMake build, the test suite or the release
tarballs.

## Layout

| Path | Contents |
|---|---|
| `cert/` | Lake package `firkode_cert`. Core Lean only, no dependencies. Exact certificates checked with `decide +kernel`. |
| `spec/` | Lake package `firkode_spec`. Depends on Mathlib and `cert`; theorems for a general number of stages. |
| `gen_firkode_tables.py` | Generates `cert/FirkodeCert/Generated.lean` from the C sources. |

Modules of `cert/FirkodeCert`:

| Module | Contents |
|---|---|
| `Basic` | Exact arithmetic in ℚ(√6) (`Q6`), dense polynomials with division and the extended Euclidean algorithm over ℚ, small matrices, Padé approximants of `e^z`. |
| `Generated` | Tables, `firk_gamma0[]`, constants and routine fingerprints extracted from the C sources. |
| `Mirror` | Transcriptions of the C routines that build and use the tables. |
| `Tables` | Certificates for the closed-form tables, s ≤ 3. |
| `NodePoly` | Certificates from the Radau node polynomial, s ≤ 9. |
| `Gamma0` | Certificates for the `firk_gamma0` literals. |
| `Api` | The stage count chosen by `FIRKodeSetOrder`. |
| `Axioms` | Audit that every theorem uses only the standard axioms. |

Modules of `spec/FirkodeSpec`:

| Module | Contents |
|---|---|
| `ListPoly` | The coefficient lists of `cert` as Mathlib polynomials. |
| `Bridge` | `Q6.toReal : ℚ(√6) → ℝ` preserves the field operations and is injective. |
| `Collocation` | Collocation methods for any s: C(s), exact quadrature modulo the node polynomial, B(2s−1) and the quadrature error in degree 2s − 1 from orthogonality, D(s−1), the failure of C(s+1) and D(s), uniqueness, `det A = Π c_i / s!`, and the coefficients under a change of field. |
| `Integral` | The algebraic antiderivative is the interval integral over the reals; the weights of a method with B(2s−1) are `b_j = ∫_0^1 l_j² ≥ 0`. |
| `AStability` | The energy identity `|Y_s|² − |y0|² = 2 Re z Σ b_j |Y_j|² + 2s |u_s|² J` for complex `z`, hence A-stability and L-stability of collocation methods with B(2s−1) and `J = ∫_0^1 τ^(s−1) M < 0`, in particular of the Radau IIA methods (s ≤ 9). |
| `ErrorEstimate` | The error-estimate weights `e_j = −γ0 l_j(0)/c_j`, the estimate identity, and its order `p = s`, for any s. |
| `DenseOutput` | The dense-output basis `L_j`: interpolation of the stages, `L′(c) A = I`, and the `FIRKodeGetDky` derivative formula, for any s. |
| `Mirror` | The transcribed C routines compute the collocation quantities, for any s. |
| `Stability` | The collocation polynomial of `y′ = zy`, Nørsett's formula `R = N/D` for collocation methods, and stiff decay `R(z) → 0` as `z → −∞`, for any s. |
| `RuleOfSigns` | Descartes' rule of signs on intervals: the sign-variation counts of `cert` are Mathlib's `Polynomial.signVariations`, and bound the roots of `p` in an interval through the transformed polynomials `mobius`, `reflect` and `shift`. |
| `RootCount` | The γ0 certificates as root counts of `μ^s Q_s(1/μ)`: for odd s exactly one real root, within half a unit in the last printed digit of `firk_gamma0[s]`; for even s none. |
| `Nodes` | The Radau nodes exist (s ≤ 9): the grid scan of `firk_radau_nodes` and the intermediate value theorem give the `s` roots `0 < c_1 < … < c_s = 1` of the node polynomial `M_s`, one in each bracket of the scan; the `Radau` results then hold for them without hypotheses. |
| `Charpoly` | For s ≤ 9, `det(I − zA)` and `det(I − z(A − 𝟙bᵀ))` are the Padé polynomials, from Nørsett's formula and their coprimality; the real eigenvalues of `A` are the real roots of `det(μI − A)`, which gives the `firk_gamma0` comment: a unique real eigenvalue of `A⁻¹` for odd s, none for even s. |
| `CheckOrder` | `FIRKodeTable_CheckOrder` in exact arithmetic returns `(2s − 1, s, s − 1, 2s − 1)` on the Radau IIA tables (s ≤ 9). |
| `Radau` | The FIRKODE tables over the reals: the closed-form tables are the collocation tables (s ≤ 3); the node recurrence gives the shifted Legendre polynomials, and B(2s−1), the B(2s) residual, `e_s` and the Padé stability function hold for the Radau nodes (s ≤ 9). |
| `Axioms` | Audit that every theorem of `spec`, and of the `cert` modules it uses, uses only the standard axioms. |

## Building

The toolchain is pinned in `cert/lean-toolchain`
(`leanprover/lean4:v4.32.0-rc1`); install it with
[elan](https://github.com/leanprover/elan). Then:

```sh
cd verification/lean/cert
lake build
```

A full build takes under 20 seconds, most of it in `NodePoly`. Build outputs
go to `.lake/`, which is ignored by git and can be deleted at any time.

`cert/FirkodeCert/Generated.lean` holds the closed-form tables, the
`firk_gamma0[]` literals and a few constants, extracted from the C sources by
the generator. After changing those sources, run

```sh
python3 verification/lean/gen_firkode_tables.py --check
```

It exits with status 1 when the generated file is out of date. The file also
records a SHA-256 fingerprint of each C routine that is transcribed by hand,
computed after removing comments and normalizing whitespace, so a code change
in one of them fails the check too. Review the Lean transcription of that
routine, then run the generator without `--check` to update the file.

### The Mathlib package

`spec` pins Mathlib to commit `4efb186f102ebfd2eea1545c151d6fbcfdff0e43`,
whose toolchain is also `leanprover/lean4:v4.32.0-rc1`. Mathlib must not be
rebuilt from source (that takes hours), so unpack its prebuilt files for the
modules `spec` imports, and check that nothing is left to build, before the
first `lake build`:

```sh
cd verification/lean/spec
MATHLIB_NO_CACHE_ON_UPDATE=1 lake update --keep-toolchain
MODS=$(grep -ho '^import Mathlib[.A-Za-z0-9_]*' FirkodeSpec/*.lean | cut -d' ' -f2 | sort -u)
lake exe cache get $MODS
lake build --no-build $(printf '+%s ' $MODS)
lake build
```

`lake update` clones Mathlib and its dependencies (about 0.7 GB) into
`spec/.lake`. `lake exe cache get` takes the prebuilt files from
`~/.cache/mathlib` when they are there and downloads the rest; the imported
part of Mathlib unpacks to a few GB. If `lake build --no-build` reports
anything to build, stop: Mathlib would be compiled from source. Do not run
`lake exe cache clean`, which empties the shared `~/.cache/mathlib`. In zsh,
write `${=MODS}` instead of `$MODS`.

## Trust model

| Kind | What |
|---|---|
| Kernel-checked | Every certificate is proved with `decide +kernel`, or `omega` in `Api`; the `spec` theorems are ordinary Mathlib proofs. The `Axioms` modules of both packages fail the build if a theorem depends on an axiom other than `propext`, `Classical.choice` and `Quot.sound`; that rules out `sorry` and `native_decide`. |
| Generated from C | The closed-form tables of `firk_radau_closed_form`, the `firk_gamma0[]` literals, and the constants `FIRK_MAX_STAGES`, `FIRK_MAX_CHECK_ORDER`, `FIRK_NODE_GRID_PER_STAGE` and `FIRK_DEFAULT_STAGES`. |
| Hand-mirrored | The Lean transcriptions in `Mirror`, checked by review; each C routine is fingerprinted by the generator. The inverse `A⁻¹` that `firk_derive_ld` computes by LU factorization is taken to be the exact inverse. |
| Cited, not formalized | Butcher's theorem (B, C, D ⇒ order). |
| Not covered | Rounding in the C code; the semantics of the C code itself. |

## Claims and where they are proved

"Proved" means kernel-checked. "Needs spec" means the certificate is proved
but the step from it to the claim is one of the `spec` items above.

| Claim | Source | Lean | Status |
|---|---|---|---|
| The nodes are the roots of `P_s(2c−1) − P_{s−1}(2c−1)` with `0 < c_1 < … < c_s = 1` | `Mathematics.rst` | `Tables.closedForm_nodes_radau`, `Tables.closedForm_nodes_increasing` (s ≤ 3); `Nodes.radauNodes_nodePoly`, `Nodes.radauNodes_strictMono`, `Nodes.radauNodes_pos`, `Nodes.radauNodes_last` (s ≤ 9); the recurrence of `firk_radau_poly` gives Mathlib's shifted Legendre polynomials, `Radau.radauPoly_eq_shiftedLegendre` (s ≤ 9) | Proved |
| The bracketing grid of `firk_radau_nodes` separates the interior nodes | `firkode_tables.c` | `NodePoly.nodes_grid`, `Nodes.radauNodes_bracket`, `Nodes.radauNodes_mem_bracket` (s ≤ 9) | Proved |
| The closed-form tables equal the computed collocation tables | `firkode_tables.c`, `firk_test_tables.c` | `Tables.closedForm_eq_collocation`; over the reals `Radau.closedForm_eq_collocA` | Proved |
| C(s), and `FIRKodeTable_CheckOrder` returns q = 2s − 1 | `Mathematics.rst`, `Tables.rst` | `Tables.closedForm_checkOrder` (s ≤ 3); `CheckOrder.radau_checkOrder` (s ≤ 9), with D(s−1) from `Colloc.collocA_D` and the failures of B(2s), C(s+1) and D(s) from `Colloc.collocA_B2s`, `Colloc.collocA_not_C_succ` and `Colloc.collocA_not_D` | Proved, in exact arithmetic; order cited (Butcher) |
| B(2s−1) | `Mathematics.rst` | `Tables.closedForm_checkOrder` (s ≤ 3); `Radau.radau_B` from `NodePoly.nodePoly_orthogonal`, giving `Nodes.radauNodes_B` (s ≤ 9) | Proved |
| The smallest violated-condition residual is about 1e−10 (B(18), s = 9) | `FIRKodeTable_CheckOrder` comment | `Radau.radau_B2s_residual` from `NodePoly.nodePoly_B2s_residual` (9.4e−11), giving `Nodes.radauNodes_B2s_residual` | Proved |
| `d = bᵀA⁻¹ = e_s`, so `y_{n+1} = y_n + Z_s`, and `R(∞) = 0` | `Mathematics.rst` | `Tables.closedForm_d` (s ≤ 3) | Proved |
| The stability function is the (s−1, s) Padé approximant | `Mathematics.rst` | `Tables.closedForm_pade` (s ≤ 3); `Colloc.norsett` with `NodePoly.nodePoly_norsett`, giving `Radau.radau_pade` and `Nodes.radauNodes_pade`; `det(I − zA)` and `det(I − z(A − 𝟙bᵀ))` are the Padé polynomials, `Charpoly.radau_charpolyRev`, `Charpoly.radau_charpolyRev_num` (s ≤ 9) | Proved |
| A-stability and L-stability | `Mathematics.rst` | `Colloc.radau_A_stable`: for `Re z ≤ 0`, `I − zA` is invertible and `|R(z)| ≤ 1`; `Colloc.radau_L_stable`: `R(z) → 0` as `|z| → ∞` (s ≤ 9), from the energy identity `Colloc.energy_identity` (any s, given B(2s−1) and a negative B(2s) error). Independent checks: `Tables.closedForm_Epoly` (s ≤ 3), `NodePoly.pade_Epoly`, `NodePoly.pade_poles_right` (s ≤ 9) | Proved |
| γ0 = 1/U1 (odd s) and 1/Re λ (even s) | `Mathematics.rst`, `firk_gamma0` comment | `Gamma0.gamma0_one`, `Gamma0.gamma0_two`; for s = 3, 5, 7, 9, `A⁻¹` has exactly one real eigenvalue `U1` and `firk_gamma0[s]` is within half a unit in its last printed digit of `1/U1`, `Charpoly.radau_gamma0_odd_inv` (from `RootCount.gamma0_odd_unique`); for s = 2, 4, 6, 8 neither `A` nor `A⁻¹` has a real eigenvalue, `Charpoly.radau_even_no_real_eigenvalue` | Proved for s ≤ 2 and for odd s; for even s ≥ 4 the value `1/Re λ` is numeric only (roadmap R6) |
| The error estimate has order p = s | `Mathematics.rst`, `Tables.rst` | `Colloc.estimate_eq`, `Colloc.bracket_lt`, `Colloc.bracket_X_pow` (any s); `Tables.closedForm_moments` (s ≤ 3) | Proved; the O(h^{s+1}) bound for smooth `f` is roadmap R8 |
| `e_s = (−1)^s γ0/s` | `firk_test_tables.c` | `Tables.closedForm_e_last` (s ≤ 3); `Radau.radau_e_last` from `NodePoly.nodePoly_interior_product`, giving `Nodes.radauNodes_e_last` (s ≤ 9) | Proved |
| For s = 3 the weights are RADAU5's (`e = γ0 · DD`) | `Mathematics.rst`, `firk_test_tables.c` | `Tables.radau3_e_radau5` | Proved |
| `e_j = −γ0 l_j(0)/c_j` | `firkode_tables.c` | `Colloc.errWeights_eq` (any s); `Tables.closedForm_e` (s ≤ 3) | Proved |
| The dense output interpolates the stages and ends at `y_{n+1}` | `Mathematics.rst` | `Colloc.dense_interp`, `Colloc.denseL_eval_one` (any s); `Tables.closedForm_dense` (s ≤ 3) | Proved |
| The `FIRKodeGetDky` derivative satisfies the collocation equations | `firkode.c` | `Colloc.denseL_deriv_collocA`, `Colloc.iterate_derivative_denseL_eval` (any s); `Tables.closedForm_dense_deriv` (s ≤ 3) | Proved |
| `firk_lagrange_coeffs` and `firk_collocation_coeffs` compute the collocation coefficients | `firkode_tables.c` | `lagrangeCoeffs_eq`, `collocationA_eq` (any s) | Proved, in exact arithmetic |
| `firk_derive_ld` computes `d = e_s`, `e_j = −γ0 l_j(0)/c_j` and the dense-output coefficients `P` | `firkode_tables.c` | `derive_d_eq`, `derive_e_eq`, `derive_P_eq` (any s) | Proved, in exact arithmetic, given the exact `A⁻¹` |
| `FIRKodeGetDky` and `firkPredict` evaluate the dense output and its derivatives | `firkode.c` | `denseWeights_eq`, `dkyWeights_eq` (any s) | Proved, in exact arithmetic |
| `FIRKodeSetOrder` picks the smallest s with 2s − 1 ≥ ord | `firkode_io.c`, `User_callable.rst` | `Api.setOrder_minimal`, `Api.setOrder_supported`, `Api.setOrder_default` | Proved |
| γ0 minimizes the stiff-limit deviation of the preconditioned matrix | `firk_gamma0` comment | — | Open (roadmap R7) |
