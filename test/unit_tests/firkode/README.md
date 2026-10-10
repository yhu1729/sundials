# FIRKODE unit tests

The tests of FIRKODE (design document `src/firkode/DESIGN.md`, Chapter 11). Each
test states its pass criterion in its header comment.

## Layout

| Directory    | Tests                                                                            | Judged by                                                                                                                                                                             |
| ------------ | -------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `C_serial`   | unit tests of the integrator, FIRKLS, the stage stack and the tables             | exit code                                                                                                                                                                             |
| `CXX_serial` | accuracy, order and mass matrix tests on model problems; not in single precision | exit code and answer file (`.out`), except the Krylov-path tests `firk_test_krylov`, `firk_test_mass` and `firk_test_mass_fem_krylov`, whose iteration counts differ across platforms, and `firk_test_radau5_steps`, whose accept/reject decisions flip at roundoff |
| `C_parallel` | the stage stack with distributed vectors, on 1 and 4 ranks                       | exit code                                                                                                                                                                             |

`firk_test_mocks.{h,c}` holds the recording mocks (`SUNLinearSolver`, vector
operations) that `firk_test_firkls` and the stack tests use. The tests link the
FIRKODE object library, so they can call internal functions such as
`firkSetNumStagesAny`, which runs the tables and paths for s ≠ 3 that the public
API of v0.1.0 does not accept.

## Running

```bash
ctest --test-dir build -R firk
```

`C_parallel` is built when MPI and `NVECTOR_PARALLEL` are enabled.
`firk_test_logging` and `firk_test_logging_suntools` are registered at
`SUNDIALS_LOGGING_LEVEL` 4 or more; the second parses the log of the first with
`suntools.logs`. Tests that count warnings (`firk_test_hnil`,
`firk_test_rootfind`) count them only when warnings are logged
(`SUNDIALS_LOGGING_LEVEL` 2 or more).

## Hazards

The hazard checklist of the design (Table 11.3) and the tests that enforce each
row. The tag is the earliest milestone that must close the row.

| Row | Tag | Requirement                                                                                | Test                                                                                  |
| --- | --- | ------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------- |
| H1  | M2  | a slow Newton rate re-evaluates J and sets up at the next step, under both policies        | `firk_test_reuse`                                                                     |
| H2  | M1  | the dense and band difference-quotient Jacobians survive a recoverable RHS failure         | `firk_test_firkls`                                                                    |
| H3  | M2  | `Reset` keeps every counter, `ReInit` zeroes them                                          | `firk_test_reset`, with dense and with Krylov solvers                                 |
| H4  | M1  | stage-solver setters never free a live object                                              | `firk_test_firkls` (the object stays attached; ASan)                                  |
| H5  | M2  | `SetUserData` reaches every callback                                                       | `firk_test_userdata`                                                                  |
| H6  | M1  | `jtsetup` once per Newton iteration, and per estimate only with a matrix-free block solver | `firk_test_firkls`                                                                    |
| H7  | M1  | the difference-quotient Jv of a zero vector is 0                                           | `firk_test_firkls`                                                                    |
| H8  | M2  | J is never re-evaluated at the same (t_n, y_n)                                             | `firk_test_reuse` (failures forced by the RHS)                                        |
| H9  | M1  | no forbidden operation on a stacked vector; the `Init` check of local reductions           | `firk_test_stack_ops`; `firk_test_stack_mpi` (the `Init` check)                       |
| H10 | M1  | lengths from `N_VGetLength` only                                                           | `firk_test_stack_ops`, `firk_test_stack_mpi`                                          |
| H11 | M1  | a zero guess before every stacked, block and mass solve                                    | `firk_test_stack_ops` (stacked), `firk_test_firkls` (block, mass)                     |
| H13 | M2  | `dsm` floored at 1e-10                                                                     | `firk_test_adapt` (and a polynomial solution, not in single precision)                |
| H14 | M3  | no clone per step; identical allocation on all ranks                                       | `firk_test_stack_ops`, `firk_test_stack_mpi`; a NULL constructor return by review     |
| H15 | M1  | multi-vector operations on stacks only inside SPFGMR                                       | `firk_test_stack_ops`, `firk_test_stack_mpi`                                          |
| H16 | M1  | the weight view is never cloned and is destroyed without ownership                         | `firk_test_stack_ops` (clones; ASan); its rebuild at an order change is M5            |
| H17 | M1  | no small-residual shortcut in the estimate solve                                           | `firk_test_firkls`, `firk_test_errest`                                                |
| H18 | M1  | `SUNLS_RES_REDUCED` accepted only at Newton iteration 1                                    | `firk_test_firkls` (the block solver; the internal stacked solver takes no injection) |
| H19 | M2  | recoverable and unrecoverable RHS failures in a stage                                      | `firk_test_recoverable`, `firk_test_newton`                                           |
| H20 | M2  | an RHS failure in the refilter keeps the first estimate                                    | `firk_test_recoverable`                                                               |
| H21 | M2  | the stop time: the clamped step, refusal behind tn, clearing, fixed steps after it         | `firk_test_tstop`                                                                     |
| H22 | M2  | rootfinding on the collocation polynomial, directions, inactive roots, a root at t0        | `firk_test_root`, `firk_test_rootfind`                                                |
| H27 | M1  | the vector budget of Table 8.7                                                             | `firk_test_workspace`                                                                 |
| H28 | M1  | mass solver pairing, `time_dep`, M set up once, mass counters                              | `firk_test_mass`, `firk_test_reset`                                                   |
| H29 | M1  | `jok` and `jcur` honored for `psetup` and `linsys` under both policies                     | `firk_test_firkls`                                                                    |
| H30 | M2  | the dead band keeps the factorization                                                      | `firk_test_reuse`                                                                     |
| H31 | M1  | balanced log records on every exit                                                         | `firk_test_logging`, `firk_test_logging_suntools`                                     |
| H32 | M2  | t + h = t warnings capped and logged                                                       | `firk_test_hnil`                                                                      |
| H35 | M1  | a zero weight is an error; weights once per step                                           | `firk_test_ewt`                                                                       |
| H36 | M2  | `EstimateStep` once per attempt, before `UpdateH`; the predictive default                  | `firk_test_adapt`                                                                     |
| H37 | M1  | the stacked solve applies only the user's `psolve`                                         | `firk_test_firkls`                                                                    |

Rows without a unit test:

| Row | Tag | Requirement                                                               | Test                                                                                                                 |
| --- | --- | ------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| H12 | M8  | `SUNStepper` glue                                                         | none yet: the `SUNStepper` interface is a later milestone                                                            |
| H23 | M5  | dense output and roots across an order change                             | none yet: order changes are a later milestone                                                                        |
| H24 | M0  | literals in every precision under `-Wconversion` and `-Wdouble-promotion` | the CI builds in single, double and extended precision with `SUNDIALS_ENABLE_ALL_WARNINGS` (`test/env/setup_env.sh`) |
| H25 | M0  | C99, no variable-length arrays                                            | the same CI builds (`-Wvla -Walloca`)                                                                                |
| H26 | M3  | answer files in the answers repository for three precisions               | open: the answer files are in this directory only                                                                    |
| H33 | M5  | an allocation failure during an order change                              | none yet: order changes are a later milestone                                                                        |
| H34 | M2  | examples return a nonzero exit code when the solver fails                 | review of `examples/firkode`; CTest judges the examples by their exit code                                           |

## Milestone criteria

| Criterion (design Chapter 10)                                      | Test                                                                                                                               |
| ------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------- |
| M0: tables, order conditions, closed forms                         | `firk_test_tables`                                                                                                                 |
| M1: one fixed step reproduces R(hL) y0; L-stability; Newton count  | `firk_test_dahlquist`                                                                                                              |
| M1: fixed-step order on Prothero–Robinson and its Kværnø variant   | `firk_test_order`                                                                                                                  |
| M1: matrix-free path against the dense path                        | `firk_test_krylov`                                                                                                                 |
| M1: M y' = f against y' = M⁻¹ f, mass counters                     | `firk_test_mass`                                                                                                                   |
| M1: `GetDky` exact to degree s, not s + 1                          | `firk_test_dense_output`                                                                                                           |
| M1: vector budget of Table 8.7                                     | `firk_test_workspace`                                                                                                              |
| M1: modes and the EXPERT gate                                      | `firk_test_mode`                                                                                                                   |
| M2: estimate order p = 3; dense and iterative estimates within 50% | `firk_test_errest`                                                                                                                 |
| M2: Robertson and van der Pol accuracy and cost                    | `firk_test_robertson`, `firk_test_vdp`                                                                                             |
| M2: steps within 20% of `radau5.f`; error within twice its own or below half the tolerance (today every cell passes on the latter, two with ratios 3.6 and 23, see the test header) | `firk_test_radau5_steps` (reads `radau5/*.ref`) |
| M2: the automatic Newton tolerance (RADAU5's rule)                 | `firk_test_nlscoef`                                                                                                                |
| M2: stop time, reset, roots, options, logging                      | `firk_test_tstop`, `firk_test_reset`, `firk_test_root`, `firk_test_rootfind`, `firk_test_options`, `firk_test_logging`             |
| M2: every recoverable-failure path                                 | `firk_test_recoverable`, `firk_test_newton`                                                                                        |
| M3: 1-rank and 4-rank agreement, no clone per step                 | `firk_test_stack_mpi` (fixed steps); the example `examples/firkode/C_parallel/firk_robertson_batched_p` (fixed and adaptive steps) |
