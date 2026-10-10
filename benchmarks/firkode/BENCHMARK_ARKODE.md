# FIRKODE versus ARKODE: implicit Runge-Kutta methods and structure-exploiting steppers

Status: pre-registered. The question, the decision rules, the methodology and the predictions
and pass and fail rules of every family are fixed before the full sweep; the generated
sections are filled by `scripts/make_report.py --report arkode` from the CSV files of
`results/<run-id>/`. Sections marked **[pending]** have no data yet; Section 9 gives the
commands. The companion report `BENCHMARK_BDF.md` compares FIRKODE with CVODE and IDA.

Evidence tags: **[measured]** = a row of a CSV file in `results/`, **[predicted]** = written
down before the runs, **[development]** = seen in a single run while the harness was written
(Section 9.2; not a result, restated as a prediction where it shaped one), **[verified]** =
read from the SUNDIALS sources or examples at the commit of this report.

## 1. Question and decision rules

### 1.1 Question

At matched accuracy, is FIRKODE (Radau IIA, s = 3, order 5, stage order 3, L-stable and
stiffly accurate) faster than ARKODE, and can every win or loss be attributed to a named
mechanism? Two comparisons are asked separately:

1. **Implicit Runge-Kutta against implicit Runge-Kutta** (families A0-A4): FIRKODE against
   ARKStep with its diagonally implicit (DIRK, ESDIRK, SDIRK) tables on the full problem.
   Both are one-step methods with an A-stable (mostly L-stable) implicit method; they differ
   in stage order, error constant, the cost of a step and the step-size control.
2. **Full implicitness against structure** (families A5-A8): FIRKODE, which treats the
   whole right-hand side implicitly, against ARKODE steppers that exploit a structure the
   user knows: ARKStep IMEX (an explicit and an implicit part), LSRKStep (explicit
   stabilized methods for a near-real spectrum with an eigenvalue bound) and MRIStep
   (multirate, a slow and a fast part).

### 1.2 Decision rules

A family **supports** its candidate if FIRKODE s = 3 is faster than the stated ARKODE class
at matched mixed digits (mescd = -log10 of the mixed error at the final time, interpolated in
log time) in at least one regime **and** the counters show the predicted mechanism; it
**refutes** the candidate if the class stays faster across the sweep. A win is called only at
>= 1.2x and when adjacent tolerance cells agree; concurrent workers add noise (rows record
`nworkers`).

The DIRK class is reported three ways, because the minimum over several noisy
configurations is biased low:

- **best default**: the fastest of ARKODE's default tables of order 3, 4 and 5;
- **best of all**: the fastest of the nine declared DIRK configurations (Section 3.1), all
  chosen before the runs;
- **tuned**: the settings of `examples/arkode/C_serial/ark_robertson.c` on ROBER, tuned on
  that problem by its authors; reported apart and never part of a verdict.

A verdict is stated for the best default and for the best of all; it is "robust" only when
both agree.

### 1.3 Cost model behind the candidates [predicted]

- **FIRKODE**: each simplified Newton iteration costs 3 right-hand sides and 3 k_lin block
  solves with one N x N matrix `M - h gamma0 J` (k_lin = stacked FGMRES iterations per Newton
  iteration, 2.5-5 on the small problems [measured in `BENCHMARK_BDF.md`, B3]), about 2
  Newton iterations per step; one factorization per setup.
- **ARKStep DIRK**: an ESDIRK of s stages solves s - 1 implicit stages in sequence, each with
  1-3 Newton iterations of one right-hand side and one solve, all with the same
  `M - h gamma J` (one diagonal coefficient), so also one factorization per setup.
  ARKODE's order-4 default (6 stages) costs about 5 x 2 = 10 right-hand sides and solves per
  step, FIRKODE about 6 right-hand sides and 15-30 block solves: comparable per-step work.
- **Steps decide.** No ARKODE DIRK table has stage order above 2 (Section 3.1); Radau IIA
  has stage order 3. In the stiff limit of a problem with time-dependent forcing a method's
  global error falls to about its stage order plus one (Prothero-Robinson analysis,
  Hairer-Wanner II, Section IV.15 [verify]), so DIRKs need more steps as the tolerance
  tightens. On the small problems ARKODE's DIRK defaults were slower than FIRKODE s = 3 in
  most cells of an older suite (DIRK4 by 1.3-20x, DIRK5 by 1.3-76x at mescd 4-8, with DIRK3
  and DIRK5 faster in one HIRES cell each) [measured in the 2026-10-05 suite on
  `work/firkode-benchmark` with an older harness and older FIRKODE defaults; not comparable
  with this sweep].
- **Structure.** IMEX removes the nonlinear and the full-Jacobian work of the explicit part
  at the price of its step-size limit; LSRK has no linear algebra at all at the price of
  about sqrt(h rho) stages per step and order 2; MRI takes few slow implicit steps while the
  cheap fast part takes many explicit ones.

## 2. Environment and provenance

<!-- begin generated: provenance -->
| item | value |
| --- | --- |
| Machine | AMD Ryzen Threadripper PRO 9985WX 64-Cores, 377 GiB, 128 logical CPUs |
| OS | Ubuntu 24.04.5 LTS (Linux-6.8.0-146-generic-x86_64-with-glibc2.39) |
| Compiler | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0; CMake C++ compiler `/usr/bin/c++` |
| CMake / generator | cmake version 3.28.3 / Ninja |
| Build type and flags | Release: C `-O3 -DNDEBUG`, C++ `-O3 -DNDEBUG` |
| Precision / index size | DOUBLE / 64 |
| Logging level / error checks / profiling | 2 / OFF / OFF |
| SUNDIALS | SUNDIALS 7.9.0 |
| Git commit | v7.9.0-23-g703a69cbb on feature/firkode |
| Python / numpy / matplotlib | 3.12.3 / 2.5.3 / 3.11.2 |
| Threads | OMP_NUM_THREADS=1, VECLIB_MAXIMUM_THREADS=1 |
| Run started (UTC) | 2026-10-10T03:20:38+00:00 |
| Binary SHA-256 | 77ed377bc248df4c9ea25b4d7a4604e616370370047e43f4883180362168329a |
<!-- end generated -->

## 3. Methodology

### 3.1 Harness and configurations

`firkode_benchmark` runs one cell per process and prints one CSV row; the driver
`scripts/run_benchmarks.py` sweeps the families with `--family arkode` (A0-A8), computes the
references, and runs serially or with `--workers N` (README). FIRKODE runs as in
`BENCHMARK_BDF.md` (BEGINNER mode, s = 3). The ARKODE integrators:

| integrator | stepper | configurations (label) |
|---|---|---|
| `arkode` | ARKStep, implicit part only (fe = 0) | DIRK3/4/5 (ARKODE's defaults: ESDIRK325L2SA_5_2_3, ESDIRK436L2SA_6_3_4, ESDIRK547L2SA2_7_4_5), DIRK4/5 optimal (`ARKStepSetOptimalParams`), SDIRK_5_3_4, ESDIRK437L2SA_7_3_4, KVAERNO_7_4_5, ARK548L2SAb_DIRK_8_4_5; DIRK4 [ark_robertson settings] |
| `arkode-imex` | ARKStep IMEX on the problem's `imex` splitting | IMEX3/4/5 (ARKODE's default pairs ARK324L2SA, ARK437L2SA, ARK548L2SAb), linearly implicit (`ARKodeSetLinear`) where the implicit part is linear |
| `arkode-lsrk` | LSRKStep | RKC2 and RKL2 with the problem's eigenvalue bound; RKC2 power with the power iteration on difference quotients |
| `arkode-mri` | MRIStep, explicit ARKStep inner stepper | MRI ESDIRK34a (MRI-GARK ESDIRK34a on the `mri` splitting), MRI ESDIRK34a htol (MRIHTol control), MRI IMEX-SR32 (IMEX-MRI-SR32 on `mri-imex`) |

Properties of the DIRK tables [verified: `src/arkode/arkode_butcher_dirk.def`,
`doc/arkode/guide/source/Butcher.rst`]: the seven tables above are A-stable, L-stable and
stiffly accurate; the ESDIRK tables (explicit first stage) have stage order at most 2 and
SDIRK_5_3_4 (implicit first stage) stage order 1. ARKStep rejects tables with
entries above the diagonal, so Radau IIA cannot run in ARKODE
(`src/arkode/arkode_arkstep.c`, `arkStep_CheckButcherTables`).

Defaults that differ between the packages [verified]: ARKODE's step-size controller is the I
controller and its error estimate is the embedded difference `h sum (b - d) f` without the
`(I - h gamma J)^-1` filter of RADAU5 that FIRKODE applies; its predictor is the trivial
one; Newton tolerance 0.1, at most 3 iterations, setup every 20 steps or when gamma changes
by 20 percent, a new Jacobian every 51 steps.

### 3.2 Tolerances and references

As in `BENCHMARK_BDF.md` Section 3.2, plus: `kprmri`, `pr`, `kpr`, `analytic` and `sector`
have exact solutions; family A4 uses the extended-precision references of B13 (FIRKODE
s = 3 at rtol 1e-16 in 80-bit precision, cross-checked by CVODE).

<!-- begin generated: references -->
| problem | params | mixed disagreement (reference check) | reference steps |
| --- | --- | --- | --- |
| diurnal | mx=10;my=10 | 1.8e-08 | 3823 |
| diurnal | mx=20;my=20 | 2.9e-08 | 4181 |
| diurnal | mx=40;my=40 | 2.0e-08 | 3857 |
| diurnal | mx=80;my=80 | 5.1e-09 | 3407 |
| heat1dfem | n=101 | 1.8e-12 | 373 |
| heat1dfem | n=401 | 1.2e-12 | 231 |
| heat2d | mx=200;my=100;velocity=0.0 | 1.4e-10 | 99 |
| heat2d | mx=200;my=100;velocity=0.5 | 4.2e-11 | 113 |
| heat2d | mx=200;my=100;velocity=50.0 | 3.8e-11 | 198 |
| heat2d | mx=200;my=100;velocity=500.0 | 3.4e-10 | 223 |
| hires |  | 2.9e-14 | 22347 |
| orego |  | 3.1e-13 | 292837 |
| rober |  | 9.6e-15 | 136141 |
| vdpol |  | 1.7e-14 | 281512 |
| wave1dfem | damping=kv;mass=lumped;n=100 | 1.2e-11 | 2510 |
| wave1dfem | damping=kv;mass=lumped;n=400 | 7.9e-11 | 3979 |
| wave1dfem | damping=mass;mass=lumped;n=100 | 9.6e-11 | 2647 |
| wave1dfem | damping=mass;mass=lumped;n=400 | 1.2e-10 | 6419 |
<!-- end generated -->

### 3.3 Timing protocol

As in `BENCHMARK_BDF.md` Section 3.3: median wall time of setup plus evolve over repetitions
of about one second, one cell at a time per worker, `OMP_NUM_THREADS=1`; the sentinel cell
runs before and after every family in every worker, its rows are in
`results/<run-id>/sentinel.csv` and its spread in `results/<run-id>/tables/sentinel.md`.
Profiled twins (`[profile]`) are untimed and never enter a ratio.

### 3.4 Fairness caveats

1. **Best of K.** See Section 1.2; the configuration list is fixed before the runs.
2. **Controllers and error estimates differ** (Section 3.1). Only matched accuracy is
   compared; the only ARKODE tuning is ARKODE's own deprecated `ARKStepSetOptimalParams`,
   whose parameter sets were tuned for older default tables.
3. **IMEX, LSRK and MRI use structure FIRKODE is not given**: the splitting, a linear implicit
   part (`ARKodeSetLinear`: one Newton iteration, a setup whenever gamma changes), exact
   block preconditioners of the implicit part, and eigenvalue bounds. Families A5-A8 answer
   whether a user who knows the structure should prefer ARKODE; the problems were chosen
   because they split.
4. **No stability function is supplied** to the IMEX runs (`ARKodeSetStabilityFn`), so the
   explicit part's stability limit is enforced by the error controller.
5. **LSRK eigenvalue bounds** are exact for `heat2d`, `sector`, `pr` and `analytic` and a
   Gershgorin bound for `bruss3d` (2 percent above the power iteration at n = 6 and the
   initial state [measured by `--check-domeig` in the validation]); an overestimate costs
   about the square root of its ratio in stages.
6. **MRI has many tuning choices** (coupling, inner method, inner tolerance, controller).
   Only ARKODE's inner default (an explicit ERK at the slow tolerances) and the two declared
   controllers run.
7. **`nfe` means different things**: the implicit right-hand side (ARKStep, MRI slow implicit),
   all right-hand sides (LSRK), FIRKODE's evaluations; `nfe_e` and `nfe_f` count the explicit
   and fast parts. Figures plot wall time only.
8. **Mass matrices** (A3): ARKStep with a fixed mass matrix does about 2 mass solves per
   attempted step and 1 per accepted step for tables that are not stiffly accurate
   (`doc/arkode/guide/source/Mathematics.rst` [verified]); the tables here are stiffly
   accurate, which saves the solve for the new state.

### 3.5 Harness validation

`--validate` replicates repository examples and tests and compares their counters with the
answer files by the criterion `test/testRunner` applies to them; for this report it adds
replicas of `ark_analytic` with a table chosen by name, the ARKStep logging tests (DIRK and
ImEx), `ark_analytic_lsrk` and three `ark_kpr_mri` runs (MRI-GARK IRK21a, IMEX-MRI-SR21,
MERK21), and the self-checks `--check-split` (the pieces of every splitting add up to f, the
implicit Jacobian matches differences) and `--check-domeig` (every eigenvalue bound against
the power iteration). The sweeps run only if every replica and check passes.

<!-- begin generated: validation -->
| replica | passed | exact | tolerance | mismatch | note |
| --- | --- | --- | --- | --- | --- |
| firk_robertson | yes | yes | answer | - | examples/firkode/C_serial/firk_robertson.out |
| firk_analytic | yes | yes | answer | - | examples/firkode/C_serial/firk_analytic.out |
| firk_diurnal_kry | yes | yes | 0.05 | - | the counts of the Linux workstation of the sweeps (Krylov path); examples/firkode/C_serial/firk_diurnal_kry.out reports nst=199 nfe=1197 nni=331 nsli=1905 nli=363 njtv=6078 nsetups=117 netf=10 npe=117 nps=6290 |
| firk_heat1D_fem_mass | yes | yes | answer | - | examples/firkode/C_serial/firk_heat1D_fem_mass.out |
| firk_robertson_root | yes | yes | answer | - | 179 g evaluations on the Linux workstation of the sweeps; examples/firkode/C_serial/firk_robertson_root.out reports 180 (the example advances the output time after a root return, the harness continues to the same output time) |
| firk_test_order_kpr_3_h0 | yes | yes | 0.0001 | - | test/unit_tests/firkode/CXX_serial/firk_test_order_kpr_3.out |
| firk_test_order_kpr_3_h1 | yes | yes | 0.0001 | - | test/unit_tests/firkode/CXX_serial/firk_test_order_kpr_3.out |
| firk_test_order_kpr_3_h2 | yes | yes | 0.0001 | - | test/unit_tests/firkode/CXX_serial/firk_test_order_kpr_3.out |
| firk_test_order_kpr_3_h3 | yes | yes | 0.0001 | - | test/unit_tests/firkode/CXX_serial/firk_test_order_kpr_3.out |
| firk_test_order_pr_3_-1e4_h0 | yes | yes | 0.0001 | - | test/unit_tests/firkode/CXX_serial/firk_test_order_pr_3_-1e4.out |
| firk_test_order_pr_3_-1e4_h1 | yes | yes | 0.0001 | - | test/unit_tests/firkode/CXX_serial/firk_test_order_pr_3_-1e4.out |
| firk_test_order_pr_3_-1e4_h2 | yes | yes | 0.0001 | - | test/unit_tests/firkode/CXX_serial/firk_test_order_pr_3_-1e4.out |
| firk_test_order_pr_3_-1e4_h3 | yes | yes | 0.0001 | - | test/unit_tests/firkode/CXX_serial/firk_test_order_pr_3_-1e4.out; it prints 'error < 1e4 u' for this step size |
| firk_test_vdp_1e-6 | yes | yes | answer | - | test/unit_tests/firkode/CXX_serial/firk_test_vdp.cpp: step counts of its former answer file firk_test_vdp.out |
| firk_test_vdp_1e-10 | yes | yes | answer | - | test/unit_tests/firkode/CXX_serial/firk_test_vdp.cpp: step counts of its former answer file firk_test_vdp.out; 9340 steps on the Linux workstation of the sweeps, 9345 in that file |
| firk_test_vdp_1e-4 | yes | yes | answer | - | test/unit_tests/firkode/CXX_serial/firk_test_vdp.cpp: step counts of its former answer file firk_test_vdp.out (the loose-tolerance run of benchmark issue 2) |
| ark_analytic | yes | yes | answer | - | examples/arkode/C_serial/ark_analytic.out (ARKodeSetLinear as in the example) |
| ark_analytic_table_name | yes | yes | answer | - | examples/arkode/C_serial/ark_analytic_arkode.scalar_tolerances_1e-6_1e-8_arkode.table_names_ARKODE_ESDIRK547L2SA_7_4_5_ARKODE_ERK_NONE.out (a table chosen by name) |
| test_logging_arkode_arkstep_1_1_1 | yes | yes | answer | - | test/unit_tests/logging/test_logging_arkode_arkstep_1_1_1.out (three one-step returns towards t = 1) |
| test_logging_arkode_arkstep_2_1_1 | yes | yes | answer | - | test/unit_tests/logging/test_logging_arkode_arkstep_2_1_1.out (ImEx, Newton, dense; the splitting of problems/kpr.hpp) |
| ark_analytic_lsrk | yes | yes | 0.1 | - | examples/arkode/C_serial/ark_analytic_lsrk.out |
| ark_kpr_mri_7_2_0.002 | yes | yes | answer | - | examples/arkode/C_serial/ark_kpr_mri_7_2_0.002.out (slow ARKODE_MRI_GARK_IRK21a, fast ARKODE_HEUN_EULER_2_1_2, fixed steps) |
| ark_kpr_mri_11_2_0.001 | yes | yes | answer | - | examples/arkode/C_serial/ark_kpr_mri_11_2_0.001.out (slow ARKODE_IMEX_MRI_SR21, fast ARKODE_HEUN_EULER_2_1_2, fixed steps) |
| ark_kpr_mri_3_2_0.001 | yes | yes | answer | - | examples/arkode/C_serial/ark_kpr_mri_3_2_0.001.out (slow ARKODE_MERK21, fast ARKODE_HEUN_EULER_2_1_2, fixed steps) |
| check_split_kpr_imex | yes | yes | answer | - | --check-split: |f - (fe + fi + ff)| relative to the pieces' norms, and the implicit Jacobian against central differences of fi |
| check_split_bruss3d_imex | yes | yes | answer | - | --check-split: |f - (fe + fi + ff)| relative to the pieces' norms, and the implicit Jacobian against central differences of fi |
| check_split_heat2d_imex | yes | yes | answer | - | --check-split: |f - (fe + fi + ff)| relative to the pieces' norms, and the implicit Jacobian against central differences of fi |
| check_split_heat2d_mri | yes | yes | answer | - | --check-split: |f - (fe + fi + ff)| relative to the pieces' norms, and the implicit Jacobian against central differences of fi |
| check_split_heat2d_mri-imex | yes | yes | answer | - | --check-split: |f - (fe + fi + ff)| relative to the pieces' norms, and the implicit Jacobian against central differences of fi |
| check_split_diurnal_imex | yes | yes | answer | - | --check-split: |f - (fe + fi + ff)| relative to the pieces' norms, and the implicit Jacobian against central differences of fi |
| check_split_kprmri_mri | yes | yes | answer | - | --check-split: |f - (fe + fi + ff)| relative to the pieces' norms, and the implicit Jacobian against central differences of fi |
| check_split_kprmri_mri-imex | yes | yes | answer | - | --check-split: |f - (fe + fi + ff)| relative to the pieces' norms, and the implicit Jacobian against central differences of fi |
| check_split_kprmri_mri-explicit | yes | yes | answer | - | --check-split: |f - (fe + fi + ff)| relative to the pieces' norms, and the implicit Jacobian against central differences of fi |
| check_domeig_analytic_lambda-1e+06 | yes | yes | answer | - | --check-domeig: the problem's eigenvalue bound over the power iteration's estimate at (t0, y0) |
| check_domeig_pr_lambda-10000 | yes | yes | answer | - | --check-domeig: the problem's eigenvalue bound over the power iteration's estimate at (t0, y0) |
| check_domeig_sector_m10_theta0_rmin1_rmax100 | yes | yes | answer | - | --check-domeig: the problem's eigenvalue bound over the power iteration's estimate at (t0, y0) |
| check_domeig_heat2d_mx12_my10 | yes | yes | answer | - | --check-domeig: the problem's eigenvalue bound over the power iteration's estimate at (t0, y0) |
| check_domeig_heat2d_mx12_my10_velocity0.5 | yes | yes | answer | - | --check-domeig: the problem's eigenvalue bound over the power iteration's estimate at (t0, y0) |
| check_domeig_bruss3d_n6 | yes | yes | answer | - | --check-domeig: the problem's eigenvalue bound over the power iteration's estimate at (t0, y0) |
| cv_kpr | yes | yes | 0.1 | - | examples/cvode/CXX_serial/cv_kpr.out |
| idaRoberts_dns | yes | yes | answer | - | examples/ida/serial/idaRoberts_dns.out (the example poses the DAE with the conservation constraint; the harness poses the ODE in residual form and reproduces every counter) |
| cvRoberts_dns | yes | yes | 0.1 | - | examples/cvode/serial/cvRoberts_dns.out computes the second component as -yd1 - yd3, a different floating-point expression; informational only |
<!-- end generated -->

## 4. Problems

| name | system | N | interval | families |
|---|---|---|---|---|
| ROBER, VDPOL, OREGO, HIRES | the small stiff problems of `BENCHMARK_BDF.md` (B1) | 3, 2, 3, 8 | README | A1, A4 |
| `kpr` | Kvaerno-Prothero-Robinson of `test/unit_tests/problems/kpr.hpp` with its IMEX splitting | 2 | [0, 0.5] | A0 |
| `pr` | Prothero-Robinson `y' = lambda (y - phi) + phi'`, exact solution | 4, 1000 | [0, 1] | A2 |
| `sector` | `y' = A (y - phi) + phi'`, eigenvalues `r_k (-cos theta +- i sin theta)` (B8) | 20, 200 | [0, 1], [0, 10] | A0, A2, A4, A6 |
| `heat1dfem`, `wave1dfem` | P1 finite elements with a mass matrix (B5, B9) | 101-401, 200-800 | [0, 1], [0, 4] | A3 |
| `bruss3d` | 3-D Brusselator; IMEX: explicit reaction, implicit linear diffusion | 2 n^3 | [0, 10] | A4, A5 |
| `heat2d` | 2-D advection-diffusion, 200 x 100; IMEX: explicit advection and forcing, implicit diffusion; MRI: slow implicit diffusion, fast advection | 20000 | [0, 1] | A6 |
| `kprmri` | multirate Kvaerno-Prothero-Robinson of `ark_kpr_mri.c` (slow stiffness G = -100, separation w, coupling e = 0.5) | 2 | [0, 5] | A0, A7 |
| `diurnal` | 2-species diurnal kinetics; IMEX: explicit transport, implicit kinetics | 2 mx my | [0, 86400] | A8 |

## 5. Results

### 5.1 A0: method selection (gate)

**Purpose.** Show that every table, coupling and method the families use runs with its
order, so that the labels can be trusted.

**Cells.** Fixed steps with a tightly converged Newton iteration and exact solutions: the
seven DIRK tables and the three IMEX pairs on `kpr` (h = 0.1 / 2^k, k = 0..5, t in [0, 0.5]);
MRI-GARK ESDIRK34a and IMEX-MRI-SR32 on the non-stiff `kprmri` (G = -1, slow step
h = 0.05 / 2^k, fast step h / 100 with Dormand-Prince 5(4)); RKC2 and RKL2 on SECTOR at
theta = 0 (r in [1, 100]). 80 cells, untimed.

**Rule [predicted].** A method passes when its observed order on the finest admissible pair
is at least its order minus 0.3; an observed order above the order plus 0.3 is marked, not
failed (error cancellation on `kpr`: ESDIRK436L2SA showed order 5.06 there in the 2026-10-05
suite [measured, older harness]). A failure voids every label that uses the method.

<!-- begin generated: A0_orders -->
| problem | params | configuration | method | expected order | errors at h0, h0/2, ... | observed orders | reported | at least expected - 0.3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| kpr |  | ARKODE DIRK3 | ARKODE_ESDIRK325L2SA_5_2_3 | 3 | 1.18e-02, 2.40e-04, 4.51e-05, 5.26e-06, 6.12e-07, 7.29e-08 | -, 2.41, 3.10, 3.10, 3.07 | 3.07 | yes |
| kpr |  | ARKODE DIRK4 | ARKODE_ESDIRK436L2SA_6_3_4 | 4 | 1.14e-02, 3.47e-04, 1.29e-05, 4.23e-07, 1.33e-08, 4.12e-10 | -, 4.75, 4.93, 4.99, 5.02 | 5.02 | yes (above by 1.02) |
| kpr |  | ARKODE DIRK5 | ARKODE_ESDIRK547L2SA2_7_4_5 | 5 | 2.84e-03, 7.58e-05, 1.27e-06, 2.93e-08, 7.93e-10, 2.32e-11 | 5.23, 5.90, 5.43, 5.21, 5.09 | 5.09 | yes |
| kpr |  | ARKODE DIRK SDIRK_5_3_4 | ARKODE_SDIRK_5_3_4 | 4 | 7.32e-03, 4.46e-05, 2.66e-06, 3.13e-07, 2.40e-08, 1.64e-09 | 7.36, 4.07, 3.09, 3.70, 3.87 | 3.87 | yes |
| kpr |  | ARKODE DIRK ESDIRK437L2SA_7_3_4 | ARKODE_ESDIRK437L2SA_7_3_4 | 4 | 4.76e-03, 1.80e-04, 4.70e-06, 1.50e-07, 4.94e-09, 1.70e-10 | 4.73, 5.26, 4.97, 4.93, 4.86 | 4.86 | yes (above by 0.86) |
| kpr |  | ARKODE DIRK KVAERNO_7_4_5 | ARKODE_KVAERNO_7_4_5 | 5 | 2.20e-02, 2.47e-04, 9.79e-06, 3.12e-07, 9.81e-09, 3.07e-10 | -, 4.66, 4.97, 4.99, 5.00 | 5.00 | yes |
| kpr |  | ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | ARKODE_ARK548L2SAb_DIRK_8_4_5 | 5 | 9.02e-04, 1.35e-04, 4.99e-06, 1.62e-07, 5.08e-09, 1.59e-10 | 2.74, 4.76, 4.95, 4.99, 5.00 | 5.00 | yes |
| kpr |  | ARKODE IMEX3 | ARKODE_ARK324L2SA_DIRK_4_2_3/ARKODE_ARK324L2SA_ERK_4_2_3 | 3 | 5.37e-02, 1.30e-02, 1.55e-03, 1.95e-04, 2.44e-05, 3.05e-06 | -, -, 3.00, 3.00, 3.00 | 3.00 | yes |
| kpr |  | ARKODE IMEX4 | ARKODE_ARK437L2SA_DIRK_7_3_4/ARKODE_ARK437L2SA_ERK_7_3_4 | 4 | 6.05e-03, 4.62e-04, 1.04e-05, 2.50e-07, 3.27e-09, 1.80e-10 | 3.71, 5.48, 5.38, 6.26, 4.18 | 4.18 | yes |
| kpr |  | ARKODE IMEX5 | ARKODE_ARK548L2SAb_DIRK_8_4_5/ARKODE_ARK548L2SAb_ERK_8_4_5 | 5 | 4.02e-03, 2.79e-04, 9.82e-06, 3.14e-07, 9.86e-09, 3.09e-10 | 3.85, 4.83, 4.97, 4.99, 5.00 | 5.00 | yes |
| kprmri | G=-1.0 | ARKODE MRI ESDIRK34a | ARKODE_MRI_GARK_ESDIRK34a/ARKODE_DORMAND_PRINCE_7_4_5 | 3 | 7.42e-07, 9.40e-08, 1.18e-08, 1.48e-09, 1.86e-10 | 2.98, 2.99, 3.00, 3.00 | 3.00 | yes |
| kprmri | G=-1.0 | ARKODE MRI IMEX-SR32 | ARKODE_IMEX_MRI_SR32/ARKODE_DORMAND_PRINCE_7_4_5 | 3 | 2.10e-07, 2.69e-08, 3.40e-09, 4.28e-10, 5.38e-11 | 2.97, 2.98, 2.99, 2.99 | 2.99 | yes |
| sector | m=10;rmax=100.0;rmin=1.0;theta=0.0 | ARKODE LSRK RKC2 | ARKODE_LSRK_RKC_2+user | 2 | 5.65e-03, 1.71e-03, 4.18e-04, 1.19e-04, 2.95e-05 | 1.73, 2.03, 1.82, 2.01 | 2.01 | yes |
| sector | m=10;rmax=100.0;rmin=1.0;theta=0.0 | ARKODE LSRK RKL2 | ARKODE_LSRK_RKL_2+user | 2 | 5.14e-03, 1.37e-03, 3.90e-04, 1.03e-04, 2.36e-05 | 1.90, 1.82, 1.92, 2.12 | 2.12 | yes |
<!-- end generated -->

![errors of the fixed-step runs](../../doc/shared/figs/firkode/benchmarks/arkode/arkode_order_selection.png)

**Verdict.** [pending]

### 5.2 A1: small stiff problems (candidate 1, steps against per-step work)

**Purpose.** On ROBER, VDPOL, OREGO and HIRES (N <= 8), is FIRKODE faster than ARKODE's
DIRK methods, and is the reason the step count?

**Sweep.** rtol 1e-2 ... 1e-10 in half decades (B1's grid), FIRKODE s = 3 and the nine DIRK
configurations; the ark_robertson settings on ROBER (tuned, reported apart). 697 cells.

**Predictions [predicted].**
- P1.1: t(best default) / t(FIRKODE) > 1.2 at mescd 6 and 8 on ROBER, VDPOL and OREGO;
  HIRES within [0.5, 1.5].
- P1.2: where P1.1 holds, the DIRK step count at rtol 1e-6 is >= 1.5x FIRKODE's while the
  time per step stays within [0.5, 2]x (the win comes from the steps).
- P1.3: the optimal parameter sets change the best default by <= 1.5x and flip no verdict.
- P1.4: the best of all changes the P1.1 ratios by <= 1.5x; SDIRK_5_3_4 (stage order 1) is
  the slowest order-4 table on ROBER and VDPOL.
- P1.5: the tuned settings beat FIRKODE on ROBER at mescd 4-6.

**Pass and fail.** Supports candidate 1 if P1.1 holds for the best of all on at least 3 of
the 4 problems; refutes it if the best of all is >= 1.2x faster than FIRKODE on at least 2
problems at mescd 6.

**Results [pending].**

<!-- begin generated: A1_ratio -->
| problem | configurations | t / t(FIRKODE) at mescd=4 | t / t(FIRKODE) at mescd=6 | t / t(FIRKODE) at mescd=8 |
| --- | --- | --- | --- | --- |
| rober | best DIRK default | 3.48 (ARKODE DIRK3) | 4.28 (ARKODE DIRK3) | 8.56 (ARKODE DIRK3) |
| rober | best DIRK of all | 3.48 (ARKODE DIRK3) | 1.17 (ARKODE DIRK KVAERNO_7_4_5) | 2.4 (ARKODE DIRK KVAERNO_7_4_5) |
| rober | ARKODE DIRK4 [ark_robertson settings] | NA | 0.727 (ARKODE DIRK4 [ark_robertson settings]) | 2.05 (ARKODE DIRK4 [ark_robertson settings]) |
| vdpol | best DIRK default | 1.53 (ARKODE DIRK3) | 1.61 (ARKODE DIRK3) | 2.29 (ARKODE DIRK4) |
| vdpol | best DIRK of all | 1.53 (ARKODE DIRK3) | 1.61 (ARKODE DIRK3) | 2.29 (ARKODE DIRK4) |
| orego | best DIRK default | 1.11 (ARKODE DIRK3) | 3.42 (ARKODE DIRK3) | 4.96 (ARKODE DIRK3) |
| orego | best DIRK of all | 1.11 (ARKODE DIRK3) | 3.41 (ARKODE DIRK KVAERNO_7_4_5) | 3.66 (ARKODE DIRK KVAERNO_7_4_5) |
| hires | best DIRK default | 0.394 (ARKODE DIRK4) | 0.922 (ARKODE DIRK3) | 1.63 (ARKODE DIRK4) |
| hires | best DIRK of all | 0.364 (ARKODE DIRK KVAERNO_7_4_5) | 0.922 (ARKODE DIRK3) | 1.48 (ARKODE DIRK5 optimal) |
<!-- end generated -->

<!-- begin generated: A1_counters -->
| problem | params | method | mescd | steps | attempts | Newton | setups | stages | Newton / step | us / step |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| hires |  | FIRKODE s=3 | 6.96 | 79 | 84 | 204 | 61 | 3 | 2.58 | 9.59 |
| hires |  | ARKODE DIRK3 | 5.32 | 102 | 124 | 1225 | 138 | 5 | 12 | 2.58 |
| hires |  | ARKODE DIRK4 | 4.88 | 86 | 114 | 1310 | 125 | 6 | 15.2 | 3.27 |
| hires |  | ARKODE DIRK4 optimal | 5 | 97 | 126 | 1441 | 122 | 6 | 14.9 | 3.16 |
| hires |  | ARKODE DIRK5 | 5.6 | 100 | 137 | 1928 | 139 | 7 | 19.3 | 4.03 |
| hires |  | ARKODE DIRK5 optimal | 4.7 | 126 | 139 | 2018 | 104 | 7 | 16 | 3.2 |
| hires |  | ARKODE DIRK SDIRK_5_3_4 | 3.77 | 178 | 215 | 2534 | 163 | 5 | 14.2 | 2.96 |
| hires |  | ARKODE DIRK ESDIRK437L2SA_7_3_4 | 5.37 | 140 | 203 | 2447 | 168 | 7 | 17.5 | 3.66 |
| hires |  | ARKODE DIRK KVAERNO_7_4_5 | 5.72 | 131 | 179 | 2215 | 223 | 7 | 16.9 | 3.57 |
| hires |  | ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 5.18 | 97 | 126 | 1868 | 150 | 8 | 19.3 | 4.02 |
| orego |  | FIRKODE s=3 | 7.3 | 916 | 947 | 2297 | 747 | 3 | 2.51 | 4.01 |
| orego |  | ARKODE DIRK3 | 5.74 | 1814 | 2365 | 25814 | 2245 | 5 | 14.2 | 1.69 |
| orego |  | ARKODE DIRK4 | 4.87 | 2300 | 3555 | 42504 | 3400 | 6 | 18.5 | 2.09 |
| orego |  | ARKODE DIRK4 optimal | 2.81 | 2811 | 3982 | 45871 | 3326 | 6 | 16.3 | 1.97 |
| orego |  | ARKODE DIRK5 | 3.3 | 4201 | 6296 | 82877 | 4509 | 7 | 19.7 | 2.35 |
| orego |  | ARKODE DIRK5 optimal | 1.4 | 2432 | 2712 | 40989 | 2110 | 7 | 16.9 | 1.97 |
| orego |  | ARKODE DIRK SDIRK_5_3_4 | 2.71 | 5251 | 7039 | 82555 | 4444 | 5 | 15.7 | 1.8 |
| orego |  | ARKODE DIRK ESDIRK437L2SA_7_3_4 | 2.27 | 4729 | 7021 | 89279 | 5185 | 7 | 18.9 | 2.26 |
| orego |  | ARKODE DIRK KVAERNO_7_4_5 | 3.78 | 2066 | 2909 | 38857 | 3394 | 7 | 18.8 | 2.24 |
| orego |  | ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 2.48 | 3096 | 4754 | 73291 | 4229 | 8 | 23.7 | 2.77 |
| rober |  | FIRKODE s=3 | 8.8 | 443 | 444 | 967 | 426 | 3 | 2.18 | 4.22 |
| rober |  | ARKODE DIRK3 | 6.48 | 2710 | 4169 | 38999 | 2951 | 5 | 14.4 | 1.7 |
| rober |  | ARKODE DIRK4 | 6.16 | 6476 | 11332 | 111762 | 7502 | 6 | 17.3 | 2.11 |
| rober |  | ARKODE DIRK4 [ark_robertson settings] | 6.31 | 271 | 285 | 9183 | 240 | 6 | 33.9 | 3.5 |
| rober |  | ARKODE DIRK4 optimal | 6.07 | 9757 | 15581 | 150579 | 11100 | 6 | 15.4 | 1.91 |
| rober |  | ARKODE DIRK5 | 6.29 | 17765 | 27983 | 307842 | 18774 | 7 | 17.3 | 2.17 |
| rober |  | ARKODE DIRK5 optimal | 6.18 | 17518 | 22047 | 238542 | 13192 | 7 | 13.6 | 1.72 |
| rober |  | ARKODE DIRK SDIRK_5_3_4 | 6.57 | 14694 | 22872 | 233079 | 15722 | 5 | 15.9 | 1.88 |
| rober |  | ARKODE DIRK ESDIRK437L2SA_7_3_4 | 6.68 | 15076 | 24319 | 272067 | 15688 | 7 | 18 | 2.2 |
| rober |  | ARKODE DIRK KVAERNO_7_4_5 | 6.63 | 526 | 724 | 9855 | 880 | 7 | 18.7 | 2.19 |
| rober |  | ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 6.15 | 15148 | 25871 | 320089 | 17869 | 8 | 21.1 | 2.63 |
| vdpol |  | FIRKODE s=3 | 6.5 | 867 | 876 | 2005 | 700 | 3 | 2.31 | 2.68 |
| vdpol |  | ARKODE DIRK3 | 6.9 | 1782 | 2382 | 26528 | 2137 | 5 | 14.9 | 1.53 |
| vdpol |  | ARKODE DIRK4 | 5.32 | 2000 | 3068 | 39655 | 2733 | 6 | 19.8 | 2.04 |
| vdpol |  | ARKODE DIRK4 optimal | 4.33 | 5464 | 8158 | 92112 | 6359 | 6 | 16.9 | 1.8 |
| vdpol |  | ARKODE DIRK5 | 4.15 | 12418 | 18878 | 232460 | 11349 | 7 | 18.7 | 1.99 |
| vdpol |  | ARKODE DIRK5 optimal | 2.86 | 13293 | 16153 | 172341 | 8781 | 7 | 13 | 1.5 |
| vdpol |  | ARKODE DIRK SDIRK_5_3_4 | 4.88 | 12147 | 19088 | 201321 | 12713 | 5 | 16.6 | 1.76 |
| vdpol |  | ARKODE DIRK ESDIRK437L2SA_7_3_4 | 4.38 | 13284 | 20342 | 245302 | 12802 | 7 | 18.5 | 1.99 |
| vdpol |  | ARKODE DIRK KVAERNO_7_4_5 | 6.01 | 1665 | 2264 | 31650 | 2599 | 7 | 19 | 1.95 |
| vdpol |  | ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 4.38 | 3948 | 6764 | 112105 | 5439 | 8 | 28.4 | 2.91 |
<!-- end generated -->

![work-precision of the small problems](../../doc/shared/figs/firkode/benchmarks/arkode/arkode_small_wp.png)

<!-- begin generated: A1_wp -->
| method | rtol | mescd | steps | attempts | Newton | Jac | setups | err fails | NLS fails | problem | time (ms) | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| FIRKODE s=3 | 1e-04 | 5.66 | 309 | 341 | 888 | 268 | 311 | 4 | 28 | orego | 1.49 | ok |
| FIRKODE s=3 | 1e-04 | 5.3 | 295 | 329 | 751 | 189 | 290 | 8 | 26 | vdpol | <1 | ok |
| FIRKODE s=3 | 1e-04 | 5.79 | 158 | 160 | 346 | 152 | 160 | 1 | 1 | rober | <1 | ok |
| FIRKODE s=3 | 1e-04 | 3.7 | 37 | 46 | 106 | 22 | 43 | 0 | 9 | hires | <1 | ok |
| FIRKODE s=3 | 1e-06 | 7.3 | 916 | 947 | 2297 | 592 | 747 | 4 | 27 | orego | 3.67 | ok |
| FIRKODE s=3 | 1e-06 | 6.5 | 867 | 876 | 2005 | 403 | 700 | 7 | 2 | vdpol | 2.32 | ok |
| FIRKODE s=3 | 1e-06 | 8.8 | 443 | 444 | 967 | 411 | 426 | 1 | 0 | rober | 1.87 | ok |
| FIRKODE s=3 | 1e-06 | 6.96 | 79 | 84 | 204 | 32 | 61 | 0 | 5 | hires | <1 | ok |
| FIRKODE s=3 | 1e-08 | 11.7* | 1393 | 1395 | 3046 | 499 | 545 | 2 | 0 | rober | 5.63 | ok |
| FIRKODE s=3 | 1e-08 | 10.9* | 2752 | 2754 | 6037 | 652 | 2065 | 2 | 0 | vdpol | 6.71 | ok |
| FIRKODE s=3 | 1e-08 | 9.48 | 236 | 237 | 558 | 70 | 105 | 1 | 0 | hires | 1.9 | ok |
| FIRKODE s=3 | 1e-08 | 9.71 | 2893 | 2915 | 6576 | 1020 | 1943 | 2 | 20 | orego | 9.99 | ok |
| FIRKODE s=3 | 1e-10 | 12.9* | 9195 | 9198 | 19689 | 1695 | 5556 | 2 | 1 | orego | 28.9 | ok |
| FIRKODE s=3 | 1e-10 | 13.2* | 4416 | 4417 | 9565 | 1016 | 1205 | 1 | 0 | rober | 17.5 | ok |
| FIRKODE s=3 | 1e-10 | 13.3* | 8774 | 8775 | 18708 | 1371 | 6364 | 0 | 1 | vdpol | 20.3 | ok |
| FIRKODE s=3 | 1e-10 | 11.7* | 728 | 728 | 1553 | 152 | 197 | 0 | 0 | hires | 4.94 | ok |
| ARKODE DIRK3 | 1e-04 | 4.62 | 567 | 830 | 7682 | 68 | 597 | 263 | 63 | rober | <1 | ok |
| ARKODE DIRK3 | 1e-04 | 4.45 | 642 | 970 | 10146 | 391 | 839 | 216 | 390 | vdpol | 1.07 | ok |
| ARKODE DIRK3 | 1e-04 | 2.87 | 640 | 928 | 9376 | 570 | 967 | 121 | 568 | orego | 1.14 | ok |
| ARKODE DIRK3 | 1e-04 | 3.36 | 32 | 40 | 337 | 23 | 47 | 1 | 22 | hires | <1 | ok |
| ARKODE DIRK3 | 1e-06 | 6.9 | 1782 | 2382 | 26528 | 1309 | 2137 | 290 | 1306 | vdpol | 2.72 | ok |
| ARKODE DIRK3 | 1e-06 | 6.48 | 2710 | 4169 | 38999 | 314 | 2951 | 1459 | 290 | rober | 4.61 | ok |
| ARKODE DIRK3 | 1e-06 | 5.32 | 102 | 124 | 1225 | 105 | 138 | 1 | 104 | hires | <1 | ok |
| ARKODE DIRK3 | 1e-06 | 5.74 | 1814 | 2365 | 25814 | 1453 | 2245 | 231 | 1451 | orego | 3.06 | ok |
| ARKODE DIRK3 | 1e-08 | 8.19 | 6519 | 8016 | 90121 | 4121 | 5840 | 418 | 4106 | vdpol | 8.98 | ok |
| ARKODE DIRK3 | 1e-08 | 7.15 | 405 | 491 | 5031 | 382 | 481 | 2 | 381 | hires | 1.06 | ok |
| ARKODE DIRK3 | 1e-08 | 9.1 | 14950 | 25994 | 246644 | 1491 | 19164 | 11043 | 1272 | rober | 29.2 | ok |
| ARKODE DIRK3 | 1e-08 | 6.85 | 6807 | 8358 | 94552 | 4753 | 6680 | 486 | 4745 | orego | 11.1 | ok |
| ARKODE DIRK3 | 1e-10 | 10.8* | 43676 | 77319 | 766537 | 7316 | 57781 | 32610 | 6678 | rober | 88.9 | ok |
| ARKODE DIRK3 | 1e-10 | 10.2* | 26537 | 30250 | 339330 | 13328 | 17330 | 4 | 13247 | vdpol | 35.6 | ok |
| ARKODE DIRK3 | 1e-10 | 9.3 | 1604 | 1850 | 21738 | 1387 | 1642 | 2 | 1386 | hires | 4.47 | ok |
| ARKODE DIRK3 | 1e-10 | 9.24 | 28210 | 33660 | 382897 | 14500 | 21625 | 2701 | 14426 | orego | 45 | ok |
| ARKODE DIRK4 | 1e-04 | 4.83 | 3365 | 5879 | 50753 | 105 | 3972 | 2513 | 43 | rober | 6.46 | ok |
| ARKODE DIRK4 | 1e-04 | 1.2 | 843 | 1215 | 13672 | 657 | 1195 | 173 | 656 | orego | 1.66 | ok |
| ARKODE DIRK4 | 1e-04 | 3.6 | 32 | 42 | 422 | 20 | 44 | 4 | 19 | hires | <1 | ok |
| ARKODE DIRK4 | 1e-04 | 3.7 | 2146 | 3303 | 33394 | 605 | 2316 | 975 | 576 | vdpol | 3.65 | ok |
| ARKODE DIRK4 | 1e-06 | 4.87 | 2300 | 3555 | 42504 | 1490 | 3400 | 890 | 1486 | orego | 4.8 | ok |
| ARKODE DIRK4 | 1e-06 | 6.16 | 6476 | 11332 | 111762 | 413 | 7502 | 4835 | 311 | rober | 13.6 | ok |
| ARKODE DIRK4 | 1e-06 | 5.32 | 2000 | 3068 | 39655 | 996 | 2733 | 859 | 995 | vdpol | 4.09 | ok |
| ARKODE DIRK4 | 1e-06 | 4.88 | 86 | 114 | 1310 | 79 | 125 | 8 | 78 | hires | <1 | ok |
| ARKODE DIRK4 | 1e-08 | 6.96 | 324 | 451 | 5403 | 351 | 508 | 28 | 350 | hires | 1.13 | ok |
| ARKODE DIRK4 | 1e-08 | 7.95 | 8056 | 15084 | 173570 | 1608 | 12000 | 6903 | 1570 | rober | 20.3 | ok |
| ARKODE DIRK4 | 1e-08 | 8.81 | 4295 | 5493 | 75572 | 3563 | 5270 | 272 | 3560 | vdpol | 7.73 | ok |
| ARKODE DIRK4 | 1e-08 | 7.52 | 5667 | 8319 | 104699 | 4072 | 8068 | 1713 | 4068 | orego | 12.3 | ok |
| ARKODE DIRK4 | 1e-10 | 9.8 | 16048 | 26640 | 332463 | 5479 | 22124 | 9840 | 5408 | rober | 38.7 | ok |
| ARKODE DIRK4 | 1e-10 | 11* | 17041 | 22691 | 303532 | 14587 | 21816 | 912 | 14580 | vdpol | 30 | ok |
| ARKODE DIRK4 | 1e-10 | 9.07 | 879 | 1157 | 13416 | 916 | 1264 | 8 | 915 | hires | 2.8 | ok |
| ARKODE DIRK4 | 1e-10 | 8.21 | 13710 | 17717 | 238178 | 11795 | 17497 | 1260 | 11794 | orego | 27.8 | ok |
| ARKODE DIRK4 [ark_robertson settings] | 1e-04 | 5.48 | 101 | 107 | 3381 | 53 | 146 | 1 | 52 | rober | <1 | ok |
| ARKODE DIRK4 [ark_robertson settings] | 1e-06 | 6.31 | 271 | 285 | 9183 | 153 | 240 | 3 | 152 | rober | <1 | ok |
| ARKODE DIRK4 [ark_robertson settings] | 1e-08 | 8.03 | 876 | 908 | 28362 | 282 | 357 | 20 | 280 | rober | 2.94 | ok |
| ARKODE DIRK4 [ark_robertson settings] | 1e-10 | 9.96 | 3592 | 3951 | 101931 | 963 | 1513 | 72 | 961 | rober | 10.9 | ok |
| ARKODE DIRK4 optimal | 1e-04 | 2.78 | 2726 | 3744 | 31859 | 263 | 2694 | 974 | 219 | vdpol | 3.78 | ok |
| ARKODE DIRK4 optimal | 1e-04 | 3.96 | 10486 | 15191 | 111053 | 218 | 11854 | 4703 | 17 | rober | 15.5 | ok |
| ARKODE DIRK4 optimal | 1e-04 | 1.96 | 32 | 39 | 383 | 9 | 37 | 6 | 8 | hires | <1 | ok |
| ARKODE DIRK4 optimal | 1e-04 | 1.83 | 881 | 1193 | 13094 | 319 | 1053 | 258 | 315 | orego | 1.63 | ok |
| ARKODE DIRK4 optimal | 1e-06 | 6.07 | 9757 | 15581 | 150579 | 343 | 11100 | 5822 | 173 | rober | 18.6 | ok |
| ARKODE DIRK4 optimal | 1e-06 | 5 | 97 | 126 | 1441 | 41 | 122 | 25 | 40 | hires | <1 | ok |
| ARKODE DIRK4 optimal | 1e-06 | 2.81 | 2811 | 3982 | 45871 | 858 | 3326 | 1071 | 839 | orego | 5.54 | ok |
| ARKODE DIRK4 optimal | 1e-06 | 4.33 | 5464 | 8158 | 92112 | 652 | 6359 | 2636 | 579 | vdpol | 9.82 | ok |
| ARKODE DIRK4 optimal | 1e-08 | 7.6 | 11715 | 19042 | 209907 | 2039 | 16036 | 7314 | 1944 | rober | 24.8 | ok |
| ARKODE DIRK4 optimal | 1e-08 | 6.81 | 311 | 415 | 5006 | 164 | 411 | 90 | 163 | hires | 1.04 | ok |
| ARKODE DIRK4 optimal | 1e-08 | 7.45 | 5021 | 6387 | 85272 | 2401 | 5869 | 1236 | 2398 | vdpol | 8.82 | ok |
| ARKODE DIRK4 optimal | 1e-08 | 6.54 | 7533 | 11129 | 133626 | 2681 | 10022 | 3354 | 2639 | orego | 16 | ok |
| ARKODE DIRK4 optimal | 1e-10 | 8.21 | 16264 | 20697 | 271141 | 8696 | 20458 | 3731 | 8695 | orego | 30.8 | ok |
| ARKODE DIRK4 optimal | 1e-10 | 8.98 | 17701 | 25687 | 309832 | 6002 | 24098 | 7900 | 6001 | rober | 36.7 | ok |
| ARKODE DIRK4 optimal | 1e-10 | 9.07 | 9719 | 10563 | 165408 | 7704 | 10485 | 173 | 7697 | vdpol | 17 | ok |
| ARKODE DIRK4 optimal | 1e-10 | 9.19 | 683 | 800 | 10957 | 555 | 933 | 59 | 554 | hires | 2.26 | ok |
| ARKODE DIRK5 | 1e-04 | 4.46 | 36382 | 70855 | 597371 | 722 | 52717 | 34470 | 17 | rober | 78.3 | ok |
| ARKODE DIRK5 | 1e-04 | 3.14 | 1570 | 2240 | 27369 | 489 | 1664 | 554 | 480 | orego | 3.34 | ok |
| ARKODE DIRK5 | 1e-04 | 3.08 | 31 | 37 | 450 | 14 | 39 | 4 | 13 | hires | <1 | ok |
| ARKODE DIRK5 | 1e-04 | 2.59 | 8867 | 14500 | 126159 | 555 | 10292 | 5516 | 392 | vdpol | 15.4 | ok |
| ARKODE DIRK5 | 1e-06 | 6.29 | 17765 | 27983 | 307842 | 664 | 18774 | 10213 | 384 | rober | 38.5 | ok |
| ARKODE DIRK5 | 1e-06 | 4.15 | 12418 | 18878 | 232460 | 909 | 11349 | 6346 | 696 | vdpol | 24.8 | ok |
| ARKODE DIRK5 | 1e-06 | 3.3 | 4201 | 6296 | 82877 | 1187 | 4509 | 1902 | 1154 | orego | 9.85 | ok |
| ARKODE DIRK5 | 1e-06 | 5.6 | 100 | 137 | 1928 | 74 | 139 | 24 | 73 | hires | <1 | ok |
| ARKODE DIRK5 | 1e-08 | 7.98 | 34913 | 58315 | 698728 | 3866 | 40513 | 23390 | 3498 | rober | 83.1 | ok |
| ARKODE DIRK5 | 1e-08 | 6.69 | 9861 | 16415 | 225400 | 3005 | 13467 | 6270 | 2945 | orego | 26 | ok |
| ARKODE DIRK5 | 1e-08 | 7.01 | 406 | 645 | 8912 | 232 | 615 | 197 | 231 | hires | 1.82 | ok |
| ARKODE DIRK5 | 1e-08 | 7.38 | 11570 | 18053 | 275638 | 2810 | 14254 | 6296 | 2809 | vdpol | 28.4 | ok |
| ARKODE DIRK5 | 1e-10 | 8.9 | 11317 | 13516 | 248107 | 7602 | 11558 | 1674 | 7601 | vdpol | 25.4 | ok |
| ARKODE DIRK5 | 1e-10 | 8.91 | 595 | 751 | 11958 | 439 | 723 | 126 | 438 | hires | 2.44 | ok |
| ARKODE DIRK5 | 1e-10 | 10.2 | 48792 | 83358 | 1151069 | 8308 | 66167 | 34513 | 8049 | rober | 133 | ok |
| ARKODE DIRK5 | 1e-10 | 7.9 | 13204 | 17180 | 276142 | 8614 | 15784 | 3337 | 8612 | orego | 31.9 | ok |
| ARKODE DIRK5 optimal | 1e-04 | 2.33 | 9965 | 11462 | 90058 | 327 | 6001 | 1481 | 142 | vdpol | 11.6 | ok |
| ARKODE DIRK5 optimal | 1e-04 | 2.12 | 126 | 136 | 1210 | 8 | 45 | 10 | 6 | hires | <1 | ok |
| ARKODE DIRK5 optimal | 1e-04 | 0.822 | 2778 | 3138 | 32926 | 261 | 1246 | 343 | 231 | orego | 4.31 | ok |
| ARKODE DIRK5 optimal | 1e-04 | 3.64 | 30416 | 36452 | 269957 | 607 | 16899 | 6035 | 28 | rober | 39 | ok |
| ARKODE DIRK5 optimal | 1e-06 | 6.18 | 17518 | 22047 | 238542 | 2345 | 13192 | 4528 | 2206 | rober | 30.2 | ok |
| ARKODE DIRK5 optimal | 1e-06 | 4.7 | 126 | 139 | 2018 | 71 | 104 | 13 | 70 | hires | <1 | ok |
| ARKODE DIRK5 optimal | 1e-06 | 2.86 | 13293 | 16153 | 172341 | 1254 | 8781 | 2856 | 1087 | vdpol | 19.9 | ok |
| ARKODE DIRK5 optimal | 1e-06 | 1.4 | 2432 | 2712 | 40989 | 1320 | 2110 | 271 | 1319 | orego | 4.78 | ok |
| ARKODE DIRK5 optimal | 1e-08 | 7.33 | 357 | 405 | 5985 | 230 | 356 | 43 | 228 | hires | 1.21 | ok |
| ARKODE DIRK5 optimal | 1e-08 | 6.25 | 7928 | 8668 | 141169 | 4829 | 7303 | 738 | 4828 | vdpol | 14.5 | ok |
| ARKODE DIRK5 optimal | 1e-08 | 8.24 | 11977 | 13642 | 200532 | 8580 | 12018 | 1663 | 8578 | rober | 23.3 | ok |
| ARKODE DIRK5 optimal | 1e-08 | 6.08 | 6118 | 6653 | 107639 | 4278 | 5622 | 517 | 4272 | orego | 12 | ok |
| ARKODE DIRK5 optimal | 1e-10 | 9.03 | 26592 | 28466 | 478185 | 14521 | 25122 | 1872 | 14518 | vdpol | 46.6 | ok |
| ARKODE DIRK5 optimal | 1e-10 | 9.08 | 649 | 666 | 11556 | 551 | 604 | 12 | 550 | hires | 2.33 | ok |
| ARKODE DIRK5 optimal | 1e-10 | 7.27 | 16982 | 18120 | 307943 | 12334 | 15614 | 1074 | 12323 | orego | 35.4 | ok |
| ARKODE DIRK5 optimal | 1e-10 | 10.4 | 35761 | 39660 | 621963 | 23370 | 37427 | 3896 | 23369 | rober | 71.8 | ok |
| ARKODE DIRK SDIRK_5_3_4 | 1e-04 | 6.25 | 71117 | 128871 | 871679 | 1394 | 92890 | 57753 | 17 | rober | 116 | ok |
| ARKODE DIRK SDIRK_5_3_4 | 1e-04 | 3.06 | 16364 | 24519 | 165423 | 479 | 15971 | 8128 | 172 | vdpol | 20.4 | ok |
| ARKODE DIRK SDIRK_5_3_4 | 1e-04 | 1.72 | 2696 | 3725 | 34692 | 354 | 2251 | 977 | 325 | orego | 4.36 | ok |
| ARKODE DIRK SDIRK_5_3_4 | 1e-04 | 3.27 | 65 | 86 | 828 | 11 | 63 | 20 | 10 | hires | <1 | ok |
| ARKODE DIRK SDIRK_5_3_4 | 1e-06 | 2.71 | 5251 | 7039 | 82555 | 1121 | 4444 | 1778 | 1098 | orego | 9.45 | ok |
| ARKODE DIRK SDIRK_5_3_4 | 1e-06 | 6.57 | 14694 | 22872 | 233079 | 1181 | 15722 | 8177 | 1056 | rober | 27.7 | ok |
| ARKODE DIRK SDIRK_5_3_4 | 1e-06 | 3.77 | 178 | 215 | 2534 | 75 | 163 | 36 | 74 | hires | <1 | ok |
| ARKODE DIRK SDIRK_5_3_4 | 1e-06 | 4.88 | 12147 | 19088 | 201321 | 928 | 12713 | 6928 | 754 | vdpol | 21.3 | ok |
| ARKODE DIRK SDIRK_5_3_4 | 1e-08 | 5.61 | 15336 | 18308 | 240019 | 3942 | 9580 | 2972 | 3908 | orego | 28.1 | ok |
| ARKODE DIRK SDIRK_5_3_4 | 1e-08 | 7.73 | 560 | 667 | 8573 | 293 | 489 | 107 | 291 | hires | 1.74 | ok |
| ARKODE DIRK SDIRK_5_3_4 | 1e-08 | 8.72 | 31824 | 46793 | 546228 | 5352 | 32139 | 14969 | 5144 | rober | 63.4 | ok |
| ARKODE DIRK SDIRK_5_3_4 | 1e-08 | 5.75 | 22649 | 24101 | 338417 | 4323 | 7580 | 1452 | 4203 | vdpol | 33.3 | ok |
| ARKODE DIRK SDIRK_5_3_4 | 1e-10 | 7.09 | 63293 | 69245 | 963995 | 14328 | 27481 | 5951 | 14110 | orego | 111 | ok |
| ARKODE DIRK SDIRK_5_3_4 | 1e-10 | 8.63 | 1992 | 1997 | 29704 | 938 | 957 | 5 | 937 | hires | 5.85 | ok |
| ARKODE DIRK SDIRK_5_3_4 | 1e-10 | 6.79 | 132502 | 137093 | 1945722 | 20434 | 34016 | 4590 | 20050 | vdpol | 201 | ok |
| ARKODE DIRK SDIRK_5_3_4 | 1e-10 | 10.1 | 148865 | 182159 | 2332714 | 39149 | 106031 | 33293 | 38457 | rober | 268 | ok |
| ARKODE DIRK ESDIRK437L2SA_7_3_4 | 1e-04 | 5.34 | 22477 | 45180 | 391290 | 452 | 33889 | 22701 | 21 | rober | 51.8 | ok |
| ARKODE DIRK ESDIRK437L2SA_7_3_4 | 1e-04 | 4.05 | 29 | 36 | 442 | 16 | 39 | 4 | 15 | hires | <1 | ok |
| ARKODE DIRK ESDIRK437L2SA_7_3_4 | 1e-04 | 3.33 | 7723 | 12227 | 108569 | 455 | 8427 | 4413 | 313 | vdpol | 13 | ok |
| ARKODE DIRK ESDIRK437L2SA_7_3_4 | 1e-04 | 3.02 | 1163 | 1658 | 21214 | 462 | 1310 | 372 | 456 | orego | 2.57 | ok |
| ARKODE DIRK ESDIRK437L2SA_7_3_4 | 1e-06 | 6.68 | 15076 | 24319 | 272067 | 490 | 15688 | 9239 | 243 | rober | 33.1 | ok |
| ARKODE DIRK ESDIRK437L2SA_7_3_4 | 1e-06 | 5.37 | 140 | 203 | 2447 | 59 | 168 | 47 | 58 | hires | <1 | ok |
| ARKODE DIRK ESDIRK437L2SA_7_3_4 | 1e-06 | 2.27 | 4729 | 7021 | 89279 | 1403 | 5185 | 1940 | 1360 | orego | 10.7 | ok |
| ARKODE DIRK ESDIRK437L2SA_7_3_4 | 1e-06 | 4.38 | 13284 | 20342 | 245302 | 1161 | 12802 | 6809 | 939 | vdpol | 26.4 | ok |
| ARKODE DIRK ESDIRK437L2SA_7_3_4 | 1e-08 | 7.53 | 371 | 605 | 8094 | 250 | 596 | 166 | 249 | hires | 1.65 | ok |
| ARKODE DIRK ESDIRK437L2SA_7_3_4 | 1e-08 | 8.02 | 38635 | 67919 | 809387 | 1639 | 47721 | 29266 | 973 | rober | 96.8 | ok |
| ARKODE DIRK ESDIRK437L2SA_7_3_4 | 1e-08 | 7.68 | 7357 | 12086 | 175279 | 2526 | 10401 | 4456 | 2525 | vdpol | 17.8 | ok |
| ARKODE DIRK ESDIRK437L2SA_7_3_4 | 1e-08 | 7.08 | 10823 | 18939 | 256511 | 3465 | 16008 | 7440 | 3366 | orego | 28.6 | ok |
| ARKODE DIRK ESDIRK437L2SA_7_3_4 | 1e-10 | 9.23 | 688 | 993 | 14619 | 800 | 1168 | 53 | 799 | hires | 3.05 | ok |
| ARKODE DIRK ESDIRK437L2SA_7_3_4 | 1e-10 | 8.93 | 27854 | 48518 | 667729 | 9544 | 41539 | 18881 | 9331 | orego | 77.6 | ok |
| ARKODE DIRK ESDIRK437L2SA_7_3_4 | 1e-10 | 10.6* | 62876 | 111510 | 1465189 | 7313 | 86400 | 48497 | 6722 | rober | 170 | ok |
| ARKODE DIRK ESDIRK437L2SA_7_3_4 | 1e-10 | 9.86 | 14874 | 20908 | 348932 | 8615 | 18903 | 5350 | 8614 | vdpol | 35.1 | ok |
| ARKODE DIRK KVAERNO_7_4_5 | 1e-04 | 1.68 | 1458 | 2040 | 21931 | 708 | 1661 | 341 | 694 | orego | 2.8 | ok |
| ARKODE DIRK KVAERNO_7_4_5 | 1e-04 | 4.95 | 5275 | 8411 | 74591 | 299 | 6181 | 3080 | 199 | rober | 9.92 | ok |
| ARKODE DIRK KVAERNO_7_4_5 | 1e-04 | 3.87 | 41 | 48 | 562 | 23 | 46 | 1 | 22 | hires | <1 | ok |
| ARKODE DIRK KVAERNO_7_4_5 | 1e-04 | 4.04 | 651 | 897 | 11806 | 577 | 972 | 36 | 576 | vdpol | 1.23 | ok |
| ARKODE DIRK KVAERNO_7_4_5 | 1e-06 | 6.63 | 526 | 724 | 9855 | 577 | 880 | 1 | 576 | rober | 1.15 | ok |
| ARKODE DIRK KVAERNO_7_4_5 | 1e-06 | 5.72 | 131 | 179 | 2215 | 153 | 223 | 1 | 152 | hires | <1 | ok |
| ARKODE DIRK KVAERNO_7_4_5 | 1e-06 | 6.01 | 1665 | 2264 | 31650 | 1644 | 2599 | 1 | 1643 | vdpol | 3.24 | ok |
| ARKODE DIRK KVAERNO_7_4_5 | 1e-06 | 3.78 | 2066 | 2909 | 38857 | 2165 | 3394 | 148 | 2163 | orego | 4.62 | ok |
| ARKODE DIRK KVAERNO_7_4_5 | 1e-08 | 7.19 | 4883 | 6622 | 92532 | 5329 | 7935 | 43 | 5328 | orego | 11 | ok |
| ARKODE DIRK KVAERNO_7_4_5 | 1e-08 | 8.2 | 4256 | 5729 | 82383 | 4195 | 6536 | 0 | 4194 | vdpol | 8.25 | ok |
| ARKODE DIRK KVAERNO_7_4_5 | 1e-08 | 7.9 | 1425 | 1837 | 27753 | 1402 | 2361 | 0 | 1401 | rober | 3.23 | ok |
| ARKODE DIRK KVAERNO_7_4_5 | 1e-08 | 7.31 | 399 | 556 | 7052 | 489 | 687 | 0 | 488 | hires | 1.48 | ok |
| ARKODE DIRK KVAERNO_7_4_5 | 1e-10 | 10 | 15432 | 20926 | 298694 | 15591 | 24325 | 0 | 15588 | vdpol | 30.8 | ok |
| ARKODE DIRK KVAERNO_7_4_5 | 1e-10 | 9.88 | 4906 | 6393 | 100750 | 4775 | 7622 | 55 | 4774 | rober | 11.2 | ok |
| ARKODE DIRK KVAERNO_7_4_5 | 1e-10 | 9.08 | 14822 | 20322 | 281737 | 15369 | 23078 | 249 | 15368 | orego | 33.1 | ok |
| ARKODE DIRK KVAERNO_7_4_5 | 1e-10 | 9.17 | 912 | 1239 | 16430 | 1060 | 1504 | 0 | 1059 | hires | 3.42 | ok |
| ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 1e-04 | 3.26 | 59 | 78 | 914 | 21 | 70 | 14 | 20 | hires | <1 | ok |
| ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 1e-04 | 2.44 | 2699 | 3714 | 43414 | 437 | 2249 | 884 | 394 | vdpol | 5.03 | ok |
| ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 1e-04 | 4.95 | 16694 | 29279 | 287005 | 360 | 20611 | 12580 | 40 | rober | 38.7 | ok |
| ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 1e-04 | 1.27 | 1420 | 1982 | 25802 | 612 | 1650 | 383 | 597 | orego | 3.07 | ok |
| ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 1e-06 | 4.38 | 3948 | 6764 | 112105 | 1147 | 5439 | 2516 | 1114 | vdpol | 11.5 | ok |
| ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 1e-06 | 5.18 | 97 | 126 | 1868 | 98 | 150 | 5 | 97 | hires | <1 | ok |
| ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 1e-06 | 6.15 | 15148 | 25871 | 320089 | 542 | 17869 | 10711 | 280 | rober | 39.8 | ok |
| ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 1e-06 | 2.48 | 3096 | 4754 | 73291 | 1475 | 4229 | 1283 | 1461 | orego | 8.59 | ok |
| ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 1e-08 | 7.98 | 29510 | 52635 | 770596 | 2164 | 38573 | 23072 | 1901 | rober | 89.4 | ok |
| ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 1e-08 | 6.76 | 264 | 359 | 5680 | 263 | 409 | 37 | 262 | hires | 1.17 | ok |
| ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 1e-08 | 7.26 | 6198 | 8994 | 152450 | 3415 | 8450 | 2144 | 3414 | vdpol | 15.4 | ok |
| ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 1e-08 | 7.24 | 8506 | 13963 | 218611 | 3510 | 12156 | 4716 | 3448 | orego | 25.3 | ok |
| ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 1e-10 | 10.2 | 72710 | 136878 | 2131479 | 6374 | 106481 | 63823 | 5572 | rober | 242 | ok |
| ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 1e-10 | 10.4* | 20390 | 31117 | 497570 | 11289 | 28692 | 8017 | 11262 | orego | 56.9 | ok |
| ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 1e-10 | 9.73 | 12380 | 16432 | 289607 | 10648 | 18209 | 1185 | 10645 | vdpol | 29.1 | ok |
| ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | 1e-10 | 8.57 | 827 | 1122 | 18071 | 672 | 1172 | 149 | 671 | hires | 3.68 | ok |
<!-- end generated -->

**Verdict.** [pending]

### 5.3 A2: stage order and order reduction (candidate 2)

**Purpose.** Measure the order reduction of the DIRK tables (stage order <= 2) against Radau
IIA (stage order 3) in the stiff limit, and whether it decides the time at tight
tolerances.

**Sweep.**
- A2_fixed (untimed): Prothero-Robinson with N = 4 and lambda in {-1, -1e2, -1e4, -1e6,
  -1e8}, SECTOR at theta = 0 with r in [1, 1e6], fixed steps h = 0.5 / 2^k (k = 0..6) on
  [0, 1]; FIRKODE and the seven DIRK tables. 336 cells. FIRKODE's stacked FGMRES gets 5
  restarts here (`--stage-maxrs 5`; Section 7), so that the cells measure the method, not
  its default linear-solver budget.
- A2_adaptive: Prothero-Robinson with N = 1000 (above the timing floor) and lambda in {-1,
  -1e4, -1e8}; SECTOR with m = 100, theta in {0, 80} degrees, r in [1, 1e3] and [1, 1e6];
  rtol 1e-4 ... 1e-10 in half decades; FIRKODE, DIRK3/4/5, SDIRK_5_3_4, KVAERNO_7_4_5.
  546 cells.

**Predictions [predicted].**
- P2.1: in the stiff limit (|lambda| h >> 1) FIRKODE's fixed-step order is >= 3.5; the
  ESDIRK tables, KVAERNO_7_4_5 and ARK548L2SAb fall into [1.8, 3.2]; SDIRK_5_3_4 into
  [1.5, 2.3]; at lambda = -1 every method is within 0.3 of its classical order (or above it
  on these linear problems).
- P2.2: log10 t(best DIRK)/t(FIRKODE) grows by >= 0.3 from mescd 4 to 8 on Prothero-Robinson
  at lambda = -1e8 and on SECTOR at theta = 0, r <= 1e6; it changes by <= 0.2 at
  lambda = -1.
- P2.3: the ratio at theta = 80 degrees is within [0.67, 1.5]x of that at theta = 0 (both
  method families are A-stable; the A(alpha) mechanism of `BENCHMARK_BDF.md` B8 does not
  apply).

**Pass and fail.** Supports candidate 2 if P2.1 and P2.2 hold; refutes it if some DIRK table
keeps its classical order minus 0.5 in the stiff limit or the ratio grows by < 0.1.

**Results [pending].**

<!-- begin generated: A2_fixed -->
| problem | params | configuration | method | expected order | errors at h0, h0/2, ... | observed orders | reported | at least expected - 0.3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pr | lambda=-1.0;n=4 | FIRKODE s=3 | RadauIIA-s3 | 5 | 5.90e-07, 1.66e-08, 7.40e-10, 2.61e-11, 8.58e-13, 2.74e-14, 5.55e-16 | 5.15, 4.49, 4.83, -, -, - | 4.83 | yes |
| pr | lambda=-1.0;n=4 | ARKODE DIRK3 | ARKODE_ESDIRK325L2SA_5_2_3 | 3 | 7.49e-05, 3.98e-06, 2.06e-07, 7.97e-09, 1.11e-10, 8.31e-11, 1.47e-11 | 4.23, 4.27, 4.69, 6.17, 0.42, 2.50 | 2.50 | NO |
| pr | lambda=-1.0;n=4 | ARKODE DIRK4 | ARKODE_ESDIRK436L2SA_6_3_4 | 4 | 3.35e-05, 2.15e-06, 1.37e-07, 8.65e-09, 5.44e-10, 3.41e-11, 2.14e-12 | 3.96, 3.97, 3.98, 3.99, 3.99, - | 3.99 | yes |
| pr | lambda=-1.0;n=4 | ARKODE DIRK5 | ARKODE_ESDIRK547L2SA2_7_4_5 | 5 | 2.18e-05, 6.24e-07, 1.94e-08, 6.08e-10, 1.90e-11, 5.96e-13, 1.87e-14 | 5.12, 5.01, 5.00, 5.00, -, - | 5.00 | yes |
| pr | lambda=-1.0;n=4 | ARKODE DIRK SDIRK_5_3_4 | ARKODE_SDIRK_5_3_4 | 4 | 7.05e-05, 3.67e-06, 2.17e-07, 1.32e-08, 8.16e-10, 5.07e-11, 3.16e-12 | 4.26, 4.08, 4.04, 4.02, 4.01, 4.00 | 4.00 | yes |
| pr | lambda=-1.0;n=4 | ARKODE DIRK ESDIRK437L2SA_7_3_4 | ARKODE_ESDIRK437L2SA_7_3_4 | 4 | 1.44e-05, 3.55e-07, 7.42e-09, 2.93e-12, 1.43e-11, 1.35e-12, 9.85e-14 | 5.34, 5.58, 11.31, -2.29, -, - | -2.29 | NO |
| pr | lambda=-1.0;n=4 | ARKODE DIRK KVAERNO_7_4_5 | ARKODE_KVAERNO_7_4_5 | 5 | 4.72e-05, 1.40e-06, 4.42e-08, 1.39e-09, 4.38e-11, 1.37e-12, 4.32e-14 | 5.08, 4.98, 4.99, 4.99, -, - | 4.99 | yes |
| pr | lambda=-1.0;n=4 | ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | ARKODE_ARK548L2SAb_DIRK_8_4_5 | 5 | 3.68e-06, 9.02e-08, 2.75e-09, 8.57e-11, 2.68e-12, 8.38e-14, 2.55e-15 | 5.35, 5.04, 5.00, 5.00, -, - | 5.00 | yes |
| pr | lambda=-100.0;n=4 | FIRKODE s=3 | RadauIIA-s3 | 5 | 1.73e-05, 7.99e-07, 2.85e-08, 7.28e-10, 1.50e-11, 3.58e-13, 1.07e-14 | 4.44, 4.81, 5.29, 5.60, -, - | 5.60 | yes (above by 0.60) |
| pr | lambda=-100.0;n=4 | ARKODE DIRK3 | ARKODE_ESDIRK325L2SA_5_2_3 | 3 | 2.64e-05, 3.96e-06, 4.58e-07, 3.65e-08, 3.21e-09, 4.95e-10, 8.63e-11 | 2.74, 3.11, 3.65, 3.50, 2.70, 2.52 | 2.52 | NO |
| pr | lambda=-100.0;n=4 | ARKODE DIRK4 | ARKODE_ESDIRK436L2SA_6_3_4 | 4 | 6.19e-06, 8.38e-06, 1.97e-06, 2.40e-07, 2.04e-08, 1.44e-09, 9.36e-11 | -0.44, 2.09, 3.04, 3.55, 3.82, 3.95 | 3.95 | yes |
| pr | lambda=-100.0;n=4 | ARKODE DIRK5 | ARKODE_ESDIRK547L2SA2_7_4_5 | 5 | 2.33e-04, 3.85e-05, 5.31e-06, 4.88e-07, 2.80e-08, 1.19e-09, 4.39e-11 | 2.60, 2.86, 3.44, 4.12, 4.55, 4.76 | 4.76 | yes |
| pr | lambda=-100.0;n=4 | ARKODE DIRK SDIRK_5_3_4 | ARKODE_SDIRK_5_3_4 | 4 | 1.01e-03, 4.01e-04, 1.06e-04, 1.84e-05, 2.20e-06, 2.05e-07, 1.63e-08 | 1.33, 1.91, 2.53, 3.06, 3.42, 3.66 | 3.66 | NO |
| pr | lambda=-100.0;n=4 | ARKODE DIRK ESDIRK437L2SA_7_3_4 | ARKODE_ESDIRK437L2SA_7_3_4 | 4 | 2.46e-04, 5.93e-05, 7.77e-06, 5.33e-07, 2.55e-08, 1.10e-09, 5.13e-11 | 2.05, 2.93, 3.86, 4.38, 4.53, 4.42 | 4.42 | yes (above by 0.42) |
| pr | lambda=-100.0;n=4 | ARKODE DIRK KVAERNO_7_4_5 | ARKODE_KVAERNO_7_4_5 | 5 | 3.12e-05, 9.76e-06, 1.61e-06, 1.60e-07, 1.04e-08, 4.98e-10, 1.98e-11 | 1.68, 2.60, 3.33, 3.95, 4.38, 4.65 | 4.65 | NO |
| pr | lambda=-100.0;n=4 | ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | ARKODE_ARK548L2SAb_DIRK_8_4_5 | 5 | 2.17e-05, 7.48e-06, 1.29e-06, 1.43e-07, 1.05e-08, 5.67e-10, 2.47e-11 | 1.53, 2.54, 3.17, 3.77, 4.21, 4.52 | 4.52 | NO |
| pr | lambda=-10000.0;n=4 | FIRKODE s=3 | RadauIIA-s3 | 5 | 1.83e-07, 8.77e-09, 4.72e-10, 2.73e-11, 1.64e-12, 9.98e-14, 6.11e-15 | 4.38, 4.22, 4.11, -, -, - | 4.11 | NO |
| pr | lambda=-10000.0;n=4 | ARKODE DIRK3 | ARKODE_ESDIRK325L2SA_5_2_3 | 3 | 2.84e-07, 6.06e-08, 1.46e-08, 3.59e-09, 8.80e-10, 2.11e-10, 4.84e-11 | 2.23, 2.05, 2.02, 2.03, 2.06, 2.13 | 2.13 | NO |
| pr | lambda=-10000.0;n=4 | ARKODE DIRK4 | ARKODE_ESDIRK436L2SA_6_3_4 | 4 | 2.32e-07, 1.47e-07, 4.10e-08, 1.05e-08, 2.61e-09, 6.46e-10, 1.56e-10 | 0.66, 1.84, 1.97, 2.00, 2.02, 2.05 | 2.05 | NO |
| pr | lambda=-10000.0;n=4 | ARKODE DIRK5 | ARKODE_ESDIRK547L2SA2_7_4_5 | 5 | 1.41e-06, 8.73e-07, 2.39e-07, 5.91e-08, 1.40e-08, 3.17e-09, 6.65e-10 | 0.70, 1.87, 2.01, 2.07, 2.15, 2.25 | 2.25 | NO |
| pr | lambda=-10000.0;n=4 | ARKODE DIRK SDIRK_5_3_4 | ARKODE_SDIRK_5_3_4 | 4 | 1.37e-05, 6.28e-06, 2.87e-06, 1.36e-06, 6.52e-07, 3.13e-07, 1.46e-07 | 1.12, 1.13, 1.08, 1.06, 1.06, 1.10 | 1.10 | NO |
| pr | lambda=-10000.0;n=4 | ARKODE DIRK ESDIRK437L2SA_7_3_4 | ARKODE_ESDIRK437L2SA_7_3_4 | 4 | 1.12e-07, 1.90e-06, 5.64e-07, 1.44e-07, 3.55e-08, 8.37e-09, 1.79e-09 | -4.09, 1.75, 1.97, 2.02, 2.08, 2.23 | 2.23 | NO |
| pr | lambda=-10000.0;n=4 | ARKODE DIRK KVAERNO_7_4_5 | ARKODE_KVAERNO_7_4_5 | 5 | 5.13e-07, 1.66e-07, 4.33e-08, 1.09e-08, 2.68e-09, 6.49e-10, 1.51e-10 | 1.62, 1.94, 2.00, 2.02, 2.05, 2.10 | 2.10 | NO |
| pr | lambda=-10000.0;n=4 | ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | ARKODE_ARK548L2SAb_DIRK_8_4_5 | 5 | 3.30e-07, 1.33e-07, 3.54e-08, 8.87e-09, 2.17e-09, 5.18e-10, 1.18e-10 | 1.32, 1.91, 2.00, 2.03, 2.07, 2.13 | 2.13 | NO |
| pr | lambda=-1000000.0;n=4 | FIRKODE s=3 | RadauIIA-s3 | 5 | 1.83e-09, 8.76e-11, 4.71e-12, 2.72e-13, 1.64e-14, 9.99e-16, 0.00e+00 | 4.38, 4.22, -, -, -, - | 4.22 | NO |
| pr | lambda=-1000000.0;n=4 | ARKODE DIRK3 | ARKODE_ESDIRK325L2SA_5_2_3 | 3 | 2.84e-09, 6.08e-10, 1.47e-10, 3.66e-11, 9.13e-12, 2.28e-12, 5.70e-13 | 2.22, 2.05, 2.01, 2.00, 2.00, - | 2.00 | NO |
| pr | lambda=-1000000.0;n=4 | ARKODE DIRK4 | ARKODE_ESDIRK436L2SA_6_3_4 | 4 | 2.37e-09, 1.48e-09, 4.11e-10, 1.05e-10, 2.64e-11, 6.61e-12, 1.65e-12 | 0.68, 1.84, 1.97, 1.99, 2.00, - | 2.00 | NO |
| pr | lambda=-1000000.0;n=4 | ARKODE DIRK5 | ARKODE_ESDIRK547L2SA2_7_4_5 | 5 | 1.38e-08, 8.80e-09, 2.45e-09, 6.26e-10, 1.57e-10, 3.93e-11, 9.80e-12 | 0.64, 1.84, 1.97, 1.99, 2.00, 2.00 | 2.00 | NO |
| pr | lambda=-1000000.0;n=4 | ARKODE DIRK SDIRK_5_3_4 | ARKODE_SDIRK_5_3_4 | 4 | 1.37e-07, 6.30e-08, 2.89e-08, 1.37e-08, 6.69e-09, 3.30e-09, 1.64e-09 | 1.12, 1.12, 1.07, 1.04, 1.02, 1.01 | 1.01 | NO |
| pr | lambda=-1000000.0;n=4 | ARKODE DIRK ESDIRK437L2SA_7_3_4 | ARKODE_ESDIRK437L2SA_7_3_4 | 4 | 5.34e-10, 1.93e-08, 5.70e-09, 1.47e-09, 3.70e-10, 9.27e-11, 2.32e-11 | -5.18, 1.76, 1.95, 1.99, 2.00, 2.00 | 2.00 | NO |
| pr | lambda=-1000000.0;n=4 | ARKODE DIRK KVAERNO_7_4_5 | ARKODE_KVAERNO_7_4_5 | 5 | 5.17e-09, 1.67e-09, 4.37e-10, 1.10e-10, 2.76e-11, 6.90e-12, 1.73e-12 | 1.63, 1.94, 1.99, 2.00, 2.00, - | 2.00 | NO |
| pr | lambda=-1000000.0;n=4 | ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | ARKODE_ARK548L2SAb_DIRK_8_4_5 | 5 | 3.31e-09, 1.33e-09, 3.58e-10, 9.09e-11, 2.28e-11, 5.70e-12, 1.42e-12 | 1.31, 1.90, 1.98, 2.00, 2.00, - | 2.00 | NO |
| pr | lambda=-100000000.0;n=4 | FIRKODE s=3 | RadauIIA-s3 | 5 | 1.83e-11, 8.76e-13, 4.71e-14, 2.66e-15, 1.11e-16, 0.00e+00, 0.00e+00 | -, -, -, -, -, - | NA | NA |
| pr | lambda=-100000000.0;n=4 | ARKODE DIRK3 | ARKODE_ESDIRK325L2SA_5_2_3 | 3 | 2.84e-11, 6.08e-12, 1.47e-12, 3.66e-13, 9.16e-14, 2.31e-14, 5.77e-15 | 2.22, -, -, -, -, - | 2.22 | NO |
| pr | lambda=-100000000.0;n=4 | ARKODE DIRK4 | ARKODE_ESDIRK436L2SA_6_3_4 | 4 | 2.37e-11, 1.48e-11, 4.11e-12, 1.05e-12, 2.64e-13, 6.61e-14, 1.64e-14 | 0.68, 1.84, -, -, -, - | 1.84 | NO |
| pr | lambda=-100000000.0;n=4 | ARKODE DIRK5 | ARKODE_ESDIRK547L2SA2_7_4_5 | 5 | 1.38e-10, 8.81e-11, 2.45e-11, 6.27e-12, 1.57e-12, 3.96e-13, 9.84e-14 | 0.64, 1.84, 1.97, -, -, - | 1.97 | NO |
| pr | lambda=-100000000.0;n=4 | ARKODE DIRK SDIRK_5_3_4 | ARKODE_SDIRK_5_3_4 | 4 | 1.37e-09, 6.30e-10, 2.89e-10, 1.37e-10, 6.69e-11, 3.30e-11, 1.64e-11 | 1.12, 1.12, 1.07, 1.04, 1.02, 1.01 | 1.01 | NO |
| pr | lambda=-100000000.0;n=4 | ARKODE DIRK ESDIRK437L2SA_7_3_4 | ARKODE_ESDIRK437L2SA_7_3_4 | 4 | 5.51e-12, 1.93e-10, 5.70e-11, 1.47e-11, 3.71e-12, 9.23e-13, 2.38e-13 | -5.13, 1.76, 1.95, 1.99, -, - | 1.99 | NO |
| pr | lambda=-100000000.0;n=4 | ARKODE DIRK KVAERNO_7_4_5 | ARKODE_KVAERNO_7_4_5 | 5 | 5.17e-11, 1.67e-11, 4.37e-12, 1.10e-12, 2.76e-13, 6.92e-14, 1.72e-14 | 1.63, 1.94, -, -, -, - | 1.94 | NO |
| pr | lambda=-100000000.0;n=4 | ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | ARKODE_ARK548L2SAb_DIRK_8_4_5 | 5 | 3.31e-11, 1.33e-11, 3.58e-12, 9.09e-13, 2.28e-13, 5.71e-14, 1.42e-14 | 1.31, 1.90, -, -, -, - | 1.90 | NO |
| sector | m=10;rmax=1000000.0;rmin=1.0;theta=0.0 | FIRKODE s=3 | RadauIIA-s3 | 5 | 7.16e-04, 9.19e-05, 8.67e-06, 5.84e-07, 2.82e-08, 1.08e-09, 3.73e-11 | 2.96, 3.41, 3.89, 4.37, 4.70, 4.86 | 4.86 | yes |
| sector | m=10;rmax=1000000.0;rmin=1.0;theta=0.0 | ARKODE DIRK3 | ARKODE_ESDIRK325L2SA_5_2_3 | 3 | 1.42e-03, 1.13e-04, 1.01e-05, 1.34e-06, 2.15e-07, 3.00e-08, 3.95e-09 | 3.66, 3.48, 2.91, 2.64, 2.84, 2.92 | 2.92 | yes |
| sector | m=10;rmax=1000000.0;rmin=1.0;theta=0.0 | ARKODE DIRK4 | ARKODE_ESDIRK436L2SA_6_3_4 | 4 | 2.46e-03, 3.29e-04, 2.69e-05, 1.51e-06, 3.42e-07, 7.57e-08, 1.40e-08 | 2.90, 3.61, 4.15, 2.15, 2.18, 2.43 | 2.43 | NO |
| sector | m=10;rmax=1000000.0;rmin=1.0;theta=0.0 | ARKODE DIRK5 | ARKODE_ESDIRK547L2SA2_7_4_5 | 5 | 8.49e-03, 6.16e-04, 3.50e-05, 6.88e-06, 1.39e-06, 2.52e-07, 3.84e-08 | 3.78, 4.14, 2.35, 2.31, 2.46, 2.71 | 2.71 | NO |
| sector | m=10;rmax=1000000.0;rmin=1.0;theta=0.0 | ARKODE DIRK SDIRK_5_3_4 | ARKODE_SDIRK_5_3_4 | 4 | 1.40e-02, 3.93e-03, 1.45e-03, 2.81e-04, 4.21e-05, 8.58e-06, 1.49e-06 | -, 1.44, 2.36, 2.74, 2.29, 2.53 | 2.53 | NO |
| sector | m=10;rmax=1000000.0;rmin=1.0;theta=0.0 | ARKODE DIRK ESDIRK437L2SA_7_3_4 | ARKODE_ESDIRK437L2SA_7_3_4 | 4 | 2.38e-02, 2.11e-03, 8.75e-05, 1.85e-05, 3.75e-06, 5.97e-07, 6.75e-08 | -, 4.59, 2.24, 2.30, 2.65, 3.14 | 3.14 | NO |
| sector | m=10;rmax=1000000.0;rmin=1.0;theta=0.0 | ARKODE DIRK KVAERNO_7_4_5 | ARKODE_KVAERNO_7_4_5 | 5 | 1.60e-03, 1.95e-04, 1.26e-05, 1.42e-06, 3.25e-07, 6.75e-08, 1.14e-08 | 3.04, 3.96, 3.14, 2.13, 2.27, 2.57 | 2.57 | NO |
| sector | m=10;rmax=1000000.0;rmin=1.0;theta=0.0 | ARKODE DIRK ARK548L2SAb_DIRK_8_4_5 | ARKODE_ARK548L2SAb_DIRK_8_4_5 | 5 | 1.61e-03, 1.66e-04, 1.20e-05, 1.14e-06, 2.54e-07, 5.18e-08, 8.96e-09 | 3.28, 3.79, 3.39, 2.17, 2.29, 2.53 | 2.53 | NO |
<!-- end generated -->

![fixed-step errors in the stiff limit](../../doc/shared/figs/firkode/benchmarks/arkode/arkode_order_fixed.png)

<!-- begin generated: A2_slope -->
| problem | params | class | t(best) / t(FIRKODE) at mescd=4 | t(best) / t(FIRKODE) at mescd=6 | t(best) / t(FIRKODE) at mescd=8 | log10 ratio(8) - log10 ratio(4) |
| --- | --- | --- | --- | --- | --- | --- |
| pr | lambda=-1.0;n=1000 | best DIRK default | NA | NA | NA | NA |
| pr | lambda=-1.0;n=1000 | best DIRK of all | NA | NA | NA | NA |
| pr | lambda=-10000.0;n=1000 | best DIRK default | NA | NA | 0.513 (ARKODE DIRK3) | NA |
| pr | lambda=-10000.0;n=1000 | best DIRK of all | NA | NA | 0.513 (ARKODE DIRK3) | NA |
| pr | lambda=-100000000.0;n=1000 | best DIRK default | NA | NA | NA | NA |
| pr | lambda=-100000000.0;n=1000 | best DIRK of all | NA | NA | NA | NA |
| sector | m=100;rmax=1000.0;rmin=1.0;theta=0.0 | best DIRK default | NA | 0.914 (ARKODE DIRK3) | 1.69 (ARKODE DIRK4) | NA |
| sector | m=100;rmax=1000.0;rmin=1.0;theta=0.0 | best DIRK of all | NA | 0.914 (ARKODE DIRK3) | 1.69 (ARKODE DIRK4) | NA |
| sector | m=100;rmax=1000.0;rmin=1.0;theta=80.0 | best DIRK default | NA | 0.47 (ARKODE DIRK3) | 0.932 (ARKODE DIRK4) | NA |
| sector | m=100;rmax=1000.0;rmin=1.0;theta=80.0 | best DIRK of all | NA | 0.47 (ARKODE DIRK3) | 0.932 (ARKODE DIRK4) | NA |
| sector | m=100;rmax=1000000.0;rmin=1.0;theta=0.0 | best DIRK default | NA | 0.919 (ARKODE DIRK3) | 1.17 (ARKODE DIRK3) | NA |
| sector | m=100;rmax=1000000.0;rmin=1.0;theta=0.0 | best DIRK of all | NA | 0.919 (ARKODE DIRK3) | 1.17 (ARKODE DIRK3) | NA |
| sector | m=100;rmax=1000000.0;rmin=1.0;theta=80.0 | best DIRK default | NA | 0.488 (ARKODE DIRK3) | 1.23 (ARKODE DIRK3) | NA |
| sector | m=100;rmax=1000000.0;rmin=1.0;theta=80.0 | best DIRK of all | NA | 0.488 (ARKODE DIRK3) | 1.23 (ARKODE DIRK3) | NA |
<!-- end generated -->

![work-precision on stiff Prothero-Robinson and SECTOR](../../doc/shared/figs/firkode/benchmarks/arkode/arkode_stage_order_wp.png)

**Verdict.** [pending]

### 5.4 A3: mass matrices (candidate 3)

**Purpose.** Both packages take a constant mass matrix; compare their cost on the
finite-element heat and wave problems.

**Sweep.** `heat1dfem` with n in {101, 401}, rtol 1e-4 ... 1e-10 in decades, FIRKODE and
DIRK3/4/5 (B5's cells); `wave1dfem` with n in {100, 400}, Kelvin-Voigt and mass damping,
lumped and consistent mass, rtol 1e-3 ... 1e-8, FIRKODE and DIRK4/5. 200 cells.

**Predictions [predicted].**
- P3.1: ARKODE's mass solves per attempted step lie in [1, 3], and its consistent-mass time
  is <= 1.5x its lumped-mass time.
- P3.2: the ratio t(best DIRK)/t(FIRKODE) with a consistent mass is within [0.67, 1.5]x of
  the ratio with the lumped mass. If FIRKODE's own consistent/lumped time ratio exceeds 2
  (its slow mass-matrix path, `BENCHMARK_BDF.md` B9), the difference is FIRKODE's and goes
  to Section 7.

**Results [pending].**

<!-- begin generated: A3_mass -->
| problem | params | class | t(best) / t(FIRKODE) at mescd=4 | t(best) / t(FIRKODE) at mescd=6 | t(best) / t(FIRKODE) at mescd=8 |
| --- | --- | --- | --- | --- | --- |
| heat1dfem | n=101 | best DIRK default | NA | 1.47 (ARKODE DIRK3) | 2.51 (ARKODE DIRK3) |
| heat1dfem | n=101 | best DIRK of all | NA | 1.47 (ARKODE DIRK3) | 2.51 (ARKODE DIRK3) |
| heat1dfem | n=401 | best DIRK default | NA | 0.876 (ARKODE DIRK3) | 2.6 (ARKODE DIRK4) |
| heat1dfem | n=401 | best DIRK of all | NA | 0.876 (ARKODE DIRK3) | 2.6 (ARKODE DIRK4) |
| wave1dfem | damping=kv;mass=consistent;n=100 | best DIRK default | NA | NA | NA |
| wave1dfem | damping=kv;mass=consistent;n=100 | best DIRK of all | NA | NA | NA |
| wave1dfem | damping=kv;mass=consistent;n=400 | best DIRK default | NA | NA | NA |
| wave1dfem | damping=kv;mass=consistent;n=400 | best DIRK of all | NA | NA | NA |
| wave1dfem | damping=kv;mass=lumped;n=100 | best DIRK default | NA | 0.569 (ARKODE DIRK5) | NA |
| wave1dfem | damping=kv;mass=lumped;n=100 | best DIRK of all | NA | 0.569 (ARKODE DIRK5) | NA |
| wave1dfem | damping=kv;mass=lumped;n=400 | best DIRK default | 0.894 (ARKODE DIRK5) | 0.707 (ARKODE DIRK5) | NA |
| wave1dfem | damping=kv;mass=lumped;n=400 | best DIRK of all | 0.894 (ARKODE DIRK5) | 0.707 (ARKODE DIRK5) | NA |
| wave1dfem | damping=mass;mass=consistent;n=100 | best DIRK default | NA | NA | NA |
| wave1dfem | damping=mass;mass=consistent;n=100 | best DIRK of all | NA | NA | NA |
| wave1dfem | damping=mass;mass=consistent;n=400 | best DIRK default | NA | NA | NA |
| wave1dfem | damping=mass;mass=consistent;n=400 | best DIRK of all | NA | NA | NA |
| wave1dfem | damping=mass;mass=lumped;n=100 | best DIRK default | 0.72 (ARKODE DIRK5) | 0.502 (ARKODE DIRK5) | NA |
| wave1dfem | damping=mass;mass=lumped;n=100 | best DIRK of all | 0.72 (ARKODE DIRK5) | 0.502 (ARKODE DIRK5) | NA |
| wave1dfem | damping=mass;mass=lumped;n=400 | best DIRK default | 0.728 (ARKODE DIRK5) | 0.524 (ARKODE DIRK5) | NA |
| wave1dfem | damping=mass;mass=lumped;n=400 | best DIRK of all | 0.728 (ARKODE DIRK5) | 0.524 (ARKODE DIRK5) | NA |
<!-- end generated -->

<!-- begin generated: A3_counters -->
| problem | params | method | mescd | steps | attempts | mass solves | mass solves / attempt | Newton | stacked FGMRES / Newton |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| heat1dfem | n=101 | FIRKODE s=3 | 7.77 | 11 | 11 | 42 | 3.82 | 13 | 3.31 |
| heat1dfem | n=101 | ARKODE DIRK3 | 6.87 | 13 | 13 | 28 | 2.15 | 142 | NA |
| heat1dfem | n=101 | ARKODE DIRK4 | 6.3 | 12 | 12 | 26 | 2.17 | 163 | NA |
| heat1dfem | n=101 | ARKODE DIRK5 | 6.49 | 11 | 11 | 24 | 2.18 | 174 | NA |
| heat1dfem | n=401 | FIRKODE s=3 | 7.6 | 10 | 10 | 39 | 3.9 | 12 | 3.5 |
| heat1dfem | n=401 | ARKODE DIRK3 | 7.05 | 13 | 13 | 28 | 2.15 | 143 | NA |
| heat1dfem | n=401 | ARKODE DIRK4 | 6.32 | 12 | 12 | 26 | 2.17 | 163 | NA |
| heat1dfem | n=401 | ARKODE DIRK5 | 6.33 | 12 | 12 | 26 | 2.17 | 189 | NA |
| wave1dfem | damping=kv;mass=consistent;n=100 | FIRKODE s=3 | 7.27 | 47361 | 75385 | 286756 | 3.8 | 95584 | 12.5 |
| wave1dfem | damping=kv;mass=lumped;n=100 | FIRKODE s=3 | 7.54 | 166 | 174 | NA | NA | 199 | 3.72 |
| wave1dfem | damping=kv;mass=lumped;n=100 | ARKODE DIRK4 | 4.81 | 80 | 86 | NA | NA | 1070 | NA |
| wave1dfem | damping=kv;mass=lumped;n=100 | ARKODE DIRK5 | 5.56 | 56 | 63 | NA | NA | 896 | NA |
| wave1dfem | damping=kv;mass=lumped;n=400 | FIRKODE s=3 | 6.74 | 174 | 180 | NA | NA | 206 | 3.85 |
| wave1dfem | damping=kv;mass=lumped;n=400 | ARKODE DIRK4 | 4.36 | 87 | 94 | NA | NA | 1171 | NA |
| wave1dfem | damping=kv;mass=lumped;n=400 | ARKODE DIRK5 | 5.57 | 61 | 73 | NA | NA | 1035 | NA |
| wave1dfem | damping=mass;mass=consistent;n=100 | FIRKODE s=3 | 6.47 | 54688 | 87332 | 331705 | 3.8 | 110567 | 12.2 |
| wave1dfem | damping=mass;mass=lumped;n=100 | FIRKODE s=3 | 7.38 | 184 | 189 | NA | NA | 216 | 4.39 |
| wave1dfem | damping=mass;mass=lumped;n=100 | ARKODE DIRK4 | 4.19 | 89 | 95 | NA | NA | 1191 | NA |
| wave1dfem | damping=mass;mass=lumped;n=100 | ARKODE DIRK5 | 5.52 | 66 | 80 | NA | NA | 1114 | NA |
| wave1dfem | damping=mass;mass=lumped;n=400 | FIRKODE s=3 | 6.75 | 183 | 187 | NA | NA | 214 | 5.51 |
| wave1dfem | damping=mass;mass=lumped;n=400 | ARKODE DIRK4 | 4.19 | 89 | 95 | NA | NA | 1189 | NA |
| wave1dfem | damping=mass;mass=lumped;n=400 | ARKODE DIRK5 | 5.51 | 66 | 80 | NA | NA | 1118 | NA |
<!-- end generated -->

![work-precision with mass matrices](../../doc/shared/figs/firkode/benchmarks/arkode/arkode_mass_wp.png)

**Verdict.** [pending]

### 5.5 A4: tight tolerances, order 5 against order 5 (candidate 4)

**Purpose.** At tight tolerances the order-5 DIRK tables meet Radau IIA of order 5; where
does each pay?

**Sweep.** ROBER, VDPOL, OREGO, HIRES (references at rtol 1e-16 in extended precision) and
SECTOR at theta = 0 and 80 degrees (exact): rtol 1e-8 ... 1e-13 in half decades, FIRKODE,
DIRK5, DIRK5 optimal, KVAERNO_7_4_5, ARK548L2SAb; `bruss3d` n = 16 at rtol 1e-8 ... 1e-11
with FIRKODE and DIRK5. 344 cells; needs `--ref-bin`.

**Predictions [predicted].**
- P4.1: on the four small problems the ratio changes by < 1.5x between mescd 8 and 12 (equal
  order).
- P4.2: on SECTOR the ratio grows by >= 2x from mescd 8 to 12 (stage order, as in A2); no
  crossover in FIRKODE's disfavour up to the roundoff floor.

**Results [pending].**

<!-- begin generated: A4_crossover -->
| problem | params | ratio at mescd=8 | ratio at mescd=10 | ratio at mescd=12 | crossover mescd | floor (ref disagreement) |
| --- | --- | --- | --- | --- | --- | --- |
| bruss3d | n=16 | NA | NA | NA | none | NA |
| hires |  | NA | 1.08 (ARKODE DIRK5) | 1.05 (ARKODE DIRK5) | 10.7 | 13.5 |
| orego |  | NA | 5.34 (ARKODE DIRK KVAERNO_7_4_5) | 4.61 (ARKODE DIRK KVAERNO_7_4_5) | none | 12.5 |
| rober |  | NA | NA | 8.62 (ARKODE DIRK KVAERNO_7_4_5) | none | 14.0 |
| sector | m=100;rmax=1000.0;rmin=1.0;theta=0.0 | NA | 2.08 (ARKODE DIRK KVAERNO_7_4_5) | 3.01 (ARKODE DIRK KVAERNO_7_4_5) | none | exact |
| sector | m=100;rmax=1000.0;rmin=1.0;theta=80.0 | NA | 1.22 (ARKODE DIRK KVAERNO_7_4_5) | NA | none | exact |
| vdpol |  | NA | NA | 11.3 (ARKODE DIRK KVAERNO_7_4_5) | none | 13.8 |
<!-- end generated -->

![work-precision at tight tolerances](../../doc/shared/figs/firkode/benchmarks/arkode/arkode_tight_wp.png)

**Verdict.** [pending]

### 5.6 A5: 3-D Brusselator, full implicitness against IMEX and LSRK (candidate 5)

**Purpose.** On a reaction-diffusion problem whose stiffness is the diffusion (linear,
constant Jacobian) and whose reaction is mild, compare FIRKODE with IMEX (implicit linear
diffusion, explicit reaction) and with explicit stabilized RKC/RKL.

**Sweep.** n in {8, 12, 16} (24 with `--include-large`), KLU, rtol 1e-3 ... 1e-8 in decades;
FIRKODE, DIRK4, IMEX3/4/5, RKC2, RKL2, RKC2 power; profiled twins of FIRKODE, DIRK4 and IMEX4
at rtol 1e-6. 153 cells.

**Predictions [predicted].**
- P5.1: IMEX4 / FIRKODE <= 1/2 at mescd 4 and 6 for n >= 12 and <= 1/1.2 at mescd 8, with
  nje <= 2 and about one Newton iteration per implicit stage.
- P5.2: best LSRK / FIRKODE <= 1/3 at mescd 4 for n >= 12 but > 1 at mescd 8; the
  crossover lies in mescd [5, 8]; LSRK's most stages per step is about
  sqrt(h rho / 0.65) within 30 percent.
- P5.3: the power iteration costs <= 5 percent of LSRK's right-hand sides, but RKC2 power is
  >= 2x slower than RKC2 with the bound: with 100 iterations and a relative tolerance of
  5e-3 the iteration underestimates the spectral radius of a diffusion operator whose top
  eigenvalues cluster, and the steps it allows fail the error test [development: 103
  error-test failures against 1 on `heat2d` 40 x 40 at rtol 1e-4].
- P5.4: DIRK4 / FIRKODE >= 1.
- Profiled twins: the factorization share is >= 50 percent for FIRKODE and <= 20 percent for
  IMEX4.

**Results [pending].**

<!-- begin generated: A5_ratio -->
| n | class | t(best) / t(FIRKODE) at mescd=4 | t(best) / t(FIRKODE) at mescd=6 | t(best) / t(FIRKODE) at mescd=8 |
| --- | --- | --- | --- | --- |
| 12 | best DIRK default | NA | NA | NA |
| 12 | best DIRK of all | NA | NA | NA |
| 12 | best IMEX | NA | NA | NA |
| 12 | best LSRK | NA | NA | NA |
| 16 | best DIRK default | NA | NA | NA |
| 16 | best DIRK of all | NA | NA | NA |
| 16 | best IMEX | NA | NA | NA |
| 16 | best LSRK | NA | NA | NA |
| 8 | best DIRK default | NA | NA | NA |
| 8 | best DIRK of all | NA | NA | NA |
| 8 | best IMEX | NA | NA | NA |
| 8 | best LSRK | NA | NA | NA |
<!-- end generated -->

<!-- begin generated: A5_lsrk_crossover -->
| n | LSRK crossover mescd |
| --- | --- |
| 12 | none |
| 16 | none |
| 8 | none |
<!-- end generated -->

<!-- begin generated: A5_counters -->
| problem | params | method | mescd | steps | RHS (implicit or all) | RHS explicit | Newton | Jac | setups | LSRK max stages | power RHS / RHS |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| bruss3d | n=12 | ARKODE LSRK RKC2 | NA | 118 | 1515 | NA | NA | NA | NA | 31 | NA |
| bruss3d | n=12 | ARKODE LSRK RKC2 power | NA | 122 | 1648 | NA | NA | NA | NA | 30 | 0.0394 |
| bruss3d | n=12 | ARKODE LSRK RKL2 | NA | 116 | 1690 | NA | NA | NA | NA | 38 | NA |
| bruss3d | n=16 | ARKODE LSRK RKC2 | NA | 116 | 2049 | NA | NA | NA | NA | 42 | NA |
| bruss3d | n=16 | ARKODE LSRK RKC2 power | NA | 126 | 2338 | NA | NA | NA | NA | 43 | 0.0278 |
| bruss3d | n=16 | ARKODE LSRK RKL2 | NA | 114 | 2280 | NA | NA | NA | NA | 52 | NA |
| bruss3d | n=8 | ARKODE LSRK RKC2 | NA | 123 | 1000 | NA | NA | NA | NA | 20 | NA |
| bruss3d | n=8 | ARKODE LSRK RKC2 power | NA | 124 | 1006 | NA | NA | NA | NA | 19 | 0.0457 |
| bruss3d | n=8 | ARKODE LSRK RKL2 | NA | 121 | 1108 | NA | NA | NA | NA | 24 | NA |
<!-- end generated -->

<!-- begin generated: A5_profile -->
(no profiled runs)
<!-- end generated -->

![work-precision of the 3-D Brusselator](../../doc/shared/figs/firkode/benchmarks/arkode/arkode_bruss3d_wp.png)

![time split of the profiled runs](../../doc/shared/figs/firkode/benchmarks/arkode/arkode_bruss3d_profile.png)

**Verdict.** [pending]

### 5.7 A6: advection-diffusion over the velocity (candidate 6)

**Purpose.** Move a problem from diffusion-dominated to advection-dominated and watch each
structure-exploiting class lose its advantage at a predicted point: IMEX at the advective
CFL limit, LSRK when the spectrum leaves the real axis (cell Peclet number above 1), MRI
not at all.

**Sweep.** `heat2d` on 200 x 100 with velocity v in {0, 0.5, 50, 500} (cell Peclet numbers
0, 0.0012, 0.12, 1.24), rtol 1e-3 ... 1e-8 in decades; FIRKODE, DIRK4, IMEX4, IMEX5, RKC2,
RKL2, RKC2 power, and for v > 0 MRI ESDIRK34a and MRI IMEX-SR32. Control: SECTOR (m = 100,
r in [1, 1e3]) at theta = 0 and 30 degrees with FIRKODE, DIRK4, RKC2 and RKL2. 252 cells.

**Predictions [predicted].**
- P6.1 (IMEX against the CFL limit): IMEX4 / FIRKODE <= 1/1.2 at v <= 0.5, >= 2 at v = 50,
  >= 20 at v = 500; IMEX's step count at v = 500 is >= 10x that at v = 0.5, with
  error-test failures >= 10 percent of the attempts.
- P6.2 (LSRK against complex spectra): best LSRK / FIRKODE <= 1/2 at mescd 4 for v <= 50,
  with a crossover in mescd [5, 7]; at v = 500 LSRK takes >= 3x its v = 50 time. On SECTOR
  at theta = 30 degrees LSRK takes >= 5x its theta = 0 time or fails, while FIRKODE changes
  by <= 1.5x.
- P6.3 (MRI): MRI's slow step count at v = 500 is within 2x of that at v = 0.5 while its fast
  step count grows >= 5x from v = 50 to 500; MRI / IMEX4 <= 1/3 at v = 500. MRI against
  FIRKODE is registered as an inconclusive band [1/3, 3].

**Results [pending].**

<!-- begin generated: A6_map -->
| velocity | class | t(best) / t(FIRKODE) at mescd=4 | t(best) / t(FIRKODE) at mescd=6 | t(best) / t(FIRKODE) at mescd=8 |
| --- | --- | --- | --- | --- |
| 0.0 | best DIRK default | NA | 0.994 (ARKODE DIRK4) | 1.2 (ARKODE DIRK4) |
| 0.0 | best DIRK of all | NA | 0.994 (ARKODE DIRK4) | 1.2 (ARKODE DIRK4) |
| 0.0 | best IMEX | NA | 3.44 (ARKODE IMEX5) | 3.6 (ARKODE IMEX5) |
| 0.0 | best LSRK | NA | 0.446 (ARKODE LSRK RKC2) | 1.43 (ARKODE LSRK RKC2 power) |
| 0.5 | best DIRK default | 0.494 (ARKODE DIRK4) | 0.57 (ARKODE DIRK4) | 0.492 (ARKODE DIRK4) |
| 0.5 | best DIRK of all | 0.494 (ARKODE DIRK4) | 0.57 (ARKODE DIRK4) | 0.492 (ARKODE DIRK4) |
| 0.5 | best IMEX | 3.55 (ARKODE IMEX4) | 6.41 (ARKODE IMEX4) | NA |
| 0.5 | best LSRK | NA | 0.189 (ARKODE LSRK RKC2) | NA |
| 0.5 | best MRI | NA | 12.1 (ARKODE MRI IMEX-SR32) | 51.6 (ARKODE MRI IMEX-SR32) |
| 50.0 | best DIRK default | 0.503 (ARKODE DIRK4) | 0.78 (ARKODE DIRK4) | NA |
| 50.0 | best DIRK of all | 0.503 (ARKODE DIRK4) | 0.78 (ARKODE DIRK4) | NA |
| 50.0 | best IMEX | 74.7 (ARKODE IMEX5) | 45.4 (ARKODE IMEX5) | NA |
| 50.0 | best LSRK | 0.484 (ARKODE LSRK RKL2) | 0.335 (ARKODE LSRK RKL2) | NA |
| 50.0 | best MRI | 165 (ARKODE MRI IMEX-SR32) | NA | NA |
| 500.0 | best DIRK default | 0.508 (ARKODE DIRK4) | 1.14 (ARKODE DIRK4) | NA |
| 500.0 | best DIRK of all | 0.508 (ARKODE DIRK4) | 1.14 (ARKODE DIRK4) | NA |
| 500.0 | best IMEX | NA | NA | NA |
| 500.0 | best LSRK | 5.97 (ARKODE LSRK RKL2) | 2.7 (ARKODE LSRK RKL2) | NA |
| 500.0 | best MRI | NA | NA | NA |
<!-- end generated -->

<!-- begin generated: A6_sector -->
| theta | class | t(best) / t(FIRKODE) at mescd=4 | t(best) / t(FIRKODE) at mescd=6 | t(best) / t(FIRKODE) at mescd=8 |
| --- | --- | --- | --- | --- |
| 0.0 | best DIRK default | 0.481 (ARKODE DIRK4) | 1.53 (ARKODE DIRK4) | NA |
| 0.0 | best DIRK of all | 0.481 (ARKODE DIRK4) | 1.53 (ARKODE DIRK4) | NA |
| 0.0 | best LSRK | 1.36 (ARKODE LSRK RKL2) | 2.61 (ARKODE LSRK RKC2) | NA |
| 30.0 | best DIRK default | 1.09 (ARKODE DIRK4) | 1.06 (ARKODE DIRK4) | NA |
| 30.0 | best DIRK of all | 1.09 (ARKODE DIRK4) | 1.06 (ARKODE DIRK4) | NA |
| 30.0 | best LSRK | NA | 2.3 (ARKODE LSRK RKL2) | NA |
<!-- end generated -->

<!-- begin generated: A6_counters -->
| problem | params | method | mescd | steps | attempts | err fails | RHS (implicit or all) | RHS explicit | fast steps | LSRK max stages |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| heat2d | mx=200;my=100;velocity=0.0 | FIRKODE s=3 | 7.48 | 10 | 11 | 1 | 51 | NA | NA | NA |
| heat2d | mx=200;my=100;velocity=0.0 | ARKODE DIRK4 | 6.5 | 13 | 13 | 0 | 246 | NA | NA | NA |
| heat2d | mx=200;my=100;velocity=0.0 | ARKODE IMEX4 | 5.39 | 27 | 27 | 0 | 353 | 191 | NA | NA |
| heat2d | mx=200;my=100;velocity=0.0 | ARKODE IMEX5 | 6.51 | 40 | 40 | 0 | 602 | 322 | NA | NA |
| heat2d | mx=200;my=100;velocity=0.0 | ARKODE LSRK RKC2 | 5.89 | 57 | 57 | 0 | 4189 | NA | NA | 83 |
| heat2d | mx=200;my=100;velocity=0.0 | ARKODE LSRK RKC2 power | 7.88 | 782 | 1135 | 353 | 33911 | NA | NA | 80 |
| heat2d | mx=200;my=100;velocity=0.0 | ARKODE LSRK RKL2 | 5.99 | 55 | 55 | 0 | 4723 | NA | NA | 95 |
| heat2d | mx=200;my=100;velocity=0.5 | FIRKODE s=3 | 7.28 | 14 | 15 | 1 | 70 | NA | NA | NA |
| heat2d | mx=200;my=100;velocity=0.5 | ARKODE DIRK4 | 6.38 | 15 | 15 | 0 | 267 | NA | NA | NA |
| heat2d | mx=200;my=100;velocity=0.5 | ARKODE IMEX4 | 4.66 | 72 | 72 | 0 | 938 | 506 | NA | NA |
| heat2d | mx=200;my=100;velocity=0.5 | ARKODE IMEX5 | 5.85 | 180 | 180 | 0 | 2702 | 1442 | NA | NA |
| heat2d | mx=200;my=100;velocity=0.5 | ARKODE LSRK RKC2 | 5.9 | 62 | 62 | 0 | 4311 | NA | NA | 78 |
| heat2d | mx=200;my=100;velocity=0.5 | ARKODE LSRK RKC2 power | 7.82 | 646 | 954 | 308 | 33739 | NA | NA | 74 |
| heat2d | mx=200;my=100;velocity=0.5 | ARKODE LSRK RKL2 | 5.93 | 59 | 59 | 0 | 4812 | NA | NA | 90 |
| heat2d | mx=200;my=100;velocity=0.5 | ARKODE MRI ESDIRK34a | 7.06 | 2385 | 2385 | 0 | 16697 | NA | 7157 | NA |
| heat2d | mx=200;my=100;velocity=0.5 | ARKODE MRI IMEX-SR32 | 7.43 | 1394 | 1394 | 0 | 11154 | 5578 | 6972 | NA |
| heat2d | mx=200;my=100;velocity=50.0 | FIRKODE s=3 | 7.37 | 20 | 20 | 0 | 91 | NA | NA | NA |
| heat2d | mx=200;my=100;velocity=50.0 | ARKODE DIRK4 | 5.08 | 20 | 20 | 0 | 354 | NA | NA | NA |
| heat2d | mx=200;my=100;velocity=50.0 | ARKODE IMEX4 | 4.83 | 2349 | 3437 | 1088 | 43595 | 22973 | NA | NA |
| heat2d | mx=200;my=100;velocity=50.0 | ARKODE IMEX5 | 5.3 | 1545 | 2307 | 762 | 33845 | 17696 | NA | NA |
| heat2d | mx=200;my=100;velocity=50.0 | ARKODE LSRK RKC2 | 6.55 | 657 | 889 | 232 | 20345 | NA | NA | 53 |
| heat2d | mx=200;my=100;velocity=50.0 | ARKODE LSRK RKC2 power | 6.68 | 971 | 1330 | 359 | 27353 | NA | NA | 72 |
| heat2d | mx=200;my=100;velocity=50.0 | ARKODE LSRK RKL2 | 6.67 | 344 | 478 | 134 | 17090 | NA | NA | 84 |
| heat2d | mx=200;my=100;velocity=500.0 | FIRKODE s=3 | 7.4 | 25 | 27 | 2 | 127 | NA | NA | NA |
| heat2d | mx=200;my=100;velocity=500.0 | ARKODE DIRK4 | 5.07 | 24 | 24 | 0 | 417 | NA | NA | NA |
| heat2d | mx=200;my=100;velocity=500.0 | ARKODE LSRK RKC2 | 9.41* | 183783 | 183783 | 0 | 367569 | NA | NA | 2 |
| heat2d | mx=200;my=100;velocity=500.0 | ARKODE LSRK RKC2 power | 5.46 | 55406 | 64329 | 8923 | 202978 | NA | NA | 6 |
| heat2d | mx=200;my=100;velocity=500.0 | ARKODE LSRK RKL2 | 5.15 | 42094 | 44567 | 2473 | 196432 | NA | NA | 22 |
| sector | m=100;rmax=1000.0;rmin=1.0;theta=0.0 | FIRKODE s=3 | 6.86 | 169 | 169 | 0 | 762 | NA | NA | NA |
| sector | m=100;rmax=1000.0;rmin=1.0;theta=0.0 | ARKODE DIRK4 | 5.58 | 207 | 207 | 0 | 3945 | NA | NA | NA |
| sector | m=100;rmax=1000.0;rmin=1.0;theta=0.0 | ARKODE LSRK RKC2 | 5.34 | 2899 | 2901 | 2 | 8714 | NA | NA | 4 |
| sector | m=100;rmax=1000.0;rmin=1.0;theta=0.0 | ARKODE LSRK RKL2 | 5.05 | 1819 | 1821 | 2 | 7280 | NA | NA | 4 |
| sector | m=100;rmax=1000.0;rmin=1.0;theta=30.0 | FIRKODE s=3 | 6.86 | 170 | 170 | 0 | 766 | NA | NA | NA |
| sector | m=100;rmax=1000.0;rmin=1.0;theta=30.0 | ARKODE DIRK4 | 5.51 | 209 | 209 | 0 | 3964 | NA | NA | NA |
| sector | m=100;rmax=1000.0;rmin=1.0;theta=30.0 | ARKODE LSRK RKC2 | 5.74 | 20447 | 20447 | 0 | 40897 | NA | NA | 2 |
| sector | m=100;rmax=1000.0;rmin=1.0;theta=30.0 | ARKODE LSRK RKL2 | 5.49 | 3628 | 3628 | 0 | 14515 | NA | NA | 8 |
<!-- end generated -->

![work-precision of the advection-diffusion problem](../../doc/shared/figs/firkode/benchmarks/arkode/arkode_heat2d_wp.png)

![LSRK on real and complex spectra](../../doc/shared/figs/firkode/benchmarks/arkode/arkode_sector_lsrk_wp.png)

**Verdict.** [pending]

### 5.8 A7: multirate (candidate 7)

**Purpose.** On a problem with a fast and a slow scale, does the single-rate FIRKODE pay for
resolving the fast scale with implicit steps where MRI resolves it explicitly?

**Sweep.** `kprmri` with G = -100, e = 0.5 and w in {10, 100, 1000}, rtol 1e-3 ... 1e-8 in
decades; FIRKODE, DIRK4, MRI ESDIRK34a with decoupled and with MRIHTol control, MRI
IMEX-SR32. 90 cells.

**Predictions [predicted].**
- P7.1: FIRKODE's step count grows >= 20x from w = 10 to 1000 while MRI's slow step count
  grows <= 3x and its fast step count >= 20x.
- P7.2: best MRI / FIRKODE <= 1/1.5 at w = 1000 (mescd 4 and 6) and >= 1 at w = 10.
- P7.3: MRIHTol control equals decoupled control unless the fast error estimate exceeds the
  slow tolerance: its fast-tolerance factor starts at 1 and is capped at 1
  (`INNER_MAX_TOLFAC`, `src/sunadaptcontroller/mrihtol/sunadaptcontroller_mrihtol.c`
  [verified]), so it can only tighten the fast tolerance.

**Pass and fail.** Supports the multirate candidate (an ARKODE win) if P7.1 and P7.2 hold;
refutes it if MRI's slow step count grows in proportion to w or MRI never wins.

**Results [pending].**

<!-- begin generated: A7_w -->
| w | class | t(best) / t(FIRKODE) at mescd=4 | t(best) / t(FIRKODE) at mescd=6 | t(best) / t(FIRKODE) at mescd=8 |
| --- | --- | --- | --- | --- |
| 10.0 | best DIRK default | NA | 1.09 (ARKODE DIRK4) | 1.88 (ARKODE DIRK4) |
| 10.0 | best DIRK of all | NA | 1.09 (ARKODE DIRK4) | 1.88 (ARKODE DIRK4) |
| 10.0 | best MRI | NA | 6.61 (ARKODE MRI IMEX-SR32) | 8.33 (ARKODE MRI ESDIRK34a) |
| 100.0 | best DIRK default | NA | 1.58 (ARKODE DIRK4) | NA |
| 100.0 | best DIRK of all | NA | 1.58 (ARKODE DIRK4) | NA |
| 100.0 | best MRI | NA | 1.38 (ARKODE MRI ESDIRK34a) | 2.34 (ARKODE MRI IMEX-SR32) |
| 1000.0 | best DIRK default | 1.99 (ARKODE DIRK4) | 2.49 (ARKODE DIRK4) | NA |
| 1000.0 | best DIRK of all | 1.99 (ARKODE DIRK4) | 2.49 (ARKODE DIRK4) | NA |
| 1000.0 | best MRI | 0.641 (ARKODE MRI ESDIRK34a htol) | 0.715 (ARKODE MRI ESDIRK34a) | NA |
<!-- end generated -->

<!-- begin generated: A7_counters -->
| problem | params | method | mescd | slow steps | fast steps | slow RHS (implicit) | slow RHS (explicit) | fast RHS | fast err fails |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| kprmri | G=-100.0;e=0.5;w=10.0 | FIRKODE s=3 | 8.3 | 313 | NA | 1452 | NA | NA | NA |
| kprmri | G=-100.0;e=0.5;w=10.0 | ARKODE DIRK4 | 6.77 | 197 | NA | 3966 | NA | NA | NA |
| kprmri | G=-100.0;e=0.5;w=10.0 | ARKODE MRI ESDIRK34a | 5.91 | 725 | 2196 | 8710 | NA | 11706 | 0 |
| kprmri | G=-100.0;e=0.5;w=10.0 | ARKODE MRI ESDIRK34a htol | 5.91 | 725 | 2196 | 8710 | NA | 11706 | 0 |
| kprmri | G=-100.0;e=0.5;w=10.0 | ARKODE MRI IMEX-SR32 | 6.25 | 412 | 2123 | 6328 | 1671 | 11003 | 0 |
| kprmri | G=-100.0;e=0.5;w=100.0 | FIRKODE s=3 | 7 | 3113 | NA | 14329 | NA | NA | NA |
| kprmri | G=-100.0;e=0.5;w=100.0 | ARKODE DIRK4 | 2.52 | 1965 | NA | 37148 | NA | NA | NA |
| kprmri | G=-100.0;e=0.5;w=100.0 | ARKODE MRI ESDIRK34a | 5.29 | 725 | 2638 | 9409 | NA | 13718 | 61 |
| kprmri | G=-100.0;e=0.5;w=100.0 | ARKODE MRI ESDIRK34a htol | 5.29 | 725 | 2638 | 9409 | NA | 13718 | 61 |
| kprmri | G=-100.0;e=0.5;w=100.0 | ARKODE MRI IMEX-SR32 | 6.27 | 417 | 9377 | 6908 | 1697 | 44594 | 1133 |
| kprmri | G=-100.0;e=0.5;w=1000.0 | FIRKODE s=3 | 6.12 | 31107 | NA | 143225 | NA | NA | NA |
| kprmri | G=-100.0;e=0.5;w=1000.0 | ARKODE DIRK4 | 0.84 | 19074 | NA | 352375 | NA | NA | NA |
| kprmri | G=-100.0;e=0.5;w=1000.0 | ARKODE MRI ESDIRK34a | 4.31 | 725 | 18356 | 9736 | NA | 83404 | 1766 |
| kprmri | G=-100.0;e=0.5;w=1000.0 | ARKODE MRI ESDIRK34a htol | 4.31 | 725 | 18356 | 9736 | NA | 83404 | 1766 |
| kprmri | G=-100.0;e=0.5;w=1000.0 | ARKODE MRI IMEX-SR32 | 4.21 | 450 | 88305 | 8082 | 1921 | 394202 | 9519 |
<!-- end generated -->

![work-precision of the multirate problem](../../doc/shared/figs/firkode/benchmarks/arkode/arkode_kprmri_wp.png)

**Verdict.** [pending]

### 5.9 A8: diurnal kinetics, IMEX with exact local preconditioning (candidate 8)

**Purpose.** With the stiff kinetics implicit (2 x 2 blocks, solved exactly by the block
preconditioner) and the transport explicit, IMEX needs no global Krylov iteration; FIRKODE
treats everything implicitly with a stacked Krylov solve. The explicit transport stiffens as
the grid refines.

**Sweep.** `diurnal` on grids 10, 20, 40, 80 (160 with `--include-large`), rtol 1e-4, 1e-6,
1e-8, SPGMR; FIRKODE, DIRK4, IMEX3/4/5. 60 cells.

**Predictions [predicted].**
- P8.1: IMEX's Krylov iterations per Newton iteration are <= 1.1, FIRKODE's stacked FGMRES
  iterations per Newton iteration >= 2.
- P8.2: IMEX4 / FIRKODE rises with the grid at rtol 1e-4 and 1e-6 and is <= 1 at grids 10
  and 20; IMEX's step count doubles or more from grid 40 to 80 at rtol 1e-4 (the explicit
  vertical diffusion: rate about 4 KV0 e^10 / dy^2), FIRKODE's grows by <= 1.3x.

**Results [pending].**

<!-- begin generated: A8_grid -->
| mx | class | t(best) / t(FIRKODE) at mescd=4 | t(best) / t(FIRKODE) at mescd=6 | t(best) / t(FIRKODE) at mescd=8 |
| --- | --- | --- | --- | --- |
| 10 | best DIRK default | NA | 4.13 (ARKODE DIRK4) | NA |
| 10 | best DIRK of all | NA | 4.13 (ARKODE DIRK4) | NA |
| 10 | best IMEX | NA | 9.68 (ARKODE IMEX3) | NA |
| 20 | best DIRK default | 2.82 (ARKODE DIRK4) | NA | NA |
| 20 | best DIRK of all | 2.82 (ARKODE DIRK4) | NA | NA |
| 20 | best IMEX | 3.99 (ARKODE IMEX4) | 6.87 (ARKODE IMEX3) | NA |
| 40 | best DIRK default | NA | NA | NA |
| 40 | best DIRK of all | NA | NA | NA |
| 40 | best IMEX | NA | 2.36 (ARKODE IMEX3) | NA |
| 80 | best DIRK default | NA | 0.866 (ARKODE DIRK4) | NA |
| 80 | best DIRK of all | NA | 0.866 (ARKODE DIRK4) | NA |
| 80 | best IMEX | NA | 0.727 (ARKODE IMEX4) | NA |
<!-- end generated -->

<!-- begin generated: A8_counters -->
| problem | params | method | mescd | steps | Newton | Krylov / Newton (ARKODE) | stacked FGMRES / Newton (FIRKODE) | prec solves |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| diurnal | mx=10;my=10 | FIRKODE s=3 | 6.96* | 359 | 600 | 1.07 | 5.47 | 10871 |
| diurnal | mx=10;my=10 | ARKODE DIRK4 | 4.33 | 449 | 3482 | 2.92 | NA | 13525 |
| diurnal | mx=10;my=10 | ARKODE IMEX3 | 4.58 | 1768 | 12464 | 1.13 | NA | 25278 |
| diurnal | mx=10;my=10 | ARKODE IMEX4 | 4.01 | 797 | 6955 | 1.73 | NA | 18191 |
| diurnal | mx=10;my=10 | ARKODE IMEX5 | 4.8 | 2133 | 18196 | 1.73 | NA | 47701 |
| diurnal | mx=20;my=20 | FIRKODE s=3 | 6.4 | 406 | 719 | 1.1 | 6.01 | 14184 |
| diurnal | mx=20;my=20 | ARKODE DIRK4 | 3.86 | 465 | 3793 | 3.1 | NA | 15433 |
| diurnal | mx=20;my=20 | ARKODE IMEX3 | 4.41 | 1864 | 13457 | 1.1 | NA | 26031 |
| diurnal | mx=20;my=20 | ARKODE IMEX4 | 3.66 | 818 | 8168 | 1.55 | NA | 19348 |
| diurnal | mx=20;my=20 | ARKODE IMEX5 | 4.28 | 2141 | 20599 | 1.57 | NA | 49210 |
| diurnal | mx=40;my=40 | FIRKODE s=3 | 5.86 | 381 | 699 | 1.49 | 8.52 | 19300 |
| diurnal | mx=40;my=40 | ARKODE DIRK4 | 5.17 | 522 | 4650 | 3.14 | NA | 18417 |
| diurnal | mx=40;my=40 | ARKODE IMEX3 | 4.29 | 947 | 6083 | 1.24 | NA | 13381 |
| diurnal | mx=40;my=40 | ARKODE IMEX4 | 3.58 | 415 | 3146 | 2.09 | NA | 9558 |
| diurnal | mx=40;my=40 | ARKODE IMEX5 | 4.03 | 1115 | 8871 | 1.84 | NA | 24911 |
| diurnal | mx=80;my=80 | FIRKODE s=3 | 7.17 | 473 | 856 | 1.49 | 8.74 | 24245 |
| diurnal | mx=80;my=80 | ARKODE DIRK4 | 5.73 | 482 | 4162 | 3.57 | NA | 18822 |
| diurnal | mx=80;my=80 | ARKODE IMEX3 | 4.79 | 842 | 5426 | 1.06 | NA | 10293 |
| diurnal | mx=80;my=80 | ARKODE IMEX4 | 4.68 | 428 | 3465 | 1.67 | NA | 8949 |
| diurnal | mx=80;my=80 | ARKODE IMEX5 | 4.83 | 932 | 7490 | 1.68 | NA | 19674 |
<!-- end generated -->

![work-precision of the diurnal problem](../../doc/shared/figs/firkode/benchmarks/arkode/arkode_diurnal_wp.png)

**Verdict.** [pending]

## 6. Summary

| family | candidate | verdict (best default / best of all) |
|---|---|---|
| A0 | method selection | [pending] |
| A1 | small stiff problems | [pending] |
| A2 | stage order | [pending] |
| A3 | mass matrices | [pending] |
| A4 | tight tolerances | [pending] |
| A5 | IMEX and LSRK on the 3-D Brusselator | [pending] |
| A6 | advection-diffusion over the velocity | [pending] |
| A7 | multirate | [pending] |
| A8 | IMEX on diurnal kinetics | [pending] |

## 7. Known issues and observations found by this suite

Nothing under `src/` or `include/` is changed for this suite.

- **ARKODE: `LSRKStepSetDomEigFrequency` with a negative argument** [verified]
  (`src/arkode/arkode_lsrkstep_io.c`, the function's two `if` blocks): for `nsteps < 0` it
  sets the default frequency, then the second block's `else` branch stores the negative
  value, so the dominant eigenvalue is recomputed every step; the guide says a negative
  value restores the default (25). The harness never passes a negative value.
- **ARKODE: `ARKStepSetOptimalParams`** is deprecated without a replacement and its parameter
  sets were tuned for older default tables [verified, `src/arkode/arkode_arkstep_io.c`].
- **FIRKODE: the default stacked FGMRES space on a fixed large step** [development]: on
  SECTOR with ten distinct real stiff modes (m = 10, theta = 0, r in [1, 1e3] or [1, 1e6])
  and fixed steps h >= 0.0078, the stacked FGMRES (6s = 18 iterations, no restart) fails to
  converge on the first step and FIRKODE stops with `FIRK_CONV_FAILURE`; the preconditioned
  stage operator has a distinct eigenvalue per mode and stage, and 15-30 iterations per
  Newton iteration are needed. Adaptive runs recover by reducing the step; the A2_fixed
  cells use 5 restarts.
- [pending: FIRKODE issues found by the sweep]

## 8. Limitations

One machine, one cell at a time per worker; parallel workers add noise. The DIRK class is a
declared list of nine configurations, not every table ARKODE has; the IMEX, LSRK and MRI
classes run ARKODE's defaults plus the stated options, with no per-problem tuning. The
splittings are the natural ones of each problem and are given to ARKODE only. `kprmri` has
N = 2 and equally expensive slow and fast parts, which understates what MRI saves on
problems whose fast part is cheap. FIRKODE's mass-matrix path and its loose-tolerance step
collapse are known weaknesses (`BENCHMARK_BDF.md`, Section 7) that can decide A3 and the
loosest cells. The extended-precision references need a build whose `long double` is wider
than double.

## 9. Reproduction

### 9.1 Commands

```bash
# builds (README: build-benchmark with KLU, build-benchmark-ref with extended precision)
PY=python3   # with numpy and matplotlib
BIN=build-benchmark/benchmarks/firkode/firkode_benchmark
REF=build-benchmark-ref/benchmarks/firkode/firkode_benchmark
RES=benchmarks/firkode/results/<date>-<machine>
# validation (replicas and splitting/eigenvalue checks)
$PY benchmarks/firkode/scripts/run_benchmarks.py --bin $BIN --results $RES --validate-only
# the families A0-A8, serial or with N workers (resumable)
$PY benchmarks/firkode/scripts/run_benchmarks.py --bin $BIN --ref-bin $REF --results $RES \
    --family arkode --workers 8
# report
$PY benchmarks/firkode/scripts/make_report.py --report arkode --results $RES \
    --figdir doc/shared/figs/firkode/benchmarks/arkode --inject benchmarks/firkode/BENCHMARK_ARKODE.md
```

### 9.2 Development runs

The predictions were written after the harness was developed on the author's workstation
(Apple silicon, macOS): single cells of every new integrator and problem, the validation, and
a smoke run of all families with a budget of 5 s per cell (`--fast 5`), stopped at about two
thirds of the cells and discarded; it found the FIRKODE fixed-step issue of Section 7, fixed
in the A2_fixed cell definition before any sweep. Statements informed by those runs carry the tag [development]. No prediction
is changed after the full sweep starts; an amendment before it is a dated commit.

## Appendix: CSV columns added for this report

`split`, `nfe_e`, `nfe_f`, `nst_f`, `nst_f_att`, `netf_f`, `ndomeig`, `nfe_domeig`,
`lsrk_max_stages`, `nst_expstab`, `nst_acc`, `split_residual`, `split_jac_err`,
`domeig_ratio`; the validation adds `exact`. The full dictionary is in README.md.
