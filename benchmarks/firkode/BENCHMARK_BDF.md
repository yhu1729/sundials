# FIRKODE versus BDF (CVODE and IDA): six candidate regimes

Status: pre-registered. The question, the decision rule, the methodology and the predictions
and pass and fail rules of every family are fixed before the full sweep; the generated
sections are filled by `scripts/make_report.py` from the CSV files of `results/<run-id>/`.
Sections marked **[pending]** have no data yet; Section 9 gives the commands.

Evidence tags: **[measured]** = a row of a CSV file in `results/`, **[predicted]** = written
down before the runs (Section 1.2 and the prediction block of every family), **[verify]** =
to be checked against a primary source before it is used.

## 1. Question and decision rule

### 1.1 Question

At matched accuracy, with symmetrically tuned baselines, is FIRKODE (Radau IIA, s = 3, order
5) faster than the BDF codes of SUNDIALS, CVODE and IDA, on six problem families, and can
every win or loss be attributed to a named mechanism? The report asks where the one-step
fully implicit method pays.

### 1.2 Decision rule

A family **supports** its candidate if FIRKODE s = 3 is faster than the best BDF configuration
(the fastest CVODE or IDA configuration of the cell group) at matched mixed digits (mescd =
-log10 of the mixed error at the final time, interpolated in log time) in at least one regime
**and** the counters show the predicted mechanism. It **refutes** the candidate if the best
BDF configuration stays faster across the whole sweep.
Timing noise is measured by the sentinel spread (Section 3.3); a win is called only at
>= 1.2x and when adjacent tolerance cells agree. Concurrent workers (`--workers N`) add noise;
rows record `nworkers`.

### 1.3 Cost model behind the candidates [predicted]

Each FIRKODE Newton iteration costs 3 right-hand sides plus 3 k_lin block solves and
Jacobian-vector products, where k_lin is the mean number of stacked FGMRES iterations per
Newton iteration; CVODE costs one right-hand side and one solve, IDA one residual (a
right-hand side and, with a mass matrix, a mass-matrix product) and one solve. All three
factorize one N x N matrix per setup. The step-count ratio BDF/FIRKODE needed for FIRKODE to
break even is about 4x when right-hand sides dominate, 10-15x when block solves dominate,
about 1x when factorizations dominate, and about 20x when the fixed per-step overhead of a
small system dominates. The six families push the ratio past the threshold (BDF limited by
stability or by order ramps after restarts) or move the cost into factorizations.

IDA [predicted]: IDA uses the same BDF formulas with the same A(alpha) sectors per order but
its own order and step selection, solves for y' as well as y, has no stability-limit
detection and evaluates the Jacobian at every setup. On the ODE problems its time at matched
accuracy stays within a factor of 2 of CVODE default's. It changes a family's verdict only
where its order selection copes with a spectrum beyond a BDF sector better than CVODE's; it
has no counterpart to STALD.

## 2. Environment and provenance

<!-- begin generated: provenance -->
[pending]
<!-- end generated -->

## 3. Methodology

### 3.1 Harness

`firkode_benchmark` (this directory) runs one cell per process and prints one CSV row; the
driver `scripts/run_benchmarks.py` sweeps the families, computes references, times the
sentinel cell before and after each family, and runs serially or with `--workers N`
processes claiming cells dynamically (README).

- FIRKODE runs in its BEGINNER mode unless an EXPERT-only option is requested; the `mode`
  column records which step-size rule ran (EXPERT uses a different one). FIRKODE implements
  s = 3 only.
- References are FIRKODE s = 3 at a tight tolerance cross-checked by CVODE (IDA for a
  mass-matrix problem); the mixed disagreement is the floor of every cell (`*` in the
  tables, hollow markers in the figures). Family B13 uses an extended-precision build
  (80-bit, `--ref-bin`) for its references; SECTOR has an exact solution; the published
  problems use the values of the Test Set for IVP Solvers and of E. Hairer's test set
  (`references/`).
- Instruments: CVODE stability-limit detection (`--stald`, `nstab`), BDF order histograms
  (`--order-hist`, untimed twins labelled `[order]`), restarts with re-initialization
  (`--restart-every K`, `nrestarts`), the KLU sparse direct solver (`--ls klu`), a profiling
  shim (`--profile`, untimed twins labelled `[profile]`), Jacobian dumps with the spectrum
  analysis `scripts/spectrum.py`, and the setup-reuse options (`--reuse-policy periodic`,
  `--msbp`, `--dgmax`).

### 3.2 Tolerances and references

| problem | atol / rtol | reference | cross-check |
|---|---|---|---|
| ROBER | 1e-6 | FIRKODE s = 3 at 1e-13 (B13: 1e-16, extended) | CVODE |
| VDPOL, OREGO, HIRES | 1 | FIRKODE s = 3 at 1e-13 (B13: 1e-16, extended) | CVODE |
| SECTOR, rc_pwm | 1 | exact | - |
| wave1dfem | 1e-4 | IDA at 1e-12 | FIRKODE (lumped mass) |
| POLLU | 1 | Test Set `solut` (nb = 1, no restarts); otherwise FIRKODE at 1e-13 with the same restart schedule | CVODE |
| RINGMOD | 1 | Test Set `solut` | - |
| BEAM | 1 | Test Set `solut` | Hairer `res_exact_pic` |
| PLATE | 1e-3 | Hairer `res_exact_pic` | - |
| bruss3d | 1e-2 | FIRKODE s = 3 with KLU at 1e-10 (1e-9 for n = 32) | CVODE |
| heat2d | 1e-4 | FIRKODE s = 3 at 1e-10 | CVODE |

<!-- begin generated: references -->
[pending]
<!-- end generated -->

### 3.3 Timing protocol

Wall-clock seconds of setup plus evolve, `std::chrono::steady_clock`, repeated in one process
(pilot run, then 5-200 repetitions for about one second, a single repetition above 10 s);
the median is reported, the minimum is in the CSV. One cell at a time per worker,
`OMP_NUM_THREADS=1`. The sentinel cell (ROBER, FIRKODE s = 3, rtol 1e-6, 20 repetitions)
runs before and after every family in every worker; its rows are in
`results/<run-id>/sentinel.csv`, and `make_report.py` writes their spread to
`results/<run-id>/tables/sentinel.md`.

### 3.4 Fairness caveats

1. Baselines run at their defaults; the only tuned variants are the ones the family sweeps
   (STALD, a limited maximum order, the setup frequency, a matched initial step), and the
   matching FIRKODE knob is swept where one exists (periodic reuse policy, initial step).
2. FIRKODE solves the stage systems inexactly with a stacked FGMRES preconditioned by block
   solves; CVODE solves one N x N system directly. Counters are compared on the work they
   measure; the bottom line is wall time at matched accuracy.
3. `nfe` means different things: FIRKODE counts the refilter evaluation and the initial-step
   heuristic, IDA counts residuals.
4. FIRKODE's periodic reuse policy is an EXPERT-mode variant and runs with the EXPERT
   step-size rule; such cells are labelled `periodic`.
5. BEAM has no published Jacobian: every integrator uses dense difference quotients there.
6. IDA solves the residual form `M y' - f(t, y) = 0` of every problem, the ODE problems
   included (`M = I`), and runs every CVODE configuration it has an option for: the maximum
   order and the initial step, but neither STALD nor a setup frequency. It evaluates the
   Jacobian at every setup, because its iteration matrix `cj M - J` changes with `cj`; CVODE
   reuses J across setups.

### 3.5 Harness validation

`--validate` replicates repository examples and tests and compares their counters with the
answer files by the criterion `test/testRunner` applies to them (integers within 10 percent,
floats to 4 digits; Krylov paths within a stated relative tolerance); the table marks the
replicas that match exactly. Where the counters of the Linux workstation that runs the sweep
differ from an answer file in the last bits (a root search, a long run at a tight tolerance),
the replica expects the Linux value, so that the direct-solver replicas match exactly on that
workstation. The sweeps run only if every replica passes.

<!-- begin generated: validation -->
[pending]
<!-- end generated -->

## 4. Problems

| name | system | N | interval | family |
|---|---|---|---|---|
| `sector` | `y' = A (y - phi) + phi'`, eigenvalues `r_k (-cos theta +- i sin theta)`, `r_k` log-spaced | 2 m | [0, 10] | B8, B11 (square forcing), B13 |
| `wave1dfem` | P1 finite elements of `M u'' + C u' + K u = g(t) b`, damping `alpha M`, `beta K` or both, consistent or lumped mass | 2 n | [0, 4] | B9 |
| `beam`, `plate`, `ringmod` | published stiff problems (Test Set for IVP Solvers, E. Hairer) | 80, 80, 15 | [0, 5], [0, 7], [0, 1e-3] | B10 |
| `pollu` | pollution kinetics, `nb` block-diagonal copies, relaxation towards the initial state at every restart | 20 nb | [0, 60] | B11 |
| `bruss3d` | 3-D Brusselator, 7-point Laplacian, CSR Jacobian for KLU | 2 n^3 | [0, 10] | B12, B13 |
| ROBER, VDPOL, OREGO, HIRES, `heat2d` | small stiff problems and 2-D advection-diffusion (README) | 3, 2, 3, 8, mx my | README | B1, B13, B12 (control) |

## 5. Results

### 5.1 B8: spectral-angle sweep (candidate 1, the A(alpha) mechanism)

**Purpose.** Isolate the A(alpha) mechanism with an exact solution and one parameter, the
angle theta of the stiff eigenvalues from the negative real axis. BDF3/4/5 are
A(alpha)-stable with alpha = 86.03, 73.35 and 51.84 degrees (computed by `spectrum.py` from
the boundary locus and matching Hairer-Wanner II, Section V.2); a spectrum spread over many
magnitudes at an angle beyond alpha_q always puts some h lambda in the unstable lobe of order q.

**Sweep.** theta in {0, 30, 45, 50, 55, 60, 70, 75, 80, 85, 88, 89.5} degrees, r in [1, 1e3]
and [1, 1e6], m = 100 (N = 200), rtol 1e-4 ... 1e-10 in half decades, atol = rtol;
configurations FIRKODE s = 3, CVODE default, CVODE with STALD, and CVODE with the maximum
order limited to 4, 3 and 2 (the A(73.35 degree), A(86.03 degree) and A-stable members of
the family, so that every angle has a BDF configuration whose stability sector covers the
spectrum), IDA default and IDA with the same order limits; order histograms of CVODE, CVODE
with STALD and IDA at six angles; N = 2000 at three angles (FIRKODE, CVODE, CVODE with STALD,
IDA); a coupled (dense) and a perturbed variant at 80 degrees.

**Predictions [predicted].** theta <= 50 degrees: CVODE 3-5x faster. 52-86 degrees: CVODE's
mean order falls as theta passes 52 and 73 degrees, error-test failures rise without STALD,
the ratio moves towards FIRKODE. theta >= 86 degrees: the best BDF configuration is
effectively order 2, FIRKODE s = 3 is faster at rtol <= 1e-6 and the gap grows as rtol
tightens, roughly as tol^(-1/2) against tol^(-1/5) in steps. FIRKODE's k_lin stays bounded
(2.5-5) for all theta. IDA default, without STALD, behaves like CVODE default beyond 73
degrees; its order-limited runs track CVODE's.

**Pass and fail.** Supports candidate 1 if for some theta >= 75 degrees and rtol <= 1e-6
t(best BDF) / t(FIRKODE s = 3) > 1 and the time-weighted BDF order there is <= 2.5. Refutes it
if the best BDF configuration stays faster at 89.5 degrees for every rtol.

**Results [pending].**

<!-- begin generated: B8_ratio_theta -->
[pending]
<!-- end generated -->

![t(best BDF) / t(FIRKODE) against the angle](../../doc/shared/figs/firkode/benchmarks/bdf/bdf_sector_ratio.png)

![work-precision at four angles](../../doc/shared/figs/firkode/benchmarks/bdf/bdf_sector_wp.png)

<!-- begin generated: B8_orders -->
[pending]
<!-- end generated -->

<!-- begin generated: B8_variants -->
[pending]
<!-- end generated -->

**Verdict.** [pending]

### 5.2 B9: damped elastodynamics with a consistent mass matrix (candidate 2)

**Purpose.** A physical version of B8 that uses both FIRKODE strengths (oscillatory spectrum,
mass matrix) with a built-in control that moves the spectrum to the real axis: mass damping
`alpha M` leaves the high modes near the imaginary axis, Kelvin-Voigt damping `beta K` makes
them real.

**Sweep.** n in {100, 400, 1600}; lumped mass (ODE form): mass and kv damping with FIRKODE,
CVODE, CVODE STALD, CVODE q <= 2, IDA and IDA q <= 2; consistent mass: mass, kv and Rayleigh
damping with FIRKODE, IDA and IDA q <= 2 (CVODE cannot take a mass matrix); rtol 1e-3 ...
1e-9 in decades, atol = 1e-4 rtol; order histograms of CVODE and IDA under mass damping.

**Predictions [predicted].** Mass damping: the picture of B8 at theta >= 86 degrees, FIRKODE
faster at rtol <= 1e-5 and more so at larger n. Kelvin-Voigt damping: the BDF codes faster by
about the B1 factor (the control). Consistent mass: IDA's advantage on the finite-element
heat equation (B5) disappears or reverses under mass damping only. On the lumped cells IDA
and CVODE give the same verdict.

**Pass and fail.** Supports candidate 2 if FIRKODE wins under mass damping and loses under
Kelvin-Voigt damping at the same n and rtol.

**Results [pending].**

<!-- begin generated: B9_ratio -->
[pending]
<!-- end generated -->

![work-precision of the wave problem](../../doc/shared/figs/firkode/benchmarks/bdf/bdf_wave_wp.png)

<!-- begin generated: B9_orders -->
[pending]
<!-- end generated -->

<!-- begin generated: B9_wp -->
[pending]
<!-- end generated -->

**Verdict.** [pending]

### 5.3 B10: published stiff-oscillatory problems (candidate 3)

**Purpose.** External validity on problems with published references.

**Gate B10.0.** Each problem is transcribed from its primary source (`scripts/fetch_testset.py`
pins the digests) and must reproduce its published reference at a tight tolerance. The
spectrum of each problem along a tight trajectory:

<!-- begin generated: B10_spectrum -->
[pending]
<!-- end generated -->

Only problems with a substantial share of the spectrum beyond 73 degrees carry the candidate;
a problem that fails the gate is reported as "not oscillatory-stiff", which is a result.

**Sweep.** FIRKODE s = 3, CVODE default, CVODE STALD, CVODE q <= 2, IDA default and IDA
q <= 2; rtol 1e-3 ... 1e-10 in half decades; order histograms of CVODE and IDA at 1e-6 and
1e-8.

**Predictions [predicted].** Problems that pass the gate behave like B8 at the matching
angle. Supports candidate 3 if at least one published problem shows a FIRKODE win with the B8
attribution (BDF order <= 2.5).

**Results [pending].**

<!-- begin generated: B10_ratio -->
[pending]
<!-- end generated -->

![work-precision of the published problems](../../doc/shared/figs/firkode/benchmarks/bdf/bdf_published_wp.png)

<!-- begin generated: B10_orders -->
[pending]
<!-- end generated -->

<!-- begin generated: B10_wp -->
[pending]
<!-- end generated -->

**Verdict.** [pending]

### 5.4 B11: restart-heavy workloads (candidate 4)

**Purpose.** Measure the BDF order ramp after each re-initialization against a one-step
method that restarts at full order.

**Workloads.** `pollu` with Lie splitting: chemistry over an interval, then the exact
relaxation `y <- y + (1 - exp(-dt/tau)) (y_bg - y)` towards the initial state (tau = 10),
then a re-initialization; K in {1, 10, 100, 1000} intervals over [0, 60]; `nb` in {1, 64,
1024} independent copies (band solver, mu = ml = 19). `rc_pwm`: the square-wave forcing of
SECTOR (theta = 0, m = 25, rates in [1, 1e6], period 10/K) with a re-initialization at every
switch, where the forcing is discontinuous; at the restart time the forcing takes its right
limit, so every interval is an initial value problem with its own forcing. FIRKODE, CVODE and
IDA run every workload with their own initial-step heuristics; on the single POLLU system a
variant passes the same initial step 1e-6 to each. Order histograms of CVODE and IDA at
K >= 100.

**Predictions [predicted].** CVODE's steps per interval approach a floor of 10-20 as the
intervals shrink (the order ramp from q = 1), FIRKODE's floor is a few steps. At N = 20 the
per-step overhead hides the gain; with nb >= 64 it is amortized and FIRKODE wins once CVODE's
steps per interval exceed about 4x FIRKODE's. IDA restarts at order 1 like CVODE and shows the
same floor.

**Pass and fail.** Supports candidate 4 if at some K and nb FIRKODE is faster at matched
final accuracy and the best BDF configuration's steps per interval sit at the ramp floor
(most intervals end before q reaches 5 in the order histograms).

**Results [pending].**

<!-- begin generated: B11_ratio -->
[pending]
<!-- end generated -->

![time against the number of intervals](../../doc/shared/figs/firkode/benchmarks/bdf/bdf_restarts.png)

<!-- begin generated: B11_restarts -->
[pending]
<!-- end generated -->

<!-- begin generated: B11_orders -->
[pending]
<!-- end generated -->

**Verdict.** [pending]

### 5.5 B12: factorization-dominated 3-D reaction-diffusion (candidate 5)

**Purpose.** The regime where setups dominate and FIRKODE's fewer setups might pay, even on a
real spectrum.

**Problem.** `bruss3d` with alpha = 1 (h rho(J) grows like n^2), n in {8, 12, 16, 24} (N up
to 27648; n = 32 with `--include-large`), KLU on the fixed CSR pattern; rtol 1e-4 ... 1e-8;
FIRKODE default (RADAU5 reuse policy), FIRKODE with the periodic policy at msbp 20 and 50,
CVODE default (msbp 20), CVODE msbp 50, IDA default (KLU on `cj I - J`; IDA has no
setup-frequency option); a profiled twin of every configuration at rtol 1e-6. Control:
`heat2d` with a band solver at 100 x 100 and 200 x 200, where a factorization costs only
about the bandwidth more than a solve, with the same configurations but FIRKODE periodic at
msbp 20.

**Predictions [predicted].** The share of time in factorization grows with n and exceeds 80
percent for n >= 16. Where it exceeds about 90 percent, the time ratio tracks the setup
ratio; FIRKODE wins by 1.1-1.5x if its setup count stays lower. Refuted if the setup counts
converge once both use the same periodic frequency. IDA's setup count is close to CVODE
default's; its Jacobian at every setup costs little next to a KLU factorization.

**Results [pending].**

<!-- begin generated: B12_ratio -->
[pending]
<!-- end generated -->

<!-- begin generated: B12_profile -->
[pending]
<!-- end generated -->

![time split of the profiled runs](../../doc/shared/figs/firkode/benchmarks/bdf/bdf_bruss3d_profile.png)

<!-- begin generated: B12_wp -->
[pending]
<!-- end generated -->

**Verdict.** [pending]

### 5.6 B13: tight tolerances (candidate 6)

**Purpose.** Where does the order-5 one-step method pay against BDF at tight tolerances?
FIRKODE implements s = 3 only, so the family compares methods of equal order with different
error constants and order selection.

**Sweep.** ROBER, VDPOL, OREGO, HIRES, SECTOR at 0 and 80 degrees, `bruss3d` n = 16; rtol
1e-8 ... 1e-13 in half decades (`bruss3d` 1e-8 ... 1e-11); FIRKODE s = 3, CVODE default,
IDA default, CVODE STALD (SECTOR 80 degrees). References: FIRKODE s = 3 at rtol 1e-16 in
80-bit extended precision cross-checked by CVODE at the same tolerance (ODE problems), exact
(SECTOR), double precision with KLU at 1e-10 (`bruss3d`, floor marked).

**Predictions [predicted].** Both methods are order 5, so the step ratio does not grow as the
tolerance tightens; the ratio stays near its mescd-8 value unless BDF's order selection falls
below 5 or STALD reductions appear. A crossover exists only where the B8 or B12 mechanisms act
(SECTOR at 80 degrees, `bruss3d`). The report gives the crossover mescd per problem, or "none
up to the roundoff floor".

**Results [pending].**

<!-- begin generated: B13_crossover -->
[pending]
<!-- end generated -->

![work-precision at tight tolerances](../../doc/shared/figs/firkode/benchmarks/bdf/bdf_tight_wp.png)

<!-- begin generated: B13_wp -->
[pending]
<!-- end generated -->

**Verdict.** [pending]

### 5.7 Baseline: the small stiff problems (B1)

The matched-accuracy ratios of the other SUNDIALS implicit integrators (CVODE, IDA and
ARKODE DIRK) against FIRKODE s = 3:

<!-- begin generated: B1_rober_speedup -->
[pending]
<!-- end generated -->

<!-- begin generated: B1_vdpol_speedup -->
[pending]
<!-- end generated -->

<!-- begin generated: B1_orego_speedup -->
[pending]
<!-- end generated -->

<!-- begin generated: B1_hires_speedup -->
[pending]
<!-- end generated -->

## 6. Summary

| family | candidate | verdict |
|---|---|---|
| B8 | A(alpha) mechanism | [pending] |
| B9 | mass matrix plus oscillatory spectrum | [pending] |
| B10 | published oscillatory-stiff problems | [pending] |
| B11 | restarts | [pending] |
| B12 | factorization-dominated | [pending] |
| B13 | tight tolerances | [pending] |

## 7. Known FIRKODE issues and observations found by this suite

Nothing under `src/firkode/` is changed for this suite; issues it finds are recorded here for
the package. [pending]

## 8. Limitations

One machine, one cell at a time per worker; parallel workers add noise (sentinel rows,
Section 3.3); the published problems are transcriptions, verified against their references
by gate B10.0 only; BEAM uses difference-quotient Jacobians everywhere; IDA runs the ODE
problems in residual form, so its times include that form's overhead; the extended-precision
references need a build whose `long double` is wider than double (x86_64, natively or under
Rosetta on Apple silicon); `ringmod` cells are capped at 900 s.

## 9. Reproduction

```bash
# builds (README: build-benchmark with KLU, build-benchmark-ref with extended precision)
PY=python3   # with numpy and matplotlib
BIN=build-benchmark/benchmarks/firkode/firkode_benchmark
REF=build-benchmark-ref/benchmarks/firkode/firkode_benchmark
RES=benchmarks/firkode/results/<date>-<machine>
# validation
$PY benchmarks/firkode/scripts/run_benchmarks.py --bin $BIN --results $RES --validate-only
# the full families, serial or with N workers (resumable)
$PY benchmarks/firkode/scripts/run_benchmarks.py --bin $BIN --ref-bin $REF --results $RES \
    --family B1 B8 B9 B10 B11 B12 B13 --workers 8
# report
$PY benchmarks/firkode/scripts/make_report.py --results $RES \
    --figdir doc/shared/figs/firkode/benchmarks/bdf --inject benchmarks/firkode/BENCHMARK_BDF.md
```

## Appendix: CSV columns added by this suite

`mode` (FIRKODE beginner/expert), `nstab`, `q_hist`, `q_mean_step`, `q_mean_time`,
`nrestarts`, `t_setup`, `t_solve`, `t_rhs`, `t_jac`, `t_profiled`, `profiled`; driver columns
`worker`, `nworkers`. The full dictionary is in README.md.
