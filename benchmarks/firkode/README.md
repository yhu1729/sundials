# FIRKODE benchmarks

This directory holds two serial benchmark programs for FIRKODE (Radau IIA
fully implicit Runge-Kutta methods):

- `firk_klin` measures the stacked Krylov iteration counts of the stage
  solve (the first section below);
- `firkode_benchmark` is the harness that compares FIRKODE with CVODE (BDF),
  ARKODE (ARKStep DIRK/ESDIRK tables and IMEX pairs, LSRKStep, MRIStep) and
  IDA (BDF, residual form) on a common set of stiff test problems (the
  second section). Its results are recorded in
  [`BENCHMARK_BDF.md`](BENCHMARK_BDF.md) (CVODE and IDA) and
  [`BENCHMARK_ARKODE.md`](BENCHMARK_ARKODE.md) (ARKODE).

## firk_klin: stacked Krylov iterations per Newton iteration

`firk_klin` measures the cost of the stacked flexible GMRES iteration of
FIRKODE (design document `src/firkode/DESIGN.md`, Section 11.4): the mean
number of FGMRES iterations per Newton iteration, `k_lin`, of the three-stage
method with its BEGINNER defaults, on

| problem  | system                                              | N          | block solver                                   |
| -------- | --------------------------------------------------- | ---------- | ---------------------------------------------- |
| `rober`  | Robertson kinetics to t = 1e11                      | 3          | dense                                          |
| `vdpol`  | van der Pol, eps = 1e-6, to t = 2                   | 2          | dense                                          |
| `orego`  | the Oregonator to t = 360                           | 3          | dense                                          |
| `heat2d` | 2-D heat equation, 5-point stencil, to t = 0.1      | nx * nx    | band (mu = ml = nx), or SPGMR with Jacobi      |

at `rtol` 1e-6 and 1e-8 (`atol` = 1e-4 `rtol` for the small problems, 1e-10
for the heat equation). With a direct block solver every block of the
stacked preconditioner is solved exactly; with a matrix-free block solver
the stacked iteration sees only the user's preconditioner (design decision
D2), here the Jacobi preconditioner of the heat equation.

### Running

Configure SUNDIALS with `-DSUNDIALS_ENABLE_BENCHMARKS=ON` (and
`-DCMAKE_BUILD_TYPE=Release`), build the target `firk_klin` and run

```bash
./benchmarks/firkode/firk_klin [nx]
```

`nx` (default 100) is the number of interior grid points per direction of
the heat equation. The output is one CSV line per run with the number of
steps, Newton iterations, stacked FGMRES iterations, their ratio `k_lin`,
the stacked convergence failures, the Jacobian evaluations and the wall-clock
time. The counts of the Krylov paths depend on the platform; the recorded
values are in the FIRKODE user guide (the chapter on the mathematics, the
section on the stage solve).

Note that the `heat2d` problem of `firk_klin` (pure diffusion to t = 0.1
with a rough initial condition) is not the `heat2d` problem of
`firkode_benchmark` (advection-diffusion with forcing to t = 1).

## firkode_benchmark: FIRKODE versus the other SUNDIALS implicit time steppers

The harness is one executable, `firkode_benchmark`, that runs a single
benchmark cell (problem, integrator, options) and prints a CSV header and
one data row, plus two Python scripts: `scripts/run_benchmarks.py` sweeps
the benchmark families and `scripts/make_report.py` turns the resulting
CSV files into the tables and figures of the report.

### Building

Configure SUNDIALS with FIRKODE and at least one of CVODE, ARKODE and IDA
enabled, `SUNDIALS_ENABLE_BENCHMARKS=ON`, double precision, and a Release
build type, for example

```bash
cmake -S . -B build-benchmark -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DSUNDIALS_ENABLE_FIRKODE=ON -DSUNDIALS_ENABLE_ARKODE=ON \
  -DSUNDIALS_ENABLE_CVODE=ON -DSUNDIALS_ENABLE_IDA=ON \
  -DSUNDIALS_ENABLE_BENCHMARKS=ON -DBENCHMARK_NVECTOR=OFF \
  -DSUNDIALS_ENABLE_C_EXAMPLES=OFF -DSUNDIALS_ENABLE_CXX_EXAMPLES=OFF \
  -DBUILD_TESTING=OFF -DSUNDIALS_ENABLE_PROFILING=OFF
cmake --build build-benchmark --target firkode_benchmark
```

The executable is `build-benchmark/benchmarks/firkode/firkode_benchmark`.

Reference solutions beyond the double-precision floor come from a second,
extended-precision build. On Apple silicon `long double` is plain double,
so that build targets x86_64 and runs under Rosetta (80-bit x87 arithmetic,
`LDBL_MANT_DIG = 64`); KLU cannot be enabled with extended precision:

```bash
cmake -S . -B build-benchmark-ref -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_ARCHITECTURES=x86_64 -DSUNDIALS_PRECISION=EXTENDED \
  -DSUNDIALS_ENABLE_FIRKODE=ON -DSUNDIALS_ENABLE_ARKODE=ON \
  -DSUNDIALS_ENABLE_CVODE=ON -DSUNDIALS_ENABLE_IDA=ON \
  -DSUNDIALS_ENABLE_BENCHMARKS=ON -DBENCHMARK_NVECTOR=OFF \
  -DSUNDIALS_ENABLE_C_EXAMPLES=OFF -DSUNDIALS_ENABLE_CXX_EXAMPLES=OFF \
  -DBUILD_TESTING=OFF -DSUNDIALS_ENABLE_PROFILING=OFF
cmake --build build-benchmark-ref --target firkode_benchmark
build-benchmark-ref/benchmarks/firkode/firkode_benchmark --version   # precision: extended, LDBL_MANT_DIG = 64
```

`--write-ref` of that binary writes 21 significant digits; the double build
reads such files rounded to double (`--ref FILE`, or `--ref-bin` of the
driver).

### Problems

| name        | system                                                             | N          | interval        | default solver      | notes                                  |
|-------------|--------------------------------------------------------------------|------------|-----------------|---------------------|----------------------------------------|
| `rober`     | Robertson kinetics                                                 | 3          | [0, 1e11]       | dense               | roots u = 1e-4, w = 1e-2 (`--roots`)   |
| `vdpol`     | Van der Pol, `y2' = ((1-y1^2) y2 - y1)/eps`                        | 2          | [0, 2]          | dense               | `--eps` (default 1e-6)                 |
| `orego`     | Oregonator                                                         | 3          | [0, 360]        | dense               |                                        |
| `hires`     | HIRES                                                              | 8          | [0, 321.8122]   | dense               |                                        |
| `analytic`  | `y' = lambda y + 1/(1+t^2) - lambda atan t`                        | 1          | [0, 10]         | dense               | `--lambda`; exact solution             |
| `pr`        | Prothero-Robinson (`test/unit_tests/problems/pr.hpp`)              | `--n`      | [0, 1]          | band (mu = ml = 0)  | `--lambda`; exact solution             |
| `kpr`       | Kvaerno-Prothero-Robinson (`test/unit_tests/problems/kpr.hpp`)     | 2          | [0, 1]          | dense               | exact solution                         |
| `diurnal`   | 2-species diurnal kinetics PDE, `--mx` x `--my` grid               | 2 mx my    | [0, 86400]      | SPGMR, left prec    | block-diagonal preconditioner, jtimes  |
| `heat2d`    | 2-D advection-diffusion, `--mx` x `--my` interior points           | mx my      | [0, 1]          | band (mu = ml = mx) | `--velocity`; matrix-free `--ls pcg`/`spgmr` with `--psolve jacobi\|innerpcg` |
| `heat1dfem` | 1-D heat, P1 finite elements with a tridiagonal mass matrix         | `--n`      | [0, 1]          | band                | mass matrix (FIRKODE, ARKODE, IDA)     |
| `wave1dfem` | damped 1-D elastodynamics, P1 finite elements: `[I 0; 0 M] y' = [v; -K u - C v + g(t) b]`, interleaved `(u_i, v_i)` | `2 n` | [0, 4] | band (`mu = ml = 3`), mass band 2 | `--n`, `--damping mass\|kv\|rayleigh`, `--alpha`, `--beta`/`--beta-omega`, `--mass consistent\|lumped` (lumped: ODE form, CVODE admissible); `scripts/wave1d_fem_exact.py` checks a reference against the modal closed form |
| `pollu`     | pollution kinetics of the Test Set for IVP Solvers (20 species, 25 reactions) | `20 nb` | [0, 60] | dense (`nb` = 1), band `mu = ml = 19` | `--nb` block-diagonal copies with perturbed initial states, `--tau` relaxation towards the initial state at every `--restart-every` boundary (family B11); published reference `references/pollu.ref` |
| `ringmod`   | ring modulator of the Test Set, ODE case (`Cs` = 2e-12) | 15 | [0, 1e-3] | dense | published reference `references/ringmod.ref`; the overflow guard of the source returns a recoverable error |
| `beam`      | elastic beam of the Test Set / E. Hairer (40 nodes) | 80 | [0, 5] | dense, difference-quotient Jacobian (no analytic one published) | published reference `references/beam.ref` (Test Set) cross-checked with Hairer's `res_exact_pic` (2.3e-11) |
| `plate`     | rectangular plate with a moving load (E. Hairer, 8 x 5 grid) | 80 | [0, 7] | dense | published reference `references/plate.ref`; atol = 1e-3 rtol as in the published driver |
| `detest_b5` | problem B5 of the stiff DETEST set (Enright, Hull and Lindberg 1975): `y1' = -10 y1 + 100 y2`, `y2' = -100 y1 - 10 y2`, `y3' = -4 y3`, `y4' = -y4`, `y5' = -0.5 y5`, `y6' = -0.1 y6`, `y(0) = 1`; eigenvalues `-10 +- 100 i` (84.29 degrees from the negative real axis), -4, -1, -0.5, -0.1 | 6 | [0, 20] | dense (band `mu = ml = 1`, KLU) | exact solution; transcribed from STDTST (ACM TOMS Algorithm 648), whose end-point values agree with it to 3.9e-14 (their y1, y2 carry that error; the exact values are about 1e-87) |
| `bruss3d`   | 3-D Brusselator on the unit cube (`a` = 1, `b` = 3.4), 7-point Laplacian with Neumann conditions, interleaved `(u, v)` | `2 n^3` | [0, 10] | KLU (CSR Jacobian); band/dense for small `n` | `--n`, `--alpha` (default 1), `--pert`; family B12 (factorization-dominated) |
| `sector`    | linear `y' = A (y - phi) + phi'` with eigenvalues `r_k (-cos theta +- i sin theta)`, `r_k` log-spaced in `[rmin, rmax]`, `A = Q blockdiag(B_k) Q` | `2 m` | [0, 10] | band (`mu = ml = 1`); dense if `--coupled 1` | exact solution; `--m --theta --rmin --rmax --coupled --delta --nonlin --forcing sine\|square --period` (the square forcing is B11's switched RC ladder) |
| `kprmri`    | the multirate Kvaerno-Prothero-Robinson problem of `examples/arkode/C_serial/ark_kpr_mri.c`: slow row `G tmp1 + e tmp2 + r'/(2u)`, fast row `e tmp1 - tmp2 + s'/(2 sqrt(2 + s))`, `r = 0.5 cos t`, `s = cos(w t)` | 2 | [0, 5] | dense | exact solution; `--G` (default -100, slow stiffness), `--w` (100, time-scale separation), `--e` (0.5, coupling); splittings `mri`, `mri-imex`, `mri-explicit` (MRIStep) |

Splittings for ARKStep IMEX and MRIStep (`--split NAME`, checked by
`--check-split`): `bruss3d` `imex` (explicit reaction, implicit linear
diffusion, its Jacobian on the KLU pattern of the full one); `heat2d` `imex`
(explicit advection and forcing, implicit diffusion), `mri` (slow implicit
diffusion and forcing, fast explicit advection) and `mri-imex` (slow explicit
forcing, slow implicit diffusion, fast advection); `diurnal` `imex`
(explicit transport, implicit kinetics with an exact block preconditioner);
`kpr` `imex` (the splitting of the ARKODE logging tests). Eigenvalue bounds
for LSRKStep (`--check-domeig`): exact for `analytic`, `pr`, `sector` and
`heat2d` (complex for a cell Peclet number above one), a Gershgorin bound
for `bruss3d`.

Errors of the PDE problems are measured against a tight-tolerance solution
of the same semi-discrete system, not against the PDE solution.

### Command line

```
firkode_benchmark --problem P --integrator I [options]
```

| option                                   | meaning                                                                    |
|------------------------------------------|----------------------------------------------------------------------------|
| `--rtol R`, `--atol A`, `--atol-vec a,b,…` | tolerances (problem defaults when omitted)                               |
| `--stages S`                             | FIRKODE stage count (this branch supports only the default, 3)             |
| `--order Q`                              | ARKODE DIRK order 1..5 (default 4); CVODE/IDA maximum BDF order            |
| `--ls dense\|band\|klu\|spgmr\|pcg`, `--prec`, `--maxl`, `--maxrs`, `--epslin` | linear solver settings; `klu` (sparse direct, CSR) needs a build with `SUNDIALS_ENABLE_KLU=ON` and a problem with a sparse Jacobian; IDA forms `cj I - J` in that pattern, which must then store the diagonal, and has no sparse mass matrix |
| `--msbp N`, `--msbj N`, `--dgmax D`, `--dgmax-jbad D` | setup and Jacobian reuse: setup frequency, Jacobian frequency, the relative change of gamma that forces a setup or a new Jacobian (CVODE; FIRKODE under `--reuse-policy periodic`; ARKODE honours `--msbp` and `--dgmax`) |
| `--reuse-policy radau5\|periodic`         | FIRKODE setup-reuse policy (EXPERT; default `radau5`, under which the four options above have no effect); the `method` column reads `RadauIIA-s3+periodic` |
| `--stage-maxl`, `--stage-maxrs`, `--stage-epslin` | FIRKODE stacked FGMRES settings                                   |
| `--maxcor`, `--nlscoef`                  | nonlinear solver settings                                                  |
| `--predictor extrapolate\|trivial`, `--no-refilter`, `--no-newton-safety`, `--controller NAME`, `--impgus-k K1 K2` | FIRKODE options (EXPERT mode) |
| `--mode beginner\|expert`                 | FIRKODE mode (default BEGINNER). Any EXPERT-only option switches to EXPERT, which also changes the step-size formula; the mode is recorded in the `mode` column |
| `--stald`                                | CVODE: BDF stability-limit detection (`CVodeSetStabLimDet`); the `method` column reads `BDF+STALD` |
| `--arkode-linear`, `--ark-predictor P`, `--ark-itable NAME`, `--ark-optimal` | ARKODE options: a linearly implicit problem, the predictor, the implicit table by name (exclusive with `--order`; the `method` column reads the table name, `+optimal`, `+linear`), and `ARKStepSetOptimalParams` (deprecated, tuned for older default tables) applied before the other options |
| `--integrator arkode-imex --split NAME`, `--ark-etable NAME` | ARKStep IMEX on an additive splitting of the problem (`fe` explicit, `fi` implicit with the splitting's Jacobian or product and preconditioner); `--order q` selects ARKODE's default pair, `--ark-itable` with `--ark-etable` a pair by name; `--arkode-linear` applies when the implicit part is linear; `nfe` counts `fi`, `nfe_e` counts `fe` |
| `--integrator arkode-lsrk`, `--lsrk-method rkc2\|rkl2`, `--dom-eig user\|power`, `--dom-eig-freq N`, `--lsrk-max-stages S`, `--dom-eig-safety F` | LSRKStep (RKC or RKL, order 2, M = I) with the problem's eigenvalue bound or the power iteration (100 iterations, relative tolerance 5e-3, 10 initial warm-ups, a deterministic start vector); the eigenvalue is computed once when the bound is constant (or the problem linear, for the power iteration) and every 25 steps otherwise; `stages` and `lsrk_max_stages` are the most stages of a step |
| `--integrator arkode-mri --split NAME`, `--mri-table NAME`, `--mri-fast-table NAME`, `--mri-fast-h H`, `--mri-controller decoupled\|htol` | MRIStep: slow scale on `fe` and `fi` of the splitting with an MRI coupling by name (default MRIStep's), fast scale `ff` with an explicit ARKStep inner stepper (an ERK table by name, default ARKStep's) at the same tolerances or a fixed fast step; the slow step is adaptive (an I controller; `htol` = MRIHTol with two I controllers, which tightens the fast tolerance only: its factor is capped at 1) or fixed (`--fixed-h`); `nst_f`, `nst_f_att`, `nfe_f`, `netf_f` count the fast scale; no `--restart-every` |
| `--maxsteps N`, `--maxnef N`, `--h0 H`, `--tf T`                | limits and final time                                       |
| `--grid K`, `--grid-geom T1 MULT`, `--no-tstop`                | output grid (used to replicate the examples); with the default stop-time handling every output time is a stop time, so `--grid K` makes the integrator land on K points exactly instead of interpolating |
| `--steps K`                              | K internal steps (one-step returns towards the final time, as the logging tests take them); no error metrics |
| `--check-split --split NAME`, `--check-domeig` | checks instead of a run: the pieces of a splitting add up to f and its implicit Jacobian matches differences of `fi` (`split_residual`, `split_jac_err`); the eigenvalue bound against the power iteration (`domeig_ratio`) |
| `--fixed-h H`                            | fixed-step run (FIRKODE, ARKODE)                                           |
| `--one-step --h H`                       | one adaptive step of size H; prints the error-estimate norms               |
| `--roots`                                | rootfinding run; prints the root times                                     |
| `--dense-k K`                            | dense output at K points, compared with an in-process reference            |
| `--dump-jac t1,t2,…`, `--dump-dir DIR`, `--dump-jac-dq` | write the dense Jacobian (N ≤ 2000) at the given times to `DIR/jac_<k>.txt` (the user Jacobian, or forward differences without one or with `--dump-jac-dq`) and the mass matrix to `DIR/mass.txt`; analysed by `scripts/spectrum.py` |
| `--restart-every K`                      | split `[t0, tf]` into K equal intervals; at each internal boundary the problem's `restart_update` hook is applied and the integrator is re-initialized (`*ReInit`, IDA with a consistent `y'`); counters are summed over the intervals and `nrestarts` = K-1; in-process references follow the same schedule |
| `--order-hist`                           | CVODE/IDA: step one internal step at a time and record the BDF order of every step (`q_hist`, `q_mean_step`, `q_mean_time`); timings of such runs are not comparable |
| `--ref FILE`, `--write-ref FILE`, `--ref-inprocess`, `--ref-check`, `--ref-integrator`, `--ref-check-integrator`, `--ref-rtol` | reference solution; `--ref-check` repeats the reference with a second integrator (CVODE, or IDA with a mass matrix) |
| `--reps R`                               | repetitions of create/evolve (min and median times are printed)            |
| `--profile`                              | wrap the linear solvers and the right-hand-side/Jacobian callbacks in timers (`t_setup`, `t_solve`, `t_rhs`, `t_jac`, `t_profiled`, `profiled` = 1); preconditioner setups and solves count as setups and solves; the wall times of such runs are not comparable with unprofiled ones |
| `--mx`, `--my`, `--n`, `--velocity`, `--lambda`, `--eps`, `--psolve`, `--psolve-rtol` | problem parameters   |
| `--KEY VALUE`                            | any further problem parameter; each problem declares the keys it reads and other keys are rejected |
| `--version`                              | SUNDIALS version and the FIRKODE default parameters                        |

Exit codes: 0 success, 1 harness or usage error, 2 unsupported combination
(for example CVODE on a mass-matrix problem), 3 the integrator returned a
failure flag (a CSV row is still printed).

### CSV columns

| column | meaning |
|--------|---------|
| `problem`, `integrator`, `method`, `order`, `stages`, `mode`, `ls`, `prec` | cell description; `method` is `RadauIIA-sS`, `BDF`, the ARKODE table name, or `IDA-BDF`; `mode` is FIRKODE's `beginner` or `expert` (`NA` otherwise) |
| `rtol`, `atol`, `N`, `fixed_h` | tolerances, system size, fixed step (NA in adaptive runs) |
| `nst`, `nst_att` | steps and step attempts (attempts: FIRKODE and ARKODE only) |
| `nfe`, `nfe_ls` | right-hand-side evaluations by the integrator and by difference-quotient Jacobian or Jacobian-vector products. FIRKODE counts the refilter evaluation and the initial-step heuristic; ARKODE counts implicit RHS evaluations; IDA counts residual evaluations (one RHS plus one mass product each) |
| `nsetups`, `nje` | linear solver setups and Jacobian evaluations |
| `nni`, `nncf`, `nsf`, `netf` | Newton iterations, nonlinear convergence failures, step attempts abandoned by the nonlinear solver, error test failures |
| `nstab` | CVODE only: order reductions by the stability-limit detection (0 unless `--stald`) |
| `nli`, `nlcf`, `npe`, `nps`, `njtv` | iterations and failures of the user's Krylov solver, preconditioner setups and solves, Jacobian-vector products |
| `nsli`, `nscf`, `nbs` | FIRKODE only: stacked FGMRES iterations, stacked solve failures, block solves with M - gamma J (including the preconditioner blocks of the stacked iteration). The mean number of FGMRES iterations per Newton iteration is `nsli/nni` |
| `nmass_setups`, `nmass_solves`, `nmass_mult` | mass-matrix setups, solves and products |
| `ngeval`, `nroots`, `root_times` | root-function evaluations, roots found, `t:dir0/dir1;…` |
| `nrestarts` | re-initializations of a `--restart-every` run (NA otherwise) |
| `err_mixed`, `err_wrms`, `err_abs_max` | errors at the final time: `max_i \|y_i - yref_i\| / (atol/rtol + \|yref_i\|)` (10 to the power of minus the IVP Test Set `mescd`), the WRMS norm with weights `1/(rtol \|yref\| + atol)`, and the maximum absolute error |
| `dense_err_mixed`, `dense_err_abs_max` | maxima over the dense-output samples |
| `hlast`, `ele_wrms`, `ele_max` | last step size and the norms of the local error estimate after one step (`--one-step`) |
| `q_hist`, `q_mean_step`, `q_mean_time` | `--order-hist` only: steps taken at BDF orders 1 to 5 (`n1;n2;n3;n4;n5`), the step-weighted and the step-size-weighted mean order |
| `nst_ref`, `ref_disagreement` | steps of the reference run and the mixed error between the reference and the cross-check integrator (`--ref-check`) |
| `time_min`, `time_median`, `reps` | wall-clock seconds of create/evolve, minimum and median over the repetitions |
| `t_setup`, `t_solve`, `t_rhs`, `t_jac`, `t_profiled`, `profiled` | `--profile` only: seconds in linear-solver setups (factorizations, preconditioner setups), solves, right-hand sides and Jacobians, and the wall time of that profiled repetition; `profiled` = 1 marks the row |
| `flag` | the integrator's return flag |
| `split` | the additive splitting of a run (`--split`, `arkode-imex` and `arkode-mri`; also the checked one of `--check-split`) |
| `nfe_e`, `nfe_f` | explicit right-hand sides (ARKStep IMEX `fe`, MRIStep slow explicit `fe`) and fast right-hand sides (MRIStep); `nfe` counts the implicit (IMEX, MRI slow implicit) or the only (LSRK) right-hand side |
| `nst_f`, `nst_f_att`, `netf_f` | fast steps, fast step attempts and fast error-test failures (MRIStep inner stepper) |
| `ndomeig`, `nfe_domeig`, `lsrk_max_stages` | LSRKStep: dominant-eigenvalue updates, right-hand sides of the power iteration (`--dom-eig power`), and the most stages of a step (last interval of a restarted run) |
| `nst_expstab`, `nst_acc` | ARKODE: steps limited by stability and by accuracy |
| `split_residual`, `split_jac_err` | `--check-split`: max over three times and two states of `\|f - (fe + fi + ff)\|` relative to the sum of the pieces' max norms, and of the implicit Jacobian (or Jacobian-vector products) against central differences of `fi`, relative to its max norm |
| `domeig_ratio` | `--check-domeig`: modulus of the problem's eigenvalue bound at `(t0, y0)` over the power iteration's estimate (>= 1 for an upper bound) |

The driver prepends `family`, `cell_key`, `label`, `params`, `driver_flag`
(`ok`, `TIMEOUT`, `UNSUPPORTED`, `SOLVER_FAIL_<flag>`, `FAIL_rc<n>`),
`timing_reliable` (0 when the run is shorter than 0.1 ms), `wall_s` and
`timestamp`.

### Running the suite

```bash
PY=python3   # needs numpy and matplotlib for the report
$PY benchmarks/firkode/scripts/run_benchmarks.py \
    --bin build-benchmark/benchmarks/firkode/firkode_benchmark \
    --results benchmarks/firkode/results/<run-id> --validate --family all
$PY benchmarks/firkode/scripts/make_report.py \
    --results benchmarks/firkode/results/<run-id> \
    --figdir doc/shared/figs/firkode/benchmarks/bdf \
    --inject benchmarks/firkode/BENCHMARK_BDF.md
$PY benchmarks/firkode/scripts/make_report.py --report arkode \
    --results benchmarks/firkode/results/<run-id> \
    --figdir doc/shared/figs/firkode/benchmarks/arkode \
    --inject benchmarks/firkode/BENCHMARK_ARKODE.md
```

`--family` takes family names (B1-B13, A0-A8) and the aliases `bdf` (the
families of `BENCHMARK_BDF.md`: B1, B8-B13), `arkode` (those of
`BENCHMARK_ARKODE.md`: A0-A8) and `all`; an unknown name is an error.
`--ref-bin PATH` computes the reference solutions with another binary (the
extended-precision build) and records it in `provenance.json`; B13 and A4
refuse to run without it.
`--fast [SECONDS]` runs the part of the suite that finishes fast: the full cell
lists are kept, cheap cells (by their time cap) run first, every cell is cut
off after SECONDS (default 10; a reference computation gets three times that)
and recorded with the driver flag `FAST_TIMEOUT`, and the repetitions are
limited to the budget. Such rows are redone by any later run of the same
results directory without `--fast`, and a reference that exceeded the budget
is recomputed then (`--workers` combines with it).
`--progress` prints the number of finished cells of the selected families
and the elapsed time, after every cell (with `--workers`, the parent reads
the worker CSV files every 60 s). `--validate` and the B3 logging run are not counted.
`--list` prints the number of cells of the selected families and exits; the
Jacobian dumps of the B10 gate land in `results/<run-id>/jac/<problem>/` (not
committed, like `reference/`, `logs/` and `claims/`).

### Parallel runs on a multi-core workstation

```bash
$PY benchmarks/firkode/scripts/run_benchmarks.py \
    --bin build-benchmark/benchmarks/firkode/firkode_benchmark \
    --results benchmarks/firkode/results/<run-id> --family all --workers 8 [--pin]
```

`--workers N` starts N worker processes (each runs one cell at a time with
`OMP_NUM_THREADS=1`; `--pin` binds worker k to a core with `taskset` where
available). The workers share the deterministic cell list of every family
and claim cells dynamically through atomic claim files in
`results/<run-id>/claims/`, expensive cells (by their time cap) first, so
the load balances without a fixed shard assignment. Each worker appends to
its own `<family>.worker<k>.csv`, `sentinel.worker<k>.csv` and
`provenance.worker<k>.json`; when all workers have finished the parent
merges them into `<family>.csv`, `sentinel.csv` and the `runs` list of
`provenance.json` (`--merge` does the same by hand, for example after
workers launched by a job scheduler with `--worker k --nworkers N`). Every
row records `worker` and `nworkers`. A reference solution is computed by
the first worker that needs it while the others wait on a lock file.
Resuming works as before: cells present in the merged file or in any
worker file are skipped; `--clear-claims` removes the claims a crashed
worker left behind (the parent does this before it starts its workers).

Concurrent workers share caches, memory bandwidth and the turbo budget of
the machine, so their timings carry more noise than a serial run; the
sentinel rows of every worker measure that, and the report marks rows
with `nworkers` > 1.
`--validate` first replicates repository examples and tests and compares
the counters with their `.out` files by the criterion `test/testRunner`
applies to answer files (integers within 10 percent, floats to 4 digits;
Krylov paths and replicas of a different program by a stated relative
tolerance). Root-search counts and long tight-tolerance runs differ in the
last bits between platforms; where the counters of the Linux workstation
that runs the full sweeps differ from an answer file, the replica expects the
Linux value (the note gives the answer file's), so the direct-solver replicas
match exactly on that workstation. `validation.csv` and the
console record whether each replica matched exactly. The validation also
runs `--check-split` on every splitting (residual below 1e-13, Jacobian
differences below 1e-6) and `--check-domeig` on every eigenvalue bound
(ratio within [0.99, 1.02] for the exact bounds, [0.99, 1.6] for the
Gershgorin bound of `bruss3d`). Sweeps are resumable: cells already present in
the family's CSV are skipped; `--retry-failed` repeats failed or timed-out
cells and `--force` repeats everything. Reference solutions are computed
once per problem and size into `results/<run-id>/reference/` (not committed).

See `results/README.md` for the layout of a results directory.
