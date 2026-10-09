# Benchmark results

One directory per run, named `<date>-<machine>`, written by
`../scripts/run_benchmarks.py` and summarised by `../scripts/make_report.py`.
Every number in `../BENCHMARK_BDF.md` comes from a row of these files.

| file | content |
|------|---------|
| `provenance.json` | machine, compiler, CMake cache, git commit, versions and arguments of every driver invocation |
| `validation.csv` | replicas of repository examples and tests with the expected counters and pass/fail |
| `sentinel.csv` | timing of a fixed cell before and after every family (drift detection) |
| `B8.csv` ... `B13.csv`, `B10_spectrum.csv` | the FIRKODE-versus-BDF families of `../BENCHMARK_BDF.md` |
| `<family>.worker<k>.csv`, `sentinel.worker<k>.csv`, `provenance.worker<k>.json`, `claims/` | per-worker files and cell claims of a parallel run (`--workers N`), folded into the files above by the merge step |
| `B1.csv` | work-precision on the small stiff problems (ROBER, VDPOL, OREGO, HIRES) |
| `B2.csv` | PDE problems at scale (diurnal kinetics, 2-D heat and advection-diffusion) |
| `B3.csv`, `B3_klin_log.csv` | FIRKODE stacked Krylov iterations versus stage count; the `_log` file holds per-solve maxima from a logging build |
| `B4_fixed.csv`, `B4_estimate.csv` | fixed-step convergence and error-estimate order |
| `B5.csv` | mass-matrix problem |
| `B6.csv` | FIRKODE option sensitivity |
| `B7_roots.csv`, `B7_dense.csv` | rootfinding and dense output |
| `tables/*.md` | generated Markdown fragments injected into the report |

The `reference/` (tight-tolerance reference solutions) and `logs/`
(stderr of failed cells, INFO logs) subdirectories are not committed.
