# radau5.f reference files

Reference runs of `radau5.f` (E. Hairer and G. Wanner) for the comparison of
FIRKODE with RADAU5 that design milestone M2 requires (`src/firkode/DESIGN.md`,
Section 10.3 and Table 6.2).

| Problem | Equations                                               | Final time | atol / rtol |
| ------- | ------------------------------------------------------- | ---------- | ----------- |
| `rober` | Robertson kinetics, y(0) = (1, 0, 0)                    | 1e11       | 1e-6        |
| `vdpol` | van der Pol, eps = 1e-6 (mu = 1000), y(0) = (2, 0)      | 2          | 1           |
| `orego` | the Oregonator, y(0) = (1, 2, 3)                        | 360        | 1           |

Each problem runs at the nominal RADAU5 tolerances 1e-4, 1e-6 and 1e-8, whose
files list the accepted steps (`step t_n h_n`, to 7 digits), and at 1e-14, whose
solution at the final time is the reference solution. RADAU5 replaces the
nominal rtol by 0.1 rtol^(2/3) and keeps atol / rtol; each file records the
nominal values (`rtol_nominal`, `atol_nominal`) and the values RADAU5 used
(`rtol`, `atol`), at which FIRKODE is compared. RADAU5 runs with its default
parameters, the analytic Jacobian and the initial step 1e-6. The 1e-14 solutions
agree with the published reference solutions (Test Set for IVP Solvers, and
y1(2) = 1.706167732170469 for van der Pol) to 2e-12 relative or better.

## Regenerating

```bash
python3 test/unit_tests/firkode/radau5/run_radau5.py
```

The script downloads `radau5.f`, `decsol.f` and `dc_decsol.f` from
<https://www.unige.ch/~hairer/prog/stiff/> into a temporary directory (or
`--workdir`), checks their SHA-256 against the digests it pins, compiles them
with `radau5_driver.f90` (gfortran, `--fc` and `--fflags` to change), and writes
the `.ref` files to this directory (or `--outdir`). The downloaded sources are
not part of SUNDIALS and are never written here. The first lines of each file
record the source digest, the compiler and the platform.

## The test

`firk_test_radau5_steps` (`test/unit_tests/firkode/CXX_serial`) reads from
each file the lines `n`, `tend`, `h0`, `rtol`, `atol`, `nstep`, `naccpt`,
`njac`, `nsol` and `y`, runs FIRKODE at the recorded `rtol` and `atol` with the
recorded `h0` and the final time as a stop time, and compares the number of
steps with `naccpt` and the error at the final time, in the WRMS norm of the
recorded tolerances, with RADAU5's error against the 1e-14 solution. The
`step` lines are not read.
