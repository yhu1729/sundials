#!/usr/bin/env python3
# -----------------------------------------------------------------
# Programmer(s): Yifan Hu @ UMBC
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------
# Sweep driver for the FIRKODE benchmark suite. Runs firkode_benchmark over
# the benchmark families, writes one CSV per family under the results
# directory, and is resumable (cells already recorded are skipped).
#
#   run_benchmarks.py --bin build-benchmark/benchmarks/firkode/firkode_benchmark \
#       --results benchmarks/firkode/results/<run-id> [--validate] \
#       [--family B1 B4 ... | --family all] [--fast [SECONDS]] [--include-large] \
#       [--log-bin build-benchmark-log/benchmarks/firkode/firkode_benchmark] \
#       [--workers N] [--progress]
# -----------------------------------------------------------------
"""Sweep driver for the FIRKODE benchmark suite."""

from __future__ import annotations

import argparse
import csv
import dataclasses
import datetime as dt
import hashlib
import io
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from benchmark_common import (
    ATOL_FACTOR,
    REF_RTOL,
    REF_RTOL_DOUBLE_FLOOR,
    TIMING_FLOOR_S,
    Cell,
    ResultsTable,
    binary_version,
    claim_cell,
    clear_claims,
    collect_provenance,
    ReferenceTimeout,
    ensure_reference,
    family_files,
    merge_results,
    read_csv,
    run_binary,
    save_json,
    to_float,
    to_int,
)

REPO = Path(__file__).resolve().parents[3]
MAXSTEPS = ["--maxsteps", "2000000"]

# Integrator configurations: key -> (integrator, label, extra args)
CONFIGS = {
    # this branch of FIRKODE supports the three-stage method only
    "firk3": ("firkode", "FIRKODE s=3", ["--stages", "3"]),
    "cvode": ("cvode", "CVODE BDF", []),
    "ark2": ("arkode", "ARKODE DIRK2", ["--order", "2"]),
    "ark3": ("arkode", "ARKODE DIRK3", ["--order", "3"]),
    "ark4": ("arkode", "ARKODE DIRK4", ["--order", "4"]),
    "ark5": ("arkode", "ARKODE DIRK5", ["--order", "5"]),
    "ida": ("ida", "IDA BDF", []),
    # baselines of the FIRKODE-versus-BDF families (BENCHMARK_BDF.md)
    "cvode_stald": ("cvode", "CVODE BDF+STALD", ["--stald"]),
    "cvode_q2": ("cvode", "CVODE BDF q<=2", ["--order", "2"]),
    "cvode_q3": ("cvode", "CVODE BDF q<=3", ["--order", "3"]),
    "cvode_q4": ("cvode", "CVODE BDF q<=4", ["--order", "4"]),
    "cvode_msbp50": ("cvode", "CVODE BDF msbp=50", ["--msbp", "50"]),
    "cvode_h0": ("cvode", "CVODE BDF h0=1e-6", ["--h0", "1e-6"]),
    "ida_q2": ("ida", "IDA BDF q<=2", ["--order", "2"]),
    "ida_q3": ("ida", "IDA BDF q<=3", ["--order", "3"]),
    "ida_q4": ("ida", "IDA BDF q<=4", ["--order", "4"]),
    "ida_h0": ("ida", "IDA BDF h0=1e-6", ["--h0", "1e-6"]),
    "firk3_h0": ("firkode", "FIRKODE s=3 h0=1e-6", ["--stages", "3", "--h0", "1e-6"]),
    "firk3_periodic": (
        "firkode",
        "FIRKODE s=3 periodic msbp=20",
        ["--stages", "3", "--reuse-policy", "periodic"],
    ),
    "firk3_periodic50": (
        "firkode",
        "FIRKODE s=3 periodic msbp=50",
        ["--stages", "3", "--reuse-policy", "periodic", "--msbp", "50"],
    ),
    # ARKODE with the settings of examples/arkode/C_serial/ark_robertson.c
    "ark4tuned": (
        "arkode",
        "ARKODE DIRK4 [ark_robertson settings]",
        [
            "--order",
            "4",
            "--h0",
            "1e-8",
            "--maxnef",
            "20",
            "--maxcor",
            "8",
            "--nlscoef",
            "1e-7",
            "--ark-predictor",
            "1",
        ],
    ),
}


def rtol_grid(lo_exp: int, hi_exp: int, half: bool) -> list[float]:
    """rtol = 10^-e for e from lo_exp to hi_exp (half-decades if half)."""
    step = 0.5 if half else 1.0
    out = []
    e = float(lo_exp)
    while e <= hi_exp + 1e-9:
        out.append(10.0**-e)
        e += step
    return out


def make_cell(
    family, problem, cfg, rtol, params=None, args=None, maxsteps=None, ref_rtol=None, **kw
) -> Cell:
    integ, label, cargs = CONFIGS[cfg]
    atol = ATOL_FACTOR[problem] * rtol if rtol is not None else None
    steps = ["--maxsteps", str(maxsteps)] if maxsteps else MAXSTEPS
    return Cell(
        family,
        problem,
        integ,
        label,
        rtol,
        atol,
        dict(params or {}),
        list(cargs) + list(args or []) + steps,
        ref_rtol=ref_rtol,
        **kw,
    )


# ---------------------------------------------------------------------------
# Families
# ---------------------------------------------------------------------------


def cells_b1() -> list[Cell]:
    cells = []
    rtols = rtol_grid(2, 10, half=True)
    for p in ["rober", "vdpol", "orego", "hires"]:
        for cfg in ["firk3", "cvode", "ark3", "ark4", "ark5", "ida"]:
            for r in rtols:
                cells.append(make_cell("B1", p, cfg, r, cap_s=120.0))
        if p == "rober":
            for r in rtols:
                cells.append(make_cell("B1", p, "ark4tuned", r, cap_s=120.0))
    return cells


def cells_b2(include_large: bool) -> list[Cell]:
    cells = []
    # diurnal kinetics, SPGMR + block-diagonal preconditioner
    for g in [10, 20, 40, 80]:
        for cfg in ["firk3", "cvode", "ark3", "ark4", "ida"]:
            for r in [1e-4, 1e-6, 1e-8]:
                cells.append(make_cell("B2", "diurnal", cfg, r, {"mx": g, "my": g}, cap_s=900.0))
    # 2-D heat, band Jacobian (direct block solver)
    for cfg in ["firk3", "cvode", "ark4", "ida"]:
        for r in [1e-4, 1e-6, 1e-8]:
            cells.append(make_cell("B2", "heat2d", cfg, r, {"mx": 100, "my": 100}, cap_s=900.0))
    # 2-D heat, matrix-free PCG with Jacobi or inner-PCG block preconditioning
    sizes = [100, 316] + ([1000] if include_large else [])
    for n in sizes:
        prm = {"mx": n, "my": n}
        cells.append(
            make_cell(
                "B2",
                "heat2d",
                "firk3",
                1e-6,
                prm,
                ["--ls", "pcg", "--psolve", "jacobi"],
                cap_s=1800.0,
            )
        )
        cells.append(
            make_cell(
                "B2",
                "heat2d",
                "firk3",
                1e-6,
                prm,
                ["--ls", "pcg", "--psolve", "innerpcg"],
                cap_s=1800.0,
            )
        )
        cells.append(make_cell("B2", "heat2d", "cvode", 1e-6, prm, ["--ls", "pcg"], cap_s=1800.0))
        cells.append(make_cell("B2", "heat2d", "ark4", 1e-6, prm, ["--ls", "pcg"], cap_s=1800.0))
    # advection-diffusion, band, increasing cell Peclet number
    for v in [0.5, 50.0, 500.0]:
        prm = {"mx": 200, "my": 100, "velocity": v}
        for cfg in ["firk3", "cvode", "ark4", "ida"]:
            cells.append(make_cell("B2", "heat2d", cfg, 1e-6, prm, cap_s=900.0))
    return cells


B3_PROBLEMS = [
    ("rober", {}),
    ("vdpol", {}),
    ("orego", {}),
    ("heat2d", {"mx": 100, "my": 100}),
    ("heat2d", {"mx": 200, "my": 100, "velocity": 500.0}),
    ("diurnal", {"mx": 20, "my": 20}),
]


def cells_b3() -> list[Cell]:
    cells = []
    for p, prm in B3_PROBLEMS:
        for r in [1e-6, 1e-8]:
            cells.append(make_cell("B3", p, "firk3", r, prm, cap_s=600.0, reps="single"))
    return cells


B4_PROBLEMS = [
    # problem, params, h0 (fixed-step), h0 (one-step estimate), tf
    ("pr", {"n": 4, "lambda": -1.0}, 0.5, 0.8, 1.0),
    ("pr", {"n": 4, "lambda": -1.0e4}, 0.5, 0.8, 1.0),
    ("pr", {"n": 4, "lambda": -1.0e6}, 0.5, 0.8, 1.0),
    ("kpr", {}, 0.1, 0.4, 0.5),
]

# Fixed-step runs: tight Newton as in firk_test_order. One-step runs: loose
# tolerances so that the step is accepted, Newton converged to ~1e-12 so that
# neither the estimate nor the one-step error is polluted by the iteration.
FIXED_ARGS = ["--rtol", "1e-12", "--atol", "1e-12", "--nlscoef", "1e-3", "--maxcor", "20"]
ONESTEP_ARGS = ["--rtol", "1", "--atol", "1", "--nlscoef", "1e-12", "--maxcor", "30"]


def cells_b4() -> list[Cell]:
    cells = []
    for p, prm, h0, _, tf in B4_PROBLEMS:
        for cfg in ["firk3", "ark2", "ark3", "ark4", "ark5"]:
            integ, label, cargs = CONFIGS[cfg]
            extra = ["--stage-epslin", "1e-3"] if integ == "firkode" else []
            for k in range(7):
                h = h0 / 2**k
                cells.append(
                    Cell(
                        "B4_fixed",
                        p,
                        integ,
                        label,
                        None,
                        None,
                        dict(prm),
                        list(cargs)
                        + ["--fixed-h", f"{h:.17g}", "--tf", f"{tf:g}"]
                        + FIXED_ARGS
                        + extra
                        + MAXSTEPS,
                        cap_s=60.0,
                        reps="single",
                        timed=False,
                        needs_ref=False,
                    )
                )
    # error-estimate order: one adaptive step of prescribed size
    for p, prm, _, h0, tf in [B4_PROBLEMS[0], B4_PROBLEMS[1], B4_PROBLEMS[3]]:
        for cfg in ["firk3", "ark3", "ark4", "ark5"]:
            integ, label, cargs = CONFIGS[cfg]
            extra = ["--stage-epslin", "1"] if integ == "firkode" else []
            for k in range(8):
                h = h0 / 2**k
                cells.append(
                    Cell(
                        "B4_estimate",
                        p,
                        integ,
                        label,
                        None,
                        None,
                        dict(prm),
                        list(cargs)
                        + ["--one-step", "--h", f"{h:.17g}", "--tf", f"{tf:g}"]
                        + ONESTEP_ARGS
                        + extra,
                        cap_s=60.0,
                        reps="single",
                        timed=False,
                        needs_ref=False,
                    )
                )
    return cells


def cells_b5() -> list[Cell]:
    cells = []
    # a step cap of 200 000 (instead of 2 000 000) makes the FIRKODE cells that
    # stall on the fine mesh fail within half a minute instead of timing out
    for n in [101, 401, 1601]:
        for cfg in ["firk3", "ark3", "ark4", "ark5", "ida"]:
            for r in rtol_grid(4, 10, half=False):
                cells.append(
                    make_cell("B5", "heat1dfem", cfg, r, {"n": n}, cap_s=300.0, maxsteps=200000)
                )
    return cells


B6_VARIANTS = [
    ("default", []),
    ("trivial-predictor", ["--predictor", "trivial"]),
    ("no-refilter", ["--no-refilter"]),
    ("no-newton-safety", ["--no-newton-safety"]),
    ("impgus-k1=k2=1", ["--impgus-k", "1", "1"]),
    ("controller-I", ["--controller", "i"]),
    ("controller-PI", ["--controller", "pi"]),
    ("controller-PID", ["--controller", "pid"]),
    ("controller-H0321", ["--controller", "h0321"]),
    ("msbp-1", ["--msbp", "1"]),
    ("msbp-50", ["--msbp", "50"]),
    ("stage-epslin-0.01", ["--stage-epslin", "0.01"]),
    ("stage-epslin-0.2", ["--stage-epslin", "0.2"]),
    ("stage-maxl-18", ["--stage-maxl", "18"]),
]


def cells_b6() -> list[Cell]:
    cells = []
    for p, prm in [("rober", {}), ("vdpol", {}), ("diurnal", {"mx": 20, "my": 20})]:
        for name, args in B6_VARIANTS:
            for r in [1e-6, 1e-8]:
                c = make_cell("B6", p, "firk3", r, prm, args, cap_s=600.0)
                c.label = f"FIRKODE s=3 [{name}]"
                cells.append(c)
    return cells


def cells_b7() -> list[Cell]:
    cells = []
    for cfg in ["firk3", "cvode", "ark4", "ida"]:
        for r in rtol_grid(4, 10, half=False):
            cells.append(
                make_cell("B7_roots", "rober", cfg, r, {}, ["--roots"], cap_s=300.0, reps="single")
            )
    # tight reference for the root times
    ref = make_cell(
        "B7_roots", "rober", "firk3", 1e-13, {}, ["--roots"], cap_s=1800.0, reps="single"
    )
    ref.label = "reference FIRKODE s=3 rtol=1e-13"
    cells.append(ref)
    for cfg in ["firk3", "cvode", "ark4", "ark5"]:
        for r in [1e-6, 1e-8]:
            cells.append(
                make_cell(
                    "B7_dense",
                    "vdpol",
                    cfg,
                    r,
                    {},
                    ["--dense-k", "1000"],
                    cap_s=600.0,
                    reps="single",
                    timed=False,
                    needs_ref=False,
                )
            )
    return cells


# ---------------------------------------------------------------------------
# Families B8-B13: where FIRKODE might beat CVODE BDF (BENCHMARK_BDF.md)
# ---------------------------------------------------------------------------


def order_twin(cell: Cell) -> Cell:
    """The untimed --order-hist twin of a CVODE or IDA cell."""
    twin = dataclasses.replace(cell, args=list(cell.args) + ["--order-hist"])
    twin.label = cell.label + " [order]"
    twin.timed = False
    twin.reps = "single"
    return twin


def profile_twin(cell: Cell) -> Cell:
    """The untimed --profile twin of a cell (time split of one run)."""
    twin = dataclasses.replace(cell, args=list(cell.args) + ["--profile"])
    twin.label = cell.label + " [profile]"
    twin.timed = False
    twin.reps = "single"
    return twin


B8_THETAS = [0.0, 30.0, 45.0, 50.0, 55.0, 60.0, 70.0, 75.0, 80.0, 85.0, 88.0, 89.5]
B8_RANGES = [(1.0, 1.0e3), (1.0, 1.0e6)]
# the maximum-order variants are the A(alpha)-stable members of the family:
# q <= 4 (73.35 degrees), q <= 3 (86.03 degrees) and q <= 2 (A-stable); IDA
# runs every CVODE variant but STALD, which it does not have
B8_CONFIGS = [
    "firk3",
    "cvode",
    "cvode_stald",
    "cvode_q4",
    "cvode_q3",
    "cvode_q2",
    "ida",
    "ida_q4",
    "ida_q3",
    "ida_q2",
]
B10_CONFIGS = ["firk3", "cvode", "cvode_stald", "cvode_q2", "ida", "ida_q2"]
# the BDF configurations with an untimed order-histogram twin
ORDER_TWIN_CONFIGS = ["cvode", "ida"]


def cells_b8() -> list[Cell]:
    """SECTOR: the spectral angle sweep (candidate 1)."""
    cells = []
    for theta in B8_THETAS:
        for rmin, rmax in B8_RANGES:
            prm = {"m": 100, "theta": theta, "rmin": rmin, "rmax": rmax}
            for cfg in B8_CONFIGS:
                for r in rtol_grid(4, 10, half=True):
                    cells.append(make_cell("B8", "sector", cfg, r, prm, cap_s=300.0))
    # BDF order histograms along the angle (untimed twins)
    for theta in [0.0, 60.0, 75.0, 80.0, 85.0, 89.5]:
        prm = {"m": 100, "theta": theta, "rmin": 1.0, "rmax": 1.0e3}
        for cfg in ["cvode", "cvode_stald", "ida"]:
            for r in [1e-6, 1e-8]:
                cells.append(order_twin(make_cell("B8", "sector", cfg, r, prm, cap_s=300.0)))
    # larger systems: per-step overhead amortized
    for theta in [0.0, 80.0, 89.5]:
        prm = {"m": 1000, "theta": theta, "rmin": 1.0, "rmax": 1.0e3}
        for cfg in ["firk3", "cvode", "cvode_stald", "ida"]:
            for r in [1e-6, 1e-8]:
                cells.append(make_cell("B8", "sector", cfg, r, prm, cap_s=600.0))
    # coupled (dense) spot check and the perturbed initial state
    for extra in [{"coupled": 1, "m": 50}, {"delta": 1e-2, "m": 100}]:
        prm = dict(theta=80.0, rmin=1.0, rmax=1.0e3, **extra)
        for cfg in B8_CONFIGS:
            for r in [1e-6, 1e-8]:
                cells.append(make_cell("B8", "sector", cfg, r, prm, cap_s=300.0))
    return cells


B9_SIZES = [100, 400, 1600]


def cells_b9() -> list[Cell]:
    """Damped finite-element wave problem (candidate 2)."""
    cells = []
    rtols = rtol_grid(3, 9, half=False)
    for n in B9_SIZES:
        lumped = [("mass", "lumped"), ("kv", "lumped")]
        consistent = [("mass", "consistent"), ("kv", "consistent"), ("rayleigh", "consistent")]
        for damping, mass in lumped:
            prm = {"n": n, "damping": damping, "mass": mass}
            for cfg in ["firk3", "cvode", "cvode_stald", "cvode_q2", "ida", "ida_q2"]:
                for r in rtols:
                    cells.append(make_cell("B9", "wave1dfem", cfg, r, prm, cap_s=600.0))
            if damping == "mass":
                for cfg in ORDER_TWIN_CONFIGS:
                    for r in [1e-6, 1e-8]:
                        cells.append(order_twin(make_cell("B9", "wave1dfem", cfg, r, prm)))
        for damping, mass in consistent:
            prm = {"n": n, "damping": damping, "mass": mass}
            for cfg in ["firk3", "ida", "ida_q2"]:
                for r in rtols:
                    # FIRKODE's mass-matrix path is slow here (known issue): a
                    # smaller step cap fails such cells in minutes
                    cells.append(
                        make_cell("B9", "wave1dfem", cfg, r, prm, cap_s=600.0, maxsteps=400000)
                    )
            if damping == "mass":
                for r in [1e-6, 1e-8]:
                    cells.append(order_twin(make_cell("B9", "wave1dfem", "ida", r, prm)))
    return cells


B10_PROBLEMS = ["beam", "plate", "ringmod"]


def cells_b10() -> list[Cell]:
    """Published stiff-oscillatory problems (candidate 3), with the spectrum
    dumps of gate B10.0 and the BDF order histograms."""
    cells = []
    for p in B10_PROBLEMS:
        for cfg in B10_CONFIGS:
            for r in rtol_grid(3, 10, half=True):
                cells.append(make_cell("B10", p, cfg, r, cap_s=900.0))
        for cfg in ORDER_TWIN_CONFIGS:
            for r in [1e-6, 1e-8]:
                cells.append(order_twin(make_cell("B10", p, cfg, r, cap_s=900.0)))
        # gate B10.0: the Jacobian at ten times along a tight trajectory
        tf = {"beam": 5.0, "plate": 7.0, "ringmod": 1e-3}[p]
        times = ",".join(f"{tf * k / 10:.17g}" for k in range(1, 11))
        gate = make_cell(
            "B10_spectrum",
            p,
            "cvode",
            1e-8,
            {},
            ["--dump-jac", times, "--dump-dir", f"jac/{p}"],
            cap_s=900.0,
            reps="single",
            timed=False,
        )
        gate.label = "CVODE BDF [spectrum]"
        cells.append(gate)
    return cells


def cells_b11() -> list[Cell]:
    """Restart-heavy workloads (candidate 4)."""
    cells = []
    rtols = rtol_grid(4, 8, half=False)
    # Lie splitting of the pollution problem with an exact relaxation step
    for K in [1, 10, 100, 1000]:
        prm = {"tau": 10.0, "restart-every": K}
        for cfg in ["firk3", "cvode", "ida", "firk3_h0", "cvode_h0", "ida_h0"]:
            for r in rtols:
                cells.append(make_cell("B11", "pollu", cfg, r, prm, cap_s=600.0))
        if K >= 100:
            for cfg in ORDER_TWIN_CONFIGS:
                cells.append(order_twin(make_cell("B11", "pollu", cfg, 1e-6, prm, cap_s=600.0)))
    # batches of independent cells amortize the per-step overhead
    for nb in [64, 1024]:
        for K in [1, 10, 100]:
            prm = {"nb": nb, "tau": 10.0, "restart-every": K}
            for cfg in ["firk3", "cvode", "ida"]:
                for r in rtols:
                    cells.append(make_cell("B11", "pollu", cfg, r, prm, cap_s=1800.0))
    # the switched RC ladder (square-wave forcing of SECTOR, theta = 0)
    for K in [10, 100, 1000]:
        prm = {
            "forcing": "square",
            "theta": 0.0,
            "m": 25,
            "rmin": 1.0,
            "rmax": 1.0e6,
            "period": 10.0 / K,
            "restart-every": K,
        }
        for cfg in ["firk3", "cvode", "ida"]:
            for r in rtols:
                cells.append(make_cell("B11", "sector", cfg, r, prm, cap_s=900.0))
        if K >= 100:
            for cfg in ORDER_TWIN_CONFIGS:
                cells.append(order_twin(make_cell("B11", "sector", cfg, 1e-6, prm, cap_s=900.0)))
    return cells


def cells_b12(include_large: bool) -> list[Cell]:
    """Factorization-dominated 3-D reaction-diffusion with KLU (candidate 5)."""
    cells = []
    sizes = [8, 12, 16, 24] + ([32] if include_large else [])
    # IDA has no setup-frequency option and runs its default only
    cfgs = ["firk3", "firk3_periodic", "firk3_periodic50", "cvode", "cvode_msbp50", "ida"]
    # the reference tolerance loosens with the size (one factorization grows
    # like the fill of a 3-D grid)
    ref_rtol = {8: 1e-10, 12: 1e-10, 16: 1e-10, 24: 1e-10, 32: 1e-9}
    caps = {8: 600.0, 12: 1200.0, 16: 3600.0, 24: 7200.0, 32: 14400.0}
    for n in sizes:
        prm = {"n": n}
        for cfg in cfgs:
            for r in rtol_grid(4, 8, half=False):
                cells.append(
                    make_cell("B12", "bruss3d", cfg, r, prm, cap_s=caps[n], ref_rtol=ref_rtol[n])
                )
        # the time split of one run per configuration
        for cfg in cfgs:
            cells.append(
                profile_twin(
                    make_cell(
                        "B12", "bruss3d", cfg, 1e-6, prm, cap_s=caps[n], ref_rtol=ref_rtol[n]
                    )
                )
            )
    # control: the 2-D heat equation with a band solver, where the
    # factorization costs only about the bandwidth more than a solve
    for mx in [100, 200]:
        prm = {"mx": mx, "my": mx}
        for cfg in ["firk3", "firk3_periodic50", "cvode", "cvode_msbp50", "ida"]:
            for r in [1e-4, 1e-6, 1e-8]:
                cells.append(make_cell("B12", "heat2d", cfg, r, prm, cap_s=1800.0))
            cells.append(profile_twin(make_cell("B12", "heat2d", cfg, 1e-6, prm, cap_s=1800.0)))
    return cells


def cells_b13() -> list[Cell]:
    """Tight tolerances against extended-precision references (candidate 6).
    Needs --ref-bin (the x86_64 extended build) for the ODE problems."""
    cells = []
    for p in ["rober", "vdpol", "orego", "hires", "sector"]:
        variants = [{}]
        if p == "sector":
            variants = [
                {"m": 100, "theta": 0.0, "rmin": 1.0, "rmax": 1.0e3},
                {"m": 100, "theta": 80.0, "rmin": 1.0, "rmax": 1.0e3},
            ]
        for prm in variants:
            cfgs = ["firk3", "cvode", "ida"]
            if prm.get("theta") == 80.0:
                cfgs.append("cvode_stald")
            for cfg in cfgs:
                for r in rtol_grid(8, 13, half=True):
                    # SECTOR has an exact solution; the others get an 80-bit reference
                    rr = None if p == "sector" else 1e-16
                    cells.append(make_cell("B13", p, cfg, r, prm, cap_s=900.0, ref_rtol=rr))
    for cfg in ["firk3", "cvode", "ida"]:
        for r in rtol_grid(8, 11, half=True):
            cells.append(
                make_cell("B13", "bruss3d", cfg, r, {"n": 16}, cap_s=3600.0, ref_rtol=1e-10)
            )
    return cells


FAMILY_BUILDERS = {
    "B1": lambda large: cells_b1(),
    "B2": lambda large: cells_b2(large),
    "B3": lambda large: cells_b3(),
    "B4": lambda large: cells_b4(),
    "B5": lambda large: cells_b5(),
    "B6": lambda large: cells_b6(),
    "B7": lambda large: cells_b7(),
    "B8": lambda large: cells_b8(),
    "B9": lambda large: cells_b9(),
    "B10": lambda large: cells_b10(),
    "B11": lambda large: cells_b11(),
    "B12": lambda large: cells_b12(large),
    "B13": lambda large: cells_b13(),
}

FAMILY_ORDER = ["B4", "B1", "B5", "B7", "B6", "B3", "B2", "B8", "B9", "B10", "B11", "B12", "B13"]

# ---------------------------------------------------------------------------
# Validation replicas of repository examples and tests
# ---------------------------------------------------------------------------

EX = "examples/firkode/C_serial"
UT = "test/unit_tests/firkode/CXX_serial"


def validation_cells() -> list[dict]:
    """Replicas of repository examples with the counters of their .out files.

    tol: "exact" for direct linear solvers; a relative tolerance where
    floating-point differences between builds are expected (Krylov paths;
    answer files generated on another platform).
    """
    out = []

    def add(name, cell, expected, tol="exact", note=""):
        out.append(dict(name=name, cell=cell, expected=expected, tol=tol, note=note))

    common = ["--grid", "10", "--no-tstop", "--maxsteps", "100000"]
    add(
        "firk_robertson",
        Cell(
            "validation",
            "rober",
            "firkode",
            "FIRKODE s=3",
            1e-4,
            1e-11,
            {},
            ["--stages", "3", "--grid", "100", "--no-tstop", "--maxsteps", "100000"],
            300,
            "single",
            False,
            False,
        ),
        dict(
            nst=170,
            nst_att=170,
            nfe=1241,
            nni=357,
            nncf=0,
            nsetups=169,
            nje=156,
            netf=0,
            nbs=4505,
            nsli=1445,
        ),
        note=f"{EX}/firk_robertson.out",
    )
    add(
        "firk_analytic",
        Cell(
            "validation",
            "analytic",
            "firkode",
            "FIRKODE s=3",
            1e-6,
            1e-10,
            {"lambda": -100.0},
            ["--stages", "3"] + common,
            300,
            "single",
            False,
            False,
        ),
        dict(nst=38, nst_att=38, nfe=175, nsetups=29, nje=1, nni=45, netf=0, nsli=93),
        note=f"{EX}/firk_analytic.out",
    )
    add(
        "firk_diurnal_kry",
        Cell(
            "validation",
            "diurnal",
            "firkode",
            "FIRKODE s=3",
            1e-5,
            1e-3,
            {"mx": 10, "my": 10},
            ["--stages", "3", "--grid", "12", "--no-tstop", "--maxsteps", "10000"],
            300,
            "single",
            False,
            False,
        ),
        dict(
            nst=199,
            nfe=1197,
            nni=331,
            nsli=1905,
            nli=363,
            njtv=6078,
            nsetups=117,
            netf=10,
            npe=117,
            nps=6290,
        ),
        tol=0.05,
        note=f"{EX}/firk_diurnal_kry.out (Krylov path)",
    )
    add(
        "firk_heat1D_fem_mass",
        Cell(
            "validation",
            "heat1dfem",
            "firkode",
            "FIRKODE s=3",
            1e-6,
            1e-10,
            {"n": 101},
            ["--stages", "3", "--grid", "10", "--no-tstop", "--maxsteps", "10000"],
            300,
            "single",
            False,
            False,
        ),
        dict(
            nst=11,
            nst_att=11,
            nfe=52,
            nsetups=3,
            nje=1,
            nni=13,
            netf=0,
            nmass_setups=1,
            nmass_solves=42,
            nmass_mult=182,
        ),
        note=f"{EX}/firk_heat1D_fem_mass.out",
    )
    roots_common = [
        "--tf",
        "4e10",
        "--grid",
        "12",
        "--grid-geom",
        "0.4",
        "10",
        "--no-tstop",
        "--roots",
        "--maxsteps",
        "100000",
    ]
    add(
        "firk_robertson_root",
        Cell(
            "validation",
            "rober",
            "firkode",
            "FIRKODE s=3",
            1e-4,
            None,
            {},
            ["--stages", "3", "--atol-vec", "1e-8,1e-11,1e-8"] + roots_common,
            300,
            "single",
            False,
            False,
        ),
        dict(
            nst=140,
            nst_att=141,
            nfe=1063,
            nsetups=141,
            nje=126,
            nni=307,
            nncf=0,
            netf=1,
            nroots=2,
            ngeval=182,
        ),
        note=f"{EX}/firk_robertson_root.out reports 180 g evaluations; the example advances the output time after "
        "a root return while the harness continues to the same output time, which costs one g evaluation per root",
    )
    order_common = FIXED_ARGS + ["--stage-epslin", "1e-3", "--maxsteps", "100000"]
    for k, (h, err) in enumerate(
        [(0.02, 9.3065e-07), (0.01, 2.9403e-08), (0.005, 9.2170e-10), (0.0025, 2.8888e-11)]
    ):
        add(
            f"firk_test_order_kpr_3_h{k}",
            Cell(
                "validation",
                "kpr",
                "firkode",
                "FIRKODE s=3",
                None,
                None,
                {},
                ["--stages", "3", "--fixed-h", f"{h:g}", "--tf", "0.5"] + order_common,
                60,
                "single",
                False,
                False,
            ),
            dict(nst=int(round(0.5 / h)), err_abs_max=err),
            tol=1e-4,
            note=f"{UT}/firk_test_order_kpr_3.out",
        )
    # the answer file lists the steps taken (11 at h = 0.1: one step was retried
    # after a Newton failure with a stale Jacobian)
    for k, (h, nst, err) in enumerate(
        [
            (0.2, 5, 3.3880e-09),
            (0.1, 11, 1.8742e-10),
            (0.05, 20, 1.0999e-11),
            (0.025, 40, 6.6480e-13),
        ]
    ):
        add(
            f"firk_test_order_pr_3_-1e4_h{k}",
            Cell(
                "validation",
                "pr",
                "firkode",
                "FIRKODE s=3",
                None,
                None,
                {"n": 4, "lambda": -1.0e4},
                ["--stages", "3", "--fixed-h", f"{h:g}", "--tf", "1"] + order_common,
                60,
                "single",
                False,
                False,
            ),
            dict(nst=nst, err_abs_max=err),
            tol=1e-4,
            note=f"{UT}/firk_test_order_pr_3_-1e4.out",
        )
    add(
        "firk_test_vdp_1e-6",
        Cell(
            "validation",
            "vdpol",
            "firkode",
            "FIRKODE s=3",
            1e-6,
            1e-8,
            {},
            ["--stages", "3", "--no-tstop", "--maxsteps", "1000000"],
            300,
            "single",
            False,
            False,
        ),
        dict(nst=922),
        note=f"{UT}/firk_test_vdp.out",
    )
    add(
        "firk_test_vdp_1e-10",
        Cell(
            "validation",
            "vdpol",
            "firkode",
            "FIRKODE s=3",
            1e-10,
            1e-12,
            {},
            ["--stages", "3", "--no-tstop", "--maxsteps", "1000000"],
            300,
            "single",
            False,
            False,
        ),
        dict(nst=9345),
        note=f"{UT}/firk_test_vdp.out",
    )
    add(
        "firk_test_vdp_1e-4",
        Cell(
            "validation",
            "vdpol",
            "firkode",
            "FIRKODE s=3",
            1e-4,
            1e-4,
            {},
            ["--stages", "3", "--no-tstop", "--maxsteps", "1000000"],
            300,
            "single",
            False,
            False,
        ),
        dict(nst=295),
        note=f"{UT}/firk_test_vdp.out (the loose-tolerance run of benchmark issue 2)",
    )
    add(
        "ark_analytic",
        Cell(
            "validation",
            "analytic",
            "arkode",
            "ARKODE DIRK4",
            1e-5,
            1e-10,
            {"lambda": -100.0},
            ["--arkode-linear"] + common,
            300,
            "single",
            False,
            False,
        ),
        dict(nst=34, nst_att=35, nfe=353, nsetups=35, nje=2, nni=175, netf=1),
        note="examples/arkode/C_serial/ark_analytic.out (ARKodeSetLinear as in the example)",
    )
    add(
        "cv_kpr",
        Cell(
            "validation",
            "kpr",
            "cvode",
            "CVODE BDF",
            1e-6,
            1e-10,
            {},
            ["--tf", "10"] + common,
            300,
            "single",
            False,
            False,
        ),
        dict(nst=1577, nfe=2197, nni=2194, nsetups=272, nje=29, netf=147),
        tol=0.10,
        note="examples/cvode/CXX_serial/cv_kpr.out was generated on another platform; the example compiled here "
        "against the same Release library gives nst=1561 nfe=2117 nni=2114 nsetups=253 nje=29 netf=133, "
        "which the harness reproduces exactly",
    )
    add(
        "idaRoberts_dns",
        Cell(
            "validation",
            "rober",
            "ida",
            "IDA BDF",
            1e-4,
            None,
            {},
            ["--atol-vec", "1e-8,1e-6,1e-6"] + roots_common,
            300,
            "single",
            False,
            False,
        ),
        dict(nst=362, nfe=537, nni=537, nsetups=60, nje=60, netf=15, nncf=5, ngeval=404, nroots=2),
        note="examples/ida/serial/idaRoberts_dns.out (the example poses the DAE with the conservation constraint; "
        "the harness poses the ODE in residual form and reproduces every counter)",
    )
    add(
        "cvRoberts_dns",
        Cell(
            "validation",
            "rober",
            "cvode",
            "CVODE BDF",
            1e-4,
            None,
            {},
            ["--atol-vec", "1e-8,1e-14,1e-6"] + roots_common,
            300,
            "single",
            False,
            False,
        ),
        dict(
            nst=542, nfe=754, nni=751, nsetups=107, nje=11, netf=22, nncf=3, ngeval=570, nroots=2
        ),
        tol=0.10,
        note="examples/cvode/serial/cvRoberts_dns.out computes the second component as -yd1 - yd3, a different "
        "floating-point expression; informational only",
    )
    return out


def compare(expected: dict, row: dict, tol) -> tuple[bool, list[str]]:
    msgs = []
    ok = True
    for k, v in expected.items():
        got = row.get(k, "NA")
        if isinstance(v, float):
            g = to_float(got)
            rel = abs(g - v) / abs(v) if v != 0 else abs(g)
            passed = rel <= (tol if tol != "exact" else 1e-12)
        else:
            g = to_int(got)
            if g is None:
                passed = False
            elif tol == "exact":
                passed = g == v
            else:
                passed = abs(g - v) <= max(1, tol * abs(v))
        if not passed:
            ok = False
            msgs.append(f"{k}: expected {v}, got {got}")
    return ok, msgs


# ---------------------------------------------------------------------------
# Running cells
# ---------------------------------------------------------------------------


def choose_reps(cell: Cell, pilot_time: float, reps_opt: str) -> int:
    if reps_opt != "auto":
        return int(reps_opt)
    if cell.reps == "single" or not cell.timed:
        return 1
    if isinstance(cell.reps, int):
        return cell.reps
    if pilot_time < TIMING_FLOOR_S:
        return 1
    if pilot_time >= 10.0:
        return 1
    # about one second of repeated work, at least 5 and at most 200 repetitions
    return int(min(200, max(5, math.ceil(1.0 / pilot_time))))


def run_cell(
    cell: Cell,
    binary: str,
    env: dict,
    logs: Path,
    ref_path,
    reps_opt: str,
    budget: float | None = None,
) -> dict:
    # the key of the cell as planned, before --dump-dir is made absolute
    key = cell.key()
    extra = []
    if ref_path is not None:
        extra += ["--ref", str(ref_path)]
    # a relative --dump-dir lives in the results directory
    if "--dump-dir" in cell.args:
        i = cell.args.index("--dump-dir") + 1
        if i < len(cell.args) and not Path(cell.args[i]).is_absolute():
            d = logs.parents[1] / cell.args[i]
            d.mkdir(parents=True, exist_ok=True)
            cell = dataclasses.replace(cell, args=cell.args[:i] + [str(d)] + cell.args[i + 1 :])
    cmd = cell.command(binary, extra)
    # --fast: the budget bounds the pilot run and the repetitions of a cell
    cap = cell.cap_s if budget is None else min(cell.cap_s, budget)
    row, flag, wall = run_binary(cmd, cap, env, logs / f"{sanitize(key)}.err")
    if flag == "TIMEOUT" and budget is not None and cap < cell.cap_s:
        flag = "FAST_TIMEOUT"
    reps = 1
    if flag == "ok" and row is not None:
        reps = choose_reps(cell, to_float(row.get("time_min")), reps_opt)
        if budget is not None and to_float(row.get("time_min")) > 0:
            reps = max(1, min(reps, int(budget / to_float(row.get("time_min")))))
        if reps > 1:
            cmd2 = cell.command(binary, extra + ["--reps", str(reps)])
            row2, flag2, wall2 = run_binary(
                cmd2, cap * reps + 60, env, logs / f"{sanitize(key)}.err"
            )
            if flag2 == "ok" and row2 is not None:
                row, flag, wall = row2, flag2, wall + wall2
    out = {
        "family": cell.family,
        "cell_key": key,
        "label": cell.label,
        "params": cell.params_str(),
        "driver_flag": flag,
        "timing_reliable": int(
            cell.timed
            and flag == "ok"
            and row is not None
            and to_float(row.get("time_min")) >= TIMING_FLOOR_S
        ),
        "wall_s": f"{wall:.3f}",
        "timestamp": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
    }
    if row is not None:
        out.update(row)
    else:
        out.update(
            {
                "problem": cell.problem,
                "integrator": cell.integrator,
                "rtol": cell.rtol,
                "atol": cell.atol,
            }
        )
    return out


def sanitize(key: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.=-]+", "_", key)[:180]


def run_sentinel(
    binary: str, env: dict, results: Path, family: str, phase: str, suffix: str = ""
) -> None:
    cell = Cell(
        "sentinel",
        "rober",
        "firkode",
        "FIRKODE s=3",
        1e-6,
        1e-12,
        {},
        ["--stages", "3", "--reps", "20"] + MAXSTEPS,
        300,
        20,
    )
    row, flag, _ = run_binary(cell.command(binary), 300, env)
    table = ResultsTable(results / f"sentinel{suffix}.csv")
    table.append(
        {
            "family": family,
            "cell_key": f"sentinel|{family}|{phase}{suffix}|{dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
            "label": phase,
            "params": "",
            "driver_flag": flag,
            "timing_reliable": 1,
            "wall_s": "",
            "timestamp": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
            "time_min": row.get("time_min") if row else "NA",
            "time_median": row.get("time_median") if row else "NA",
            "nst": row.get("nst") if row else "NA",
        }
    )


def pending_cells(cells: list[Cell], args) -> list[Cell]:
    """The cells of a family still to run, in the order they run."""
    results = Path(args.results)
    # cells already present in the merged file or in any worker file
    done = set()
    for fam in sorted({c.family for c in cells}):
        for path in family_files(results, fam):
            for r in read_csv(path):
                # cells cut short by a --fast budget are redone by any run
                # without --fast
                if r.get("driver_flag") == "FAST_TIMEOUT" and args.fast is None:
                    continue
                if not args.retry_failed or r.get("driver_flag") == "ok":
                    done.add(r.get("cell_key", ""))
    if args.force:
        done = set()
    todo = [c for c in cells if c.key() not in done]
    # expensive cells first, so that the parallel tail is short (the cap is
    # the only cost estimate available before a cell has run)
    # expensive cells first for a short parallel tail; cheap cells first in a
    # --fast run, whose point is the quickly finishing part of the suite
    todo.sort(key=lambda c: c.cap_s if args.fast is not None else -c.cap_s)
    return todo


def run_family(name: str, cells: list[Cell], args, env: dict, progress=None) -> None:
    results = Path(args.results)
    parallel = args.nworkers > 1
    suffix = f".worker{args.worker}" if parallel else ""
    tag = f"[w{args.worker}] " if parallel else ""
    # sub-families (B4_fixed, B4_estimate, B7_roots, B7_dense) get their own
    # CSV; a worker appends to its own file, merged by --merge afterwards
    subfamilies = sorted({c.family for c in cells})
    tables = {fam: ResultsTable(results / f"{fam}{suffix}.csv") for fam in subfamilies}
    logs = results / "logs" / name
    refdir = results / "reference"
    todo = pending_cells(cells, args)
    print(f"{tag}== {name}: {len(cells)} cells, {len(todo)} to run", flush=True)
    if not todo:
        return
    if args.retry_failed:
        for t in tables.values():
            t.remove_keys({c.key() for c in todo})
    elif args.fast is None:
        # rows cut by an earlier --fast run are replaced, not duplicated (a
        # worker's file only; the merge keeps the latest row per key)
        redo = {c.key() for c in todo}
        for t in tables.values():
            stale = {r["cell_key"] for r in t.rows if r.get("driver_flag") == "FAST_TIMEOUT"}
            if stale & redo:
                t.remove_keys(stale & redo)
    run_sentinel(args.bin, env, results, name, "start", suffix)
    t_family = dt.datetime.now(dt.UTC)
    failed_refs: set = set()
    n_run = 0
    for i, cell in enumerate(todo, 1):
        if parallel and not claim_cell(results, cell.family, cell.key(), args.worker):
            continue
        n_run += 1
        ref_path = None
        ref_key = (cell.problem, cell.params_str())
        if cell.needs_ref and cell.problem in REF_RTOL and ref_key not in failed_refs:
            try:
                # the extended-precision binary only for references that need it
                want = cell.ref_rtol or REF_RTOL.get(cell.problem, 1.0)
                ref_binary = (
                    args.ref_bin if (args.ref_bin and want < REF_RTOL_DOUBLE_FLOOR) else args.bin
                )
                ref_path, _ = ensure_reference(
                    ref_binary,
                    cell.problem,
                    cell.params,
                    cell.rtol,
                    cell.atol,
                    refdir,
                    env,
                    cell.ref_rtol,
                    bool(args.ref_bin),
                    None if args.fast is None else 3.0 * args.fast,
                )
            except ReferenceTimeout as exc:
                print(f"  {exc}; its cells get no error metrics in this run", flush=True)
                failed_refs.add(ref_key)
            except RuntimeError as exc:
                # keep going without an error column for this problem size and
                # record the failure (reported in the benchmark document)
                print(f"  WARNING: {exc}; cells of {ref_key} get no error metrics", flush=True)
                failed_refs.add(ref_key)
        row = run_cell(cell, args.bin, env, logs, ref_path, args.reps, args.fast)
        if ref_key in failed_refs:
            row["driver_note"] = (
                "reference skipped (fast budget)" if args.fast is not None else "reference failed"
            )
        if args.fast is not None:
            row["fast_budget_s"] = f"{args.fast:g}"
        row["worker"] = args.worker if parallel else 0
        row["nworkers"] = args.nworkers
        tables[cell.family].append(row)
        print(
            f"{tag}  [{i}/{len(todo)}] {cell.problem} {cell.params_str()} {cell.label} rtol={cell.rtol} -> "
            f"{row['driver_flag']} nst={row.get('nst', 'NA')} err={row.get('err_mixed', 'NA')} "
            f"t={row.get('time_median', 'NA')}",
            flush=True,
        )
        if progress is not None:
            progress.record(row)
            print(f"{tag}  {progress.report()}", flush=True)
    run_sentinel(args.bin, env, results, name, "end", suffix)
    print(
        f"{tag}== {name} finished in {dt.datetime.now(dt.UTC) - t_family} ({n_run} cells run here)",
        flush=True,
    )


# ---------------------------------------------------------------------------
# Progress (--progress)
# ---------------------------------------------------------------------------

# how often the parent of --workers reads the worker CSV files, and the
# longest silence between two progress lines while no cell finishes
PROGRESS_POLL_S = 30.0
PROGRESS_QUIET_S = 60.0


def read_complete_rows(path: Path) -> list[dict]:
    """Rows of a CSV file another process may be writing; a trailing partial
    line is dropped (it is read complete at the next poll)."""
    try:
        text = path.read_text()
    except FileNotFoundError:
        return []
    end = text.rfind("\n")
    if end < 0:
        return []
    return list(csv.DictReader(io.StringIO(text[: end + 1], newline="")))


def row_ident(row: dict) -> tuple:
    return (row.get("family"), row.get("cell_key"), row.get("timestamp"), row.get("worker"))


def format_duration(seconds: float) -> str:
    seconds = max(0, round(seconds))
    return f"{seconds // 3600}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}"


class Progress:
    """Counts the finished cells of a sweep."""

    def __init__(self, plan: list[Cell], results: Path):
        self.t0 = dt.datetime.now(dt.UTC)
        # cell keys repeat across families (B1 and B8 share cells)
        self.pending = {(c.family, c.key()) for c in plan}
        self.total = len(self.pending)
        self.n_done = 0
        self.families = sorted({c.family for c in plan})
        # rows present before this run, not counted as finished by poll()
        self.seen: set = set()
        for fam in self.families:
            for path in family_files(results, fam):
                for r in read_csv(path):
                    self.seen.add(row_ident(r))

    def record(self, row: dict) -> None:
        """Count a row written in this run (ignored unless it is a pending cell)."""
        key = (row.get("family"), row.get("cell_key"))
        if key in self.pending:
            self.pending.remove(key)
            self.n_done += 1

    def poll(self, results: Path) -> int:
        """Count the new rows of the worker files; returns their number."""
        before = self.n_done
        for fam in self.families:
            for path in sorted(results.glob(f"{fam}.worker*.csv")):
                for r in read_complete_rows(path):
                    ident = row_ident(r)
                    if ident not in self.seen:
                        self.seen.add(ident)
                        self.record(r)
        return self.n_done - before

    def report(self, now: dt.datetime | None = None) -> str:
        now = now or dt.datetime.now(dt.UTC)
        elapsed = (now - self.t0).total_seconds()
        pct = 100.0 * self.n_done / self.total if self.total else 100.0
        return (
            f"progress: {self.n_done}/{self.total} cells ({pct:.1f}%), "
            f"elapsed {format_duration(elapsed)}"
        )

def spawn_workers(args, progress: Progress | None = None) -> int:
    """Run the sweep with --workers N worker processes claiming cells dynamically."""
    results = Path(args.results)
    logs = results / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    n = clear_claims(results)
    if n:
        print(f"removed {n} stale claims", flush=True)
    # the worker command line: this invocation without the parent-only flags
    argv = sys.argv[1:]
    worker_argv = []
    skip = 0
    for i, a in enumerate(argv):
        if skip:
            skip -= 1
            continue
        if a in ("--workers",):
            skip = 1
            continue
        if a.startswith("--workers="):
            continue
        if a in (
            "--pin",
            "--validate",
            "--validate-only",
            "--merge",
            "--clear-claims",
            "--progress",
        ):
            continue
        worker_argv.append(a)
    cores = None
    if args.pin:
        if shutil.which("taskset") and hasattr(os, "sched_getaffinity"):
            cores = sorted(os.sched_getaffinity(0))
        else:
            print("--pin: no taskset/sched_getaffinity on this platform; not pinning", flush=True)
    procs = []
    for k in range(args.workers):
        cmd = [sys.executable, str(Path(__file__).resolve())]
        cmd += worker_argv + ["--worker", str(k), "--nworkers", str(args.workers)]
        if cores:
            cmd = ["taskset", "-c", str(cores[k % len(cores)])] + cmd
        log = (logs / f"worker{k}.log").open("a")
        procs.append((k, subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT), log))
        print(f"started worker {k}: {' '.join(cmd)}", flush=True)
    rc = 0
    last_report = time.monotonic()
    for k, proc, log in procs:
        while True:
            try:
                code = proc.wait(timeout=PROGRESS_POLL_S if progress else None)
                break
            except subprocess.TimeoutExpired:
                if progress.poll(results) or time.monotonic() - last_report >= PROGRESS_QUIET_S:
                    print(progress.report(), flush=True)
                    last_report = time.monotonic()
        log.close()
        print(f"worker {k} finished with exit code {code}", flush=True)
        rc = rc or code
    if progress is not None:
        progress.poll(results)
        print(progress.report(), flush=True)
    merged = merge_results(results)
    for base, count in merged.items():
        print(f"merged {base}.csv: {count} rows", flush=True)
    return rc


def run_validation(args, env: dict) -> bool:
    results = Path(args.results)
    table = ResultsTable(results / "validation.csv")
    table.remove_keys(table.done_keys())
    logs = results / "logs" / "validation"
    all_ok = True
    print("== validation", flush=True)
    for item in validation_cells():
        row = run_cell(item["cell"], args.bin, env, logs, None, "1")
        ok, msgs = compare(item["expected"], row, item["tol"])
        all_ok &= ok
        row.update(
            {
                "family": "validation",
                "label": item["name"],
                "expected": ";".join(f"{k}={v}" for k, v in item["expected"].items()),
                "tolerance": str(item["tol"]),
                "passed": int(ok),
                "mismatch": "; ".join(msgs),
                "note": item["note"],
            }
        )
        table.append(row)
        print(
            f"  {'PASS' if ok else 'FAIL'} {item['name']}"
            + (f": {'; '.join(msgs)}" if msgs else ""),
            flush=True,
        )
    # self-consistency of FIRKODE direct-solver counters (adaptive mode):
    # block solves = s per stacked FGMRES iteration + one error estimate per
    # attempt whose Newton iteration converged + one per refilter
    for r in table.rows:
        if (
            r.get("integrator") == "firkode"
            and r.get("ls") in ("dense", "band")
            and r.get("driver_flag") == "ok"
        ):
            if r.get("fixed_h", "NA") != "NA":
                continue
            s, nsli, nbs, nst_att, nst, nsf = (
                to_int(r.get(k)) for k in ("stages", "nsli", "nbs", "nst_att", "nst", "nsf")
            )
            if None not in (s, nsli, nbs, nst_att, nst, nsf):
                n_refilter = nbs - s * nsli - (nst_att - nsf)
                if n_refilter < 0 or nst_att < nst:
                    all_ok = False
                    print(
                        f"  FAIL self-consistency for {r['label']}: n_refilter={n_refilter}",
                        flush=True,
                    )
    print(f"== validation {'passed' if all_ok else 'FAILED'}", flush=True)
    return all_ok


# ---------------------------------------------------------------------------
# k_lin per stacked solve from the INFO log (needs SUNDIALS_LOGGING_LEVEL >= 3)
# ---------------------------------------------------------------------------

LOG_RE = re.compile(
    r"\[end-stage-linear-solve\] status = (success small rhs|success|failed)(?:, retval = (-?\d+))?(?:, iters = (\d+))?"
)


def parse_stage_log(path: Path) -> dict:
    iters, n_small, n_failed = [], 0, 0
    with path.open() as fp:
        for line in fp:
            m = LOG_RE.search(line)
            if not m:
                continue
            if m.group(1) == "success small rhs":
                n_small += 1
                continue
            if m.group(1) == "failed":
                n_failed += 1
            if m.group(3) is not None:
                iters.append(int(m.group(3)))
    out = {"n_solves": len(iters) + n_small, "n_small_rhs": n_small, "n_failed": n_failed}
    if iters:
        iters.sort()
        out.update(
            {
                "klin_max": iters[-1],
                "klin_p95": iters[min(len(iters) - 1, int(0.95 * len(iters)))],
                "klin_mean_log": f"{statistics.fmean(iters):.4f}",
            }
        )
    return out


def run_b3_log(args, env: dict) -> None:
    results = Path(args.results)
    table = ResultsTable(results / "B3_klin_log.csv")
    done = table.done_keys()
    logs = results / "logs" / "B3_log"
    logs.mkdir(parents=True, exist_ok=True)
    cells = [c for c in cells_b3() if c.rtol == 1e-6]
    print(f"== B3 logging run: {len(cells)} cells", flush=True)
    for cell in cells:
        if cell.key() in done and not args.force:
            continue
        logfile = logs / f"{sanitize(cell.key())}.log"
        if logfile.exists():
            logfile.unlink()
        env2 = dict(env)
        env2["SUNLOGGER_INFO_FILENAME"] = str(logfile)
        row = run_cell(cell, args.log_bin, env2, logs, None, "1")
        row["timing_reliable"] = 0
        if logfile.exists():
            row.update(parse_stage_log(logfile))
            s = to_int(row.get("stages")) or 0
            maxl = max(5, 3 * s)
            with logfile.open() as fp:
                row["n_restarted"] = sum(
                    1
                    for line in fp
                    if (m := LOG_RE.search(line))
                    and m.group(3) is not None
                    and int(m.group(3)) > maxl
                )
        table.append(row)
        print(
            f"  {cell.problem} {cell.params_str()} s={row.get('stages')} -> {row['driver_flag']} "
            f"klin_max={row.get('klin_max', 'NA')} n_solves={row.get('n_solves', 'NA')}",
            flush=True,
        )


# ---------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bin", required=True, help="path to firkode_benchmark (Release build)")
    ap.add_argument(
        "--ref-bin",
        help="firkode_benchmark used to compute the reference solutions (e.g. the "
        "extended-precision build); the files it writes are read by --bin",
    )
    ap.add_argument("--results", required=True, help="results directory for this run")
    ap.add_argument("--family", nargs="*", default=[], help="families to run (B1..B13 or all)")
    ap.add_argument("--list", action="store_true", help="print the cell counts and exit")
    ap.add_argument(
        "--fast",
        nargs="?",
        const=10.0,
        type=float,
        metavar="SECONDS",
        default=None,
        help="partial run of the cells that finish fast: every cell (and its reference, with "
        "three times the budget) is cut off after SECONDS (default 10) and recorded as "
        "FAST_TIMEOUT, cheap cells run first; a later run without --fast redoes the cut cells",
    )
    ap.add_argument("--validate", action="store_true", help="run the example replicas first")
    ap.add_argument("--validate-only", action="store_true")
    ap.add_argument("--include-large", action="store_true", help="include the N=1e6 heat cells")
    ap.add_argument(
        "--log-bin", default=None, help="firkode_benchmark from a SUNDIALS_LOGGING_LEVEL>=3 build"
    )
    ap.add_argument("--retry-failed", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--reps", default="auto", help="auto or a fixed repetition count")
    ap.add_argument(
        "--build-dir", default=None, help="build directory for provenance (default: from --bin)"
    )
    ap.add_argument(
        "--workers",
        type=int,
        default=0,
        help="run the sweep with this many worker processes (one cell at a time each), "
        "then merge their CSV files",
    )
    ap.add_argument("--pin", action="store_true", help="pin each worker to a core (taskset)")
    ap.add_argument(
        "--progress",
        action="store_true",
        help="print the finished cells of the sweep and the elapsed time "
        f"(after every cell, or every {PROGRESS_POLL_S:g} s with --workers)",
    )
    ap.add_argument("--worker", type=int, default=None, help="(internal) this worker's index")
    ap.add_argument("--nworkers", type=int, default=1, help="(internal) number of workers")
    ap.add_argument(
        "--merge", action="store_true", help="merge the worker CSV files of --results and exit"
    )
    ap.add_argument(
        "--clear-claims",
        action="store_true",
        help="remove stale cell claims before running (after a worker died)",
    )
    ap.add_argument(
        "--machine-note",
        default="",
        help="hardware description recorded in provenance.json when sysctl is unavailable (sandbox)",
    )
    args = ap.parse_args()

    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("VECLIB_MAXIMUM_THREADS", "1")

    binary = Path(args.bin).resolve()
    args.bin = str(binary)
    if args.ref_bin:
        args.ref_bin = str(Path(args.ref_bin).resolve())
    if not binary.exists():
        print(f"binary not found: {binary}", file=sys.stderr)
        return 1
    results = Path(args.results)
    results.mkdir(parents=True, exist_ok=True)
    if args.merge:
        for base, count in merge_results(results).items():
            print(f"merged {base}.csv: {count} rows")
        return 0
    if args.clear_claims:
        print(f"removed {clear_claims(results)} claims")
    if args.worker is not None and args.nworkers < 1:
        print("--nworkers must be >= 1", file=sys.stderr)
        return 1
    build_dir = Path(args.build_dir) if args.build_dir else binary.parents[2]
    worker_mode = args.worker is not None
    prov_path = results / (
        f"provenance.worker{args.worker}.json" if worker_mode else "provenance.json"
    )
    prov = collect_provenance(binary, build_dir, REPO, env)
    prov["workers"] = args.nworkers if worker_mode else max(args.workers, 1)
    if worker_mode:
        prov["worker"] = args.worker
    if args.ref_bin:
        prov["ref_binary"] = args.ref_bin
        prov["ref_binary_sha256"] = hashlib.sha256(Path(args.ref_bin).read_bytes()).hexdigest()
        prov["ref_benchmark_version_output"] = binary_version(Path(args.ref_bin))
    prov["arguments"] = sys.argv[1:]
    if args.machine_note:
        prov["machine_note"] = args.machine_note
    if worker_mode:
        prov["runs"] = []
        save_json(prov_path, prov)
    elif prov_path.exists():
        import json

        old = json.loads(prov_path.read_text())
        old.setdefault("runs", []).append(prov)
        save_json(prov_path, old)
    else:
        prov["runs"] = []
        save_json(prov_path, prov)

    if args.validate or args.validate_only:
        ok = run_validation(args, env)
        if args.validate_only:
            return 0 if ok else 1
        if not ok:
            print(
                "validation failed; refusing to run the sweeps (use --family without --validate to override)"
            )
            return 1

    families = args.family
    if families == ["all"]:
        families = FAMILY_ORDER
    families = [f for f in FAMILY_ORDER if f in families]
    if args.list:
        total = 0
        for name in families:
            cells = FAMILY_BUILDERS[name](args.include_large)
            timed = sum(1 for c in cells if c.timed)
            print(f"{name}: {len(cells)} cells ({timed} timed)")
            total += len(cells)
        print(f"total: {total}")
        return 0
    progress = None
    if args.progress and not worker_mode:
        plan = []
        for name in families:
            plan += pending_cells(FAMILY_BUILDERS[name](args.include_large), args)
        progress = Progress(plan, results)
        print(progress.report(), flush=True)
    if args.workers > 1 and not worker_mode:
        return spawn_workers(args, progress)
    for name in families:
        cells = FAMILY_BUILDERS[name](args.include_large)
        run_family(name, cells, args, env, progress)
    if args.log_bin and "B3" in families:
        run_b3_log(args, env)
    return 0


if __name__ == "__main__":
    sys.exit(main())
