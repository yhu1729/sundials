#!/usr/bin/env python3
# ---------------------------------------------------------------
# Programmer(s): Yifan Hu @ UMBC
# ---------------------------------------------------------------
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
# ---------------------------------------------------------------
# Spectrum analysis of the Jacobians written by firkode_benchmark
# --dump-jac: the angle of every eigenvalue from the negative real
# axis and the share of the spectrum outside the A(alpha) sectors of
# the BDF methods of orders 3, 4 and 5.
#
# The A(alpha) angles are computed here from the boundary locus of the
# BDF methods rather than copied from a table: for BDF of order q the
# boundary is z(zeta) = sum_{j=1}^{q} (1 - 1/zeta)^j / j with zeta on
# the unit circle, and alpha is the smallest angle between a boundary
# point in the left half-plane and the negative real axis. The values
# 86.03, 73.35 and 51.84 degrees (Hairer and Wanner, Solving Ordinary
# Differential Equations II, Section V.2) serve as a self-check.
#
# Usage: spectrum.py JAC... [--mass FILE] [--json FILE]
#   JAC: jac_<k>.txt files (header "# t=... N=...", then N rows of N
#        values); --mass: the mass matrix of a generalized problem
#        (M y' = f, eigenvalues of M^{-1} J).
# ---------------------------------------------------------------

import argparse
import json
import math
import sys

import numpy as np

HW_ANGLES = {1: 90.0, 2: 90.0, 3: 86.03, 4: 73.35, 5: 51.84, 6: 17.84}


def bdf_alpha(q: int, n: int = 200000) -> float:
    """A(alpha) angle of the BDF method of order q, in degrees, from its boundary locus."""
    theta = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)[1:]
    w = 1.0 - np.exp(-1j * theta)
    z = sum(w**j / j for j in range(1, q + 1))
    left = z[z.real < 0.0]
    if left.size == 0:
        return 90.0
    ang = np.degrees(np.arctan2(np.abs(left.imag), -left.real))
    return float(ang.min())


def read_matrix(path: str) -> tuple[float, np.ndarray]:
    t = math.nan
    with open(path) as fp:
        header = fp.readline()
    for tok in header.lstrip("# ").split():
        if tok.startswith("t="):
            t = float(tok[2:])
    a = np.loadtxt(path, comments="#")
    a = np.atleast_2d(a)
    if a.shape[0] != a.shape[1]:
        raise ValueError(f"{path}: matrix is {a.shape}, not square")
    return t, a


def analyze(jac: np.ndarray, mass: np.ndarray | None, alphas: dict) -> dict:
    if mass is None:
        lam = np.linalg.eigvals(jac)
    else:
        import scipy.linalg

        lam = scipy.linalg.eigvals(jac, mass)
    lam = lam[np.isfinite(lam)]
    radius = float(np.max(np.abs(lam))) if lam.size else 0.0
    left = lam[lam.real < 0.0]
    angle = np.degrees(np.arctan2(np.abs(left.imag), -left.real))
    weight = np.abs(left)
    out = {
        "n": int(lam.size),
        "spectral_radius": radius,
        "n_right_half_plane": int(np.sum(lam.real >= 0.0)),
        "max_angle_deg": float(angle.max()) if angle.size else 0.0,
        "max_abs_imag": float(np.max(np.abs(lam.imag))) if lam.size else 0.0,
    }
    for q, alpha in alphas.items():
        mask = angle > alpha
        out[f"share_count_beyond_bdf{q}"] = float(mask.mean()) if angle.size else 0.0
        out[f"share_weight_beyond_bdf{q}"] = (
            float(weight[mask].sum() / weight.sum()) if weight.sum() > 0 else 0.0
        )
        out[f"max_abs_beyond_bdf{q}"] = float(weight[mask].max()) if mask.any() else 0.0
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("jac", nargs="+", help="Jacobian files written by --dump-jac")
    ap.add_argument("--mass", help="mass matrix file (generalized eigenvalues)")
    ap.add_argument("--json", help="write the per-file summary to this JSON file")
    args = ap.parse_args()

    alphas = {q: bdf_alpha(q) for q in (3, 4, 5)}
    for q, a in alphas.items():
        if abs(a - HW_ANGLES[q]) > 0.05:
            sys.exit(f"BDF{q}: computed alpha {a:.3f} differs from {HW_ANGLES[q]}")
    print("BDF A(alpha) angles (deg): " + ", ".join(f"q={q}: {a:.2f}" for q, a in alphas.items()))

    mass = read_matrix(args.mass)[1] if args.mass else None
    summary = []
    print(
        f"{'file':30s} {'t':>12s} {'N':>6s} {'rho(J)':>10s} {'max angle':>9s} "
        f"{'>BDF3 (#, |.|)':>16s} {'>BDF4 (#, |.|)':>16s} {'>BDF5 (#, |.|)':>16s}"
    )
    for path in args.jac:
        t, jac = read_matrix(path)
        res = analyze(jac, mass, alphas)
        res["file"] = path
        res["t"] = t
        summary.append(res)
        cols = " ".join(
            f"{res[f'share_count_beyond_bdf{q}']:7.3f} {res[f'share_weight_beyond_bdf{q}']:7.3f} "
            for q in (3, 4, 5)
        )
        print(
            f"{path:30s} {t:12.5g} {res['n']:6d} {res['spectral_radius']:10.3e} "
            f"{res['max_angle_deg']:9.2f} {cols}"
        )
    if args.json:
        with open(args.json, "w") as fp:
            json.dump({"bdf_alpha_deg": alphas, "files": summary}, fp, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
