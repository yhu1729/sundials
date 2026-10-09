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
# Independent check of the wave1dfem reference solutions (family B9,
# mass-proportional damping): the semi-discrete system
#   M u'' + alpha M u' + K u = g(t) b,  u(0) = u'(0) = 0,
# with the P1 matrices M = (dx/6) tridiag(1,4,1), K = (1/dx)
# tridiag(-1,2,-1) on n interior nodes, b = M phi, phi = sin(pi x) +
# 0.5 sin(2 pi x), g = sin^2(pi t/2) = (1 - cos(pi t))/2 on [0, 1] and 1
# afterwards, is diagonalized by the generalized eigenvectors of (K, M)
# (Phi^T M Phi = I). Every mode is a damped oscillator
#   xi'' + alpha xi' + lambda xi = f(t) (Phi^T b)_k
# whose solution is written in closed form: on [0, 1] the particular
# solution of the constant plus cosine forcing, matched to zero initial
# data; on [1, tf] the constant-forcing solution matched at t = 1.
#
# Usage: wave1d_fem_exact.py REF_FILE [--n N] [--alpha A] [--mass consistent|lumped]
#   REF_FILE: a reference written by firkode_benchmark --write-ref for
#   wave1dfem with --damping mass; the script prints the mixed error
#   (max |y - yref| / (1 + |yref|)) and the maximum absolute error.
# ---------------------------------------------------------------

import argparse
import sys

import numpy as np
import scipy.linalg


def matrices(n: int, lumped: bool):
    dx = 1.0 / (n + 1)
    K = (np.diag(2.0 * np.ones(n)) - np.diag(np.ones(n - 1), 1) - np.diag(np.ones(n - 1), -1)) / dx
    if lumped:
        M = dx * np.eye(n)
    else:
        M = dx * (4.0 * np.eye(n) + np.diag(np.ones(n - 1), 1) + np.diag(np.ones(n - 1), -1)) / 6.0
    x = dx * np.arange(1, n + 1)
    phi = np.sin(np.pi * x) + 0.5 * np.sin(2.0 * np.pi * x)
    b = M @ phi
    return dx, M, K, b


def mode_solution(lam: float, alpha: float, f: float, t: float):
    """xi(t), xi'(t) of xi'' + alpha xi' + lam xi = f g(t), xi(0) = xi'(0) = 0."""
    wd = np.sqrt(lam - 0.25 * alpha * alpha)  # lam >> alpha^2/4 for every mode
    w = np.pi

    def homogeneous(c1, c2, tau):
        e = np.exp(-0.5 * alpha * tau)
        xi = e * (c1 * np.cos(wd * tau) + c2 * np.sin(wd * tau))
        dxi = e * (
            (-0.5 * alpha * c1 + wd * c2) * np.cos(wd * tau)
            + (-0.5 * alpha * c2 - wd * c1) * np.sin(wd * tau)
        )
        return xi, dxi

    # phase 1: forcing f/2 - (f/2) cos(w t); particular solution
    a0 = 0.5 * f / lam
    h = -0.5 * f / complex(lam - w * w, alpha * w)  # amplitude of Re[h e^{i w t}]

    def particular1(tau):
        xi = a0 + (h * np.exp(1j * w * tau)).real
        dxi = (1j * w * h * np.exp(1j * w * tau)).real
        return xi, dxi

    p0, dp0 = particular1(0.0)
    c1 = -p0
    c2 = (-dp0 + 0.5 * alpha * c1) / wd
    if t <= 1.0:
        xh, dxh = homogeneous(c1, c2, t)
        xp, dxp = particular1(t)
        return xh + xp, dxh + dxp
    xh, dxh = homogeneous(c1, c2, 1.0)
    xp, dxp = particular1(1.0)
    x1, dx1 = xh + xp, dxh + dxp
    # phase 2: constant forcing f, matched at t = 1
    a1 = f / lam
    c1 = x1 - a1
    c2 = (dx1 + 0.5 * alpha * c1) / wd
    xh, dxh = homogeneous(c1, c2, t - 1.0)
    return xh + a1, dxh


def exact(n: int, alpha: float, lumped: bool, t: float):
    dx, M, K, b = matrices(n, lumped)
    lam, Phi = scipy.linalg.eigh(K, M)  # K Phi = M Phi diag(lam), Phi^T M Phi = I
    fb = Phi.T @ b
    xi = np.zeros(n)
    dxi = np.zeros(n)
    for k in range(n):
        xi[k], dxi[k] = mode_solution(lam[k], alpha, fb[k], t)
    u = Phi @ xi
    v = Phi @ dxi
    y = np.empty(2 * n)
    y[0::2] = u
    y[1::2] = v
    return y


def read_reference(path: str):
    with open(path) as fp:
        lines = [l for l in fp if not l.startswith("#")]
    N, tf = lines[0].split()
    return float(tf), np.array([float(v) for v in lines[1 : 1 + int(N)]])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("ref")
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--alpha", type=float, default=0.1)
    ap.add_argument("--mass", default="consistent", choices=["consistent", "lumped"])
    args = ap.parse_args()
    tf, yref = read_reference(args.ref)
    if yref.size != 2 * args.n:
        sys.exit(f"reference has {yref.size} values, expected {2 * args.n}")
    y = exact(args.n, args.alpha, args.mass == "lumped", tf)
    mixed = np.max(np.abs(y - yref) / (1.0 + np.abs(yref)))
    print(
        f"tf = {tf:g}  n = {args.n}  max|u| = {np.max(np.abs(y[0::2])):.4e}  max|v| = {np.max(np.abs(y[1::2])):.4e}"
    )
    print(f"mixed error of the reference against the modal closed form: {mixed:.3e}")
    print(f"max absolute error: {np.max(np.abs(y - yref)):.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
