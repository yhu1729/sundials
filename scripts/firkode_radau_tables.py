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
# Derive the Radau IIA constants used by FIRKODE (src/firkode/
# firkode_tables.c and test/unit_tests/firkode) in high precision.
#
# For each stage count s this prints the nodes c, the real shift
# gamma0 (1/U1 with U1 the real eigenvalue of A^{-1} for odd s, and
# 1/Re(lambda) with lambda the smallest-modulus eigenvalue of A^{-1}
# for even s), the error-estimate weights e/gamma0, and R(inf).
#
# Requires mpmath:  python3 scripts/firkode_radau_tables.py [smax]
# ---------------------------------------------------------------

import sys

from mpmath import eig, findroot, inverse, legendre, matrix, mp, mpf, nstr

mp.dps = 40


def radau_nodes(s):
    """Interior roots of P_s(x) - P_{s-1}(x) on (-1, 1), mapped to (0, 1),
    followed by the endpoint node c_s = 1."""
    if s == 1:
        return [mpf(1)]
    f = lambda x: legendre(s, x) - legendre(s - 1, x)
    n = 400 * s
    xs = [mpf(-1) + mpf(2) * j / n for j in range(n)]
    roots = []
    for a, b in zip(xs[:-1], xs[1:]):
        if f(a) * f(b) < 0:
            roots.append(
                findroot(f, (a, b), solver="bisect", tol=mpf(10) ** (-38), maxsteps=400)
            )
    assert len(roots) == s - 1, (s, len(roots))
    return [(r + 1) / 2 for r in roots] + [mpf(1)]


def lagrange_coeffs(c, j):
    """Monomial coefficients of the Lagrange basis polynomial l_j on nodes c."""
    coeffs = [mpf(1)]
    for k, ck in enumerate(c):
        if k == j:
            continue
        new = [mpf(0)] * (len(coeffs) + 1)
        for m, cm in enumerate(coeffs):
            new[m + 1] += cm
            new[m] -= cm * ck
        coeffs = [x / (c[j] - ck) for x in new]
    return coeffs


def radau_table(s):
    c = radau_nodes(s)
    A = matrix(s, s)
    for j in range(s):
        r = lagrange_coeffs(c, j)
        for i in range(s):
            A[i, j] = sum(rm * c[i] ** (m + 1) / (m + 1) for m, rm in enumerate(r))
    return c, A


def main(smax):
    for s in range(1, smax + 1):
        c, A = radau_table(s)
        Ainv = inverse(A)
        ev = eig(Ainv, left=False, right=False)
        reals = [e.real for e in ev if abs(e.imag) < mpf(10) ** -25]
        if reals:
            gamma0, how = 1 / reals[0], "1/U1, U1 real eigenvalue of A^{-1}"
        else:
            lam = min(ev, key=abs)
            gamma0, how = (
                1 / lam.real,
                "1/Re(lambda), lambda min-modulus eigenvalue of A^{-1}",
            )
        l0 = [lagrange_coeffs(c, i)[0] for i in range(s)]
        e_over_g = [-sum(l0[i] * Ainv[i, j] for i in range(s)) for j in range(s)]
        one = matrix([1] * s)
        b = [A[s - 1, j] for j in range(s)]
        rinf = 1 - sum(b[i] * (Ainv * one)[i] for i in range(s))
        print(f"s = {s}")
        print(f"  gamma0   = {nstr(gamma0, 25)}L  ({how})")
        print("  c        = " + ", ".join(nstr(x, 25) + "L" for x in c))
        print("  e/gamma0 = " + ", ".join(nstr(x, 20) for x in e_over_g))
        print(f"  R(inf)   = {nstr(rinf, 5)}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 9)
