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
# Generator of the Radau IIA tables of FIRKODE (design document,
# Section 5.6 and Milestone M0).
#
# For 1 <= s <= 9 the script computes, with mpmath at 50 digits, the
# nodes c, the coefficients A, the weights b, A^{-1}, d = b^T A^{-1},
# the error-estimate weights e, the real shift gamma0 and the
# dense-output coefficients P, checks the order conditions and the
# identities the integrator relies on, and writes the tables as
# SUN_RCONST literals with 40 significant digits into the generated
# region of src/firkode/firkode_tables.c. Nothing is computed at run
# time.
#
# Usage (requires mpmath):
#   firkode_radau_tables.py [--smax N]      print a summary
#   firkode_radau_tables.py --emit-c [PATH] rewrite the generated region
#   firkode_radau_tables.py --check [PATH]  exit 1 if the region is stale
#   firkode_radau_tables.py --table1        print Table 4.1 of the design
#                                           and check the gamma0 literals
# PATH defaults to src/firkode/firkode_tables.c. Exit status: 0 on
# success, 1 on a failed check, 2 on a usage error or missing mpmath.
# ---------------------------------------------------------------

import argparse
import difflib
import itertools
import math
import os
import re
import sys

try:
    from mpmath import eig, findroot, inverse, legendre, matrix, mp, mpf, nstr, polyroots
except ImportError:  # pragma: no cover
    sys.stderr.write("firkode_radau_tables.py requires mpmath (pip install mpmath)\n")
    sys.exit(2)

GENERATOR_VERSION = "1.0"
WORKING_DPS = 50  # digits of the computation
CHECK_DPS = 60  # digits of the reproducibility check
PRINTED_DIGITS = 40  # significant digits of the literals
FIRK_MAX_STAGES = 9
BEGIN_MARK = "/* BEGIN GENERATED RADAU IIA TABLES */"
END_MARK = "/* END GENERATED RADAU IIA TABLES */"

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_C_FILE = os.path.join(REPO, "src", "firkode", "firkode_tables.c")


class CheckError(Exception):
    """A self-check or consistency check failed."""


# ---------------------------------------------------------------
# Construction of the tables (at the current mp.dps)
# ---------------------------------------------------------------


def radau_nodes(s):
    """Interior roots of P_s(x) - P_{s-1}(x) on (-1, 1), mapped to (0, 1),
    followed by the endpoint node c_s = 1."""
    if s == 1:
        return [mpf(1)]
    f = lambda x: legendre(s, x) - legendre(s - 1, x)
    n = 400 * s
    xs = [mpf(-1) + mpf(2) * j / n for j in range(n)]
    tol = mpf(10) ** (-(mp.dps - 2))
    roots = []
    for a, b in itertools.pairwise(xs):
        if f(a) * f(b) < 0:
            roots.append(findroot(f, (a, b), solver="bisect", tol=tol, maxsteps=4 * mp.prec))
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
    """Nodes c and coefficients A of the s-stage Radau IIA method."""
    c = radau_nodes(s)
    A = matrix(s, s)
    for j in range(s):
        r = lagrange_coeffs(c, j)
        for i in range(s):
            A[i, j] = sum(rm * c[i] ** (m + 1) / (m + 1) for m, rm in enumerate(r))
    return c, A


def pade_denominator(s):
    """Coefficients q_0..q_s of Q(z) = det(I - zA), the denominator of the
    (s-1, s) Pade approximant of e^z (Appendix A of the design)."""
    f = math.factorial
    return [
        mpf(f(2 * s - 1 - j) * f(s)) / (f(2 * s - 1) * f(j) * f(s - j)) * (-1) ** j
        for j in range(s + 1)
    ]


def select_gamma0(s, eigs):
    """The shift rule of Section 4.5 applied to the eigenvalues of A^{-1}:
    1/U1 for the unique real eigenvalue U1 (odd s), 1/Re(lambda) for the
    eigenvalue lambda of smallest modulus (even s)."""
    tol = mpf(10) ** (-(mp.dps // 2))
    reals = [e.real for e in eigs if abs(e.imag) <= tol * abs(e)]
    if s % 2 == 1:
        if len(reals) != 1:
            raise CheckError(f"s = {s}: {len(reals)} real eigenvalues of A^-1, expected 1")
        return 1 / reals[0]
    if reals:
        raise CheckError(f"s = {s}: A^-1 has a real eigenvalue")
    return 1 / min(eigs, key=abs).real


def tables(s):
    """All table data of the s-stage method at the current mp.dps, after the
    self-checks of the design (Milestone M0)."""
    c, A = radau_table(s)
    b = [A[s - 1, j] for j in range(s)]
    Ainv = inverse(A)
    d = [sum(b[i] * Ainv[i, j] for i in range(s)) for j in range(s)]

    # gamma0 by two independent routes: the eigenvalues of A^{-1}, and the
    # roots of the Pade denominator, which are the same numbers
    g_eig = select_gamma0(s, list(eig(Ainv, left=False, right=False)))
    q = pade_denominator(s)
    g_pade = select_gamma0(s, list(polyroots(q[::-1], maxsteps=200, extraprec=4 * mp.prec)))

    l = [lagrange_coeffs(c, i) for i in range(s)]
    l0 = [l[i][0] for i in range(s)]
    e = [-g_eig * sum(l0[i] * Ainv[i, j] for i in range(s)) for j in range(s)]
    P = [[l[i][k] / c[i] for i in range(s)] for k in range(s)]

    tol = mpf(10) ** (-(mp.dps - 5))

    def check(name, residual):
        if abs(residual) > tol:
            raise CheckError(f"s = {s}: {name} residual {nstr(residual, 5)}")

    check("gamma0 routes", (g_eig - g_pade) / g_eig)
    for k in range(1, 2 * s):  # B(2s - 1)
        check(f"B({k})", sum(b[j] * c[j] ** (k - 1) for j in range(s)) - mpf(1) / k)
    for i in range(s):  # C(s)
        for k in range(1, s + 1):
            check(
                f"C({k}), row {i}",
                sum(A[i, j] * c[j] ** (k - 1) for j in range(s)) - c[i] ** k / k,
            )
    for k in range(1, s):  # D(s - 1)
        for j in range(s):
            lhs = sum(b[i] * c[i] ** (k - 1) * A[i, j] for i in range(s))
            check(f"D({k}), column {j}", lhs - b[j] * (1 - c[j] ** k) / k)
    for i in range(s):
        for j in range(s):
            check(f"A A^-1 ({i}, {j})", sum(A[i, m] * Ainv[m, j] for m in range(s)) - (i == j))
            check(f"L_{j}(c_{i})", sum(P[k][j] * c[i] ** (k + 1) for k in range(s)) - (i == j))
    for j in range(s):
        check(f"d_{j} = e_s", d[j] - (j == s - 1))
        check(f"e_{j} = -gamma0 l_{j}(0)/c_{j}", e[j] + g_eig * l0[j] / c[j])
    check("e_s = (-1)^s gamma0/s", e[s - 1] - (-1) ** s * g_eig / s)
    if s % 2 == 1:
        Ag = Ainv - (1 / g_eig) * mp.eye(s)
        check("det(A^-1 - gamma0^-1 I)", mp.det(Ag) / (1 + mp.norm(Ainv)) ** s)

    return {
        "s": s,
        "c": c,
        "A": [A[i, j] for i in range(s) for j in range(s)],
        "b": b,
        "Ainv": [Ainv[i, j] for i in range(s) for j in range(s)],
        "d": [mpf(1) if j == s - 1 else mpf(0) for j in range(s)],
        "e": e,
        "gamma0": g_eig,
        "P": [P[k][i] for k in range(s) for i in range(s)],
    }


# ---------------------------------------------------------------
# C output
# ---------------------------------------------------------------


def literal(x):
    """x as a C literal: the shortest exact decimal when x has one with at
    most 17 digits, otherwise 40 significant digits."""
    if x == 0:
        return "0.0"
    tol = abs(x) * mpf(10) ** (-(mp.dps - 8))
    for k in range(1, 18):
        t = nstr(x, k, min_fixed=-30, max_fixed=30)
        if abs(mpf(t) - x) <= tol:
            return t if ("." in t or "e" in t) else t + ".0"
    return nstr(x, PRINTED_DIGITS, strip_zeros=False, min_fixed=-30, max_fixed=30)


def emit_array(name, values, s, rows):
    """Initializer of an array field, one literal per line; rows of a
    row-major s x s array are labelled."""
    out = [f"  .{name} = {{"]
    for n, v in enumerate(values):
        if rows is not None and n % s == 0:
            out.append(f"    /* {rows}[{n // s}][:] */")
        out.append(f"    SUN_RCONST({literal(v)}),")
    out[-1] = out[-1].rstrip(",")
    out.append("  },")
    return out


def emit_region(dps):
    """The text of the generated region of firkode_tables.c."""
    with mp.workdps(dps):
        data = [tables(s) for s in range(1, FIRK_MAX_STAGES + 1)]
        lines = [
            BEGIN_MARK,
            "/* clang-format off */",
            f"/* Generated by scripts/firkode_radau_tables.py, generator version {GENERATOR_VERSION},",
            f"   with mpmath at {WORKING_DPS} digits; literals carry {PRINTED_DIGITS} significant digits.",
            "   Do not edit: rerun the generator with --emit-c. CI runs it with --check.",
            "   A and A^{-1} are row-major, A[i * s + j]; P is P[k * s + i], the",
            "   coefficient of theta^(k+1) in the dense-output basis polynomial L_i. */",
            "",
        ]
        for t in data:
            s = t["s"]
            lines.append(f"static const struct FIRKodeTableMem firk_radau_iia_{s} = {{")
            lines.append(f"  .s = {s},")
            lines.append(f"  .q = {2 * s - 1},")
            lines.append(f"  .p = {s},")
            lines += emit_array("c", t["c"], s, None)
            lines += emit_array("A", t["A"], s, "A")
            lines += emit_array("b", t["b"], s, None)
            lines += emit_array("Ainv", t["Ainv"], s, "Ainv")
            lines += emit_array("d", t["d"], s, None)
            lines += emit_array("e", t["e"], s, None)
            lines.append(f"  .gamma0 = SUN_RCONST({literal(t['gamma0'])}),")
            lines += emit_array("P", t["P"], s, "P")
            lines.append("  .stiffly_accurate = SUNTRUE};")
            lines.append("")
        lines.append(
            "static const struct FIRKodeTableMem* const firk_radau_iia[FIRK_MAX_STAGES + 1] = {"
        )
        lines.append("  NULL,")
        for s in range(1, FIRK_MAX_STAGES + 1):
            lines.append(f"  &firk_radau_iia_{s}" + ("," if s < FIRK_MAX_STAGES else "};"))
        # clang-format separates the last declaration from the next comment
        lines.append("")
        lines.append("/* clang-format on */")
        lines.append(END_MARK)
    return "\n".join(lines) + "\n"


def generated_region():
    """The generated region, after checking that a computation with more
    digits prints the same literals."""
    text = emit_region(WORKING_DPS)
    if emit_region(CHECK_DPS) != text:
        raise CheckError(
            f"the literals printed at {WORKING_DPS} and {CHECK_DPS} digits differ;"
            " increase WORKING_DPS"
        )
    return text


def split_file(path):
    """(head, region, tail) of a C file with exactly one generated region."""
    with open(path, encoding="utf-8") as f:
        text = f.read()
    if text.count(BEGIN_MARK) != 1 or text.count(END_MARK) != 1:
        raise CheckError(f"{path}: expected one '{BEGIN_MARK}' and one '{END_MARK}' line")
    i = text.index(BEGIN_MARK)
    j = text.index(END_MARK) + len(END_MARK) + 1
    if j < i:
        raise CheckError(f"{path}: the end marker precedes the begin marker")
    return text[:i], text[i:j], text[j:]


def emit_c(path):
    head, _, tail = split_file(path)
    with open(path, "w", encoding="utf-8") as f:
        f.write(head + generated_region() + tail)
    print(f"wrote the generated region of {path}")


def check(path):
    _, region, _ = split_file(path)
    expected = generated_region()
    if region == expected:
        print(f"{path}: generated region is current")
        return 0
    sys.stdout.writelines(
        difflib.unified_diff(
            region.splitlines(keepends=True),
            expected.splitlines(keepends=True),
            fromfile=path + " (committed)",
            tofile=path + " (generated)",
        )
    )
    print(f"{path}: generated region is stale; run with --emit-c", file=sys.stderr)
    return 1


# ---------------------------------------------------------------
# Table 4.1 (Appendix B of the design)
# ---------------------------------------------------------------


def contraction(s, n=4000):
    """rho_*(gamma0) and the gamma0 of the s-stage method, by the procedure
    of Appendix B in double precision."""
    with mp.workdps(WORKING_DPS):
        lam = list(polyroots(pade_denominator(s)[::-1], maxsteps=200, extraprec=4 * mp.prec))
        g0 = select_gamma0(s, lam)
        mu = [complex(1 / x) for x in lam]
    gamma0 = float(g0)
    if s == 1:
        return 0.0, g0
    pts = []
    for k in range(n):
        z = 1j * math.tan(math.pi * (k + 0.5) / n - math.pi / 2)
        pts += [(1 - z * m) / (1 - gamma0 * z) for m in mu]
    pts += [m / gamma0 for m in mu]

    def worst(cc):
        return max(abs(1 - w / cc) for w in pts)

    # bracketing scan, then golden-section search around the best point
    grid = [0.05 * (100.0 ** (i / 199)) for i in range(200)]
    vals = [worst(cc) for cc in grid]
    i = min(range(len(grid)), key=vals.__getitem__)
    lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    r = (math.sqrt(5) - 1) / 2
    x1, x2 = hi - r * (hi - lo), lo + r * (hi - lo)
    f1, f2 = worst(x1), worst(x2)
    for _ in range(80):
        if f1 <= f2:
            hi, x2, f2 = x2, x1, f1
            x1 = hi - r * (hi - lo)
            f1 = worst(x1)
        else:
            lo, x1, f1 = x1, x2, f2
            x2 = lo + r * (hi - lo)
            f2 = worst(x2)
    rho = min(f1, f2, vals[i])
    return rho, g0


def literal_gamma0(path):
    """The gamma0 literals of the generated region, keyed by s."""
    _, region, _ = split_file(path)
    out = {}
    for m in re.finditer(
        r"firk_radau_iia_(\d+) = \{.*?\.gamma0 = SUN_RCONST\(([-+0-9.eE]+)\)", region, re.DOTALL
    ):
        out[int(m.group(1))] = m.group(2)
    return out


def table1(path):
    literals = literal_gamma0(path)
    status = 0
    rows = []
    for s in range(1, FIRK_MAX_STAGES + 1):
        rho, g0 = contraction(s)
        k = 0 if rho == 0 else math.ceil(math.log(0.05) / math.log(rho))
        rows.append((s, rho, k))
        with mp.workdps(WORKING_DPS):
            lit = literals.get(s)
            if lit is None or abs(mpf(lit) - g0) > abs(g0) * mpf(10) ** (1 - PRINTED_DIGITS) / 2:
                print(f"s = {s}: gamma0 literal {lit} does not match {nstr(g0, PRINTED_DIGITS)}")
                status = 1
    print(
        "| $s$                          | " + " | ".join(f"{s:<7}" for s, _, _ in rows[1:]) + " |"
    )
    print("| ---------------------------- | " + " | ".join("-------" for _ in rows[1:]) + " |")
    print(
        "| $\\rho_\\ast(\\gamma_0)$        | "
        + " | ".join(f"{rho:<7.3f}" for _, rho, _ in rows[1:])
        + " |"
    )
    print(
        "| $k$ for a 20-fold reduction  | " + " | ".join(f"{k:<7}" for _, _, k in rows[1:]) + " |"
    )
    if status == 0:
        print(f"gamma0 literals of {path} agree with the Pade roots to all printed digits")
    return status


# ---------------------------------------------------------------
# Summary
# ---------------------------------------------------------------


def summary(smax):
    with mp.workdps(WORKING_DPS):
        for s in range(1, smax + 1):
            t = tables(s)
            print(f"s = {s}")
            print(f"  gamma0 = {nstr(t['gamma0'], 25)}")
            print("  c      = " + ", ".join(nstr(x, 25) for x in t["c"]))
            print("  e      = " + ", ".join(nstr(x, 20) for x in t["e"]))


def main(argv=None):
    parser = argparse.ArgumentParser(description="Radau IIA tables of FIRKODE")
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--emit-c",
        nargs="?",
        const=DEFAULT_C_FILE,
        metavar="PATH",
        help="rewrite the generated region of PATH",
    )
    group.add_argument(
        "--check",
        nargs="?",
        const=DEFAULT_C_FILE,
        metavar="PATH",
        help="exit 1 if the generated region of PATH is stale",
    )
    group.add_argument(
        "--table1",
        nargs="?",
        const=DEFAULT_C_FILE,
        metavar="PATH",
        help="print Table 4.1 and check the gamma0 literals of PATH",
    )
    parser.add_argument("--smax", type=int, default=FIRK_MAX_STAGES, help="summary up to s")
    args = parser.parse_args(argv)
    try:
        if args.emit_c:
            emit_c(args.emit_c)
            return 0
        if args.check:
            return check(args.check)
        if args.table1:
            return table1(args.table1)
        if not 1 <= args.smax <= FIRK_MAX_STAGES:
            parser.error(f"--smax must be between 1 and {FIRK_MAX_STAGES}")
        summary(args.smax)
        return 0
    except CheckError as err:
        print(f"error: {err}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
