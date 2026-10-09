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
# Stability plots of the Radau IIA methods used by FIRKODE
# (Appendix A of src/firkode/DESIGN.md).
#
# The stability function R(z) of the s-stage method is the (s-1, s)
# Pade approximant of exp(z). It is evaluated from exact rational
# coefficients and checked against R(z) = 1 + z b^T (I - zA)^{-1} 1
# for the tables of scripts/firkode_radau_tables.py. The figures
# show the stability regions, the order stars, |R(z)| along the
# imaginary and the negative real axis, and the filtered error
# estimate for the Dahlquist test equation.
#
# Requires numpy, matplotlib and mpmath:
#   python3 src/firkode/firkode_stability_plot.py [--outdir DIR]
#       [--smax S] [--grid N] [--dpi D]
# ---------------------------------------------------------------

import argparse
import math
import os
import platform
import sys
from dataclasses import dataclass
from fractions import Fraction

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import mpmath  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from numpy.polynomial import polynomial as npoly  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))
from firkode_radau_tables import lagrange_coeffs, radau_table  # noqa: E402

mpmath.mp.dps = 50
EPS = np.finfo(float).eps


def _mp(x):
    """Fraction -> mpf."""
    return mpmath.mpf(x.numerator) / x.denominator


def pade_coefficients(s):
    """Ascending exact coefficients of P and Q in R = P/Q, the (s-1, s) Pade approximant of exp."""
    f = math.factorial
    k, m = s - 1, s

    def coeff(n, j):  # coefficient of z^j in the degree-n Pade polynomial, sign aside
        return Fraction(f(k + m - j) * f(n), f(k + m) * f(j) * f(n - j))

    P = [coeff(k, j) for j in range(k + 1)]
    Q = [(-1) ** j * coeff(m, j) for j in range(m + 1)]
    return P, Q


def table_stability_function(A, z):
    """R(z) = 1 + z b^T (I - zA)^{-1} 1 from the Butcher table; b is the last row of A."""
    s = A.rows
    Z = mpmath.lu_solve(mpmath.eye(s) - z * A, mpmath.ones(s, 1))
    return 1 + z * sum(A[s - 1, j] * Z[j] for j in range(s))


def real_shift(s, poles):
    """gamma0 of DESIGN.md: 1/U1 with U1 the real eigenvalue of A^{-1} for odd s, and
    1/Re(lambda) with lambda the smallest-modulus eigenvalue for even s."""
    if s % 2:
        real = [p for p in poles if abs(p.imag) < mpmath.mpf(10) ** -30 * abs(p)]
        if len(real) != 1:
            raise RuntimeError(f"s = {s}: expected one real eigenvalue of A^-1, found {len(real)}")
        return 1 / real[0].real
    return 1 / min(poles, key=abs).real


def estimate_weights(c, A, gamma0):
    """e = -gamma0 l(0)^T A^{-1}, with l_i the Lagrange basis on the nodes c."""
    s = len(c)
    Ainv = mpmath.inverse(A)
    l0 = [lagrange_coeffs(c, i)[0] for i in range(s)]
    return [-gamma0 * sum(l0[i] * Ainv[i, j] for i in range(s)) for j in range(s)]


def error_estimate(m, z, filtered=True):
    """err / y_n for y' = lambda y and z = h lambda: the stage increments are
    Z = ((I - zA)^{-1} - I) 1 and err = [gamma0 z + e^T Z] / (1 - gamma0 z);
    the unfiltered estimate is the bracket alone."""
    one = mpmath.ones(m.s, 1)
    Z = mpmath.lu_solve(mpmath.eye(m.s) - z * m.A, one) - one
    num = m.gamma0 * z + sum(m.e[j] * Z[j] for j in range(m.s))
    return num / (1 - m.gamma0 * z) if filtered else num


@dataclass
class Method:
    s: int
    P: np.ndarray  # ascending float coefficients of the numerator of R
    Q: np.ndarray  # ascending float coefficients of the denominator of R
    Pmp: list  # the same, descending mpf, for mpmath.polyval
    Qmp: list
    poles: np.ndarray  # eigenvalues of A^{-1}
    zeros: np.ndarray  # roots of P
    gamma0: mpmath.mpf
    C: float  # exp(z) - R(z) = C z^{2s} + ...
    A: mpmath.matrix
    e: list

    def R(self, z):
        """R(z) on a complex array."""
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            return npoly.polyval(z, self.P) / npoly.polyval(z, self.Q)

    def R_mp(self, z):
        return mpmath.polyval(self.Pmp, z) / mpmath.polyval(self.Qmp, z)

    @property
    def unresolved_radius(self):
        """Below this |z|, double precision cannot tell |R(z)| from |exp(z)|."""
        return 2.0 * (EPS / self.C) ** (1.0 / (2 * self.s))


def build_method(s):
    P, Q = pade_coefficients(s)
    c, A = radau_table(s)
    Pmp, Qmp = [_mp(x) for x in P[::-1]], [_mp(x) for x in Q[::-1]]
    for z in (mpmath.mpf(-1), mpmath.mpc(1, 2), mpmath.mpc(10, -3), mpmath.mpf(-50)):
        pade = mpmath.polyval(Pmp, z) / mpmath.polyval(Qmp, z)
        table = table_stability_function(A, z)
        if abs(pade - table) > mpmath.mpf(10) ** -30 * abs(table):
            raise RuntimeError(f"s = {s}: Pade and table stability functions differ at z = {z}")
    poles = mpmath.eig(mpmath.inverse(A), left=False, right=False)
    zeros = mpmath.polyroots(Pmp, maxsteps=200, extraprec=100) if s > 1 else []
    gamma0 = real_shift(s, poles)
    f = math.factorial
    C = Fraction(f(s - 1) * f(s), f(2 * s - 1) * f(2 * s))
    return Method(
        s,
        np.array([float(x) for x in P]),
        np.array([float(x) for x in Q]),
        Pmp,
        Qmp,
        np.array([complex(p) for p in poles]),
        np.array([complex(x) for x in zeros], dtype=complex),
        gamma0,
        float(C),
        A,
        estimate_weights(c, A, gamma0),
    )


def colors(methods):
    return {m.s: plt.cm.viridis(0.9 * i / max(len(methods) - 1, 1)) for i, m in enumerate(methods)}


def plot_stability_region(methods, outfile, grid, dpi):
    """|R(z)| = 1 for every s on one axes; the regions |R(z)| > 1 are shaded."""
    X, Y = np.meshgrid(np.linspace(-10, 150, 641), np.linspace(-100, 100, 801))
    unstable = np.any([np.abs(m.R(X + 1j * Y)) > 1 for m in methods], axis=0)
    if unstable[0].any() or unstable[-1].any() or unstable[:, -1].any():
        raise RuntimeError("an instability region leaves the scan box")
    xmax, ymax = 1.08 * X[unstable].max(), 1.08 * np.abs(Y[unstable]).max()
    xmin = -0.15 * xmax
    X, Y = np.meshgrid(np.linspace(xmin, xmax, grid), np.linspace(-ymax, ymax, grid))
    col = colors(methods)
    fig, ax = plt.subplots(figsize=(7.5, 7.5 * 2 * ymax / (xmax - xmin)))
    for m in methods:
        Rabs = np.nan_to_num(np.abs(m.R(X + 1j * Y)), nan=0.0, posinf=1e300)
        box = Rabs > 1
        print(
            f"s = {m.s}: |R(z)| > 1 reaches Re z = {X[box].max():.1f}, "
            f"|Im z| = {np.abs(Y[box]).max():.1f}"
        )
        ax.contourf(X, Y, Rabs, levels=[1.0, 1e300], colors=[col[m.s]], alpha=0.08)
        ax.contour(X, Y, Rabs, levels=[1.0], colors=[col[m.s]], linewidths=1.3)
    ax.axhline(0.0, color="0.4", lw=0.7)
    ax.axvline(0.0, color="0.4", lw=0.7)
    ax.set_aspect("equal")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-ymax, ymax)
    ax.set_xlabel(r"$\mathrm{Re}\,z$")
    ax.set_ylabel(r"$\mathrm{Im}\,z$")
    ax.set_title(r"Radau IIA: $|R(z)| = 1$, shaded where $|R(z)| > 1$")
    handles = [Line2D([], [], color=col[m.s], label=f"$s = {m.s}$") for m in methods]
    ax.legend(handles=handles, loc="upper left", fontsize=9)
    fig.savefig(outfile, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def order_star(m, X, Y):
    """log|R(z) exp(-z)| on the grid; points inside the unresolved disk are redone in mpmath."""
    Z = X + 1j * Y
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        F = np.log(np.abs(m.R(Z))) - X
    for i, j in zip(*np.nonzero(np.abs(Z) < m.unresolved_radius)):
        z = mpmath.mpc(X[i, j], Y[i, j])  # subtract in mpmath: the difference is O(C z^{2s})
        F[i, j] = float(mpmath.log(abs(m.R_mp(z))) - z.real)
    return np.clip(np.nan_to_num(F, nan=0.0, posinf=1e3, neginf=-1e3), -1e3, 1e3)


def plot_order_star(methods, outfile, grid, dpi):
    """{z : |R(z)| > |exp(z)|} shaded, with the poles (o) and zeros (x) of R."""
    ncols = 3
    nrows = math.ceil(len(methods) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 3.6 * nrows), squeeze=False)
    col = colors(methods)
    n = max(grid // 2, 200)
    for ax, m in zip(axes.flat, methods):
        L = 1.35 * max(np.abs(np.concatenate([m.poles, m.zeros])).max(), 2.0)
        X, Y = np.meshgrid(np.linspace(-L, L, n), np.linspace(-L, L, n))
        F = order_star(m, X, Y)
        ax.contourf(X, Y, F, levels=[0.0, 1e3], colors=[col[m.s]], alpha=0.45)
        ax.contour(X, Y, F, levels=[0.0], colors="k", linewidths=0.5)
        ax.plot(m.poles.real, m.poles.imag, "o", mfc="none", mec="#b03030", ms=5, mew=1.2)
        ax.plot(m.zeros.real, m.zeros.imag, "kx", ms=5, mew=1.2)
        ax.axhline(0.0, color="0.4", lw=0.6)
        ax.axvline(0.0, color="0.4", lw=0.6)
        ax.set_aspect("equal")
        ax.set_xlim(-L, L)
        ax.set_ylim(-L, L)
        ax.set_title(f"$s = {m.s}$, order {2 * m.s - 1}", fontsize=10)
        ax.tick_params(labelsize=8)
    for ax in axes.flat[len(methods) :]:
        ax.set_axis_off()
    fig.suptitle(r"Order stars $\{z : |R(z)| > |e^z|\}$; $\circ$ poles and $\times$ zeros of $R$")
    fig.tight_layout()
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def plot_axis_profile(methods, outfile, dpi):
    """|R(iy)| and |R(x)| for x < 0 on log-log axes."""
    t = np.logspace(-2, 4, 600)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    col = colors(methods)
    for m in methods:
        Ri = np.abs(m.R(1j * t))
        if Ri.max() > 1 + 1e-12:
            raise RuntimeError(f"s = {m.s}: |R(iy)| = {Ri.max()} > 1")
        ax1.loglog(t, Ri, color=col[m.s], label=f"$s = {m.s}$")
        ax2.loglog(t, np.abs(m.R(-t)), color=col[m.s])
    ax2.loglog(t, 1 / t, "k--", lw=0.8, label="$1/|z|$")
    for ax in (ax1, ax2):
        ax.axhline(1.0, color="0.5", lw=0.7)
        ax.grid(True, alpha=0.3)
    ax1.set_xlabel("$y$")
    ax1.set_ylabel("$|R(iy)|$")
    ax1.set_title("imaginary axis: A-stability")
    ax1.legend(fontsize=8, ncol=3)
    ax2.set_xlabel("$-x$")
    ax2.set_ylabel("$|R(x)|$")
    ax2.set_title(r"negative real axis: L-stability, $R(x) \sim s / x$")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def plot_error_estimate(methods, outfile, dpi):
    """|err / y_n| of the filtered estimate along the negative real and the imaginary axis."""
    t = np.logspace(-2, 6, 240)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    col = colors(methods)
    for m in methods:
        tail = error_estimate(m, mpmath.mpf(-1e10))
        if abs(tail + 1) > 1e-6:
            raise RuntimeError(f"s = {m.s}: filtered estimate at z = -1e10 is {tail}, not -1")
        Ex = [float(abs(error_estimate(m, mpmath.mpf(-x)))) for x in t]
        Ey = [float(abs(error_estimate(m, mpmath.mpc(0, y)))) for y in t]
        ax1.loglog(t, Ex, color=col[m.s], label=f"$s = {m.s}$")
        ax2.loglog(t, Ey, color=col[m.s])
        if m.s == 3:
            raw = [float(abs(error_estimate(m, mpmath.mpf(-x), filtered=False))) for x in t]
            ax1.loglog(t, raw, "--", color=col[m.s], lw=0.9, label="$s = 3$, unfiltered")
    for ax in (ax1, ax2):
        ax.axhline(1.0, color="0.5", lw=0.7)
        ax.set_ylim(1e-16, 1e2)
        ax.grid(True, alpha=0.3)
    ax1.set_xlabel("$-x$")
    ax1.set_ylabel(r"$|\mathrm{err} / y_n|$")
    ax1.set_title(r"negative real axis: $\mathrm{err} \to -y_n$ as $x \to -\infty$")
    ax1.legend(fontsize=8, ncol=2)
    ax2.set_xlabel("$y$")
    ax2.set_ylabel(r"$|\mathrm{err} / y_n|$")
    ax2.set_title(r"imaginary axis, $z = iy$")
    fig.tight_layout()
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Stability plots of the FIRKODE Radau IIA methods.")
    ap.add_argument("--outdir", default=os.path.join(REPO, "doc", "shared", "figs", "firkode"))
    ap.add_argument("--smax", type=int, default=9, help="largest stage count (default: 9)")
    ap.add_argument("--grid", type=int, default=800, help="grid points per axis (default: 800)")
    ap.add_argument("--dpi", type=int, default=200, help="PNG resolution (default: 200)")
    args = ap.parse_args(argv)

    print(
        f"python {platform.python_version()}, numpy {np.__version__}, "
        f"matplotlib {matplotlib.__version__}, mpmath {mpmath.__version__}"
    )
    methods = [build_method(s) for s in range(1, args.smax + 1)]
    for m in methods:
        # err = O(z^{s+1}) for z -> 0: the log-log slope between z = -1e-3 and -1e-2 is s + 1.
        e1, e2 = (abs(error_estimate(m, mpmath.mpf(-x))) for x in (1e-3, 1e-2))
        print(
            f"s = {m.s}: gamma0 = {mpmath.nstr(m.gamma0, 20)}  max|pole| = {abs(m.poles).max():.3f}"
            f"  unresolved disk r = {m.unresolved_radius:.1e}"
            f"  estimate slope = {float(mpmath.log10(e2 / e1)):.3f} (s + 1 = {m.s + 1})"
        )

    os.makedirs(args.outdir, exist_ok=True)
    out = lambda name: os.path.join(args.outdir, f"RADAU_IIA_{name}.png")  # noqa: E731
    plot_stability_region(methods, out("stability_region"), args.grid, args.dpi)
    plot_order_star(methods, out("order_star"), args.grid, args.dpi)
    plot_axis_profile(methods, out("axis_profile"), args.dpi)
    plot_error_estimate(methods, out("error_estimate"), args.dpi)
    print(f"wrote 4 figures to {args.outdir}")


if __name__ == "__main__":
    main()
