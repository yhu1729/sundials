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
"""Tables and figures of the FIRKODE-versus-BDF families B8-B13 (BENCHMARK_BDF.md).

Every table is a Markdown fragment keyed by the name of its generated
block; make_report.py injects them. The matched-accuracy ratios use the
same interpolation in (mescd, log time) as the B1 speed-up tables; "best
BDF" is the fastest of the BDF configurations of a cell group at the
target accuracy. Rows of order-histogram and profiled twins (labels
ending in "[order]" and "[profile]") never enter a timing comparison.
"""

from __future__ import annotations

import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from benchmark_common import style_for, to_float  # noqa: E402
from make_report import (  # noqa: E402
    at_floor,
    fig_grid_wp,
    fmt,
    fmt_ms,
    label_sort_key,
    md_table,
    mescd,
    ok_rows,
    speedup_table,
    status_of,
    wp_table,
)
from report_common import (  # noqa: E402
    BASE,
    best_ratio,
    fig_profile,
    groups_by_params,
    is_bdf,
    is_twin,
    params_of,
    pnum,
    profile_table,
    timed_rows,
)
from report_common import crossover_mescd as class_crossover_mescd  # noqa: E402
from report_common import ratio_cells as class_ratio_cells  # noqa: E402

TARGETS = (4.0, 6.0, 8.0)


def best_bdf_ratio(rows, target, base=BASE):
    """(ratio, label) of the fastest BDF configuration (CVODE or IDA)
    against the base at the target mescd; NaN when either side lacks data."""
    return best_ratio(rows, target, is_bdf, base)


def ratio_cells(rows):
    return class_ratio_cells(rows, is_bdf, TARGETS)


def order_summary(rows, rtol, label):
    """time-weighted mean BDF order of the [order] twin of the given
    configuration (e.g. "CVODE BDF") at the given rtol."""
    for r in rows:
        if (
            r.get("label") == f"{label} [order]"
            and abs(to_float(r.get("rtol")) - rtol) <= 1e-3 * rtol
        ):
            q = to_float(r.get("q_mean_time"))
            if not math.isnan(q):
                return q, r.get("q_hist", "NA"), r.get("label", "")
    return math.nan, "NA", ""


def order_table(rows):
    header = [
        "method",
        "params",
        "rtol",
        "steps",
        "q_hist (q=1..5)",
        "mean q (steps)",
        "mean q (time)",
    ]
    out = []
    for r in sorted(
        [r for r in rows if "[order]" in r.get("label", "")],
        key=lambda r: (r.get("params", ""), label_sort_key(r["label"]), -to_float(r.get("rtol"))),
    ):
        out.append(
            [
                r["label"].replace(" [order]", ""),
                r.get("params", ""),
                f"{to_float(r.get('rtol')):.0e}",
                r.get("nst", "NA"),
                r.get("q_hist", "NA"),
                fmt(to_float(r.get("q_mean_step")), 3),
                fmt(to_float(r.get("q_mean_time")), 3),
            ]
        )
    return md_table(header, out) if out else "(no order histograms)\n"


# ---------------------------------------------------------------------------
# B8 SECTOR
# ---------------------------------------------------------------------------


def b8_tables(rows, figdir, plt):
    tables = {}
    main = [
        r
        for r in rows
        if params_of(r).get("m") == "100"
        and "coupled" not in r.get("params", "")
        and "delta" not in r.get("params", "")
    ]
    header = (
        ["theta (deg)", "r range"]
        + [f"t(best BDF) / t(FIRKODE) at mescd={t:g}" for t in TARGETS]
        + [f"{code} order (time-weighted) at rtol 1e-6" for code in ("CVODE", "IDA")]
        + ["k_lin at rtol 1e-6"]
    )
    out = []
    series = defaultdict(list)  # (rmax, target) -> [(theta, ratio)]
    for (theta, rmax), g in sorted(
        groups_by_params(main, ["theta", "rmax"]).items(),
        key=lambda kv: (float(kv[0][1]), float(kv[0][0])),
    ):
        cells = [theta, f"[1, {float(rmax):g}]"] + ratio_cells(g)
        for label in ("CVODE BDF", "IDA BDF"):
            q, _, _ = order_summary(g, 1e-6, label)
            cells.append(fmt(q, 3))
        klin = math.nan
        for r in ok_rows(g):
            if r.get("label") == BASE and abs(to_float(r.get("rtol")) - 1e-6) <= 1e-9:
                nni = to_float(r.get("nni"))
                klin = to_float(r.get("nsli")) / nni if nni > 0 else math.nan
        cells.append(fmt(klin, 3))
        out.append(cells)
        for t in TARGETS:
            ratio, _ = best_bdf_ratio(g, t)
            if not math.isnan(ratio):
                series[(float(rmax), t)].append((float(theta), ratio))
    tables["B8_ratio_theta"] = md_table(header, out) if out else "(no B8 data)\n"
    tables["B8_orders"] = order_table(rows)
    other = [r for r in rows if r not in main and not is_twin(r.get("label", ""))]
    if other:
        tables["B8_variants"] = wp_table(other, [1e-6, 1e-8], extra_cols=("params", "nsli"))
    if plt is not None and series:
        fig, ax = plt.subplots(figsize=(6, 4))
        styles = {4.0: "o-", 6.0: "s-", 8.0: "^-"}
        for (rmax, t), pts in sorted(series.items()):
            pts.sort()
            ax.plot(
                [p[0] for p in pts],
                [p[1] for p in pts],
                styles[t],
                ms=4,
                lw=1.2,
                label=f"mescd {t:g}, r in [1, {rmax:g}]",
                alpha=0.9 if rmax < 1e4 else 0.5,
            )
        ax.axhline(1.0, color="k", lw=0.8)
        for q, a in ((3, 86.03), (4, 73.35), (5, 51.84)):
            ax.axvline(a, color="grey", lw=0.6, ls=":")
            ax.text(a, ax.get_ylim()[1], f"BDF{q}", fontsize=7, ha="center", va="bottom")
        ax.set_yscale("log")
        ax.set_xlabel("angle of the stiff eigenvalues from the negative real axis (deg)")
        ax.set_ylabel("t(best BDF) / t(FIRKODE s=3) at matched accuracy")
        ax.legend(loc="best")
        fig.tight_layout()
        fig.savefig(Path(figdir) / "bdf_sector_ratio.png", dpi=150)
        plt.close(fig)
        panels = []
        for theta in ["0", "60", "80", "89.5"]:
            g = [r for r in main if pnum(r, "theta") == float(theta) and pnum(r, "rmax") == 1000.0]
            if g:
                panels.append((f"SECTOR theta = {theta} deg, r in [1, 1e3]", g))
        if panels:
            fig_grid_wp(plt, panels, Path(figdir) / "bdf_sector_wp.png")
    return tables


# ---------------------------------------------------------------------------
# generic family tables (B9, B10, B12, B13)
# ---------------------------------------------------------------------------


def family_ratio_table(rows, keys, title_keys=None):
    header = [k for k in keys] + [f"t(best BDF) / t(FIRKODE) at mescd={t:g}" for t in TARGETS]
    out = []
    for key, g in sorted(groups_by_params(rows, keys).items()):
        out.append(list(key) + ratio_cells(g))
    return md_table(header, out) if out else "(no data)\n"


def problem_ratio_table(rows):
    header = ["problem", "params"] + [f"t(best BDF) / t(FIRKODE) at mescd={t:g}" for t in TARGETS]
    out = []
    groups = defaultdict(list)
    for r in rows:
        groups[(r.get("problem", ""), r.get("params", ""))].append(r)
    for (p, prm), g in sorted(groups.items()):
        out.append([p, prm] + ratio_cells(g))
    return md_table(header, out) if out else "(no data)\n"


def b9_tables(rows, figdir, plt):
    tables = {"B9_ratio": family_ratio_table(timed_rows(rows), ["n", "damping", "mass"])}
    tables["B9_wp"] = wp_table(
        timed_rows(rows), [1e-4, 1e-6, 1e-8], extra_cols=("params", "nsli", "nmass_solves")
    )
    tables["B9_orders"] = order_table(rows)
    if plt is not None:
        panels = []
        for n in ["100", "400"]:
            for damping, mass in (
                ("mass", "lumped"),
                ("kv", "lumped"),
                ("mass", "consistent"),
                ("kv", "consistent"),
            ):
                g = [
                    r
                    for r in timed_rows(rows)
                    if params_of(r) == {"n": n, "damping": damping, "mass": mass}
                ]
                if g:
                    panels.append((f"wave1dfem n={n}, {damping} damping, {mass} mass", g))
        if panels:
            fig_grid_wp(plt, panels, Path(figdir) / "bdf_wave_wp.png")
    return tables


def spectrum_table(results: Path):
    """Gate B10.0: the share of each dumped spectrum beyond the BDF sectors."""
    jacdir = results / "jac"
    if not jacdir.exists():
        return "(no Jacobian dumps)\n"
    try:
        import numpy as np  # noqa: F401

        from spectrum import analyze, bdf_alpha, read_matrix
    except ImportError as exc:
        return f"(spectrum analysis unavailable: {exc})\n"
    alphas = {q: bdf_alpha(q) for q in (3, 4, 5)}
    header = [
        "problem",
        "t",
        "rho(J)",
        "max angle (deg)",
        "share > 86 deg (count / modulus)",
        "share > 73 deg",
        "share > 52 deg",
    ]
    out = []
    for pdir in sorted(jacdir.iterdir()):
        for path in sorted(pdir.glob("jac_*.txt"), key=lambda p: int(p.stem.split("_")[1])):
            t, jac = read_matrix(str(path))
            res = analyze(jac, None, alphas)
            out.append(
                [
                    pdir.name,
                    f"{t:.4g}",
                    f"{res['spectral_radius']:.3g}",
                    f"{res['max_angle_deg']:.2f}",
                ]
                + [
                    f"{res[f'share_count_beyond_bdf{q}']:.2f} / {res[f'share_weight_beyond_bdf{q}']:.2f}"
                    for q in (3, 4, 5)
                ]
            )
    return md_table(header, out) if out else "(no Jacobian dumps)\n"


def b10_tables(rows, results, figdir, plt):
    tables = {"B10_ratio": problem_ratio_table(timed_rows(rows))}
    tables["B10_wp"] = wp_table(timed_rows(rows), [1e-4, 1e-6, 1e-8, 1e-10], extra_cols=("nsli",))
    tables["B10_orders"] = order_table(rows)
    tables["B10_spectrum"] = spectrum_table(results)
    if plt is not None:
        panels = []
        for p in ("beam", "plate", "ringmod"):
            g = [r for r in timed_rows(rows) if r.get("problem") == p]
            if g:
                panels.append((p, g))
        if panels:
            fig_grid_wp(plt, panels, Path(figdir) / "bdf_published_wp.png")
    return tables


def b11_tables(rows, figdir, plt):
    header = [
        "problem",
        "params",
        "method",
        "rtol",
        "steps",
        "steps / interval",
        "Newton",
        "setups",
        "mescd",
        "time (ms)",
        "status",
    ]
    out = []
    series = defaultdict(list)
    for r in sorted(
        timed_rows(rows),
        key=lambda r: (
            r.get("problem", ""),
            r.get("params", ""),
            label_sort_key(r["label"]),
            -to_float(r.get("rtol")),
        ),
    ):
        K = pnum(r, "restart-every", 1.0)
        nst = to_float(r.get("nst"))
        m = mescd(r)
        out.append(
            [
                r.get("problem", ""),
                r.get("params", ""),
                r["label"],
                f"{to_float(r.get('rtol')):.0e}",
                r.get("nst", "NA"),
                fmt(nst / K, 3) if nst > 0 else "NA",
                r.get("nni", "NA"),
                r.get("nsetups", "NA"),
                (fmt(m, 3) + ("*" if at_floor(r) else "")) if not math.isnan(m) else "NA",
                fmt_ms(to_float(r.get("time_median")))
                if r.get("timing_reliable") == "1"
                else "NA",
                status_of(r),
            ]
        )
        if r.get("timing_reliable") == "1" and abs(to_float(r.get("rtol")) - 1e-6) <= 1e-9:
            variant = r.get("problem", "") + (
                " nb=" + params_of(r)["nb"] if "nb" in params_of(r) else ""
            )
            series[(variant, r["label"])].append((K, to_float(r.get("time_median"))))
    tables = {"B11_restarts": md_table(header, out) if out else "(no B11 data)\n"}
    tables["B11_ratio"] = problem_ratio_table(timed_rows(rows))
    tables["B11_orders"] = order_table(rows)
    if plt is not None and series:
        fig, ax = plt.subplots(figsize=(6, 4))
        for (variant, label), pts in sorted(series.items()):
            pts.sort()
            st = style_for(label)
            ax.plot(
                [p[0] for p in pts],
                [p[1] for p in pts],
                color=st["color"],
                ls=st["ls"],
                marker=st["marker"],
                mfc=st["mfc"],
                ms=5,
                lw=1.2,
                label=f"{variant}: {label}",
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("number of intervals K (re-initializations K - 1)")
        ax.set_ylabel("wall time at rtol 1e-6 (s)")
        ax.legend(loc="best", fontsize=7)
        fig.tight_layout()
        fig.savefig(Path(figdir) / "bdf_restarts.png", dpi=150)
        plt.close(fig)
    return tables


def b12_tables(rows, figdir, plt):
    tables = {"B12_ratio": problem_ratio_table(timed_rows(rows))}
    tables["B12_wp"] = wp_table(timed_rows(rows), None, extra_cols=("params", "nsli"))
    tables["B12_profile"] = profile_table(rows)
    if plt is not None:
        fig_profile(plt, rows, Path(figdir) / "bdf_bruss3d_profile.png")
    return tables


def crossover_mescd(rows, base=BASE):
    """Accuracy at which the best BDF time equals the FIRKODE time, by linear
    interpolation of log(ratio) over the targets 8..13; NaN when the ratio
    never crosses 1."""
    return class_crossover_mescd(rows, is_bdf, base=base)


def b13_tables(rows, figdir, plt):
    header = (
        ["problem", "params"]
        + [f"ratio at mescd={t:g}" for t in (8.0, 10.0, 12.0)]
        + ["crossover mescd", "floor (ref disagreement)"]
    )
    out = []
    groups = defaultdict(list)
    for r in timed_rows(rows):
        groups[(r.get("problem", ""), r.get("params", ""))].append(r)
    for (p, prm), g in sorted(groups.items()):
        cells = [p, prm]
        for t in (8.0, 10.0, 12.0):
            ratio, label = best_bdf_ratio(g, t)
            cells.append("NA" if math.isnan(ratio) else fmt(ratio, 3) + f" ({label})")
        c = crossover_mescd(g)
        cells.append(fmt(c, 3) if not math.isnan(c) else "none")
        floors = [to_float(r.get("ref_disagreement")) for r in g]
        floors = [f for f in floors if not math.isnan(f) and f > 0]
        cells.append(
            f"{-math.log10(max(floors)):.1f}" if floors else "exact" if p == "sector" else "NA"
        )
        out.append(cells)
    tables = {"B13_crossover": md_table(header, out) if out else "(no B13 data)\n"}
    tables["B13_wp"] = wp_table(
        timed_rows(rows), [1e-8, 1e-10, 1e-12, 1e-13], extra_cols=("params",)
    )
    if plt is not None:
        panels = [(f"{p} {prm}".strip(), g) for (p, prm), g in sorted(groups.items())]
        if panels:
            fig_grid_wp(plt, panels, Path(figdir) / "bdf_tight_wp.png")
    return tables


def bdf_tables(data: dict, results: Path, figdir: Path, plt) -> dict:
    tables = {}
    if data.get("B8"):
        tables.update(b8_tables(data["B8"], figdir, plt))
    if data.get("B9"):
        tables.update(b9_tables(data["B9"], figdir, plt))
    if data.get("B10") or data.get("B10_spectrum"):
        tables.update(b10_tables(data.get("B10", []), results, figdir, plt))
    if data.get("B11"):
        tables.update(b11_tables(data["B11"], figdir, plt))
    if data.get("B12"):
        tables.update(b12_tables(data["B12"], figdir, plt))
    if data.get("B13"):
        tables.update(b13_tables(data["B13"], figdir, plt))
    # the baseline families in the same report
    for name in ("B1",):
        for p in ["rober", "vdpol", "orego", "hires"]:
            rows = [r for r in data.get(name, []) if r.get("problem") == p]
            if rows:
                tables[f"{name}_{p}_speedup"] = speedup_table(rows)
    return tables
