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
"""Tables and figures of the FIRKODE-versus-ARKODE families A0-A8
(BENCHMARK_ARKODE.md).

Every table is a Markdown fragment keyed by the name of its generated block;
make_report.py --report arkode injects them. Ratios are t(best of a class) /
t(FIRKODE s=3) at matched mixed digits (mescd), interpolated in log time as
in the B1 speed-up tables; the classes are the ARKODE DIRK defaults of order
3-5 ("best default"), every declared DIRK configuration ("best of K"), IMEX,
LSRK and MRI. Figures plot wall time only: nfe counts different right-hand
sides for different steppers.
"""

from __future__ import annotations

import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from benchmark_common import to_float, to_int  # noqa: E402
from make_report import (  # noqa: E402
    at_floor,
    fig_fixed_step,
    fig_grid_wp,
    fmt,
    label_sort_key,
    md_table,
    mescd,
    observed_orders,
    ok_rows,
    wp_table,
)
from report_common import (  # noqa: E402
    best_ratio,
    crossover_mescd,
    fig_profile,
    groups_by_params,
    is_dirk,
    is_dirk_default,
    is_imex,
    is_lsrk,
    is_mri,
    profile_table,
    ratio_cells,
    timed_rows,
)

TARGETS = (4.0, 6.0, 8.0)
# family A0: a method passes when its observed fixed-step order is at least
# the expected order minus this; an order above the expected one plus this is
# marked (error cancellation on a particular problem, e.g. ESDIRK436L2SA on
# KPR), not failed
ORDER_TOL = 0.3


def order_verdict(observed: float, expected) -> str:
    if math.isnan(observed) or expected is None:
        return "NA"
    if observed < expected - ORDER_TOL:
        return "NO"
    if observed > expected + ORDER_TOL:
        return f"yes (above by {observed - expected:.2f})"
    return "yes"


CLASSES = [
    ("best DIRK default", is_dirk_default),
    ("best DIRK of all", is_dirk),
    ("best IMEX", is_imex),
    ("best LSRK", is_lsrk),
    ("best MRI", is_mri),
]


def rows_at(rows, label, rtol):
    return [
        r
        for r in ok_rows(rows)
        if r.get("label") == label and abs(to_float(r.get("rtol")) - rtol) <= 1e-3 * rtol
    ]


def classes_present(rows):
    labels = {r.get("label", "") for r in timed_rows(rows)}
    return [(name, pred) for name, pred in CLASSES if any(pred(lb) for lb in labels)]


def class_ratio_table(rows, keys, targets=TARGETS):
    """One row per parameter group and class present: the ratio of the
    class's fastest configuration to FIRKODE at each target."""
    header = list(keys) + ["class"] + [f"t(best) / t(FIRKODE) at mescd={t:g}" for t in targets]
    out = []
    for key, g in sorted(groups_by_params(timed_rows(rows), keys).items()):
        for name, pred in classes_present(g):
            out.append(list(key) + [name] + ratio_cells(g, pred, targets))
    return md_table(header, out) if out else "(no data)\n"


def problem_class_ratio_table(rows, targets=TARGETS):
    header = ["problem", "params", "class"] + [
        f"t(best) / t(FIRKODE) at mescd={t:g}" for t in targets
    ]
    groups = defaultdict(list)
    for r in timed_rows(rows):
        groups[(r.get("problem", ""), r.get("params", ""))].append(r)
    out = []
    for (p, prm), g in sorted(groups.items()):
        for name, pred in classes_present(g):
            out.append([p, prm, name] + ratio_cells(g, pred, targets))
    return md_table(header, out) if out else "(no data)\n"


def counters_table(rows, rtol, cols, keys=("problem", "params")):
    """Counters of every configuration at one rtol."""
    header = list(keys) + ["method", "mescd"] + [c for c, _ in cols]
    out = []
    for r in sorted(
        [
            r
            for r in ok_rows(timed_rows(rows))
            if abs(to_float(r.get("rtol")) - rtol) <= 1e-3 * rtol
        ],
        key=lambda r: (tuple(r.get(k, "") for k in keys), label_sort_key(r["label"])),
    ):
        m = mescd(r)
        out.append(
            [r.get(k, "") for k in keys]
            + [
                r["label"],
                (fmt(m, 3) + ("*" if at_floor(r) else "")) if not math.isnan(m) else "NA",
            ]
            + [f(r) for _, f in cols]
        )
    return md_table(header, out) if out else "(no data)\n"


def per(a, b):
    def f(r):
        x, y = to_float(r.get(a)), to_float(r.get(b))
        return fmt(x / y, 3) if y > 0 and not math.isnan(x) else "NA"

    return f


def col(key):
    return lambda r: r.get(key, "NA")


def time_per_step(r):
    t, n = to_float(r.get("time_median")), to_float(r.get("nst"))
    return fmt(1e6 * t / n, 3) if n > 0 and r.get("timing_reliable") == "1" else "NA"


def grid_figure(plt, rows, keys, title, path):
    panels = []
    for key, g in sorted(groups_by_params(timed_rows(rows), keys).items()):
        panels.append((title + " " + ", ".join(f"{k}={v}" for k, v in zip(keys, key)), g))
    if panels:
        fig_grid_wp(plt, panels, path)


# ---------------------------------------------------------------------------
# A0 method selection
# ---------------------------------------------------------------------------


def a0_orders(rows):
    header = [
        "problem",
        "params",
        "configuration",
        "method",
        "expected order",
        "errors at h0, h0/2, ...",
        "observed orders",
        "reported",
        "at least expected - 0.3",
    ]
    out = []
    groups = defaultdict(list)
    for r in ok_rows(rows):
        groups[(r["problem"], r.get("params", ""), r["label"])].append(r)
    for (p, prm, label), g in sorted(
        groups.items(), key=lambda kv: (kv[0][0], kv[0][1], label_sort_key(kv[0][2]))
    ):
        g.sort(key=lambda r: -to_float(r.get("fixed_h")))
        hs = [to_float(r.get("fixed_h")) for r in g]
        errs = [to_float(r.get("err_abs_max")) for r in g]
        orders, rep = observed_orders(hs, errs)
        q = to_int(g[0].get("order"))
        ok = order_verdict(rep, q)
        out.append(
            [
                p,
                prm,
                label,
                g[0].get("method", ""),
                q,
                ", ".join(f"{e:.2e}" for e in errs),
                ", ".join("-" if math.isnan(o) else f"{o:.2f}" for o in orders),
                "NA" if math.isnan(rep) else f"{rep:.2f}",
                ok,
            ]
        )
    return md_table(header, out) if out else "(no A0 data)\n"


# ---------------------------------------------------------------------------
# A1 small stiff ODEs
# ---------------------------------------------------------------------------


def a1_tables(rows, figdir, plt):
    header = ["problem", "configurations"] + [f"t / t(FIRKODE) at mescd={t:g}" for t in TARGETS]
    out = []
    for p in ["rober", "vdpol", "orego", "hires"]:
        g = [r for r in timed_rows(rows) if r.get("problem") == p]
        if not g:
            continue
        out.append([p, "best DIRK default"] + ratio_cells(g, is_dirk_default, TARGETS))
        out.append([p, "best DIRK of all"] + ratio_cells(g, is_dirk, TARGETS))
        tuned = [lb for lb in {r["label"] for r in g} if lb.startswith("ARKODE DIRK4 [")]
        for lb in tuned:
            out.append([p, lb] + ratio_cells(g, lambda x, lb=lb: x == lb, TARGETS))
    tables = {"A1_ratio": md_table(header, out) if out else "(no A1 data)\n"}
    tables["A1_counters"] = counters_table(
        rows,
        1e-6,
        [
            ("steps", col("nst")),
            ("attempts", col("nst_att")),
            ("Newton", col("nni")),
            ("setups", col("nsetups")),
            ("stages", col("stages")),
            ("Newton / step", per("nni", "nst")),
            ("us / step", time_per_step),
        ],
    )
    tables["A1_wp"] = wp_table(
        timed_rows(rows), [1e-4, 1e-6, 1e-8, 1e-10], extra_cols=("problem",)
    )
    if plt is not None:
        panels = []
        for p in ["rober", "vdpol", "orego", "hires"]:
            g = [r for r in timed_rows(rows) if r.get("problem") == p]
            if g:
                panels.append((p, g))
        if panels:
            fig_grid_wp(plt, panels, Path(figdir) / "arkode_small_wp.png")
    return tables


# ---------------------------------------------------------------------------
# A2 stage order
# ---------------------------------------------------------------------------


def a2_slope(rows):
    """Ratio of the best DIRK to FIRKODE at mescd 4, 6, 8 and the change of
    its logarithm from 4 to 8 (a ratio that grows with the accuracy)."""
    header = (
        ["problem", "params", "class"]
        + [f"t(best) / t(FIRKODE) at mescd={t:g}" for t in TARGETS]
        + ["log10 ratio(8) - log10 ratio(4)"]
    )
    groups = defaultdict(list)
    for r in timed_rows(rows):
        groups[(r.get("problem", ""), r.get("params", ""))].append(r)
    out = []
    for (p, prm), g in sorted(groups.items()):
        for name, pred in (("best DIRK default", is_dirk_default), ("best DIRK of all", is_dirk)):
            r4, _ = best_ratio(g, 4.0, pred)
            r8, _ = best_ratio(g, 8.0, pred)
            d = (
                fmt(math.log10(r8) - math.log10(r4), 3)
                if not (math.isnan(r4) or math.isnan(r8))
                else "NA"
            )
            out.append([p, prm, name] + ratio_cells(g, pred, TARGETS) + [d])
    return md_table(header, out) if out else "(no A2 adaptive data)\n"


def a2_tables(fixed, adaptive, figdir, plt):
    tables = {"A2_fixed": a0_orders(fixed) if fixed else "(no A2 fixed-step data)\n"}
    tables["A2_slope"] = a2_slope(adaptive)
    if plt is not None:
        if fixed:
            fig_fixed_step(plt, fixed, Path(figdir) / "arkode_order_fixed.png")
        grid_figure(
            plt,
            adaptive,
            ["n", "lambda", "theta", "rmax"],
            "",
            Path(figdir) / "arkode_stage_order_wp.png",
        )
    return tables


# ---------------------------------------------------------------------------
# A3 mass matrices, A4 tight tolerances
# ---------------------------------------------------------------------------


def a3_tables(rows, figdir, plt):
    tables = {"A3_mass": problem_class_ratio_table(rows)}
    tables["A3_counters"] = counters_table(
        rows,
        1e-6,
        [
            ("steps", col("nst")),
            ("attempts", col("nst_att")),
            ("mass solves", col("nmass_solves")),
            ("mass solves / attempt", per("nmass_solves", "nst_att")),
            ("Newton", col("nni")),
            ("stacked FGMRES / Newton", per("nsli", "nni")),
        ],
    )
    if plt is not None:
        grid_figure(plt, rows, ["n", "damping", "mass"], "", Path(figdir) / "arkode_mass_wp.png")
    return tables


def a4_tables(rows, figdir, plt):
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
        cells = [p, prm] + ratio_cells(g, is_dirk, (8.0, 10.0, 12.0))
        c = crossover_mescd(g, is_dirk)
        cells.append(fmt(c, 3) if not math.isnan(c) else "none")
        floors = [to_float(r.get("ref_disagreement")) for r in g]
        floors = [f for f in floors if not math.isnan(f) and f > 0]
        cells.append(
            f"{-math.log10(max(floors)):.1f}" if floors else "exact" if p == "sector" else "NA"
        )
        out.append(cells)
    tables = {"A4_crossover": md_table(header, out) if out else "(no A4 data)\n"}
    if plt is not None:
        panels = [(f"{p} {prm}".strip(), g) for (p, prm), g in sorted(groups.items())]
        if panels:
            fig_grid_wp(plt, panels, Path(figdir) / "arkode_tight_wp.png")
    return tables


# ---------------------------------------------------------------------------
# A5-A8 the other ARKODE steppers
# ---------------------------------------------------------------------------


def a5_tables(rows, figdir, plt):
    tables = {"A5_ratio": class_ratio_table(rows, ["n"])}
    tables["A5_counters"] = counters_table(
        rows,
        1e-6,
        [
            ("steps", col("nst")),
            ("RHS (implicit or all)", col("nfe")),
            ("RHS explicit", col("nfe_e")),
            ("Newton", col("nni")),
            ("Jac", col("nje")),
            ("setups", col("nsetups")),
            ("LSRK max stages", col("lsrk_max_stages")),
            ("power RHS / RHS", per("nfe_domeig", "nfe")),
        ],
    )
    crosses = []
    for (n,), g in sorted(groups_by_params(timed_rows(rows), ["n"]).items()):
        c = crossover_mescd(g, is_lsrk, (3.0, 4.0, 5.0, 6.0, 7.0, 8.0))
        crosses.append([n, fmt(c, 3) if not math.isnan(c) else "none"])
    tables["A5_lsrk_crossover"] = md_table(["n", "LSRK crossover mescd"], crosses)
    tables["A5_profile"] = profile_table(rows)
    if plt is not None:
        grid_figure(plt, rows, ["n"], "bruss3d", Path(figdir) / "arkode_bruss3d_wp.png")
        fig_profile(plt, rows, Path(figdir) / "arkode_bruss3d_profile.png")
    return tables


def a6_tables(rows, figdir, plt):
    heat = [r for r in rows if r.get("problem") == "heat2d"]
    sector = [r for r in rows if r.get("problem") == "sector"]
    tables = {
        "A6_map": class_ratio_table(heat, ["velocity"]),
        "A6_sector": class_ratio_table(sector, ["theta"]),
    }
    tables["A6_counters"] = counters_table(
        rows,
        1e-6,
        [
            ("steps", col("nst")),
            ("attempts", col("nst_att")),
            ("err fails", col("netf")),
            ("RHS (implicit or all)", col("nfe")),
            ("RHS explicit", col("nfe_e")),
            ("fast steps", col("nst_f")),
            ("LSRK max stages", col("lsrk_max_stages")),
        ],
    )
    if plt is not None:
        grid_figure(plt, heat, ["velocity"], "heat2d", Path(figdir) / "arkode_heat2d_wp.png")
        grid_figure(plt, sector, ["theta"], "SECTOR", Path(figdir) / "arkode_sector_lsrk_wp.png")
    return tables


def a7_tables(rows, figdir, plt):
    tables = {"A7_w": class_ratio_table(rows, ["w"])}
    tables["A7_counters"] = counters_table(
        rows,
        1e-6,
        [
            ("slow steps", col("nst")),
            ("fast steps", col("nst_f")),
            ("slow RHS (implicit)", col("nfe")),
            ("slow RHS (explicit)", col("nfe_e")),
            ("fast RHS", col("nfe_f")),
            ("fast err fails", col("netf_f")),
        ],
    )
    if plt is not None:
        grid_figure(plt, rows, ["w"], "kprmri", Path(figdir) / "arkode_kprmri_wp.png")
    return tables


def a8_tables(rows, figdir, plt):
    tables = {"A8_grid": class_ratio_table(rows, ["mx"])}
    tables["A8_counters"] = counters_table(
        rows,
        1e-6,
        [
            ("steps", col("nst")),
            ("Newton", col("nni")),
            ("Krylov / Newton (ARKODE)", per("nli", "nni")),
            ("stacked FGMRES / Newton (FIRKODE)", per("nsli", "nni")),
            ("prec solves", col("nps")),
        ],
    )
    if plt is not None:
        grid_figure(plt, rows, ["mx"], "diurnal", Path(figdir) / "arkode_diurnal_wp.png")
    return tables


def arkode_tables(data: dict, results: Path, figdir: Path, plt) -> dict:
    tables = {}
    if data.get("A0"):
        tables["A0_orders"] = a0_orders(data["A0"])
        if plt is not None:
            fig_fixed_step(plt, data["A0"], Path(figdir) / "arkode_order_selection.png")
    if data.get("A1"):
        tables.update(a1_tables(data["A1"], figdir, plt))
    if data.get("A2_fixed") or data.get("A2_adaptive"):
        tables.update(
            a2_tables(data.get("A2_fixed", []), data.get("A2_adaptive", []), figdir, plt)
        )
    if data.get("A3"):
        tables.update(a3_tables(data["A3"], figdir, plt))
    if data.get("A4"):
        tables.update(a4_tables(data["A4"], figdir, plt))
    if data.get("A5"):
        tables.update(a5_tables(data["A5"], figdir, plt))
    if data.get("A6"):
        tables.update(a6_tables(data["A6"], figdir, plt))
    if data.get("A7"):
        tables.update(a7_tables(data["A7"], figdir, plt))
    if data.get("A8"):
        tables.update(a8_tables(data["A8"], figdir, plt))
    refs = {}
    for name in ("A1", "A2_adaptive", "A3", "A4", "A5", "A6", "A8"):
        for r in data.get(name, []):
            d = to_float(r.get("ref_disagreement"))
            if not math.isnan(d):
                refs[(r["problem"], r.get("params", ""))] = (d, r.get("nst_ref"))
    tables["references"] = md_table(
        ["problem", "params", "mixed disagreement (reference check)", "reference steps"],
        [[p, prm, f"{d:.1e}", n] for (p, prm), (d, n) in sorted(refs.items())],
    )
    return tables
