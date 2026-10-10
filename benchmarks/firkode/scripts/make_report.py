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
# Report generator for the FIRKODE benchmark suite: reads the CSV files
# written by run_benchmarks.py, derives the quantities quoted in
# benchmarks/firkode/BENCHMARK_BDF.md (--report bdf, the default) or
# BENCHMARK_ARKODE.md (--report arkode), and writes Markdown table fragments
# and PNG figures. The fragments are injected into (or checked against) the
# <!-- begin generated: NAME --> ... <!-- end generated --> blocks of the
# report so that every number in it is traceable to a CSV row.
#
#   make_report.py --results benchmarks/firkode/results/<run-id> \
#       --figdir doc/shared/figs/firkode/benchmarks/bdf \
#       [--inject benchmarks/firkode/BENCHMARK_BDF.md | --check benchmarks/firkode/BENCHMARK_BDF.md]
#   make_report.py --report arkode --results benchmarks/firkode/results/<run-id> \
#       --figdir doc/shared/figs/firkode/benchmarks/arkode \
#       [--inject benchmarks/firkode/BENCHMARK_ARKODE.md]
# -----------------------------------------------------------------
"""Report generator for the FIRKODE benchmark suite."""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from benchmark_common import read_csv, style_for, to_float, to_int

EPS = 2.220446049250313e-16
ORDER_WINDOW = (1.0e4 * EPS, 1.0e-2)  # admissible errors for order estimates

LABEL_ORDER = [
    "FIRKODE s=1",
    "FIRKODE s=2",
    "FIRKODE s=3",
    "FIRKODE s=4",
    "FIRKODE s=5",
    "FIRKODE s=6",
    "FIRKODE s=7",
    "FIRKODE s=8",
    "FIRKODE s=9",
    "CVODE BDF+STALD",
    "CVODE BDF q<=4",
    "CVODE BDF q<=3",
    "CVODE BDF q<=2",
    "CVODE BDF msbp=50",
    "CVODE BDF h0=1e-6",
    "CVODE BDF",
    "ARKODE DIRK2",
    "ARKODE DIRK3",
    "ARKODE DIRK4",
    "ARKODE DIRK5",
    "ARKODE DIRK SDIRK_5_3_4",
    "ARKODE DIRK ESDIRK437L2SA_7_3_4",
    "ARKODE DIRK KVAERNO_7_4_5",
    "ARKODE DIRK ARK548L2SAb_DIRK_8_4_5",
    "ARKODE IMEX3",
    "ARKODE IMEX4",
    "ARKODE IMEX5",
    "ARKODE LSRK RKC2",
    "ARKODE LSRK RKL2",
    "ARKODE MRI ESDIRK34a",
    "ARKODE MRI IMEX-SR32",
    "IDA BDF q<=2",
    "IDA BDF",
]


def label_sort_key(label: str):
    for i, k in enumerate(LABEL_ORDER):
        if label.startswith(k):
            return (i, label)
    return (len(LABEL_ORDER), label)


def fmt(v, digits=3) -> str:
    if v is None:
        return "NA"
    if isinstance(v, str):
        return v
    if isinstance(v, int):
        return str(v)
    if math.isnan(v):
        return "NA"
    if abs(v) >= 1e4 or (abs(v) < 1e-3 and v != 0):
        return f"{v:.{digits - 1}e}"
    return f"{v:.{digits}g}"


def fmt_ms(seconds) -> str:
    if seconds is None or math.isnan(seconds):
        return "NA"
    if seconds < 1e-3:
        return "<1"
    return f"{1e3 * seconds:.3g}"


def mescd(row) -> float:
    e = to_float(row.get("err_mixed"))
    return -math.log10(e) if e > 0 else math.nan


def md_table(header: list[str], rows: list[list], align=None) -> str:
    out = ["| " + " | ".join(header) + " |"]
    align = align or ["---"] * len(header)
    out.append("| " + " | ".join(align) + " |")
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out) + "\n"


def status_of(row) -> str:
    flag = row.get("driver_flag", "")
    if flag == "ok":
        return "ok"
    if flag == "TIMEOUT":
        return "timeout"
    if flag == "FAST_TIMEOUT":
        return "skipped (fast)"
    if flag == "UNSUPPORTED":
        return "n/a"
    return f"fail ({row.get('flag', '?')})"


def ok_rows(rows):
    return [r for r in rows if r.get("driver_flag") == "ok"]


def by_label(rows) -> dict:
    out = defaultdict(list)
    for r in rows:
        out[r["label"]].append(r)
    return out


def at_floor(row) -> bool:
    d = to_float(row.get("ref_disagreement"))
    e = to_float(row.get("err_mixed"))
    return not math.isnan(d) and d > 0 and e < 10.0 * d


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def setup_matplotlib():
    os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "mpl"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "legend.fontsize": 8}
    )
    return plt


def plot_wp_axis(ax, rows, xkey, xlabel, time_only_reliable=True):
    for label, group in sorted(
        by_label(ok_rows(rows)).items(), key=lambda kv: label_sort_key(kv[0])
    ):
        pts = []
        for r in group:
            if time_only_reliable and xkey.startswith("time") and r.get("timing_reliable") != "1":
                continue
            x = to_float(r.get(xkey))
            e = to_float(r.get("err_mixed"))
            if x > 0 and e > 0:
                pts.append((x, e, at_floor(r)))
        if not pts:
            continue
        pts.sort(key=lambda p: -p[1])
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
            label=label,
        )
        floor_pts = [p for p in pts if p[2]]
        if floor_pts:
            ax.plot(
                [p[0] for p in floor_pts],
                [p[1] for p in floor_pts],
                ls="none",
                marker=st["marker"],
                mfc="none",
                mec=st["color"],
                ms=8,
            )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("mixed error at $t_f$")


def fig_work_precision(plt, rows, title, path):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
    plot_wp_axis(axes[0], rows, "time_median", "wall time (s, median of reps)")
    plot_wp_axis(axes[1], rows, "nfe", "RHS evaluations")
    axes[0].set_title(title)
    axes[1].legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_grid_wp(
    plt, groups: list[tuple[str, list]], path, xkey="time_median", xlabel="wall time (s)"
):
    n = len(groups)
    cols = 2 if n > 1 else 1
    nrows = math.ceil(n / cols)
    fig, axes = plt.subplots(nrows, cols, figsize=(4.6 * cols, 3.6 * nrows), squeeze=False)
    for ax, (title, rows) in zip(axes.flat, groups):
        plot_wp_axis(ax, rows, xkey, xlabel)
        ax.set_title(title)
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc="lower center",
            ncol=min(5, len(labels)),
            bbox_to_anchor=(0.5, -0.02),
        )
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# B1 / B5 work-precision tables
# ---------------------------------------------------------------------------


def interp_log_time(group, target_mescd):
    pts = sorted(
        (mescd(r), to_float(r.get("time_median")))
        for r in group
        if r.get("timing_reliable") == "1" and not math.isnan(mescd(r))
    )
    pts = [(m, t) for m, t in pts if t > 0]
    if len(pts) < 2:
        return math.nan
    # monotone envelope: cheapest time achieving at least this accuracy
    for i in range(len(pts)):
        m, t = pts[i]
        if m >= target_mescd:
            if i == 0:
                return math.nan
            m0, t0 = pts[i - 1]
            if m == m0:
                return t
            w = (target_mescd - m0) / (m - m0)
            return 10 ** ((1 - w) * math.log10(t0) + w * math.log10(t))
    return math.nan


def wp_table(rows, rtols, with_stages=False, extra_cols=()):
    header = [
        "method",
        "rtol",
        "mescd",
        "steps",
        "attempts",
        "Newton",
        "Jac",
        "setups",
        "err fails",
        "NLS fails",
    ]
    header += list(extra_cols) + ["time (ms)", "status"]
    out = []
    for label, group in sorted(by_label(rows).items(), key=lambda kv: label_sort_key(kv[0])):
        for r in sorted(group, key=lambda r: -to_float(r.get("rtol"))):
            if rtols is not None and not any(
                abs(to_float(r.get("rtol")) - t) <= 1e-3 * t for t in rtols
            ):
                continue
            m = mescd(r)
            cells = [
                label,
                f"{to_float(r.get('rtol')):.0e}",
                (fmt(m, 3) + ("*" if at_floor(r) else "")) if not math.isnan(m) else "NA",
                r.get("nst", "NA"),
                r.get("nst_att", "NA"),
                r.get("nni", "NA"),
                r.get("nje", "NA"),
                r.get("nsetups", "NA"),
                r.get("netf", "NA"),
                r.get("nncf", "NA"),
            ]
            for c in extra_cols:
                cells.append(r.get(c, "NA"))
            cells.append(
                fmt_ms(to_float(r.get("time_median")))
                if r.get("timing_reliable") == "1"
                else "<1"
                if r.get("driver_flag") == "ok"
                else "NA"
            )
            cells.append(status_of(r))
            out.append(cells)
    return md_table(header, out)


def speedup_table(rows, targets=(4.0, 6.0, 8.0), base="FIRKODE s=3"):
    groups = by_label(ok_rows(rows))
    if base not in groups:
        return "(no FIRKODE s=3 data)\n"
    header = ["method"] + [f"t / t(FIRKODE s=3) at mescd={t:g}" for t in targets]
    out = []
    base_t = [interp_log_time(groups[base], t) for t in targets]
    for label, group in sorted(groups.items(), key=lambda kv: label_sort_key(kv[0])):
        cells = [label]
        for t, bt in zip(targets, base_t):
            v = interp_log_time(group, t)
            cells.append(fmt(v / bt, 3) if not (math.isnan(v) or math.isnan(bt)) else "NA")
        out.append(cells)
    return md_table(header, out)


# ---------------------------------------------------------------------------
# B3 Krylov iterations vs stages
# ---------------------------------------------------------------------------


def klin_table(rows, log_rows):
    log_by_key = {r["cell_key"]: r for r in log_rows}
    header = [
        "problem",
        "params",
        "rtol",
        "s",
        "k_lin mean",
        "k_lin p95",
        "k_lin max",
        "6s cap",
        "restarted",
        "stacked fails",
        "block solves / Newton",
        "steps",
        "Newton",
        "RHS",
        "time (ms)",
        "status",
    ]
    out = []
    for r in sorted(
        rows,
        key=lambda r: (
            r["problem"],
            r.get("params", ""),
            to_float(r.get("rtol")),
            to_int(r.get("stages")) or 0,
        ),
    ):
        s = to_int(r.get("stages")) or 0
        nni = to_float(r.get("nni"))
        nsli = to_float(r.get("nsli"))
        nbs = to_float(r.get("nbs"))
        lg = log_by_key.get(r["cell_key"], {})
        out.append(
            [
                r["problem"],
                r.get("params", ""),
                f"{to_float(r.get('rtol')):.0e}",
                s,
                fmt(nsli / nni, 3) if nni > 0 else "NA",
                lg.get("klin_p95", "NA"),
                lg.get("klin_max", "NA"),
                6 * s if s > 1 else "-",
                lg.get("n_restarted", "NA"),
                r.get("nscf", "NA"),
                fmt(nbs / nni, 3) if nni > 0 else "NA",
                r.get("nst", "NA"),
                r.get("nni", "NA"),
                r.get("nfe", "NA"),
                fmt_ms(to_float(r.get("time_median"))) if r.get("driver_flag") == "ok" else "NA",
                status_of(r),
            ]
        )
    return md_table(header, out)


def fig_klin(plt, rows, log_rows, path):
    log_by_key = {r["cell_key"]: r for r in log_rows}
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    groups = defaultdict(list)
    for r in ok_rows(rows):
        if abs(to_float(r.get("rtol")) - 1e-6) < 1e-9:
            groups[(r["problem"], r.get("params", ""))].append(r)
    markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
    for i, (key, group) in enumerate(sorted(groups.items())):
        group.sort(key=lambda r: to_int(r.get("stages")) or 0)
        ss = [to_int(r.get("stages")) for r in group]
        nni = [to_float(r.get("nni")) for r in group]
        klin = [to_float(r.get("nsli")) / n if n > 0 else math.nan for r, n in zip(group, nni)]
        name = key[0] + (f" ({key[1]})" if key[1] else "")
        mk = markers[i % len(markers)]
        axes[0].plot(ss, klin, marker=mk, label=name)
        kmax = [to_float(log_by_key.get(r["cell_key"], {}).get("klin_max", "nan")) for r in group]
        if any(not math.isnan(k) for k in kmax):
            axes[1].plot(ss, kmax, marker=mk, label=name)
        axes[2].plot(
            ss,
            [to_float(r.get("nbs")) / n if n > 0 else math.nan for r, n in zip(group, nni)],
            marker=mk,
            label=name,
        )
    s_axis = list(range(1, 10))
    axes[0].plot(s_axis, [6 * s if s > 1 else 0 for s in s_axis], "k:", lw=1, label="6s cap")
    axes[1].plot(s_axis, [6 * s if s > 1 else 0 for s in s_axis], "k:", lw=1, label="6s cap")
    axes[0].set_ylabel("mean FGMRES iterations per Newton iteration")
    axes[1].set_ylabel("max FGMRES iterations per solve (log build)")
    axes[2].set_ylabel("block solves per Newton iteration")
    for ax in axes:
        ax.set_xlabel("stages s")
        ax.set_xticks(s_axis)
    axes[0].legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# B4 convergence order
# ---------------------------------------------------------------------------


def observed_orders(hs, errs):
    """Orders between consecutive admissible refinements and the reported one."""
    pairs = sorted(zip(hs, errs), key=lambda p: -p[0])
    orders = []
    for (h0, e0), (h1, e1) in zip(pairs, pairs[1:]):
        if (
            ORDER_WINDOW[0] < e0 < ORDER_WINDOW[1]
            and ORDER_WINDOW[0] < e1 < ORDER_WINDOW[1]
            and h0 != h1
        ):
            orders.append(math.log(e0 / e1) / math.log(h0 / h1))
        else:
            orders.append(math.nan)
    valid = [o for o in orders if not math.isnan(o)]
    return orders, (valid[-1] if valid else math.nan)


def fixed_step_table(rows):
    header = [
        "problem",
        "params",
        "method",
        "expected",
        "errors at h0, h0/2, ...",
        "observed orders",
        "reported",
    ]
    out = []
    groups = defaultdict(list)
    for r in ok_rows(rows):
        groups[(r["problem"], r.get("params", ""), r["label"])].append(r)
    for (p, prm, label), group in sorted(
        groups.items(), key=lambda kv: (kv[0][0], kv[0][1], label_sort_key(kv[0][2]))
    ):
        group.sort(key=lambda r: -to_float(r.get("fixed_h")))
        hs = [to_float(r.get("fixed_h")) for r in group]
        errs = [to_float(r.get("err_abs_max")) for r in group]
        orders, rep = observed_orders(hs, errs)
        expected = to_int(group[0].get("order"))
        out.append(
            [
                p,
                prm,
                label,
                expected,
                ", ".join(f"{e:.2e}" for e in errs),
                ", ".join("-" if math.isnan(o) else f"{o:.2f}" for o in orders),
                "NA" if math.isnan(rep) else f"{rep:.2f}",
            ]
        )
    return md_table(header, out)


def fig_fixed_step(plt, rows, path):
    groups = defaultdict(list)
    for r in ok_rows(rows):
        groups[(r["problem"], r.get("params", ""))].append(r)
    keys = sorted(groups)
    cols = 2
    nrows = math.ceil(len(keys) / cols)
    fig, axes = plt.subplots(nrows, cols, figsize=(9.5, 3.8 * nrows), squeeze=False)
    for ax, key in zip(axes.flat, keys):
        for label, g in sorted(
            by_label(groups[key]).items(), key=lambda kv: label_sort_key(kv[0])
        ):
            g.sort(key=lambda r: to_float(r.get("fixed_h")))
            st = style_for(label)
            ax.plot(
                [to_float(r.get("fixed_h")) for r in g],
                [to_float(r.get("err_abs_max")) for r in g],
                color=st["color"],
                ls=st["ls"],
                marker=st["marker"],
                mfc=st["mfc"],
                ms=5,
                label=label,
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("step size h")
        ax.set_ylabel("max error at $t_f$")
        ax.set_title(key[0] + (f" ({key[1]})" if key[1] else ""))
        ax.axhline(ORDER_WINDOW[0], color="grey", lw=0.8, ls="--")
    for ax in list(axes.flat)[len(keys) :]:
        ax.axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="lower center", ncol=min(5, len(labels)), bbox_to_anchor=(0.5, -0.02)
    )
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def lsq_slope(xs, ys):
    pts = [(math.log(x), math.log(y)) for x, y in zip(xs, ys) if x > 0 and y > 0]
    if len(pts) < 2:
        return math.nan, 0
    n = len(pts)
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    sxy = sum((p[0] - mx) * (p[1] - my) for p in pts)
    return (sxy / sxx if sxx > 0 else math.nan), n


def estimate_table(rows):
    header = [
        "problem",
        "params",
        "method",
        "s",
        "steps accepted at h",
        "slope of ||ele||_inf vs h",
        "expected s+1",
        "slope of true one-step error",
        "expected 2s",
    ]
    out = []
    groups = defaultdict(list)
    for r in ok_rows(rows):
        groups[(r["problem"], r.get("params", ""), r["label"])].append(r)
    for (p, prm, label), group in sorted(
        groups.items(), key=lambda kv: (kv[0][0], kv[0][1], label_sort_key(kv[0][2]))
    ):
        # keep the cells where exactly one step of the prescribed size was taken
        good = []
        for r in group:
            h = None
            args = r["cell_key"].split("|")[-1].split()
            if "--h" in args:
                h = to_float(args[args.index("--h") + 1])
            hl = to_float(r.get("hlast"))
            if (
                h
                and to_int(r.get("nst")) == 1
                and to_int(r.get("nst_att")) in (1, None)
                and abs(hl - h) <= 1e-12 * h
            ):
                good.append((h, r))
        if not good:
            continue
        good.sort(key=lambda p: p[0])
        hs = [h for h, _ in good]
        ele = [to_float(r.get("ele_max")) for _, r in good]
        true_err = [to_float(r.get("err_abs_max")) for _, r in good]
        adm = [(h, e) for h, e in zip(hs, ele) if ORDER_WINDOW[0] < e < ORDER_WINDOW[1]]
        admt = [(h, e) for h, e in zip(hs, true_err) if ORDER_WINDOW[0] < e < ORDER_WINDOW[1]]
        s_ele, n_ele = lsq_slope([a[0] for a in adm], [a[1] for a in adm])
        s_true, n_true = lsq_slope([a[0] for a in admt], [a[1] for a in admt])
        s = to_int(group[0].get("stages")) or 0
        is_firk = label.startswith("FIRKODE")
        out.append(
            [
                p,
                prm,
                label,
                s if is_firk else "-",
                f"{len(good)}/{len(group)}",
                f"{s_ele:.2f} ({n_ele} pts)" if not math.isnan(s_ele) else "NA",
                s + 1 if is_firk else "-",
                f"{s_true:.2f} ({n_true} pts)" if not math.isnan(s_true) else "NA",
                2 * s if is_firk else "-",
            ]
        )
    return md_table(header, out)


def fig_estimate(plt, rows, path):
    groups = defaultdict(list)
    for r in ok_rows(rows):
        if to_int(r.get("nst")) == 1:
            groups[(r["problem"], r.get("params", ""))].append(r)
    keys = sorted(groups)
    fig, axes = plt.subplots(
        1, max(1, len(keys)), figsize=(4.6 * max(1, len(keys)), 3.8), squeeze=False
    )
    for ax, key in zip(axes.flat, keys):
        for label, g in sorted(
            by_label(groups[key]).items(), key=lambda kv: label_sort_key(kv[0])
        ):
            g.sort(key=lambda r: to_float(r.get("hlast")))
            st = style_for(label)
            ax.plot(
                [to_float(r.get("hlast")) for r in g],
                [to_float(r.get("ele_max")) for r in g],
                color=st["color"],
                ls=st["ls"],
                marker=st["marker"],
                mfc=st["mfc"],
                ms=4,
                label=label,
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("step size h")
        ax.set_ylabel("$||$ele$||_\\infty$ after one step")
        ax.set_title(key[0] + (f" ({key[1]})" if key[1] else ""))
        ax.axhline(ORDER_WINDOW[0], color="grey", lw=0.8, ls="--")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="lower center", ncol=min(6, len(labels)), bbox_to_anchor=(0.5, -0.05)
    )
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# B6 ablation, B7 roots and dense output, sentinel, validation
# ---------------------------------------------------------------------------


def ablation_table(rows):
    header = [
        "problem",
        "rtol",
        "variant",
        "mescd",
        "steps",
        "RHS",
        "Newton",
        "FGMRES iters",
        "Jac",
        "setups",
        "time / default",
    ]
    out = []
    groups = defaultdict(dict)
    for r in ok_rows(rows):
        groups[(r["problem"], r.get("params", ""), to_float(r.get("rtol")))][r["label"]] = r
    for (p, prm, rtol), variants in sorted(
        groups.items(), key=lambda kv: (kv[0][0], kv[0][1], -kv[0][2])
    ):
        base = variants.get("FIRKODE s=3 [default]")
        bt = to_float(base.get("time_median")) if base else math.nan
        for label, r in sorted(
            variants.items(), key=lambda kv: (kv[0] != "FIRKODE s=3 [default]", kv[0])
        ):
            t = to_float(r.get("time_median"))
            out.append(
                [
                    p + (f" ({prm})" if prm else ""),
                    f"{rtol:.0e}",
                    label.replace("FIRKODE s=3 ", ""),
                    fmt(mescd(r), 3),
                    r.get("nst"),
                    r.get("nfe"),
                    r.get("nni"),
                    r.get("nsli"),
                    r.get("nje"),
                    r.get("nsetups"),
                    fmt(t / bt, 3) if bt > 0 and r.get("timing_reliable") == "1" else "NA",
                ]
            )
    return md_table(header, out)


def parse_roots(text: str):
    out = []
    if not text or text == "NA":
        return out
    for item in text.split(";"):
        t, _, dirs = item.partition(":")
        out.append((float(t), dirs))
    return out


def roots_table(rows):
    rows = ok_rows(rows)
    ref = [r for r in rows if r["label"].startswith("reference")]
    ref_roots = parse_roots(ref[0].get("root_times", "")) if ref else []
    header = ["method", "rtol", "steps", "g evals", "roots found"] + [
        f"|t_root{k + 1} - t_ref| / t_ref" for k in range(len(ref_roots))
    ]
    out = []
    if ref:
        out.append(
            [
                "reference FIRKODE s=5",
                f"{to_float(ref[0].get('rtol')):.0e}",
                ref[0].get("nst"),
                ref[0].get("ngeval"),
                "; ".join(f"{t:.10e} ({d})" for t, d in ref_roots),
            ]
            + ["-"] * len(ref_roots)
        )
    for label, group in sorted(
        by_label([r for r in rows if not r["label"].startswith("reference")]).items(),
        key=lambda kv: label_sort_key(kv[0]),
    ):
        for r in sorted(group, key=lambda r: -to_float(r.get("rtol"))):
            roots = parse_roots(r.get("root_times", ""))
            cells = [
                label,
                f"{to_float(r.get('rtol')):.0e}",
                r.get("nst"),
                r.get("ngeval"),
                "; ".join(f"{t:.6e} ({d})" for t, d in roots),
            ]
            for k, (tr, _) in enumerate(ref_roots):
                cells.append(f"{abs(roots[k][0] - tr) / tr:.1e}" if k < len(roots) else "missed")
            out.append(cells)
    return md_table(header, out)


def dense_table(rows):
    header = [
        "method",
        "rtol",
        "steps",
        "max mixed error over 1000 samples",
        "mixed error at t_f",
        "ratio",
    ]
    out = []
    for label, group in sorted(
        by_label(ok_rows(rows)).items(), key=lambda kv: label_sort_key(kv[0])
    ):
        for r in sorted(group, key=lambda r: -to_float(r.get("rtol"))):
            d = to_float(r.get("dense_err_mixed"))
            e = to_float(r.get("err_mixed"))
            out.append(
                [
                    label,
                    f"{to_float(r.get('rtol')):.0e}",
                    r.get("nst"),
                    fmt(d, 3),
                    fmt(e, 3),
                    fmt(d / e, 3) if e > 0 else "NA",
                ]
            )
    return md_table(header, out)


def fig_roots_dense(plt, roots_rows, dense_rows, path):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
    rows = ok_rows(roots_rows)
    ref = [r for r in rows if r["label"].startswith("reference")]
    ref_roots = parse_roots(ref[0].get("root_times", "")) if ref else []
    for label, group in sorted(
        by_label([r for r in rows if not r["label"].startswith("reference")]).items(),
        key=lambda kv: label_sort_key(kv[0]),
    ):
        group.sort(key=lambda r: -to_float(r.get("rtol")))
        xs, ys = [], []
        for r in group:
            roots = parse_roots(r.get("root_times", ""))
            if ref_roots and len(roots) >= len(ref_roots):
                xs.append(to_float(r.get("rtol")))
                ys.append(max(abs(roots[k][0] - tr) / tr for k, (tr, _) in enumerate(ref_roots)))
        st = style_for(label)
        axes[0].plot(
            xs, ys, color=st["color"], ls=st["ls"], marker=st["marker"], mfc=st["mfc"], label=label
        )
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].invert_xaxis()
    axes[0].set_xlabel("rtol")
    axes[0].set_ylabel("max relative root-time error")
    axes[0].set_title("Robertson rootfinding")
    axes[0].legend()
    for label, group in sorted(
        by_label(ok_rows(dense_rows)).items(), key=lambda kv: label_sort_key(kv[0])
    ):
        st = style_for(label)
        axes[1].plot(
            [to_float(r.get("err_mixed")) for r in group],
            [to_float(r.get("dense_err_mixed")) for r in group],
            color=st["color"],
            ls="none",
            marker=st["marker"],
            mfc=st["mfc"],
            ms=7,
            label=label,
        )
    lim = [1e-10, 1e-2]
    axes[1].plot(lim, lim, "k:", lw=0.8)
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("mixed error at $t_f$")
    axes[1].set_ylabel("max mixed error of dense output")
    axes[1].set_title("Van der Pol dense output (1000 samples)")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def sentinel_table(rows):
    header = ["family", "phase", "time_min (ms)", "time_median (ms)", "timestamp (UTC)"]
    out = [
        [
            r["family"],
            r["label"],
            fmt_ms(to_float(r.get("time_min"))),
            fmt_ms(to_float(r.get("time_median"))),
            r.get("timestamp", ""),
        ]
        for r in rows
    ]
    vals = [to_float(r.get("time_median")) for r in rows if to_float(r.get("time_median")) > 0]
    spread = (max(vals) - min(vals)) / min(vals) if vals else math.nan
    return (
        md_table(header, out)
        + f"\nSentinel spread (max - min) / min of the median: {fmt(100 * spread, 3)} %.\n"
    )


def validation_table(rows):
    # exact: every counter equal to the answer file; "answer" in the tolerance
    # column is the answer-file criterion of test/testRunner (run_benchmarks.py)
    header = ["replica", "passed", "exact", "tolerance", "mismatch", "note"]
    out = [
        [
            r["label"],
            "yes" if r.get("passed") == "1" else "NO",
            {"1": "yes", "0": "no"}.get(r.get("exact", ""), "NA"),
            r.get("tolerance", ""),
            r.get("mismatch", "") or "-",
            r.get("note", ""),
        ]
        for r in rows
    ]
    return md_table(header, out)


def known(value) -> bool:
    """A provenance value that was actually recorded (collect_provenance
    writes "unavailable" where a lookup failed, e.g. sysctl on Linux)."""
    return value is not None and str(value).strip() not in ("", "unavailable", "NA", "None")


def provenance_os(prov: dict) -> str:
    """The operating system of a run: recorded by newer drivers, otherwise
    derived from platform.platform() ("macOS-15.1-arm64-arm-64bit",
    "Linux-6.8.0-45-generic-x86_64-with-glibc2.39")."""
    if known(prov.get("os")):
        return {"Darwin": "macOS"}.get(prov["os"], prov["os"])
    plat = str(prov.get("platform", ""))
    if plat.startswith(("macOS", "Darwin")):
        return "macOS"
    if plat.startswith("Linux"):
        return "Linux"
    return plat or "NA"


def provenance_machine(prov: dict) -> str:
    if known(prov.get("machine_note")):
        return prov["machine_note"]
    parts = []
    cpu = prov.get("cpu_brand") if known(prov.get("cpu_brand")) else prov.get("machine")
    if known(cpu):
        perf = prov.get("perf_levels") or {}
        p, e = perf.get("performance_cores"), perf.get("efficiency_cores")
        parts.append(f"{cpu} ({p}P+{e}E cores)" if known(p) and known(e) else str(cpu))
    mem = prov.get("memory_bytes")
    if known(mem) and str(mem).isdigit():
        parts.append(f"{int(mem) // 2**30} GiB")
    if known(prov.get("cpu_count")):
        parts.append(f"{prov['cpu_count']} logical CPUs")
    return ", ".join(parts) if parts else "NA"


def provenance_table(prov: dict):
    cache = prov.get("cmake_cache", {})
    git = prov.get("git", {})
    os_name = provenance_os(prov)
    plat = prov.get("platform", "NA")
    version = prov.get("os_version")
    if os_name == "macOS":
        os_row = f"macOS {version} ({plat})" if known(version) else f"macOS ({plat})"
    else:
        # Linux: the distribution name (newer drivers) and the platform string,
        # which carries the kernel release and the C library
        os_row = f"{version} ({plat})" if known(version) else str(plat)
    compiler = prov.get("compiler") if known(prov.get("compiler")) else "NA"
    if known(cache.get("CMAKE_CXX_COMPILER")):
        compiler = f"{compiler}; CMake C++ compiler `{cache['CMAKE_CXX_COMPILER']}`"
    threads = ", ".join(f"{k}={v}" for k, v in (prov.get("env") or {}).items() if known(v))
    items = [
        ("Machine", provenance_machine(prov)),
        ("OS", os_row),
        ("Compiler", compiler),
        ("CMake / generator", f"{prov.get('cmake')} / {cache.get('CMAKE_GENERATOR')}"),
        (
            "Build type and flags",
            f"{cache.get('CMAKE_BUILD_TYPE')}: C `{cache.get('CMAKE_C_FLAGS_RELEASE')}`, C++ `{cache.get('CMAKE_CXX_FLAGS_RELEASE')}`",
        ),
        (
            "Precision / index size",
            f"{cache.get('SUNDIALS_PRECISION')} / {cache.get('SUNDIALS_INDEX_SIZE')}",
        ),
        (
            "Logging level / error checks / profiling",
            f"{cache.get('SUNDIALS_LOGGING_LEVEL')} / {cache.get('SUNDIALS_ENABLE_ERROR_CHECKS')} / {cache.get('SUNDIALS_ENABLE_PROFILING')}",
        ),
        (
            "SUNDIALS",
            prov.get("benchmark_version_output", "").splitlines()[0]
            if prov.get("benchmark_version_output")
            else "NA",
        ),
        ("Git commit", f"{git.get('describe')} on {git.get('branch')}"),
        (
            "Python / numpy / matplotlib",
            f"{prov.get('python')} / {prov.get('numpy')} / {prov.get('matplotlib')}",
        ),
        ("Threads", threads or "NA"),
        ("Run started (UTC)", prov.get("started_utc")),
        ("Binary SHA-256", prov.get("binary_sha256", "")),
    ]
    return md_table(["item", "value"], [[k, v] for k, v in items])


# ---------------------------------------------------------------------------
# Inject / check
# ---------------------------------------------------------------------------

BLOCK_RE = re.compile(
    r"(<!-- begin generated: ([A-Za-z0-9_.-]+) -->\n)(.*?)(<!-- end generated -->)", re.DOTALL
)


def inject(report: Path, tables: dict, check_only: bool) -> int:
    text = report.read_text()
    missing, mismatched = [], []

    def repl(m):
        name = m.group(2)
        if name not in tables:
            missing.append(name)
            return m.group(0)
        new = tables[name].rstrip("\n") + "\n"
        if m.group(3).rstrip("\n") + "\n" != new:
            mismatched.append(name)
        return m.group(1) + new + m.group(4)

    new_text = BLOCK_RE.sub(repl, text)
    if missing:
        print("blocks without a generated table:", ", ".join(missing))
    if check_only:
        if mismatched:
            print("stale generated blocks:", ", ".join(mismatched))
            return 1
        print("all generated blocks are current")
        return 0
    report.write_text(new_text)
    print(f"updated {len(mismatched)} block(s) in {report}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--inject", default=None, help="BENCHMARK.md to update in place")
    ap.add_argument("--check", default=None, help="BENCHMARK.md to verify")
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument(
        "--report",
        choices=("bdf", "arkode"),
        default="bdf",
        help="the report whose tables and figures to build (BENCHMARK_BDF.md or BENCHMARK_ARKODE.md)",
    )
    args = ap.parse_args()

    results = Path(args.results)
    figdir = Path(args.figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    tables_dir = results / "tables"
    tables_dir.mkdir(exist_ok=True)

    data = {
        name: read_csv(results / f"{name}.csv")
        for name in [
            "validation",
            "sentinel",
            "B1",
            "B2",
            "B3",
            "B3_klin_log",
            "B4_fixed",
            "B4_estimate",
            "B5",
            "B6",
            "B7_roots",
            "B7_dense",
            "B8",
            "B9",
            "B10",
            "B10_spectrum",
            "B11",
            "B12",
            "B13",
            "A0",
            "A1",
            "A2_fixed",
            "A2_adaptive",
            "A3",
            "A4",
            "A5",
            "A6",
            "A7",
            "A8",
        ]
    }
    prov = (
        json.loads((results / "provenance.json").read_text())
        if (results / "provenance.json").exists()
        else {}
    )

    tables = {}
    tables["provenance"] = provenance_table(prov)
    shas = {prov.get("binary_sha256")} | {
        run.get("binary_sha256") for run in prov.get("runs", []) if isinstance(run, dict)
    }
    shas.discard(None)
    if len(shas) > 1:
        print(f"warning: the results come from {len(shas)} different binaries (provenance.json)")
    if args.report == "arkode":
        return arkode_report(args, data, results, figdir, tables_dir, tables)
    plt_bdf = None if args.no_figures else setup_matplotlib()
    from report_bdf import bdf_tables

    tables.update(bdf_tables(data, results, figdir, plt_bdf))
    tables["validation"] = validation_table(data["validation"])
    tables["sentinel"] = sentinel_table(data["sentinel"])

    # B1
    b1 = data["B1"]
    for p in ["rober", "vdpol", "orego", "hires"]:
        rows = [r for r in b1 if r["problem"] == p]
        if rows:
            tables[f"B1_{p}"] = wp_table(rows, [1e-4, 1e-6, 1e-8, 1e-10])
            tables[f"B1_{p}_speedup"] = speedup_table(rows)
    refs = {}
    for r in b1 + data["B2"] + data["B5"]:
        d = to_float(r.get("ref_disagreement"))
        if not math.isnan(d):
            refs[(r["problem"], r.get("params", ""))] = (d, r.get("nst_ref"))
    tables["references"] = md_table(
        ["problem", "params", "mixed disagreement s=5 vs s=7", "reference steps"],
        [[p, prm, f"{d:.1e}", n] for (p, prm), (d, n) in sorted(refs.items())],
    )

    # B2
    b2 = data["B2"]
    diurnal = [r for r in b2 if r["problem"] == "diurnal"]
    if diurnal:
        tables["B2_diurnal"] = wp_table(
            diurnal, None, extra_cols=("params", "nli", "nsli", "nps", "nfe")
        )
    heat_band = [
        r
        for r in b2
        if r["problem"] == "heat2d"
        and r.get("ls") == "band"
        and "velocity" not in r.get("params", "")
    ]
    if heat_band:
        tables["B2_heat_band"] = wp_table(heat_band, None, extra_cols=("params", "nfe", "nsli"))
    heat_mf = [r for r in b2 if r["problem"] == "heat2d" and r.get("ls") == "pcg"]
    if heat_mf:
        tables["B2_heat_matrixfree"] = wp_table(
            heat_mf, None, extra_cols=("params", "cell_key", "nfe", "nli", "nsli", "nps")
        )
    adr = [r for r in b2 if r["problem"] == "heat2d" and "velocity" in r.get("params", "")]
    if adr:
        tables["B2_adr"] = wp_table(adr, None, extra_cols=("params", "nfe", "nsli", "nscf"))

    # B3
    if data["B3"]:
        tables["B3_klin"] = klin_table(data["B3"], data["B3_klin_log"])

    # B4
    if data["B4_fixed"]:
        tables["B4_fixed"] = fixed_step_table(data["B4_fixed"])
    if data["B4_estimate"]:
        tables["B4_estimate"] = estimate_table(data["B4_estimate"])

    # B5
    if data["B5"]:
        tables["B5_heat1dfem"] = wp_table(
            data["B5"], None, extra_cols=("params", "nmass_setups", "nmass_solves", "nmass_mult")
        )
        for n in ["101", "401", "1601"]:
            rows = [r for r in data["B5"] if r.get("params") == f"n={n}"]
            if rows:
                tables[f"B5_speedup_n{n}"] = speedup_table(rows)

    # B6, B7
    if data["B6"]:
        tables["B6_ablation"] = ablation_table(data["B6"])
    if data["B7_roots"]:
        tables["B7_roots"] = roots_table(data["B7_roots"])
    if data["B7_dense"]:
        tables["B7_dense"] = dense_table(data["B7_dense"])

    for name, text in tables.items():
        (tables_dir / f"{name}.md").write_text(
            f"<!-- generated by make_report.py from {results.name} -->\n" + text
        )

    if not args.no_figures:
        plt = setup_matplotlib()
        for p in ["rober", "vdpol", "orego", "hires"]:
            rows = [r for r in b1 if r["problem"] == p]
            if rows:
                fig_work_precision(plt, rows, p, figdir / f"wp_{p}.png")
        groups = [
            (p, [r for r in b1 if r["problem"] == p]) for p in ["rober", "vdpol", "orego", "hires"]
        ]
        groups = [g for g in groups if g[1]]
        if groups:
            fig_grid_wp(plt, groups, figdir / "wp_summary_time.png")
        if diurnal:
            grids = sorted(
                {r.get("params", "") for r in diurnal},
                key=lambda s: to_int(s.split("=")[1].split(";")[0]) or 0,
            )
            fig_grid_wp(
                plt,
                [(f"diurnal {g}", [r for r in diurnal if r.get("params") == g]) for g in grids],
                figdir / "pde_diurnal_wp.png",
            )
        if heat_band or adr:
            groups = []
            if heat_band:
                groups.append(("heat2d 100x100, band", heat_band))
            for v in sorted({r.get("params", "") for r in adr}):
                groups.append(
                    (f"advection-diffusion {v}", [r for r in adr if r.get("params") == v])
                )
            fig_grid_wp(plt, groups, figdir / "pde_heat_wp.png")
        if data["B3"]:
            fig_klin(plt, data["B3"], data["B3_klin_log"], figdir / "klin_vs_s.png")
        if data["B4_fixed"]:
            fig_fixed_step(plt, data["B4_fixed"], figdir / "order_fixed_step.png")
        if data["B4_estimate"]:
            fig_estimate(plt, data["B4_estimate"], figdir / "order_estimate.png")
        if data["B5"]:
            sizes = sorted(
                {r.get("params", "") for r in data["B5"]},
                key=lambda s: to_int(s.split("=")[1]) or 0,
            )
            fig_grid_wp(
                plt,
                [
                    (f"heat1dfem {n}", [r for r in data["B5"] if r.get("params") == n])
                    for n in sizes
                ],
                figdir / "mass_heat1d_wp.png",
            )
        if data["B7_roots"] or data["B7_dense"]:
            fig_roots_dense(
                plt, data["B7_roots"], data["B7_dense"], figdir / "roots_dense_output.png"
            )

    if args.inject:
        return inject(Path(args.inject), tables, check_only=False)
    if args.check:
        return inject(Path(args.check), tables, check_only=True)
    return 0


def arkode_report(args, data, results, figdir, tables_dir, tables) -> int:
    """The tables and figures of BENCHMARK_ARKODE.md (families A0-A8)."""
    from report_arkode import arkode_tables

    plt = None if args.no_figures else setup_matplotlib()
    tables.update(arkode_tables(data, results, figdir, plt))
    tables["validation"] = validation_table(data["validation"])
    tables["sentinel"] = sentinel_table(data["sentinel"])
    for name, text in tables.items():
        (tables_dir / f"{name}.md").write_text(
            f"<!-- generated by make_report.py --report arkode from {results.name} -->\n" + text
        )
    if args.inject:
        return inject(Path(args.inject), tables, check_only=False)
    if args.check:
        return inject(Path(args.check), tables, check_only=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
