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
"""Helpers shared by the family reports (report_bdf.py, report_arkode.py).

Row selection (the untimed twins never enter a timing comparison),
problem parameters, the classes of configurations by label, and the
matched-accuracy ratio of the fastest configuration of a class against a
base configuration (the same interpolation in (mescd, log time) as the B1
speed-up tables).
"""

from __future__ import annotations

import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from benchmark_common import to_float  # noqa: E402
from make_report import (  # noqa: E402
    by_label,
    fmt,
    interp_log_time,
    label_sort_key,
    md_table,
    ok_rows,
)

BASE = "FIRKODE s=3"

# labels of the untimed twins of a cell (order histogram, profile, spectrum)
TWIN_SUFFIXES = (" [order]", " [profile]", " [spectrum]")


def is_twin(label: str) -> bool:
    return label.endswith(TWIN_SUFFIXES)


def timed_rows(rows):
    """Rows that may enter a timing comparison (no untimed twins; a label
    such as "ARKODE DIRK4 [ark_robertson settings]" is a timed cell)."""
    return [r for r in rows if not is_twin(r.get("label", ""))]


def params_of(r) -> dict:
    out = {}
    for kv in r.get("params", "").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def pnum(r, key, default=math.nan) -> float:
    v = params_of(r).get(key)
    return to_float(v) if v is not None else default


def groups_by_params(rows, keys):
    out = defaultdict(list)
    for r in rows:
        p = params_of(r)
        out[tuple(p.get(k, "") for k in keys)].append(r)
    return out


# ---------------------------------------------------------------------------
# Classes of configurations, by label
# ---------------------------------------------------------------------------


def is_firkode(label: str) -> bool:
    return label.startswith("FIRKODE")


def is_bdf(label: str) -> bool:
    return label.startswith(("CVODE", "IDA"))


# ARKODE's default DIRK tables of order 3-5 ("best default")
DIRK_DEFAULT_LABELS = ("ARKODE DIRK3", "ARKODE DIRK4", "ARKODE DIRK5")


def is_dirk_default(label: str) -> bool:
    return label in DIRK_DEFAULT_LABELS


def is_dirk(label: str) -> bool:
    """Every declared ARKODE DIRK configuration ("best of K"); the
    configuration tuned on ROBER (a "[...]" label) is reported apart."""
    return label.startswith("ARKODE DIRK") and "[" not in label


def is_imex(label: str) -> bool:
    return label.startswith("ARKODE IMEX")


def is_lsrk(label: str) -> bool:
    return label.startswith("ARKODE LSRK")


def is_mri(label: str) -> bool:
    return label.startswith("ARKODE MRI")


def is_arkode(label: str) -> bool:
    return label.startswith("ARKODE") and "[" not in label


# ---------------------------------------------------------------------------
# Matched-accuracy ratios
# ---------------------------------------------------------------------------


def best_ratio(rows, target, is_candidate, base=BASE):
    """(ratio, label) of the fastest candidate configuration against the
    base at the target mescd; NaN when either side lacks data."""
    groups = by_label(ok_rows(timed_rows(rows)))
    if base not in groups:
        return math.nan, ""
    tb = interp_log_time(groups[base], target)
    best, best_label = math.nan, ""
    for label, g in groups.items():
        if label == base or not is_candidate(label):
            continue
        t = interp_log_time(g, target)
        if not math.isnan(t) and (math.isnan(best) or t < best):
            best, best_label = t, label
    if math.isnan(best) or math.isnan(tb):
        return math.nan, best_label
    return best / tb, best_label


def ratio_cells(rows, is_candidate, targets, base=BASE):
    cells = []
    for t in targets:
        ratio, label = best_ratio(rows, t, is_candidate, base)
        cells.append("NA" if math.isnan(ratio) else fmt(ratio, 3) + f" ({label})")
    return cells


def crossover_mescd(rows, is_candidate, targets=(8.0, 9.0, 10.0, 11.0, 12.0, 13.0), base=BASE):
    """Accuracy at which the best candidate time equals the base time, by
    linear interpolation of log(ratio) over the targets; NaN when the
    ratio never crosses 1."""
    pts = []
    for t in targets:
        ratio, _ = best_ratio(rows, t, is_candidate, base)
        if not math.isnan(ratio):
            pts.append((t, math.log10(ratio)))
    for (t0, r0), (t1, r1) in zip(pts, pts[1:]):
        if (r0 < 0) != (r1 < 0) or r0 == 0:
            if r1 == r0:
                return t0
            return t0 + (0 - r0) / (r1 - r0) * (t1 - t0)
    return math.nan


# ---------------------------------------------------------------------------
# Profiled twins
# ---------------------------------------------------------------------------


def profile_rows_and_bars(rows):
    """Table rows and stacked-bar data of the "[profile]" twins."""
    out = []
    bars = []
    for r in sorted(
        [r for r in rows if "[profile]" in r.get("label", "") and r.get("driver_flag") == "ok"],
        key=lambda r: (
            r.get("problem", ""),
            to_float(params_of(r).get("n", params_of(r).get("mx", "0"))),
            label_sort_key(r["label"]),
        ),
    ):
        ts, tso, tr, tj, tp = (
            to_float(r.get(k)) for k in ("t_setup", "t_solve", "t_rhs", "t_jac", "t_profiled")
        )
        out.append(
            [
                r.get("problem", ""),
                r.get("params", ""),
                r["label"].replace(" [profile]", ""),
                f"{to_float(r.get('rtol')):.0e}",
                r.get("nst", "NA"),
                r.get("nsetups", "NA"),
                r.get("nje", "NA"),
                fmt(ts, 3),
                fmt(tso, 3),
                fmt(tr, 3),
                fmt(tj, 3),
                fmt(tp, 3),
                fmt(ts / tp, 3) if tp > 0 else "NA",
            ]
        )
        bars.append(
            (
                f"{r.get('problem', '')} {r.get('params', '')}\n{r['label'].replace(' [profile]', '')}",
                ts,
                tso,
                tr,
                tj,
                max(tp - ts - tso - tr - tj, 0.0),
            )
        )
    return out, bars


PROFILE_HEADER = [
    "problem",
    "params",
    "method",
    "rtol",
    "steps",
    "setups",
    "Jac",
    "t_setup (s)",
    "t_solve (s)",
    "t_rhs (s)",
    "t_jac (s)",
    "wall (s)",
    "setup share",
]


def profile_table(rows) -> str:
    out, _ = profile_rows_and_bars(rows)
    return md_table(PROFILE_HEADER, out) if out else "(no profiled runs)\n"


def fig_profile(plt, rows, path):
    _, bars = profile_rows_and_bars(rows)
    if not bars:
        return
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(bars)), 4.5))
    x = list(range(len(bars)))
    bottom = [0.0] * len(bars)
    for idx, name in ((1, "setup"), (2, "solve"), (3, "rhs"), (4, "jac"), (5, "other")):
        vals = [b[idx] for b in bars]
        ax.bar(x, vals, bottom=bottom, label=name)
        bottom = [b + v for b, v in zip(bottom, vals)]
    ax.set_xticks(x)
    ax.set_xticklabels([b[0] for b in bars], rotation=90, fontsize=6)
    ax.set_ylabel("seconds (profiled run, rtol 1e-6)")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
