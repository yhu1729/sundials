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
# Shared pieces of the FIRKODE benchmark driver and report generator:
# benchmark cells, running firkode_benchmark, CSV files, references,
# provenance and the plotting style.
# -----------------------------------------------------------------
"""Shared utilities for the FIRKODE benchmark scripts."""

from __future__ import annotations

import csv
import dataclasses
import datetime as dt
import hashlib
import json
import math
import os
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

# Driver-side columns prepended to the columns printed by firkode_benchmark
DRIVER_COLUMNS = [
    "family",
    "cell_key",
    "label",
    "params",
    "driver_flag",
    "timing_reliable",
    "wall_s",
    "timestamp",
]

# Absolute-tolerance factor atol = ATOL_FACTOR[problem] * rtol
ATOL_FACTOR = {
    "rober": 1.0e-6,
    "vdpol": 1.0,
    "orego": 1.0,
    "hires": 1.0,
    "analytic": 1.0e-4,
    "pr": 1.0e-4,
    "kpr": 1.0e-4,
    "diurnal": 100.0,
    "heat2d": 1.0e-4,
    "heat1dfem": 1.0e-4,
    "sector": 1.0,
    "wave1dfem": 1.0e-4,
    "pollu": 1.0,
    "ringmod": 1.0,
    "beam": 1.0,
    "plate": 1.0e-3,
    "detest_b5": 1.0,
    "bruss3d": 1.0e-2,
    "kprmri": 1.0e-4,
}

# Reference-solution tolerance per problem (FIRKODE s=3, cross-checked
# with CVODE, or IDA for a mass-matrix problem, by --ref-check)
# cap of one reference computation (seconds)
REF_CAP_S = 2 * 3600.0

REF_RTOL = {
    "rober": 1.0e-13,
    "vdpol": 1.0e-13,
    "orego": 1.0e-13,
    "hires": 1.0e-13,
    "diurnal": 1.0e-10,
    "heat2d": 1.0e-10,
    "heat1dfem": 1.0e-12,
    "wave1dfem": 1.0e-12,
    "pollu": 1.0e-13,
    "ringmod": 1.0e-11,
    "beam": 1.0e-13,
    "plate": 1.0e-13,
    "bruss3d": 1.0e-10,
}

# Published reference solutions (benchmarks/firkode/references), used for a
# problem run with its published parameters (no params, no restarts)
PUBLISHED_REF = {
    "pollu": "pollu.ref",
    "ringmod": "ringmod.ref",
    "beam": "beam.ref",
    "plate": "plate.ref",
}


# Cross-check integrator of a computed reference: "default" lets the harness
# choose (CVODE, or IDA for a mass-matrix problem; FIRKODE when the reference
# integrator is not FIRKODE), a name selects it, None skips the check
def ref_check_integrator(problem: str, params: dict) -> str | None:
    if problem == "wave1dfem":
        # FIRKODE cannot cross-check the consistent-mass reference at the
        # reference tolerance (see the known issues); the lumped ODE form it can
        return "firkode" if params.get("mass") == "lumped" else None
    return "default"


# Floor on the reference's absolute tolerance. The finite-element heat problem
# has components that are exactly zero (Dirichlet rows); an absolute tolerance
# below the roundoff level of those components stalls FIRKODE (see the known
# issues in src/firkode/BENCHMARK.md), so the reference is not run that tight.
REF_ATOL_FLOOR = {"heat1dfem": 1.0e-13}

# Reference integrator per problem (FIRKODE s=5 unless listed). FIRKODE cannot
# integrate the finite-element heat problem on fine meshes at reference
# tolerances (its step size collapses, see the known issues), so IDA computes
# that reference; it is not cross-checked.
REF_INTEGRATOR = {"heat1dfem": "ida", "wave1dfem": "ida"}

# below this reference tolerance only the extended-precision build (--ref-bin)
# can compute a meaningful reference
REF_RTOL_DOUBLE_FLOOR = 1.0e-14

TIMING_FLOOR_S = 1.0e-4  # below this, timings are reported as unreliable


@dataclasses.dataclass
class Cell:
    """One benchmark run of firkode_benchmark."""

    family: str
    problem: str
    integrator: str
    label: str
    rtol: float | None = None
    atol: float | None = None
    params: dict = dataclasses.field(default_factory=dict)  # --mx 20 etc.
    args: list = dataclasses.field(default_factory=list)  # extra CLI flags
    cap_s: float = 300.0
    reps: str | int = "auto"  # "auto", "single" or an integer
    timed: bool = True
    needs_ref: bool = True
    ref_rtol: float | None = None  # override of REF_RTOL for this cell

    def params_str(self) -> str:
        return ";".join(f"{k}={v}" for k, v in sorted(self.params.items()))

    def key(self) -> str:
        parts = [
            self.problem,
            self.params_str(),
            self.integrator,
            self.label,
            "NA" if self.rtol is None else f"{self.rtol:g}",
            "NA" if self.atol is None else f"{self.atol:g}",
            " ".join(self.args),
        ]
        return "|".join(parts)

    def command(self, binary: str, extra: list | None = None) -> list:
        cmd = [binary, "--problem", self.problem, "--integrator", self.integrator]
        for k, v in sorted(self.params.items()):
            cmd += [f"--{k}", str(v)]
        if self.rtol is not None:
            cmd += ["--rtol", f"{self.rtol:.17g}"]
        if self.atol is not None:
            cmd += ["--atol", f"{self.atol:.17g}"]
        cmd += list(self.args)
        if extra:
            cmd += list(extra)
        return cmd


def problem_hash(problem: str, params: dict, ref_rtol: float | None = None) -> str:
    text = problem + "|" + ";".join(f"{k}={v}" for k, v in sorted(params.items()))
    if ref_rtol is not None:
        text += f"|ref_rtol={ref_rtol:g}"
    return hashlib.sha1(text.encode()).hexdigest()[:12]


def parse_csv_output(text: str) -> dict | None:
    """Parse the header+row printed by firkode_benchmark."""
    lines = [ln for ln in text.strip().splitlines() if ln.strip()]
    if len(lines) < 2:
        return None
    header = lines[-2].split(",")
    row = lines[-1].split(",")
    if len(header) != len(row) or header[0] != "problem":
        return None
    return dict(zip(header, row))


def run_binary(cmd: list, cap_s: float, env: dict, log_path: Path | None = None):
    """Run firkode_benchmark; returns (row dict or None, driver_flag, wall_s)."""
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=cap_s, env=env, check=False
        )
    except subprocess.TimeoutExpired:
        return None, "TIMEOUT", time.perf_counter() - t0
    wall = time.perf_counter() - t0
    row = parse_csv_output(proc.stdout)
    if log_path is not None and (proc.returncode != 0 or row is None):
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(
            "$ " + " ".join(cmd) + f"\nrc={proc.returncode}\n" + proc.stderr + proc.stdout
        )
    if proc.returncode == 0 and row is not None:
        return row, "ok", wall
    if proc.returncode == 2:
        return row, "UNSUPPORTED", wall
    if proc.returncode == 3 and row is not None:
        return row, f"SOLVER_FAIL_{row.get('flag', '?')}", wall
    return row, f"FAIL_rc{proc.returncode}", wall


def to_float(value) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def to_int(value) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


class ResultsTable:
    """Append-only CSV of benchmark rows, resumable by cell key."""

    def __init__(self, path: Path):
        self.path = path
        self.rows: list[dict] = []
        self.columns: list[str] = []
        if path.exists():
            with path.open(newline="") as fp:
                reader = csv.DictReader(fp)
                self.columns = list(reader.fieldnames or [])
                self.rows = list(reader)

    def done_keys(self, include_failed: bool = True) -> set:
        keys = set()
        for r in self.rows:
            if include_failed or r.get("driver_flag") == "ok":
                keys.add(r["cell_key"])
        return keys

    def append(self, row: dict) -> None:
        new_cols = [c for c in row if c not in self.columns]
        if new_cols:
            self.columns += new_cols
            self._rewrite()
        self.rows.append(row)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        write_header = not self.path.exists() or self.path.stat().st_size == 0
        with self.path.open("a", newline="") as fp:
            writer = csv.DictWriter(fp, fieldnames=self.columns, extrasaction="ignore")
            if write_header:
                writer.writeheader()
            writer.writerow({c: row.get(c, "NA") for c in self.columns})
            fp.flush()

    def remove_keys(self, keys: set) -> None:
        self.rows = [r for r in self.rows if r["cell_key"] not in keys]
        self._rewrite()

    def _rewrite(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("w", newline="") as fp:
            writer = csv.DictWriter(fp, fieldnames=self.columns, extrasaction="ignore")
            writer.writeheader()
            for r in self.rows:
                writer.writerow({c: r.get(c, "NA") for c in self.columns})


def read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as fp:
        return list(csv.DictReader(fp))


def read_reference_comments(path: Path) -> dict:
    out = {}
    if not path.exists():
        return out
    with path.open() as fp:
        for line in fp:
            if not line.startswith("#"):
                break
            m = re.match(r"#\s*([a-z_]+)=(.*)", line.strip())
            if m:
                out[m.group(1)] = m.group(2)
    return out


# ---------------------------------------------------------------------------
# Parallel workers: cell claims, the reference lock and the merge of the
# per-worker CSV files
# ---------------------------------------------------------------------------


def claim_cell(results: Path, family: str, key: str, worker: int) -> bool:
    """Atomically claim a cell for this worker (one claim file per cell).

    Returns False when another worker holds the claim. Claims are only a
    guard against concurrent duplicates: completed cells are recognized
    from the CSV rows, and stale claims (a worker that died) are removed
    with --clear-claims or by the parent of --workers before it starts.
    """
    cdir = results / "claims" / family
    cdir.mkdir(parents=True, exist_ok=True)
    path = cdir / (hashlib.sha1(key.encode()).hexdigest() + ".claim")
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False
    with os.fdopen(fd, "w") as fp:
        fp.write(f"worker={worker}\nkey={key}\ntime={time.time():.0f}\n")
    return True


def clear_claims(results: Path) -> int:
    """Remove every claim file; returns the number removed."""
    cdir = results / "claims"
    n = 0
    if cdir.exists():
        for path in cdir.rglob("*.claim"):
            path.unlink()
            n += 1
    return n


def family_files(results: Path, family: str) -> list[Path]:
    """The merged CSV of a (sub-)family and every worker file of it."""
    files = [results / f"{family}.csv"]
    files += sorted(results.glob(f"{family}.worker*.csv"))
    return [f for f in files if f.exists()]


def merge_results(results: Path, remove_worker_files: bool = True) -> dict:
    """Fold the per-worker CSV files into the family files (and the worker
    provenance into provenance.json); later rows win for a repeated key."""
    merged = {}
    worker_files = sorted(results.glob("*.worker*.csv"))
    bases = sorted({re.sub(r"\.worker\d+\.csv$", "", f.name) for f in worker_files})
    for base in bases:
        rows: dict[str, dict] = {}
        columns: list[str] = []
        for path in family_files(results, base):
            for r in read_csv(path):
                for c in r:
                    if c not in columns:
                        columns.append(c)
                k = r.get("cell_key", "")
                if k in rows and r.get("timestamp", "") < rows[k].get("timestamp", ""):
                    continue
                rows[k] = r
        target = results / f"{base}.csv"
        with target.open("w", newline="") as fp:
            writer = csv.DictWriter(fp, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            for r in rows.values():
                writer.writerow({c: r.get(c, "NA") for c in columns})
        merged[base] = len(rows)
        if remove_worker_files:
            for path in results.glob(f"{base}.worker*.csv"):
                path.unlink()
    prov_path = results / "provenance.json"
    worker_provs = sorted(results.glob("provenance.worker*.json"))
    if worker_provs and prov_path.exists():
        prov = json.loads(prov_path.read_text())
        for path in worker_provs:
            prov.setdefault("runs", []).append(json.loads(path.read_text()))
            if remove_worker_files:
                path.unlink()
        save_json(prov_path, prov)
    return merged


def ensure_reference(
    binary: str,
    problem: str,
    params: dict,
    rtol: float,
    atol: float,
    refdir: Path,
    env: dict,
    ref_rtol: float | None = None,
    extended: bool = False,
    budget_s: float | None = None,
) -> tuple[Path | None, dict]:
    """Compute (once) and return the reference file for a problem/size.

    The reference is FIRKODE s=3 at REF_RTOL[problem] (or REF_INTEGRATOR),
    cross-checked by a second integrator (the mixed error between the two
    is stored in the file). With --ref-bin the extended-precision build
    computes it and the double build reads the 21-digit file.
    """
    if problem not in REF_RTOL:
        return None, {}
    refdir.mkdir(parents=True, exist_ok=True)
    path = refdir / f"ref_{problem}_{problem_hash(problem, params, ref_rtol)}.dat"
    if path.exists():
        return path, read_reference_comments(path)
    # the published value of a problem run as published
    if ref_rtol is None and problem in PUBLISHED_REF and not params:
        src = Path(__file__).resolve().parents[1] / "references" / PUBLISHED_REF[problem]
        path.write_text(src.read_text())
        return path, read_reference_comments(path)
    if (ref_rtol or REF_RTOL[problem]) < REF_RTOL_DOUBLE_FLOOR and not extended:
        raise RuntimeError(
            f"a reference at rtol {ref_rtol or REF_RTOL[problem]:g} needs the extended-precision "
            "build (--ref-bin)"
        )
    failed = path.with_suffix(".failed")
    if failed.exists():
        raise RuntimeError(f"reference computation failed earlier for {problem} ({failed})")
    # a --fast run that timed out leaves a marker only other --fast runs respect
    fast_marker = path.with_suffix(".fast_timeout")
    if budget_s is not None and fast_marker.exists():
        raise ReferenceTimeout(f"the reference of {problem} exceeded the fast budget earlier")
    if budget_s is None and fast_marker.exists():
        fast_marker.unlink()
    # one worker computes, the others wait for the file (or for the failure
    # marker); the lock is a file created atomically
    lock = path.with_suffix(".lock")
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.close(fd)
    except FileExistsError:
        print(
            f"  waiting for the reference of {problem} computed by another worker ...", flush=True
        )
        t0 = time.time()
        while not path.exists() and not failed.exists():
            if not lock.exists():
                break
            if time.time() - t0 > REF_WAIT_S:
                raise RuntimeError(f"timed out waiting for the reference of {problem}")
            time.sleep(2.0)
        if path.exists():
            return path, read_reference_comments(path)
        raise RuntimeError(f"reference computation failed in another worker for {problem}")
    try:
        return compute_reference(
            binary, problem, params, rtol, atol, path, env, ref_rtol, budget_s
        )
    except ReferenceTimeout:
        # the harness writes the file before its own timed run; a complete
        # file is usable even if the run was cut afterwards
        if path.exists():
            print(f"  (the reference of {problem} was written before the cut)", flush=True)
            return path, read_reference_comments(path)
        fast_marker.write_text(f"budget {budget_s} s\n")
        raise
    except RuntimeError:
        failed.write_text("see the .err file\n")
        raise
    finally:
        if lock.exists():
            lock.unlink()


# seconds a worker waits for another worker's reference computation
REF_WAIT_S = 4 * 3600.0


class ReferenceTimeout(RuntimeError):
    """A reference computation exceeded the time budget of a --fast run."""


def compute_reference(
    binary: str,
    problem: str,
    params: dict,
    rtol: float,
    atol: float,
    path: Path,
    env: dict,
    ref_rtol: float | None = None,
    budget_s: float | None = None,
) -> tuple[Path, dict]:
    refdir = path.parent
    # --write-ref is a side effect of an ordinary run, so run that ordinary
    # run at the reference tolerance as well (a loose-tolerance run with the
    # large step cap of the reference could take hours if the integrator
    # stalls, see the known issues in BENCHMARK_BDF.md)
    ref_rtol = ref_rtol or REF_RTOL[problem]
    ref_atol = max(ref_rtol * (atol / rtol), REF_ATOL_FLOOR.get(problem, 0.0))
    integ = REF_INTEGRATOR.get(problem, "firkode")
    cell = Cell("ref", problem, integ, "ref", ref_rtol, ref_atol, params, [])
    extra = ["--maxsteps", "200000000", "--write-ref", str(path), "--ref-rtol", f"{ref_rtol:g}"]
    extra += ["--ref-integrator", integ]
    check = ref_check_integrator(problem, params)
    if check == "default" and integ != "firkode":
        check = "firkode"
    if check is not None:
        extra += ["--ref-check"]
        if check != "default":
            extra += ["--ref-check-integrator", check]
    cmd = cell.command(binary, extra)
    print(f"  computing reference for {problem} {cell.params_str()} ...", flush=True)
    cap = min(REF_CAP_S, budget_s) if budget_s is not None else REF_CAP_S
    row, flag, wall = run_binary(cmd, cap, env, refdir / f"{path.stem}.err")
    if flag == "TIMEOUT" and budget_s is not None and cap < REF_CAP_S:
        raise ReferenceTimeout(
            f"the reference of {problem} exceeded the fast budget of {budget_s:g} s"
        )
    if flag != "ok" or not path.exists():
        raise RuntimeError(f"reference computation failed ({flag}) for {problem}")
    info = read_reference_comments(path)
    print(
        f"    done in {wall:.1f} s: nst_ref={info.get('nst_ref')} "
        f"ref_disagreement={info.get('ref_disagreement')}",
        flush=True,
    )
    return path, info


def git_info(repo: Path) -> dict:
    def run(args):
        try:
            return subprocess.run(
                ["git", "-C", str(repo)] + args, capture_output=True, text=True, check=True
            ).stdout.strip()
        except Exception:
            return "unavailable"

    return {
        "commit": run(["rev-parse", "HEAD"]),
        "describe": run(["describe", "--always", "--dirty"]),
        "branch": run(["rev-parse", "--abbrev-ref", "HEAD"]),
        "dirty_files": run(["status", "--porcelain"]).splitlines(),
    }


def cmake_cache_info(build_dir: Path) -> dict:
    keys = [
        "CMAKE_BUILD_TYPE",
        "CMAKE_C_COMPILER",
        "CMAKE_CXX_COMPILER",
        "CMAKE_C_FLAGS",
        "CMAKE_C_FLAGS_RELEASE",
        "CMAKE_CXX_FLAGS_RELEASE",
        "SUNDIALS_PRECISION",
        "SUNDIALS_INDEX_SIZE",
        "SUNDIALS_LOGGING_LEVEL",
        "SUNDIALS_ENABLE_ERROR_CHECKS",
        "SUNDIALS_ENABLE_PROFILING",
        "SUNDIALS_ENABLE_ALL_WARNINGS",
        "BUILD_SHARED_LIBS",
        "CMAKE_GENERATOR",
    ]
    out = {}
    cache = build_dir / "CMakeCache.txt"
    if not cache.exists():
        return {"error": f"{cache} not found"}
    for line in cache.read_text().splitlines():
        for k in keys:
            if line.startswith(k + ":") or line.startswith(k + "="):
                out[k] = line.split("=", 1)[1]
    return out


def sysctl(name: str) -> str:
    try:
        return subprocess.run(
            ["sysctl", "-n", name], capture_output=True, text=True, check=True, timeout=5
        ).stdout.strip()
    except Exception:
        return "unavailable"


def tool_version(cmd: list) -> str:
    try:
        return (
            subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=20)
            .stdout.strip()
            .splitlines()[0]
        )
    except Exception:
        return "unavailable"


def linux_cpu_model(cpuinfo: str | None = None) -> str:
    """The CPU model on Linux: the first "model name" of /proc/cpuinfo (x86),
    else the "Model name:" line of lscpu (aarch64 kernels omit the former)."""
    try:
        text = cpuinfo if cpuinfo is not None else Path("/proc/cpuinfo").read_text()
    except OSError:
        text = ""
    for line in text.splitlines():
        key, sep, value = line.partition(":")
        if sep and key.strip() == "model name" and value.strip():
            return value.strip()
    if cpuinfo is not None:
        return "unavailable"
    try:
        out = subprocess.run(
            ["lscpu"], capture_output=True, text=True, check=True, timeout=5
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return "unavailable"
    for line in out.splitlines():
        key, sep, value = line.partition(":")
        if sep and key.strip() == "Model name" and value.strip():
            return value.strip()
    return "unavailable"


def linux_memory_bytes() -> str:
    try:
        return str(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES"))
    except (ValueError, OSError):
        return "unavailable"


def linux_os_release() -> str:
    try:
        return platform.freedesktop_os_release().get("PRETTY_NAME", "unavailable")
    except OSError:
        return "unavailable"


def machine_facts() -> dict:
    """CPU model, memory and OS version by platform (sysctl and sw_vers on
    macOS, /proc and os-release on Linux); a failed lookup is "unavailable"."""
    system = platform.system()
    if system == "Linux":
        return {
            "os": system,
            "cpu_brand": linux_cpu_model(),
            "memory_bytes": linux_memory_bytes(),
            "os_version": linux_os_release(),
        }
    return {
        "os": system,
        "cpu_brand": sysctl("machdep.cpu.brand_string"),
        "perf_levels": {
            "performance_cores": sysctl("hw.perflevel0.physicalcpu"),
            "efficiency_cores": sysctl("hw.perflevel1.physicalcpu"),
        },
        "memory_bytes": sysctl("hw.memsize"),
        "os_version": tool_version(["sw_vers", "-productVersion"]),
    }


def binary_version(binary: Path) -> str:
    """The --version output of a firkode_benchmark binary (or why it is unavailable)."""
    try:
        return subprocess.run(
            [str(binary), "--version"], capture_output=True, text=True, check=True, timeout=30
        ).stdout
    except (OSError, subprocess.SubprocessError) as exc:
        return f"unavailable: {exc}"


def collect_provenance(binary: Path, build_dir: Path, repo: Path, env: dict) -> dict:
    try:
        benchmark_version = subprocess.run(
            [str(binary), "--version"], capture_output=True, text=True, check=True, timeout=30
        ).stdout
    except Exception as exc:
        benchmark_version = f"unavailable: {exc}"
    try:
        import numpy

        numpy_version = numpy.__version__
    except ImportError:
        numpy_version = "unavailable"
    try:
        import matplotlib

        mpl_version = matplotlib.__version__
    except ImportError:
        mpl_version = "unavailable"
    return {
        "started_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
        "hostname": platform.node(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        **machine_facts(),
        "compiler": tool_version(["cc", "--version"]),
        "cmake": tool_version(["cmake", "--version"]),
        "git": git_info(repo),
        "cmake_cache": cmake_cache_info(build_dir),
        "binary": str(binary),
        "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "benchmark_version_output": benchmark_version,
        "python": sys.version.split()[0],
        "numpy": numpy_version,
        "matplotlib": mpl_version,
        "env": {k: env.get(k) for k in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
        "loadavg_start": os.getloadavg(),
    }


def save_json(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, default=str) + "\n")


# ---------------------------------------------------------------------------
# Plotting style shared by every figure (never color-only: marker and line
# style also encode the integrator; FIRKODE is filled, baselines open)
# ---------------------------------------------------------------------------

OKABE_ITO = {
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermilion": "#D55E00",
    "purple": "#CC79A7",
    "black": "#000000",
}

STYLE = {
    "FIRKODE s=1": dict(color=OKABE_ITO["black"], marker="o", ls="-", mfc="full"),
    "FIRKODE s=2": dict(color=OKABE_ITO["sky"], marker="s", ls="-", mfc="full"),
    "FIRKODE s=3": dict(color=OKABE_ITO["blue"], marker="^", ls="-", mfc="full"),
    "FIRKODE s=4": dict(color=OKABE_ITO["green"], marker="<", ls="-", mfc="full"),
    "FIRKODE s=5": dict(color=OKABE_ITO["green"], marker="D", ls="-", mfc="full"),
    "FIRKODE s=7": dict(color=OKABE_ITO["purple"], marker="v", ls="-", mfc="full"),
    "FIRKODE s=9": dict(color=OKABE_ITO["yellow"], marker="P", ls="-", mfc="full"),
    "CVODE BDF": dict(color=OKABE_ITO["vermilion"], marker="x", ls="--", mfc="none"),
    "CVODE BDF+STALD": dict(color=OKABE_ITO["vermilion"], marker="+", ls=":", mfc="none"),
    "CVODE BDF q<=4": dict(color=OKABE_ITO["orange"], marker=">", ls="--", mfc="none"),
    "CVODE BDF q<=3": dict(color=OKABE_ITO["yellow"], marker="<", ls="--", mfc="none"),
    "CVODE BDF q<=2": dict(color=OKABE_ITO["purple"], marker="v", ls="--", mfc="none"),
    "CVODE BDF msbp=50": dict(color=OKABE_ITO["vermilion"], marker="d", ls="-.", mfc="none"),
    "CVODE BDF h0=1e-6": dict(color=OKABE_ITO["vermilion"], marker="*", ls="--", mfc="none"),
    "ARKODE DIRK2": dict(color=OKABE_ITO["orange"], marker="1", ls=":", mfc="none"),
    "ARKODE DIRK3": dict(color=OKABE_ITO["orange"], marker="+", ls=":", mfc="none"),
    "ARKODE DIRK4": dict(color=OKABE_ITO["orange"], marker="*", ls=":", mfc="none"),
    "ARKODE DIRK5": dict(color=OKABE_ITO["vermilion"], marker="h", ls=":", mfc="none"),
    # IDA: green, dash-dot, the markers of the matching CVODE variants
    "IDA BDF": dict(color=OKABE_ITO["green"], marker="X", ls="-.", mfc="none"),
    "IDA BDF q<=4": dict(color=OKABE_ITO["green"], marker=">", ls="-.", mfc="none"),
    "IDA BDF q<=3": dict(color=OKABE_ITO["green"], marker="<", ls="-.", mfc="none"),
    "IDA BDF q<=2": dict(color=OKABE_ITO["green"], marker="v", ls="-.", mfc="none"),
    "IDA BDF h0=1e-6": dict(color=OKABE_ITO["green"], marker="*", ls="-.", mfc="none"),
    "FIRKODE s=3 periodic msbp=20": dict(color=OKABE_ITO["sky"], marker="^", ls="-.", mfc="none"),
    "FIRKODE s=3 periodic msbp=50": dict(color=OKABE_ITO["sky"], marker="D", ls="-.", mfc="none"),
    "FIRKODE s=3 h0=1e-6": dict(color=OKABE_ITO["blue"], marker="*", ls="-", mfc="none"),
    # the ARKODE configurations of BENCHMARK_ARKODE.md: DIRK orange and
    # vermilion (dotted), IMEX sky (dash-dot), LSRK purple (dashed), MRI
    # yellow (solid)
    "ARKODE DIRK4 optimal": dict(color=OKABE_ITO["orange"], marker="*", ls="-.", mfc="full"),
    "ARKODE DIRK5 optimal": dict(color=OKABE_ITO["vermilion"], marker="h", ls="-.", mfc="full"),
    "ARKODE DIRK SDIRK_5_3_4": dict(color=OKABE_ITO["orange"], marker="s", ls=":", mfc="none"),
    "ARKODE DIRK ESDIRK437L2SA_7_3_4": dict(
        color=OKABE_ITO["orange"], marker="D", ls=":", mfc="none"
    ),
    "ARKODE DIRK KVAERNO_7_4_5": dict(
        color=OKABE_ITO["vermilion"], marker="p", ls=":", mfc="none"
    ),
    "ARKODE DIRK ARK548L2SAb_DIRK_8_4_5": dict(
        color=OKABE_ITO["vermilion"], marker="8", ls=":", mfc="none"
    ),
    "ARKODE IMEX3": dict(color=OKABE_ITO["sky"], marker="+", ls="-.", mfc="none"),
    "ARKODE IMEX4": dict(color=OKABE_ITO["sky"], marker="*", ls="-.", mfc="none"),
    "ARKODE IMEX5": dict(color=OKABE_ITO["sky"], marker="h", ls="-.", mfc="none"),
    "ARKODE LSRK RKC2": dict(color=OKABE_ITO["purple"], marker="o", ls="--", mfc="none"),
    "ARKODE LSRK RKL2": dict(color=OKABE_ITO["purple"], marker="s", ls="--", mfc="none"),
    "ARKODE LSRK RKC2 power": dict(color=OKABE_ITO["purple"], marker="o", ls=":", mfc="full"),
    "ARKODE MRI ESDIRK34a": dict(color=OKABE_ITO["yellow"], marker="^", ls="-", mfc="none"),
    "ARKODE MRI ESDIRK34a htol": dict(color=OKABE_ITO["yellow"], marker="v", ls="-", mfc="full"),
    "ARKODE MRI IMEX-SR32": dict(color=OKABE_ITO["yellow"], marker="D", ls="-", mfc="none"),
}


def style_for(label: str) -> dict:
    base = STYLE.get(label)
    if base is None:
        for k, v in STYLE.items():
            if label.startswith(k):
                base = v
                break
    if base is None:
        base = dict(color="#777777", marker=".", ls="-", mfc="full")
    out = dict(color=base["color"], marker=base["marker"], ls=base["ls"])
    out["mfc"] = base["color"] if base["mfc"] == "full" else "none"
    return out
