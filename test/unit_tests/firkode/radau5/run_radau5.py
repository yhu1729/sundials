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
# Generator of the radau5.f reference files of firk_test_radau5_steps
# (FIRKODE design milestone M2).
#
# The script downloads radau5.f, decsol.f and dc_decsol.f (E. Hairer
# and G. Wanner) into a work directory, checks them against pinned
# SHA-256 digests, compiles them with radau5_driver.f90, and runs
# Robertson, van der Pol (mu = 1000) and the Oregonator at the nominal
# RADAU5 tolerances 1e-4, 1e-6 and 1e-8, recording the accepted steps,
# and at 1e-14 for the reference solutions. RADAU5 replaces the
# nominal rtol by 0.1 rtol^(2/3) and keeps atol / rtol; the .ref files
# record both, and FIRKODE runs at the transformed values (design
# Table 6.2). The downloaded sources stay in the work directory; only
# the .ref files are written to the output directory, this one by
# default.
#
# Usage: run_radau5.py [--workdir DIR] [--outdir DIR] [--fc FC]
#                      [--fflags FLAGS]
# ---------------------------------------------------------------

import argparse
import hashlib
import pathlib
import platform
import subprocess
import sys
import tempfile
import urllib.request

URL = "https://www.unige.ch/~hairer/prog/stiff/"

# the sources this script was written against (downloaded 2026-10-07)
SOURCES = {
    "radau5.f": "a3f1be7dbf719e232b710e1c7a50ad2079bef27c0a5f3b87b012f64ea65f39d3",
    "decsol.f": "da948678a939bf4f98774edb1ca02380decda3ae97f212e38704ae92deeae01f",
    "dc_decsol.f": "4c22deafa9619766b90489c63513c096630ccf63d0e418dce213a3443b7aaafc",
}

# atol / rtol of each problem, as in the FIRKODE benchmarks
PROBLEMS = {"rober": "1e-6", "vdpol": "1", "orego": "1"}

# nominal RADAU5 tolerances: the compared runs record their steps, the
# last gives the reference solutions
RTOLS = [("1e-4", True), ("1e-6", True), ("1e-8", True), ("1e-14", False)]


def sha256(path):
    """The SHA-256 digest of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch(workdir):
    """Download the sources unless present, and check their digests."""
    for name, digest in SOURCES.items():
        path = workdir / name
        if not path.exists():
            with urllib.request.urlopen(URL + name, timeout=60) as response:
                path.write_bytes(response.read())
        if sha256(path) != digest:
            sys.exit(f"{name}: SHA-256 {sha256(path)} differs from {digest}")


def build(workdir, fc, fflags, driver):
    """Compile radau5.f, its linear algebra and the driver."""
    exe = workdir / "radau5_driver"
    flags = fflags.split()
    sources = [str(workdir / name) for name in SOURCES]
    subprocess.run([fc, *flags, "-std=legacy", "-w", "-c", *sources], cwd=workdir, check=True)
    subprocess.run([fc, *flags, "-c", str(driver)], cwd=workdir, check=True)
    objects = [pathlib.Path(name).stem + ".o" for name in SOURCES]
    subprocess.run(
        [fc, *flags, "-o", str(exe), "radau5_driver.o", *objects], cwd=workdir, check=True
    )
    return exe


def main():
    here = pathlib.Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=pathlib.Path, default=None)
    parser.add_argument("--outdir", type=pathlib.Path, default=here)
    parser.add_argument("--fc", default="gfortran")
    parser.add_argument("--fflags", default="-O2")
    args = parser.parse_args()

    workdir = args.workdir or pathlib.Path(tempfile.mkdtemp(prefix="firkode_radau5_"))
    workdir.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)

    fetch(workdir)
    exe = build(workdir, args.fc, args.fflags, here / "radau5_driver.f90")
    version = subprocess.run(
        [args.fc, "--version"], capture_output=True, text=True, check=True
    ).stdout.splitlines()[0]
    provenance = [
        "# radau5.f reference for firk_test_radau5_steps, written by run_radau5.py",
        f"# sources {URL} radau5.f sha256 {SOURCES['radau5.f']}",
        f"# compiler {version}, flags {args.fflags}, platform {platform.platform()}",
    ]

    print(f"{'problem':8} {'rtol':>6} {'accepted':>8} {'rejected':>8} {'jac':>6}")
    for problem, quot in PROBLEMS.items():
        for rtol, steps in RTOLS:
            out = workdir / f"{problem}_{rtol}.ref"
            subprocess.run(
                [str(exe), problem, rtol, quot, str(out), "1" if steps else "0"], check=True
            )
            lines = out.read_text().splitlines()
            stats = dict(l.split(None, 1) for l in lines if len(l.split()) == 2 and l[0] != "#")
            print(
                f"{problem:8} {rtol:>6} {stats['naccpt']:>8} {stats['nrejct']:>8}"
                f" {stats['njac']:>6}"
            )
            (args.outdir / out.name).write_text("\n".join(provenance + lines) + "\n")


if __name__ == "__main__":
    main()
