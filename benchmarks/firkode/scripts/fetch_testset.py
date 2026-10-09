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
# Download the primary sources of the published test problems of
# problems_testset.cpp (families B10 and B11, and DETEST B5) into a work
# directory and check them against the pinned SHA-256 digests: the
# Fortran 77 sources of the Test Set for IVP Solvers (INdAM Bari) for
# BEAM, POLLU and RINGMOD, the Fortran sources and reference values of
# E. Hairer's stiff test set for BEAM and PLATE, and the stiff DETEST
# package STDTST (ACM TOMS Algorithm 648, Enright and Pryce), whose FCN,
# PDERV, IVALU and EVALU define problem B5 and its end-point values. The
# C++ problems were transcribed from these files; the script lets anyone
# re-verify the transcription and the reference values in
# benchmarks/firkode/references/. Nothing downloaded is part of SUNDIALS
# or is written into the repository.
#
# Usage: fetch_testset.py [--workdir DIR] [--no-check]
# ---------------------------------------------------------------

import argparse
import hashlib
import http.client
import pathlib
import sys
import tempfile
import time
import urllib.request

BARI = "https://archimede.uniba.it/~testset/src/problems/"
HAIRER = "https://www.unige.ch/~hairer/testset/stiff/"
NETLIB_TOMS = "https://www.netlib.org/toms/"

# file name in the work directory -> (URL, SHA-256 as downloaded on 2026-10-09)
SOURCES = {
    "bari_beam.f": (
        BARI + "beam.f",
        "89363648bcf04e73e83d019f7b5a9ea8c10748ff49b6308621f374e6629d16fe",
    ),
    "bari_pollu.f": (
        BARI + "pollu.f",
        "2aba777ee6de34e0ee074951375e029ad5171e937dabb7ab4c6461c0736e6c20",
    ),
    "bari_ringmod.f": (
        BARI + "ringmod.f",
        "613c980563551c9146e2b23a8e03d18d1243d348d238bd3cf6419d5c59702d32",
    ),
    "hairer_beam_equation.f": (
        HAIRER + "beam/equation.f",
        "e391ab503740462ea8ab415b122a4a614f67f4472d46a54d6d0a2618c3a45941",
    ),
    "hairer_beam_driver_radau5.f": (
        HAIRER + "beam/driver_radau5.f",
        "3f6b0f7f37f5f9ae086e0562e00788671f75bc42ba945f0f4c012aeb75a12a61",
    ),
    "hairer_beam_res_exact_pic": (
        HAIRER + "beam/res_exact_pic",
        "63ef0c49745bf0aeab9bf786536cd0fe2c62015ec4ba5bd38dfd0d135aac50b1",
    ),
    "hairer_plate_equation.f": (
        HAIRER + "plate/equation.f",
        "00bd7f7b903c2c39d1f52ab71e2d8ad86d5569b09b8feb2b21a5271d795278fe",
    ),
    "hairer_plate_driver_radau5.f": (
        HAIRER + "plate/driver_radau5.f",
        "4f9f4080e384a3f3c3af1dd43f9249d8b845069db7a944bd79ca01252abcb436",
    ),
    "hairer_plate_res_exact_pic": (
        HAIRER + "plate/res_exact_pic",
        "57205c8da530558692ae0c140b09e4bbd5b7c798b2982a63639284f721b7ddcb",
    ),
    # NSDTST and STDTST in one compressed file; STDTST is the second package
    "netlib_toms648.gz": (
        NETLIB_TOMS + "648.gz",
        "d472de5f6e43484b9a3616723414a040d20bbcbbf9039f07c7224281dcc34bf6",
    ),
}


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", help="where the files go (default: a temporary directory)")
    ap.add_argument(
        "--no-check", action="store_true", help="print the digests instead of checking them"
    )
    args = ap.parse_args()
    workdir = pathlib.Path(args.workdir) if args.workdir else pathlib.Path(tempfile.mkdtemp())
    workdir.mkdir(parents=True, exist_ok=True)
    bad = 0
    for name, (url, digest) in SOURCES.items():
        path = workdir / name
        if not path.exists():
            for attempt in range(5):
                try:
                    with urllib.request.urlopen(url, timeout=60) as response:
                        path.write_bytes(response.read())
                    break
                except (OSError, http.client.HTTPException) as exc:
                    if attempt == 4:
                        raise
                    print(f"retrying {name} after: {exc}")
                    time.sleep(2.0 * (attempt + 1))
        got = sha256(path)
        if args.no_check or digest is None:
            print(f"{got}  {name}  ({url})")
        elif got != digest:
            print(f"MISMATCH {name}: {got} differs from the pinned {digest}")
            bad += 1
        else:
            print(f"ok  {name}")
    print(f"files in {workdir}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
