#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): Yifan Hu @ UMBC
# -----------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------
# Parses the log written by firk_test_logging with suntools.logs (design
# hazard H31): the parser merges the records of a region into a dictionary
# and fails on an unbalanced region or a key repeated within one.
#
# Usage: firk_test_logging.py <log file> <suntools/src directory>
# Pass criterion: the log parses into step attempts.
# -----------------------------------------------------------------------------

import sys


def main():
    if len(sys.argv) != 3:
        print("usage: firk_test_logging.py <log file> <suntools/src>")
        return 2
    sys.path.insert(0, sys.argv[2])
    from suntools import logs

    try:
        steps = logs.log_file_to_list(sys.argv[1])
    except (KeyError, IndexError, ValueError) as err:  # inconsistent records
        print(f"FAIL: suntools.logs could not parse {sys.argv[1]}: {err!r}")
        return 1
    if len(steps) == 0:
        print("FAIL: the log has no step attempts")
        return 1
    print(f"suntools.logs parsed {len(steps)} step attempts")
    print("All tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
