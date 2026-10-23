#!/bin/bash
# ---------------------------------------------------------------------------------
# Programmer(s): Yifan Hu @ UMBC
# ---------------------------------------------------------------------------------
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
# ---------------------------------------------------------------------------------
# This script uses the Frama-C WP plugin to prove the ACSL contracts annotated in
# the SUNDIALS sources, including the absence of runtime errors (-wp-rte).
#
# Usage:
#    ./frama-c.sh <path to a configured SUNDIALS build directory>
#
# The build directory provides the generated sundials_config.h and
# sundials_export.h headers, so the proofs hold for the sunrealtype and
# sunindextype selected in that configuration.
#
# The following environment variables may be set to override the defaults:
#    FRAMAC          Frama-C executable (default: frama-c)
#    FRAMAC_PROVERS  comma separated list of WP provers (default: alt-ergo,z3)
#    FRAMAC_TIMEOUT  per-goal prover timeout in seconds (default: 10)
#
# The script exits with a non-zero status if any proof obligation is not proved.
# ---------------------------------------------------------------------------------

if [ $# -ne 1 ]; then
    echo "ERROR: Path to a configured SUNDIALS build directory required"
    exit 1
fi

srcdir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
builddir="$(cd "$1" && pwd)" || exit 1

framac="${FRAMAC:-frama-c}"
provers="${FRAMAC_PROVERS:-alt-ergo,z3}"
timeout="${FRAMAC_TIMEOUT:-10}"

if ! command -v "${framac}" > /dev/null 2>&1; then
    echo "ERROR: Frama-C executable '${framac}' not found"
    exit 1
fi

for header in sundials_config.h sundials_export.h; do
    if [ ! -f "${builddir}/include/sundials/${header}" ]; then
        echo "ERROR: ${builddir}/include/sundials/${header} not found"
        echo "       Configure SUNDIALS in ${builddir} first"
        exit 1
    fi
done

# Verification targets as "<source file relative to srcdir>:<function,...>"
targets=(
    "src/nvector/serial/nvector_serial.c:N_VConst_Serial"
)

logfile="$(mktemp)" || exit 1
trap 'rm -f "${logfile}"' EXIT

status=0

for target in "${targets[@]}"; do
    file="${target%%:*}"
    functions="${target#*:}"

    echo "Verifying ${functions} in ${file}"

    "${framac}" \
        -pp-annot \
        -cpp-extra-args="-I${srcdir}/include -I${builddir}/include -I${srcdir}/src/sundials" \
        -wp -wp-rte \
        -wp-model=Typed+cast \
        -wp-fct="${functions}" \
        -wp-prover="${provers}" \
        -wp-timeout="${timeout}" \
        "${srcdir}/${file}" \
        -then -report-classify \
        -report-unclassified-unknown=ERROR \
        -report-exit 2>&1 | tee "${logfile}"
    framac_status=${PIPESTATUS[0]}

    proved="$(sed -n 's/^\[wp\] Proved goals: *\([0-9]*\) *\/ *\([0-9]*\)$/\1 \2/p' "${logfile}")"
    read -r num_proved num_goals <<< "${proved}"

    if [ ${framac_status} -ne 0 ] || [ -z "${num_goals}" ] ||
       [ "${num_goals}" -eq 0 ] || [ "${num_proved}" -ne "${num_goals}" ]; then
        echo "FAILED: ${functions} in ${file}"
        status=1
    fi
done

exit ${status}
