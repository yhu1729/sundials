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
# Generate cert/FirkodeCert/Generated.lean from the FIRKODE sources: the
# generated Radau IIA literals of src/firkode/firkode_tables.c (the nodes
# c, the coefficients A and the shifts gamma0 of every built-in table, 40
# significant digits), the closed forms of the tables for s <= 3, a few
# constants, and fingerprints of the C routines and of the functions of
# scripts/firkode_radau_tables.py (which computes the literals) that are
# transcribed by hand in the Lean files.
#
# Usage:
#   python3 gen_firkode_tables.py          write Generated.lean
#   python3 gen_firkode_tables.py --check  exit 1 if Generated.lean is
#                                          stale (including a changed
#                                          fingerprint) or the cert and
#                                          spec toolchains differ
# ---------------------------------------------------------------

import argparse
import ast
import difflib
import hashlib
import re
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "cert" / "FirkodeCert" / "Generated.lean"

TABLES_C = "src/firkode/firkode_tables.c"
TABLES_H = "include/firkode/firkode_tables.h"
IMPL_H = "src/firkode/firkode_impl.h"
FIRKODE_H = "include/firkode/firkode.h"
FIRKODE_C = "src/firkode/firkode.c"
DENSE_C = "src/firkode/firkode_dense.c"
STEP_C = "src/firkode/firkode_step.c"
FIRKODE_IO_C = "src/firkode/firkode_io.c"
GENERATOR_PY = "scripts/firkode_radau_tables.py"

# C routines transcribed by hand in the Lean files. A fingerprint change
# means the transcription must be reviewed before regenerating.
MIRRORED_C = [
    (TABLES_C, "FIRKodeTable_CheckOrder"),
    (FIRKODE_C, "FIRKodeGetDky"),
    (DENSE_C, "firkDenseEval"),
    (DENSE_C, "firkPredict"),
    (STEP_C, "firkErrorEstimateBody"),
    (FIRKODE_IO_C, "firkOrderToStages"),
    (FIRKODE_IO_C, "FIRKodeSetOrder"),
]

# Functions of the table generator transcribed by hand in Mirror.lean (the
# first implementation computed the tables in C at run time; the script
# computes the same quantities with mpmath and writes them as literals).
MIRRORED_PY = [
    (GENERATOR_PY, "radau_nodes"),
    (GENERATOR_PY, "lagrange_coeffs"),
    (GENERATOR_PY, "radau_table"),
    (GENERATOR_PY, "tables"),
]

# The grid density of the bracket scan of Mirror.radauGridScan, a device of
# the Lean proofs (the generator script scans 400 points per stage with the
# same grid formula; the brackets of either scan are valid).
NODE_GRID_PER_STAGE = 40

# The closed forms of the Radau IIA tables for s = 1, 2, 3 (Hairer and
# Wanner, Solving ODEs II, Table IV.5.6), the mathematical reference that
# the literals of those tables are certified against; r6 is sqrt(6).
CLOSED_FORM_C = {1: ["1"], 2: ["1 / 3", "1"], 3: ["(4 - r6) / 10", "(4 + r6) / 10", "1"]}
CLOSED_FORM_A = {
    1: [["1"]],
    2: [["5 / 12", "-1 / 12"], ["3 / 4", "1 / 4"]],
    3: [
        ["(88 - 7 * r6) / 360", "(296 - 169 * r6) / 1800", "(-2 + 3 * r6) / 225"],
        ["(296 + 169 * r6) / 1800", "(88 + 7 * r6) / 360", "(-2 - 3 * r6) / 225"],
        ["(16 - r6) / 36", "(16 + r6) / 36", "1 / 9"],
    ],
}

HEADER = """\
/- -----------------------------------------------------------------
   Programmer(s): Yifan Hu @ UMBC
   -----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025-2026, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------- -/
"""


class GenError(Exception):
    pass


def strip_comments(text):
    """Remove C comments, keeping string and character literals intact."""
    out = []
    i, n = 0, len(text)
    while i < n:
        ch = text[i]
        if ch in "\"'":
            j = i + 1
            while j < n and text[j] != ch:
                j += 2 if text[j] == "\\" else 1
            out.append(text[i : j + 1])
            i = j + 1
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            if j < 0:
                raise GenError("unterminated comment")
            out.append(" ")
            i = j + 2
        elif text.startswith("//", i):
            j = text.find("\n", i)
            i = n if j < 0 else j
        else:
            out.append(ch)
            i += 1
    return "".join(out)


def read_source(rel):
    path = ROOT / rel
    if not path.is_file():
        raise GenError(f"missing source file {rel}")
    return strip_comments(path.read_text())


def function_definition(code, name):
    """Return the definition of a top-level C function (signature and body)."""
    pattern = re.compile(
        r"^[A-Za-z_][\w \t*]*\b" + re.escape(name) + r"\s*\([^;{]*\)\s*\{", re.MULTILINE
    )
    matches = list(pattern.finditer(code))
    if len(matches) != 1:
        raise GenError(f"expected one definition of {name}, found {len(matches)}")
    start = matches[0].start()
    depth = 0
    for i in range(matches[0].end() - 1, len(code)):
        if code[i] == "{":
            depth += 1
        elif code[i] == "}":
            depth -= 1
            if depth == 0:
                return code[start : i + 1]
    raise GenError(f"unbalanced braces in {name}")


def fingerprint(code, name):
    normalized = " ".join(function_definition(code, name).split())
    return hashlib.sha256(normalized.encode()).hexdigest()


def eval_int_expr(expr, macros):
    """Evaluate an integer macro expression built from + - * and parentheses."""

    def ev(node):
        if isinstance(node, ast.Expression):
            return ev(node.body)
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return node.value
        if isinstance(node, ast.Name) and node.id in macros:
            return macros[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -ev(node.operand)
        if isinstance(node, ast.BinOp) and type(node.op) in (ast.Add, ast.Sub, ast.Mult):
            a, b = ev(node.left), ev(node.right)
            return {ast.Add: a + b, ast.Sub: a - b, ast.Mult: a * b}[type(node.op)]
        raise GenError(f"unsupported macro expression {expr!r}")

    return ev(ast.parse(expr, mode="eval"))


def define(code, name, macros):
    m = re.search(r"^#define\s+" + name + r"\s+(.+)$", code, re.MULTILINE)
    if m is None:
        raise GenError(f"missing #define {name}")
    return eval_int_expr(m.group(1).strip(), macros)


def python_function_source(rel, name):
    """The source of a top-level function of a Python file."""
    path = ROOT / rel
    if not path.is_file():
        raise GenError(f"missing source file {rel}")
    source = path.read_text()
    tree = ast.parse(source)
    defs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name]
    if len(defs) != 1:
        raise GenError(f"expected one definition of {name} in {rel}, found {len(defs)}")
    return ast.get_source_segment(source, defs[0])


def fingerprint_text(text):
    return hashlib.sha256(" ".join(text.split()).encode()).hexdigest()


LITERAL = re.compile(r"SUN_RCONST\(([-+0-9.eE]+)\)")
TABLE = re.compile(
    r"static const struct FIRKodeTableMem firk_radau_iia_(\d+) = \{(.*?)\.owned = SUNFALSE\};",
    re.DOTALL,
)


def field(block, name):
    """The SUN_RCONST literals of the initializer of a field."""
    m = re.search(r"\." + name + r"\s*=\s*(\{[^}]*\}|SUN_RCONST\([^)]*\))", block)
    if m is None:
        raise GenError(f"missing field {name}")
    return LITERAL.findall(m.group(1))


def parse_literals(code, nstages):
    """(c, A, gamma0) of every built-in table, as decimal strings."""
    tables = {}
    for m in TABLE.finditer(code):
        s, block = int(m.group(1)), m.group(2)
        c = field(block, "c")
        A = field(block, "A")
        g = field(block, "gamma0")
        if len(c) != s or len(A) != s * s or len(g) != 1:
            raise GenError(f"table {s}: {len(c)} nodes, {len(A)} coefficients, {len(g)} shifts")
        tables[s] = (c, [A[i * s : (i + 1) * s] for i in range(s)], g[0])
    if sorted(tables) != list(range(1, nstages + 1)):
        raise GenError(f"tables found for s = {sorted(tables)}, expected 1..{nstages}")
    return tables


def lean_rat(token):
    """A decimal literal of the C sources as a Lean rational literal."""
    value = Fraction(token)
    if value.denominator == 1:
        return str(value.numerator)
    if "e" in token.lower():
        raise GenError(f"unexpected exponent in {token}")
    sign = "-" if token.startswith("-") else ""
    digits = token.lstrip("+-")
    return sign + digits


def digits_after_point(token):
    digits = token.lstrip("+-")
    return len(digits.split(".")[1]) if "." in digits else 0


def generate():
    tables_c = read_source(TABLES_C)
    tables_h = read_source(TABLES_H)
    impl_h = read_source(IMPL_H)
    firkode_h = read_source(FIRKODE_H)
    sources = {rel: read_source(rel) for rel, _ in MIRRORED_C}

    macros = {}
    macros["FIRK_MAX_STAGES"] = define(tables_h, "FIRK_MAX_STAGES", macros)
    macros["FIRK_MAX_CHECK_ORDER"] = define(tables_c, "FIRK_MAX_CHECK_ORDER", macros)
    macros["FIRK_DEFAULT_STAGES"] = define(firkode_h, "FIRK_DEFAULT_STAGES", macros)
    macros["FIRK_SUPPORTED_STAGES"] = define(impl_h, "FIRK_SUPPORTED_STAGES", macros)
    nstages = macros["FIRK_MAX_STAGES"]

    literals = parse_literals(tables_c, nstages)
    prints = [(rel, name, fingerprint(sources[rel], name)) for rel, name in MIRRORED_C]
    prints += [
        (rel, name, fingerprint_text(python_function_source(rel, name)))
        for rel, name in MIRRORED_PY
    ]

    lines = [HEADER, "import FirkodeCert.Basic", ""]
    lines += [
        "/-!",
        "# Data generated from the FIRKODE sources",
        "",
        "GENERATED by `verification/lean/gen_firkode_tables.py`. Do not edit;",
        "rerun the generator instead.",
        "-/",
        "",
        "namespace Firkode.Generated",
        "",
        "open Q6",
        "",
    ]
    for name, where in [
        ("FIRK_MAX_STAGES", TABLES_H),
        ("FIRK_MAX_CHECK_ORDER", TABLES_C),
        ("FIRK_DEFAULT_STAGES", FIRKODE_H),
        ("FIRK_SUPPORTED_STAGES", IMPL_H),
    ]:
        lean_name = "".join(w.capitalize() for w in name.split("_")[1:])
        lean_name = lean_name[0].lower() + lean_name[1:]
        lines += [f"/-- `{name}` ({where}) -/", f"def {lean_name} : Nat := {macros[name]}", ""]
    lines += [
        "/-- Grid density of the bracket scan `Mirror.radauGridScan`, a device of the",
        "    proofs (`radau_nodes` of the generator script scans 400 points per stage",
        "    on the same grid). -/",
        f"def nodeGridPerStage : Nat := {NODE_GRID_PER_STAGE}",
        "",
    ]
    lines += [
        "/-- Nodes `c` of the closed-form Radau IIA tables for `s` = 1, 2, 3 stages",
        "    (Hairer and Wanner, Table IV.5.6), the reference for the literals. -/",
    ]
    lines += ["def closedFormC : Nat → List Q6"]
    for s, c in CLOSED_FORM_C.items():
        lines += [f"  | {s} => [{', '.join(c)}]"]
    lines += ["  | _ => []", ""]
    lines += ["/-- Coefficients `A` (row major) of the closed-form tables. -/"]
    lines += ["def closedFormA : Nat → List (List Q6)"]
    for s, A in CLOSED_FORM_A.items():
        prefix = f"  | {s} => ["
        rows = (",\n" + " " * len(prefix)).join("[" + ", ".join(row) + "]" for row in A)
        lines += [f"{prefix}{rows}]"]
    lines += ["  | _ => []", ""]
    lines += [
        "/-- The nodes `c` of the built-in table `firk_radau_iia_s` of",
        f"    `{TABLES_C}`, as exact decimals. -/",
    ]
    lines += ["def literalC : Nat → List Rat"]
    for s in range(1, nstages + 1):
        lines += [f"  | {s} => [{', '.join(lean_rat(t) for t in literals[s][0])}]"]
    lines += ["  | _ => []", ""]
    lines += [
        "/-- The coefficients `A` (row major) of the built-in tables, as exact",
        "    decimals. -/",
    ]
    lines += ["def literalA : Nat → List (List Rat)"]
    for s in range(1, nstages + 1):
        prefix = f"  | {s} => ["
        rows = (",\n" + " " * len(prefix)).join(
            "[" + ", ".join(lean_rat(t) for t in row) + "]" for row in literals[s][1]
        )
        lines += [f"{prefix}{rows}]"]
    lines += ["  | _ => []", ""]
    lines += [
        "/-- The shift `gamma0` of the built-in tables for s = 0, …, FIRK_MAX_STAGES",
        "    (0 at s = 0), as exact decimals. -/",
    ]
    lines += ["def gamma0Literal : List Rat :="]
    lines += [
        "  [0,\n   "
        + ",\n   ".join(lean_rat(literals[s][2]) for s in range(1, nstages + 1))
        + "]",
        "",
    ]
    lines += ["/-- Number of digits printed after the decimal point in each `gamma0`. -/"]
    lines += ["def gamma0Digits : List Nat :="]
    lines += [
        "  [0, "
        + ", ".join(str(digits_after_point(literals[s][2])) for s in range(1, nstages + 1))
        + "]",
        "",
    ]
    lines += [
        "/-- SHA-256 of the comment-stripped, whitespace-normalized definitions",
        "    (C routines, and functions of the table generator) transcribed by hand",
        "    in the Lean files. -/",
        "def mirroredFingerprints : List (String × String × String) :=",
    ]
    entries = [f'("{rel}", "{name}",\n    "{h}")' for rel, name, h in prints]
    lines += ["  [" + ",\n   ".join(entries) + "]", ""]
    lines += ["end Firkode.Generated", ""]
    return "\n".join(lines)


def check_toolchains():
    cert = HERE / "cert" / "lean-toolchain"
    spec = HERE / "spec" / "lean-toolchain"
    if spec.is_file() and spec.read_text() != cert.read_text():
        raise GenError("cert/lean-toolchain and spec/lean-toolchain differ")


def main():
    parser = argparse.ArgumentParser(
        description="Generate cert/FirkodeCert/Generated.lean from the FIRKODE sources."
    )
    parser.add_argument("--check", action="store_true", help="check instead of writing")
    args = parser.parse_args()
    try:
        text = generate()
        check_toolchains()
    except GenError as err:
        print(f"gen_firkode_tables.py: error: {err}", file=sys.stderr)
        return 1
    if not args.check:
        OUT.write_text(text)
        return 0
    current = OUT.read_text() if OUT.is_file() else ""
    if current == text:
        return 0
    diff = difflib.unified_diff(
        current.splitlines(), text.splitlines(), "committed", "generated", lineterm=""
    )
    print("\n".join(diff))
    print(
        "gen_firkode_tables.py: Generated.lean is out of date. If a fingerprint changed,\n"
        "review the Lean transcription of that routine before regenerating.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
