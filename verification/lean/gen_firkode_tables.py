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
# Generate cert/FirkodeCert/Generated.lean from the FIRKODE C sources:
# the closed-form Radau IIA tables (firk_radau_closed_form), the
# gamma0 literals (firk_gamma0), a few constants, and fingerprints of
# the C routines that are transcribed by hand in the Lean files.
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
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "cert" / "FirkodeCert" / "Generated.lean"

TABLES_C = "src/firkode/firkode_tables.c"
TABLES_H = "include/firkode/firkode_tables.h"
IMPL_H = "src/firkode/firkode_impl.h"
FIRKODE_C = "src/firkode/firkode.c"
FIRKODE_IO_C = "src/firkode/firkode_io.c"

# C routines transcribed by hand in the Lean files. A fingerprint change
# means the transcription must be reviewed before regenerating.
MIRRORED = [
    (TABLES_C, "firk_radau_poly"),
    (TABLES_C, "firk_radau_nodes"),
    (TABLES_C, "firk_lagrange_coeffs"),
    (TABLES_C, "firk_collocation_coeffs"),
    (TABLES_C, "firk_derive_ld"),
    (TABLES_C, "firk_radau_build"),
    (TABLES_C, "FIRKodeTable_CheckOrder"),
    (FIRKODE_C, "FIRKodeGetDky"),
    (FIRKODE_C, "firkPredict"),
    (FIRKODE_C, "firkErrorEstimate"),
    (FIRKODE_IO_C, "FIRKodeSetOrder"),
]

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


LD_CONSTANTS = {"LD_ZERO": "0", "LD_ONE": "1", "LD_TWO": "2"}
DECIMAL = re.compile(r"(\d+)\.(\d*)L")


def ld_literal(token):
    """Translate a long double literal or LD_* constant to an exact Lean numeral."""
    if token in LD_CONSTANTS:
        return LD_CONSTANTS[token]
    m = DECIMAL.fullmatch(token)
    if m is None:
        raise GenError(f"unsupported literal {token!r}")
    whole, frac = m.group(1), m.group(2).rstrip("0")
    return whole if frac == "" else f"{whole}.{frac}"


def parse_gamma0(code, nstages):
    m = re.search(
        r"static const firk_ld firk_gamma0\[FIRK_MAX_STAGES \+ 1\]\s*=\s*\{([^}]*)\};", code
    )
    if m is None:
        raise GenError("missing firk_gamma0 initializer")
    tokens = [tok.strip() for tok in m.group(1).split(",")]
    if len(tokens) != nstages + 1:
        raise GenError(f"firk_gamma0 has {len(tokens)} entries, expected {nstages + 1}")
    values = [ld_literal(tok) for tok in tokens]
    # digits printed after the decimal point, which fix the unit in the last place
    digits = [
        len(DECIMAL.fullmatch(tok).group(2)) if tok not in LD_CONSTANTS else 0 for tok in tokens
    ]
    return values, digits


TOKEN = re.compile(r"\s*(\d+\.\d*L|LD_[A-Z]+|r6|[-+*/()])")


def translate_expr(expr):
    """Translate a closed-form coefficient expression to a Lean Q6 expression.

    The C and Lean expressions have the same structure; unary minus binds
    tighter than the binary operators in both languages."""
    out, pos, prev, prev_unary = "", 0, None, False
    expr = expr.strip()
    while pos < len(expr):
        m = TOKEN.match(expr, pos)
        if m is None:
            raise GenError(f"unsupported token in {expr!r} at {expr[pos:]!r}")
        tok = m.group(1)
        pos = m.end()
        if tok[0].isdigit() or tok.startswith("LD_"):
            lit = ld_literal(tok)
            tok = lit if "." not in lit else f"(ofRat {lit})"
        unary = tok == "-" and prev in (None, "(", "+", "-", "*", "/")
        if out and tok != ")" and prev != "(" and not (prev == "-" and prev_unary):
            out += " "
        out += tok
        prev, prev_unary = tok, unary
    return out


def parse_closed_form(code):
    body = function_definition(code, "firk_radau_closed_form")
    build = " ".join(function_definition(code, "firk_radau_build").split())
    if "if (closed_form && s <= 3) { firk_radau_closed_form(s, c, A); }" not in build:
        raise GenError("firk_radau_build no longer calls firk_radau_closed_form for s <= 3 only")
    switch = body[body.index("switch") :].replace("{", " ").replace("}", " ")
    parts = re.split(r"\b(case\s+\d+\s*:|default\s*:)", switch)
    labels = [" ".join(p.split()) for p in parts[1::2]]
    if labels != ["case 1:", "case 2:", "default:"]:
        raise GenError(f"unexpected switch labels {labels} in firk_radau_closed_form")
    tables = {}
    for s, section in zip([1, 2, 3], parts[2::2]):
        c, A, have_r6 = {}, {}, False
        for stmt in section.split(";"):
            stmt = stmt.strip()
            if stmt in ("", "break"):
                continue
            lhs, sep, rhs = stmt.partition("=")
            if not sep:
                raise GenError(f"unexpected statement {stmt!r} for s = {s}")
            lhs = lhs.strip()
            if lhs == "r6":
                if " ".join(rhs.split()) != "sqrtl(6.0L)" or s != 3:
                    raise GenError(f"unexpected r6 assignment {stmt!r}")
                have_r6 = True
                continue
            mc = re.fullmatch(r"c\[(\d+)\]", lhs)
            ma = re.fullmatch(r"A\[(\d+)\]\[(\d+)\]", lhs)
            if mc and int(mc.group(1)) not in c:
                c[int(mc.group(1))] = translate_expr(rhs)
            elif ma and (int(ma.group(1)), int(ma.group(2))) not in A:
                A[(int(ma.group(1)), int(ma.group(2)))] = translate_expr(rhs)
            else:
                raise GenError(f"unexpected assignment {stmt!r} for s = {s}")
        uses_r6 = any("r6" in e for e in list(c.values()) + list(A.values()))
        if uses_r6 != have_r6:
            raise GenError(f"r6 used without assignment (or vice versa) for s = {s}")
        if sorted(c) != list(range(s)) or sorted(A) != [
            (i, j) for i in range(s) for j in range(s)
        ]:
            raise GenError(f"incomplete closed-form table for s = {s}")
        tables[s] = ([c[i] for i in range(s)], [[A[(i, j)] for j in range(s)] for i in range(s)])
    return tables


def generate():
    tables_c = read_source(TABLES_C)
    tables_h = read_source(TABLES_H)
    impl_h = read_source(IMPL_H)
    sources = {
        TABLES_C: tables_c,
        FIRKODE_C: read_source(FIRKODE_C),
        FIRKODE_IO_C: read_source(FIRKODE_IO_C),
    }

    macros = {}
    macros["FIRK_MAX_STAGES"] = define(tables_h, "FIRK_MAX_STAGES", macros)
    macros["FIRK_MAX_CHECK_ORDER"] = define(tables_c, "FIRK_MAX_CHECK_ORDER", macros)
    macros["FIRK_NODE_GRID_PER_STAGE"] = define(tables_c, "FIRK_NODE_GRID_PER_STAGE", macros)
    macros["FIRK_DEFAULT_STAGES"] = define(impl_h, "FIRK_DEFAULT_STAGES", macros)

    gamma0, gamma0_digits = parse_gamma0(tables_c, macros["FIRK_MAX_STAGES"])
    tables = parse_closed_form(tables_c)
    prints = [(rel, name, fingerprint(sources[rel], name)) for rel, name in MIRRORED]

    lines = [HEADER, "import FirkodeCert.Basic", ""]
    lines += [
        "/-!",
        "# Data generated from the FIRKODE C sources",
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
        ("FIRK_NODE_GRID_PER_STAGE", TABLES_C),
        ("FIRK_DEFAULT_STAGES", IMPL_H),
    ]:
        lean_name = "".join(w.capitalize() for w in name.split("_")[1:])
        lean_name = lean_name[0].lower() + lean_name[1:]
        lines += [f"/-- `{name}` ({where}) -/", f"def {lean_name} : Nat := {macros[name]}", ""]

    lines += ["/-- Nodes `c` of `firk_radau_closed_form` for `s` = 1, 2, 3 stages. -/"]
    lines += ["def closedFormC : Nat → List Q6"]
    for s, (c, _) in tables.items():
        lines += [f"  | {s} => [{', '.join(c)}]"]
    lines += ["  | _ => []", ""]
    lines += ["/-- Coefficients `A` (row major) of `firk_radau_closed_form`. -/"]
    lines += ["def closedFormA : Nat → List (List Q6)"]
    for s, (_, A) in tables.items():
        prefix = f"  | {s} => ["
        rows = (",\n" + " " * len(prefix)).join("[" + ", ".join(row) + "]" for row in A)
        lines += [f"{prefix}{rows}]"]
    lines += ["  | _ => []", ""]
    lines += ["/-- `firk_gamma0[s]` for s = 0, …, FIRK_MAX_STAGES, as exact decimals. -/"]
    lines += ["def gamma0Literal : List Rat :="]
    lines += ["  [" + ",\n   ".join(gamma0) + "]", ""]
    lines += ["/-- Number of digits printed after the decimal point in `firk_gamma0[s]`. -/"]
    lines += ["def gamma0Digits : List Nat :="]
    lines += ["  [" + ", ".join(str(d) for d in gamma0_digits) + "]", ""]
    lines += [
        "/-- SHA-256 of the comment-stripped, whitespace-normalized C definitions",
        "    transcribed by hand in the Lean files. -/",
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
        description="Generate cert/FirkodeCert/Generated.lean from the FIRKODE C sources."
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
        "review the Lean transcription of that C routine before regenerating.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
