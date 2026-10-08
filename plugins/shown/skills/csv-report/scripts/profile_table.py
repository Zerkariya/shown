#!/usr/bin/env python3
"""Profile one or more tables so an AI agent can decide the column mapping.

Usage:
    python3 profile_table.py FILE [FILE ...] [--sheet NAME] [--header-row N] [--json]

For every file (and every sheet of an .xlsx) it prints: detected header row,
row/column counts, and per column the inferred type, missing count, unique
count, a value summary, and *hints* about which canonical role the column may
play (patient id, chromosome, start, nMajor, SV type, purity, ...). Hints are
heuristics only - the agent makes the final call when writing mapping.json.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import statistics
import sys
from collections import Counter

sys.dont_write_bytecode = True  # keep the installed skill folder clean (no __pycache__)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tableio  # noqa: E402

CHROMS = {str(i) for i in range(1, 26)} | {"X", "Y", "M", "MT", "XY"}
SVTYPES = {"DEL", "DUP", "INV", "BND", "INS", "TRA", "CNV", "CTX", "ITX", "DELETION",
           "DUPLICATION", "INVERSION", "TRANSLOCATION", "INSERTION", "TANDEM", "TD"}
ISO_DATE = re.compile(r"^\d{4}[-/.]\d{1,2}[-/.]\d{1,2}([ T]\d{1,2}:\d{2}(:\d{2})?)?$")

# (role, name regex) - matched case-insensitively against the column name.
NAME_HINTS = [
    ("patient", r"patient|^pt|^pid$|pid$|case|subject|sample|^id$|_id$|barcode|编号|患者|病人|病例|样本|样品|受试者"),
    ("chrom", r"^#?chr(om(osome)?)?$|^chrom|^chr_?1$|^seqnames?$|^contig$|染色体"),
    ("chrom2", r"chr2|chrom2|mate_?chr|partner_?chr|^chr_?b$"),
    ("start", r"^start|startpos|^pos(ition)?([_.\s]|$)|^begin|^loc|起始|开始|位置|坐标"),
    ("end", r"^end|endpos|^stop|终止|结束|截止"),
    ("pos2", r"pos2|mate_?pos|partner_?pos"),
    ("major_cn", r"nmajor|^major|^n_?a$|^nA$|major_?cn|cn_?major|主等位"),
    ("minor_cn", r"nminor|^minor|^n_?b$|^nB$|minor_?cn|cn_?minor|次等位"),
    ("major_raw", r"naraw|major_?raw"),
    ("minor_raw", r"nbraw|minor_?raw"),
    ("total_cn", r"^cn$|total_?cn|^tcn$|copy_?number|copynumber|^rdcn|_cn$|拷贝数"),
    ("logr", r"logr|log2|^lrr$|seg\.?mean|log_?ratio|logratio"),
    ("baf", r"^baf|b_?allele"),
    ("svtype", r"svtype|sv_?type|^type$|variant_?type|event_?type|变异类型|类型"),
    ("svlen", r"svlen|sv_?len|^length$|^size$|长度|大小"),
    ("qual", r"^qual"),
    ("filter", r"^filter$|^ft$|过滤"),
    ("pe", r"^pe$|paired|discordant|^dv$"),
    ("sr", r"^sr$|split|^rv$"),
    ("genotype", r"^gt$|genotype|基因型"),
    ("purity", r"purity|cellularity|^acf$|aberrant|tumou?r_?content|纯度|肿瘤含量"),
    ("ploidy", r"ploidy|^psi([_.\s]|$)|倍性|倍体"),
    ("goodness_of_fit", r"goodness|^gof([_.\s]|$)|fit_?(quality|score|pct|percent)|拟合优度"),
    ("sex", r"^sex\b|gender|性别"),
]


# direct identifiers: never show as attributes (SKILL.md privacy rule)
PII_RX = (r"name|姓名|名字|^mrn$|medical.?record|病历号|住院号|门诊号|phone|mobile|\btel\b|电话|手机|"
          r"身份证|id.?card|national.?id|passport|^ssn$|address|地址|住址|e-?mail|邮箱|birth.?date|^dob$|出生日期")
PII_EXEMPT = r"gene|sample|file|drug|cancer|disease|tumou?r|histolog|variant|chrom|cell|pathway|signature"
ARM_RX = re.compile(r"^(chr)?(\d{1,2}|X|Y)[pq](_.*)?$", re.I)


def classify(values):
    """Return (type, info) for a list of string cells."""
    present = [v for v in values if not tableio.is_missing(v)]
    info = {"missing": len(values) - len(present), "unique": len(set(present))}
    if not present:
        return "empty", info
    nums = [tableio.parse_number(v) for v in present]
    numeric = [x for x in nums if x is not None]
    if len(numeric) >= 0.95 * len(present):
        ints = all(float(x).is_integer() for x in numeric)
        info.update(min=min(numeric), max=max(numeric), median=statistics.median(numeric))
        if any(str(v).strip().endswith("%") for v in present[:50]):
            info["percent_sign"] = True
        return ("integer" if ints else "float"), info
    if sum(1 for v in present if ISO_DATE.match(v)) >= 0.9 * len(present):
        info.update(min=min(present), max=max(present))
        return "date", info
    low = {v.strip().upper() for v in present}
    if low <= {"TRUE", "FALSE", "YES", "NO", "Y", "N", "0", "1", "是", "否"}:
        return "boolean", info
    counts = Counter(present)
    info["top"] = counts.most_common(8)
    if info["unique"] <= 30 or info["unique"] <= 0.5 * len(present):
        return "category", info
    return "text", info


def hints_for(name, typ, info, values):
    hints = []
    lname = name.lower()
    if typ != "empty" and not re.match(r"^(unnamed|column_)", name, re.I) \
            and re.search(PII_RX, name, re.I) and not re.search(PII_EXEMPT, name, re.I):
        return ["DIRECT IDENTIFIER? (name / MRN / phone / ID number ...) - exclude, never display"]
    for role, rx in NAME_HINTS:
        if re.search(rx, lname, re.I) or re.search(rx, name, re.I):
            hints.append(role)
    present = [v.strip() for v in values if not tableio.is_missing(v)]
    if present:
        normalized = {re.sub(r"^chr", "", v, flags=re.I).upper() for v in present[:2000]}
        chrom_named = "chrom" in hints or any(re.match(r"^chr", v, re.I) or v.upper() in ("X", "Y") for v in present[:5000])
        if chrom_named and typ in ("integer", "category", "text") and 2 <= len(normalized) <= 30 \
                and sum(1 for v in normalized if v in CHROMS) >= 0.9 * len(normalized):
            hints.append("looks like chromosome")
        upper = {v.upper() for v in present[:2000]}
        if upper and sum(1 for v in upper if v.strip("<>") in SVTYPES) >= 0.8 * len(upper):
            hints.append("looks like SV type")
        if all(re.match(r"^<[A-Z:]+>$|^[ACGTN]*[\[\]]", v) for v in present[:200]):
            hints.append("looks like VCF ALT")
        if re.match(r"^(DEL|DUP|INV|BND|INS|CNV)\d{5,}", present[0]):
            hints.append("Delly variant ID (type is the prefix)")
        variant_ids = sum(1 for v in present[:500] if re.match(r"^(sv|del|dup|inv|bnd|ins|tra|cnv|manta|gridss|rs)[\w:.-]*\d", v, re.I))
        if variant_ids >= 0.8 * min(len(present), 500):
            hints = [h for h in hints if h != "patient"] + ["variant / record ID (not a patient ID)"]
        if all(re.match(r"^(chr)?\w+:\d[\d,]*(-\d[\d,]*)?$", v, re.I) for v in present[:200]):
            hints.append("combined locus chr:start-end (split with transforms.regex)")
        if all(re.match(r"^[0-9.][/|][0-9.]$", v) for v in present[:200]):
            hints.append("looks like genotype")
        elif all(re.match(r"^\d+\s*[+:]\s*\d+$", v) for v in present[:200]):
            hints.append("combined allele CN like 2+1 (split with transforms.regex)")
    if typ in ("integer", "float"):
        mx, mn = info["max"], info["min"]
        if typ == "integer" and 1e5 <= mx <= 3.5e8 and mn >= 0:
            hints.append("genomic-position range")
        if 0 <= mn and mx <= 1.05 and typ == "float":
            hints.append("fraction 0-1 (purity/BAF/frequency?)")
        if 1 < mx <= 100 and mn >= 0 and ("purity" in hints or info.get("percent_sign")):
            hints.append("percent 0-100 (divide by 100 for purity)")
        if typ == "float" and -6 <= mn and mx <= 6 and mn < 0:
            hints.append("signed small values (logR?)")
        if typ == "integer" and mn >= -2 and mx <= 2 and mn < 0:
            hints.append("relative call -1/0/+1 (loss/neutral/gain), not a copy number")
    if info["unique"] == len(values) - info["missing"] and len(values) > 1 and typ in ("text", "category", "integer"):
        hints.append("unique per row")
    return hints


def per_id_constant(table, id_col):
    """Columns whose value never changes within one id -> patient-level attributes."""
    idx = table.columns.index(id_col)
    groups = {}
    for r in table.rows:
        groups.setdefault(r[idx], []).append(r)
    if len(groups) == len(table.rows):
        return None
    const = []
    for j, c in enumerate(table.columns):
        if j == idx:
            continue
        if all(len({g[j] for g in rows}) == 1 for rows in groups.values()):
            const.append(c)
    return {"id_column": id_col, "n_ids": len(groups), "constant_columns": const}


def fmt_num(x):
    if x is None:
        return "NA"
    if abs(x) >= 1e5:
        return "%.3g" % x
    if float(x).is_integer():
        return str(int(x))
    return "%.4g" % x


def profile(table, max_values=8):
    cols = []
    for j, name in enumerate(table.columns):
        values = [r[j] for r in table.rows]
        typ, info = classify(values)
        if typ in ("integer", "float"):
            summary = "min %s, median %s, max %s" % (fmt_num(info["min"]), fmt_num(info["median"]), fmt_num(info["max"]))
        elif typ == "date":
            summary = "%s .. %s" % (info["min"], info["max"])
        elif "top" in info:
            summary = ", ".join("%s (%d)" % (v, n) for v, n in info["top"][:max_values])
        else:
            summary = ""
        examples = [v for v in values if not tableio.is_missing(v)][:3]
        cols.append({
            "index": j + 1, "name": name, "type": typ, "missing": info["missing"],
            "unique": info["unique"], "summary": summary, "examples": examples,
            "hints": hints_for(name, typ, info, values),
        })
    id_like = sorted((c for c in cols if "patient" in c["hints"] and c["unique"] >= 2),
                     key=lambda c: (c["unique"], c["index"]))
    id_like = [c["name"] for c in id_like]
    grouping = None
    for cand in id_like:
        g = per_id_constant(table, cand)
        if g:
            grouping = g
            break
    layout = []
    arm_cols = [c for c in table.columns if ARM_RX.match(c)]
    if len(arm_cols) >= 10:
        layout.append("wide layout: %d columns are chromosome arms (%s ...) - one row per patient, one column per arm; "
                      "see SKILL.md 'When the shape does not fit'" % (len(arm_cols), ", ".join(arm_cols[:4])))
    return {
        "file": table.path, "sheet": table.sheet, "format": table.format, "layout_hints": layout,
        "header_row": table.header_row, "skipped_rows_above_header": table.skipped_rows,
        "rows": len(table.rows), "columns": len(table.columns), "column_profiles": cols,
        "grouping": grouping, "meta": table.meta,
        "first_rows": [dict(zip(table.columns, r)) for r in table.rows[:3]],
    }


def print_profile(p):
    title = p["file"] + (" [sheet: %s]" % p["sheet"] if p["sheet"] else "")
    print("=" * 78)
    print(title)
    print("format: %s   rows: %d   columns: %d   header row: %s%s" % (
        p["format"], p["rows"], p["columns"], p["header_row"],
        "   (skipped %d row(s) above header)" % p["skipped_rows_above_header"] if p["skipped_rows_above_header"] else ""))
    if p["meta"].get("samples"):
        print("VCF samples:", ", ".join(p["meta"]["samples"]))
    for hint in p.get("layout_hints", []):
        print("layout: " + hint)
    if p["meta"].get("decimal_comma_columns"):
        print("note: decimal commas (0,55) converted to dots in: " + ", ".join(p["meta"]["decimal_comma_columns"]))
    if p["meta"].get("footer_rows_dropped"):
        print("note: dropped note row(s) at the end of the table: " + " | ".join(p["meta"]["footer_rows_dropped"])[:300])
    print("-" * 78)
    for c in p["column_profiles"]:
        print("%3d. %s  [%s]  missing=%d unique=%d" % (c["index"], c["name"], c["type"], c["missing"], c["unique"]))
        if c["summary"]:
            print("      values: " + c["summary"][:300])
        elif c["examples"]:
            print("      e.g.: " + ", ".join(c["examples"])[:300])
        if c["hints"]:
            print("      hints: " + "; ".join(c["hints"]))
    g = p["grouping"]
    if g:
        print("-" * 78)
        print("rows group by %r -> %d distinct ids; columns constant within an id (patient-level):" % (g["id_column"], g["n_ids"]))
        print("      " + (", ".join(g["constant_columns"]) or "(none)"))
    print()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--sheet", help="xlsx sheet name or 1-based index (default: every sheet)")
    ap.add_argument("--header-row", type=int, help="1-based header row (default: auto-detect)")
    ap.add_argument("--json", action="store_true", help="machine-readable output")
    args = ap.parse_args(argv)

    results = []
    for path in args.files:
        if not os.path.exists(path):
            print("error: %s does not exist" % path, file=sys.stderr)
            return 2
        sheets = [args.sheet] if args.sheet else (tableio.list_sheets(path) or [None])
        for sh in sheets:
            try:
                t = tableio.read_table(path, sheet=sh, header_row=args.header_row)
            except Exception as exc:  # report and continue with other files
                print("error reading %s%s: %s" % (path, " [%s]" % sh if sh else "", exc), file=sys.stderr)
                continue
            results.append(profile(t))
    if args.json:
        json.dump(results, sys.stdout, ensure_ascii=False, indent=1, default=str)
        print()
    else:
        for p in results:
            print_profile(p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
