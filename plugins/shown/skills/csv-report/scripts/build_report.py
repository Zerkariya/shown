#!/usr/bin/env python3
"""Build a self-contained, offline HTML report from tables + a column mapping.

Usage:
    python3 build_report.py MAPPING.json [-o report.html] [--check] [--stats-json FILE]

MAPPING.json tells the builder which file holds what and which source column
plays which canonical role (see references/mapping.md). Relative file paths are
resolved against the directory of MAPPING.json.

--check      validate the mapping and print the summary, but do not write HTML
--stats-json write the cohort summary as JSON (useful for writing "insights")

Exit codes: 0 ok, 1 mapping/data error (message says what to fix), 2 usage error.
"""
from __future__ import annotations

import argparse
import datetime as dt
import difflib
import glob
import json
import math
import os
import re
import sys
from collections import Counter, OrderedDict, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tableio  # noqa: E402

TEMPLATE = os.path.join(HERE, "..", "assets", "report_template.html")
VERSION = "1.0.0"

GENOMES = {
    # 1..22, X, Y lengths; X/Y lengths for hg19/hg38/CHM13 (ASCAT ships 1..22,X for all three)
    "hg19": [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431,
             135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
             59128983, 63025520, 48129895, 51304566, 155270560, 59373566],
    "hg38": [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
             133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
             58617616, 64444167, 46709983, 50818468, 156040895, 57227415],
    "chm13": [248387328, 242696752, 201105948, 193574945, 182045439, 172126628, 160567428, 146259331, 150617247,
              134758134, 135127769, 133324548, 113566686, 101161492, 99753195, 96330374, 84276897, 80542538,
              61707364, 66210255, 45090682, 51324926, 154259566, 62460029],
}
GENOME_ALIASES = {"grch37": "hg19", "b37": "hg19", "hg19": "hg19", "grch38": "hg38", "hg38": "hg38",
                  "chm13": "chm13", "t2t": "chm13", "t2t-chm13": "chm13", "hs1": "chm13"}
HUMAN_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
SEX_CHROMS = {"X", "Y"}

FIELDS = {
    "patients": {"required": ["patient"],
                 "optional": ["purity", "ploidy", "goodness_of_fit", "sex"]},
    "segments": {"required": ["patient", "chrom", "start", "end"],
                 "optional": ["major_cn", "minor_cn", "total_cn", "major_raw", "minor_raw", "total_raw",
                              "logr", "baf", "n_markers"]},
    "variants": {"required": ["chrom", "pos"],
                 "optional": ["patient", "end", "chrom2", "pos2", "svtype", "svlen", "qual", "filter", "pe", "sr",
                              "genotype", "id", "precise", "alt", "cn"]},
    "bins": {"required": ["chrom"],
             "optional": ["patient", "pos", "start", "end", "logr", "baf", "cn"]},
}
NUMERIC_FIELDS = {"purity", "ploidy", "goodness_of_fit", "start", "end", "pos", "pos2", "major_cn", "minor_cn",
                  "total_cn", "major_raw", "minor_raw", "total_raw", "logr", "baf", "n_markers", "svlen", "qual",
                  "pe", "sr", "cn"}
SVTYPE_SYNONYMS = {
    "DEL": "DEL", "DELETION": "DEL", "LOSS": "DEL", "缺失": "DEL",
    "DUP": "DUP", "DUPLICATION": "DUP", "TANDEM": "DUP", "TD": "DUP", "DUP:TANDEM": "DUP", "GAIN": "DUP", "重复": "DUP",
    "INV": "INV", "INVERSION": "INV", "倒位": "INV",
    "BND": "BND", "TRA": "BND", "CTX": "BND", "ITX": "BND", "TRANSLOCATION": "BND", "BREAKEND": "BND", "易位": "BND",
    "INS": "INS", "INSERTION": "INS", "INS:ME": "INS", "插入": "INS",
    "CNV": "CNV", "COPY_NUMBER": "CNV",
}
SEX_SYNONYMS = {"XX": "XX", "F": "XX", "FEMALE": "XX", "女": "XX", "女性": "XX", "W": "XX",
                "XY": "XY", "M": "XY", "MALE": "XY", "男": "XY", "男性": "XY"}
ISO_DATE = re.compile(r"^(\d{4})[-/.](\d{1,2})[-/.](\d{1,2})")


class MappingError(Exception):
    pass


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(s))]


def norm_chrom(value):
    if value is None:
        return None
    s = str(value).strip()
    if not s or s.lower() in tableio.MISSING:
        return None
    s = re.sub(r"^chr(om(osome)?)?[_ ]?", "", s, flags=re.I)
    if re.match(r"^\d+\.0+$", s):
        s = s.split(".")[0]
    u = s.upper()
    if u == "23":
        return "X"
    if u == "24":
        return "Y"
    if u in ("M", "MT", "25"):
        return "MT"
    if u in ("X", "Y"):
        return u
    if s.isdigit():
        return s.lstrip("0") or "0"
    return s


def norm_svtype(value):
    if value is None:
        return None
    s = str(value).strip().strip("<>").upper()
    if not s or s.lower() in tableio.MISSING:
        return None
    if s in SVTYPE_SYNONYMS:
        return SVTYPE_SYNONYMS[s]
    head = re.split(r"[:_\s]", s)[0]
    return SVTYPE_SYNONYMS.get(head, s)


def r4(x):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return None
    if isinstance(x, float):
        if x.is_integer() and abs(x) < 1e15:
            return int(x)
        return round(x, 4)
    return x


def apply_ops(value, ops, field):
    """Apply a transform pipeline to one cell. Returns a string or number (or None)."""
    if ops is None:
        return value
    if isinstance(ops, dict):
        ops = [ops]
    for op in ops:
        if value is None:
            return None
        kind = op.get("op")
        if kind == "regex":
            m = re.search(op["pattern"], str(value))
            value = m.group(op.get("group", 1)) if m else None
        elif kind == "map":
            table = op.get("values", {})
            key = str(value).strip()
            if key in table:
                value = table[key]
            elif key.lower() in {k.lower(): k for k in table}:
                value = table[{k.lower(): k for k in table}[key.lower()]]
            elif "default" in op:
                value = op["default"]
        elif kind == "replace":
            value = str(value).replace(op.get("old", ""), op.get("new", ""))
        elif kind == "strip_prefix":
            v = str(value)
            pre = op.get("value", "")
            value = v[len(pre):] if v.lower().startswith(pre.lower()) else v
        elif kind in ("multiply", "divide", "add", "log2", "round", "abs", "negate"):
            x = value if isinstance(value, (int, float)) else tableio.parse_number(value)
            if x is None:
                return None
            if kind == "multiply":
                x = x * float(op["value"])
            elif kind == "divide":
                x = x / float(op["value"])
            elif kind == "add":
                x = x + float(op["value"])
            elif kind == "log2":
                x = math.log2(x) if x > 0 else None
            elif kind == "round":
                x = round(x, int(op.get("digits", 0)))
            elif kind == "abs":
                x = abs(x)
            elif kind == "negate":
                x = -x
            value = x
        else:
            raise MappingError("unknown transform op %r for field %r (allowed: regex, map, replace, strip_prefix, "
                               "multiply, divide, add, log2, round, abs, negate)" % (kind, field))
    return value


# --------------------------------------------------------------------------- #
# loading
# --------------------------------------------------------------------------- #
def resolve_files(spec, base):
    files = spec if isinstance(spec, list) else [spec]
    out = []
    for f in files:
        p = f if os.path.isabs(f) else os.path.join(base, f)
        hits = sorted(glob.glob(p), key=natural_key) if any(ch in p for ch in "*?[") else [p]
        for h in hits:
            if not os.path.exists(h):
                raise MappingError("file not found: %s" % h)
            out.append(h)
    if not out:
        raise MappingError("no files match %r" % spec)
    return out


def column_spec(spec):
    """Normalise a 'columns' entry to (column, constant, filename_regex)."""
    if spec is None:
        return None, None, None
    if isinstance(spec, str):
        return spec, None, None
    if isinstance(spec, dict):
        return spec.get("column"), spec.get("value"), spec.get("from_filename")
    raise MappingError("column spec must be a column name or an object, got %r" % (spec,))


def check_columns(table, needed, label):
    missing = [c for c in needed if c not in table.columns]
    if not missing:
        return
    lines = []
    for c in missing:
        close = difflib.get_close_matches(c, table.columns, n=3, cutoff=0.4)
        lines.append("  - %r not found%s" % (c, (" (did you mean %s?)" % ", ".join(repr(x) for x in close)) if close else ""))
    raise MappingError("%s: unknown column(s):\n%s\n  available columns: %s" % (label, "\n".join(lines), ", ".join(table.columns)))


def row_passes(row, idx, filters):
    for f in filters:
        v = row[idx[f["column"]]]
        op = f.get("op", "in")
        val = f.get("value", f.get("values"))
        num = tableio.parse_number(v)
        if op == "in" and str(v) not in [str(x) for x in val]:
            return False
        if op == "not_in" and str(v) in [str(x) for x in val]:
            return False
        if op == "equals" and str(v) != str(val):
            return False
        if op == "not_equals" and str(v) == str(val):
            return False
        if op == "nonempty" and tableio.is_missing(v):
            return False
        if op in ("gt", "ge", "lt", "le"):
            if num is None:
                return False
            x = float(val)
            if (op == "gt" and not num > x) or (op == "ge" and not num >= x) or \
               (op == "lt" and not num < x) or (op == "le" and not num <= x):
                return False
    return True


def load_table_spec(spec, base, ti):
    kind = spec.get("kind")
    if kind not in FIELDS:
        raise MappingError("tables[%d]: kind must be one of %s (got %r)" % (ti, ", ".join(FIELDS), kind))
    columns = spec.get("columns") or {}
    allowed = set(FIELDS[kind]["required"]) | set(FIELDS[kind]["optional"])
    unknown = [k for k in columns if k not in allowed]
    if unknown:
        raise MappingError("tables[%d] (%s): unknown canonical field(s) %s; allowed: %s"
                           % (ti, kind, unknown, ", ".join(sorted(allowed))))
    transforms = spec.get("transforms") or {}
    out_rows, sources = [], []
    for path in resolve_files(spec.get("file") or spec.get("files"), base):
        t = tableio.read_table(path, sheet=spec.get("sheet"), header_row=spec.get("header_row"),
                               delimiter=spec.get("delimiter"))
        label = "tables[%d] %s%s" % (ti, os.path.basename(path), " [%s]" % t.sheet if t.sheet else "")
        specs = {f: column_spec(s) for f, s in columns.items()}
        needed = [c for c, _, _ in specs.values() if c]
        needed += [f["column"] for f in spec.get("row_filter", [])]
        needed += list(spec.get("attributes") or [])
        needed += list(spec.get("exclude") or [])
        check_columns(t, needed, label)
        idx = {c: i for i, c in enumerate(t.columns)}
        for f in FIELDS[kind]["required"]:
            if f not in specs:
                raise MappingError("%s: required field %r is not mapped (kind %s needs %s)"
                                   % (label, f, kind, ", ".join(FIELDS[kind]["required"])))
        consts = {}
        for f, (col, const, fn_rx) in specs.items():
            if col is None and const is None and fn_rx:
                m = re.search(fn_rx, os.path.basename(path))
                if not m:
                    raise MappingError("%s: from_filename regex %r does not match the file name" % (label, fn_rx))
                consts[f] = m.group(1) if m.groups() else m.group(0)
            elif col is None and const is not None:
                consts[f] = const
            elif col is None:
                raise MappingError("%s: field %r needs 'column', 'value' or 'from_filename'" % (label, f))
        used = {c for c, _, _ in specs.values() if c}
        attr_cols = list(spec.get("attributes") or [])
        if kind == "patients" and "attributes" not in spec:
            skip = set(spec.get("exclude") or [])
            skip |= {specs[f][0] for f in ("patient", "purity", "ploidy", "goodness_of_fit") if f in specs and specs[f][0]}
            attr_cols = [c for c in t.columns if c not in skip]
        extra_cols = []
        if kind == "variants":
            skip = used | set(attr_cols) | set(spec.get("exclude") or [])
            extra_cols = [c for c in t.columns if c not in skip][:12]
        filters = spec.get("row_filter") or []
        n_in, n_filtered = len(t.rows), 0
        for row in t.rows:
            if filters and not row_passes(row, idx, filters):
                n_filtered += 1
                continue
            rec = {}
            for f, (col, _, _) in specs.items():
                raw = consts[f] if f in consts else row[idx[col]]
                if not isinstance(raw, (int, float)) and tableio.is_missing(raw):
                    raw = None
                v = apply_ops(raw, transforms.get(f), f)
                if f in NUMERIC_FIELDS and v is not None and not isinstance(v, (int, float)):
                    v = tableio.parse_number(v)
                rec[f] = v
            rec["_attrs"] = {c: row[idx[c]] for c in attr_cols}
            if extra_cols:
                rec["_extra"] = [row[idx[c]] for c in extra_cols]
            out_rows.append(rec)
        sources.append({"table": ti, "kind": kind, "file": os.path.relpath(path, base), "sheet": t.sheet,
                        "format": t.format, "header_row": t.header_row, "rows": n_in, "filtered_out": n_filtered,
                        "columns": t.columns,
                        "mapping": {f: (specs[f][0] if specs[f][0] else ("=" + str(consts.get(f))))
                                    for f in specs},
                        "transforms": {f: transforms[f] for f in transforms if f in specs},
                        "attributes": attr_cols, "extra_columns": extra_cols,
                        "unused_columns": [c for c in t.columns if c not in used and c not in attr_cols
                                           and c not in extra_cols]})
    return kind, out_rows, sources, (extra_cols if kind == "variants" else [])


# --------------------------------------------------------------------------- #
# attribute typing
# --------------------------------------------------------------------------- #
def infer_attr_type(values):
    present = [v for v in values if v is not None and not tableio.is_missing(v)]
    if not present:
        return "category"
    nums = [tableio.parse_number(v) for v in present]
    if all(x is not None for x in nums):
        # small integer codes (0/1 flags, grade 1-4) read better as categories
        codes = all(float(x).is_integer() for x in nums) and len(set(nums)) <= 6 and max(nums) - min(nums) <= 10
        return "category" if codes else "numeric"
    if sum(1 for v in present if ISO_DATE.match(str(v))) >= 0.9 * len(present):
        return "date"
    uniq = len(set(present))
    if uniq <= 30 or uniq <= 0.5 * len(present):
        return "category"
    return "text"


def coerce_attr(value, typ):
    if value is None or tableio.is_missing(value):
        return None
    if typ == "numeric":
        return r4(tableio.parse_number(value))
    if typ == "date":
        m = ISO_DATE.match(str(value))
        if m:
            return "%04d-%02d-%02d" % (int(m.group(1)), int(m.group(2)), int(m.group(3)))
        return None
    if typ in ("category", "ordinal"):
        s = str(value).strip()
        if re.match(r"^-?\d+\.0+$", s):
            s = s.split(".")[0]
        return s
    return str(value).strip()


# --------------------------------------------------------------------------- #
# genome
# --------------------------------------------------------------------------- #
def choose_genome(setting, maxpos, warnings):
    present = set(maxpos)
    if isinstance(setting, dict) and "chroms" in setting:
        chroms = [[str(c), int(l)] for c, l in setting["chroms"]]
        return setting.get("name", "custom"), chroms
    name = GENOME_ALIASES.get(str(setting or "auto").lower())
    human = [c for c in present if c in HUMAN_CHROMS]
    if name is None and present and len(human) < 0.5 * len(present):
        chroms = [[c, int(maxpos[c])] for c in sorted(present, key=natural_key)]
        warnings.append("chromosome names are not human (1-22, X, Y); using the largest coordinate per chromosome "
                        "as its length")
        return "from-data", chroms
    if name is None:
        fits = []
        for build in ("hg38", "hg19", "chm13"):
            lens = dict(zip(HUMAN_CHROMS, GENOMES[build]))
            if all(maxpos[c] <= lens[c] * 1.0005 for c in human):
                fits.append(build)
        name = fits[0] if fits else "hg38"
        if len(fits) != 1:
            warnings.append("genome build not given; positions fit %s - assumed %s (set \"genome\" in the mapping "
                            "to be explicit)" % (", ".join(fits) or "no human build", name))
    lens = GENOMES[name]
    chroms = [[c, l] for c, l in zip(HUMAN_CHROMS, lens) if c != "Y" or "Y" in present]
    over = [c for c in human if maxpos[c] > dict(chroms).get(c, 0) * 1.0005]
    if over:
        warnings.append("positions exceed %s chromosome lengths on %s - wrong genome build or position units?"
                        % (name, ", ".join(sorted(over, key=natural_key))))
    return name, chroms


# --------------------------------------------------------------------------- #
# ASCAT-style metrics (ascat.metrics definitions, computed from segments)
# --------------------------------------------------------------------------- #
def ascat_metrics(segs):
    """segs: list of dicts with chrom, start, end, major, minor (ints) or total."""
    out = {}
    if not segs:
        return out
    size = lambda s: s["end"] - s["start"] + 1
    allele = all(s["major"] is not None and s["minor"] is not None for s in segs)
    total_bp = sum(size(s) for s in segs)
    out["n_segs"] = len(segs)
    tot = lambda s: (s["major"] + s["minor"]) if allele else s["total"]
    with_tot = [s for s in segs if tot(s) is not None]
    if with_tot:
        out["ploidy_from_segments"] = round(sum(size(s) * tot(s) for s in with_tot) / sum(size(s) for s in with_tot), 4)
        hd = [s for s in with_tot if tot(s) == 0]
        out["homdel_segs"] = len(hd)
        out["homdel_largest"] = max((size(s) for s in hd), default=0)
        out["homdel_size"] = sum(size(s) for s in hd)
        out["homdel_fraction"] = round(out["homdel_size"] / total_bp, 4) if total_bp else None
    auto = [s for s in with_tot if s["chrom"] not in SEX_CHROMS]
    auto_bp = sum(size(s) for s in auto)
    if not auto or not auto_bp:
        return out

    def mode(vals):
        acc = defaultdict(float)
        for v, s in vals:
            acc[min(5, int(round(v)))] += size(s) / 1e6
        return max(sorted(acc), key=lambda k: acc[k])

    if allele:
        out["LOH"] = round(sum(size(s) for s in auto if s["minor"] == 0) / auto_bp, 4)
        out["mode_minA"] = mode([(s["minor"], s) for s in auto])
        mj = mode([(s["major"], s) for s in auto])
        out["mode_majA"] = mj
        wgd = None
        if mj == 1:
            wgd = 0
        elif mj == 2:
            wgd = 1
        elif mj in (3, 4, 5):
            wgd = "1+"
        out["WGD"] = wgd
        if wgd is not None:
            base = 1 if wgd == 0 else 2
            same = sum(size(s) for s in auto if s["major"] == base and s["minor"] == base)
            out["GI"] = round(1 - same / auto_bp, 4)
    totals = defaultdict(float)
    for s in auto:
        totals[int(round(tot(s)))] += size(s)
    base_tot = max(sorted(totals), key=lambda k: totals[k])
    out["FGA"] = round(1 - totals[base_tot] / auto_bp, 4)
    return out


# --------------------------------------------------------------------------- #
# build
# --------------------------------------------------------------------------- #
def build(mapping, base):
    warnings = []
    tables = mapping.get("tables")
    if not tables:
        raise MappingError("mapping has no 'tables' list")
    loaded = defaultdict(list)
    sources = []
    sv_extra_names = []
    for ti, spec in enumerate(tables):
        kind, rows, srcs, extra = load_table_spec(spec, base, ti)
        loaded[kind].append((ti, rows, spec))
        sources += srcs
        if extra and not sv_extra_names:
            sv_extra_names = extra
        elif extra and extra != sv_extra_names:
            warnings.append("tables[%d]: variant extra columns differ from the first variants table; only the "
                            "first table's extra columns are shown" % ti)

    # ---------------- patients registry ----------------
    patients = OrderedDict()

    def pid_of(rec):
        p = rec.get("patient")
        if p is None:
            return None
        if isinstance(p, float) and p.is_integer():
            p = int(p)
        return str(p).strip() or None

    def ensure(pid):
        if pid not in patients:
            patients[pid] = {"id": pid, "attrs_raw": {}, "purity": None, "ploidy": None,
                             "goodness_of_fit": None, "sex": None, "sources": set()}
        return patients[pid]

    attr_order = []
    for ti, rows, spec in loaded["patients"]:
        seen = set()
        for rec in rows:
            pid = pid_of(rec)
            if pid is None:
                continue
            if pid in seen:
                warnings.append("tables[%d]: patient %r appears more than once in a patients table; first row "
                                "kept" % (ti, pid))
                continue
            seen.add(pid)
            p = ensure(pid)
            p["sources"].add("patients")
            for f in ("purity", "ploidy", "goodness_of_fit", "sex"):
                if rec.get(f) is not None and p[f] is None:
                    p[f] = rec[f]
            for k, v in rec["_attrs"].items():
                if k not in attr_order:
                    attr_order.append(k)
                if k not in p["attrs_raw"] and not tableio.is_missing(v):
                    p["attrs_raw"][k] = v

    def lift_attrs(kind_rows):
        conflicts = Counter()
        for ti, rows, spec in kind_rows:
            for rec in rows:
                pid = pid_of(rec)
                if pid is None:
                    continue
                p = ensure(pid)
                for k, v in rec["_attrs"].items():
                    if k not in attr_order:
                        attr_order.append(k)
                    if tableio.is_missing(v):
                        continue
                    if k not in p["attrs_raw"]:
                        p["attrs_raw"][k] = v
                    elif p["attrs_raw"][k] != v:
                        conflicts[k] += 1
                for f in ("purity", "ploidy", "goodness_of_fit", "sex"):
                    if rec.get(f) is not None and p[f] is None:
                        p[f] = rec[f]
        for k, n in conflicts.items():
            warnings.append("attribute %r has different values within the same patient (%d rows); the first value "
                            "is shown" % (k, n))

    # ---------------- genomic records ----------------
    unknown_contigs = Counter()
    maxpos = defaultdict(int)

    def track_pos(c, *ps):
        for p_ in ps:
            if p_ is not None:
                maxpos[c] = max(maxpos[c], p_)

    segs_by = defaultdict(list)
    seg_missing = 0
    for ti, rows, spec in loaded["segments"]:
        for rec in rows:
            pid, c = pid_of(rec), norm_chrom(rec.get("chrom"))
            s, e = rec.get("start"), rec.get("end")
            if pid is None or c is None or s is None or e is None:
                seg_missing += 1
                continue
            s, e = int(round(s)), int(round(e))
            if e < s:
                s, e = e, s
            major, minor = rec.get("major_cn"), rec.get("minor_cn")
            if major is not None and minor is not None and minor > major:
                major, minor = minor, major
            majr, minr = rec.get("major_raw"), rec.get("minor_raw")
            if majr is not None and minr is not None and minr > majr:
                majr, minr = minr, majr
            total = rec.get("total_cn")
            if total is None and major is not None and minor is not None:
                total = major + minor
            total_raw = rec.get("total_raw")
            if total_raw is None and majr is not None and minr is not None:
                total_raw = majr + minr
            if total is None and total_raw is not None:
                total = total_raw  # e.g. Delly seg.bed: only a float CN estimate
            seg = {"chrom": c, "start": s, "end": e,
                   "major": None if major is None else int(round(major)),
                   "minor": None if minor is None else int(round(minor)),
                   "major_raw": majr, "minor_raw": minr, "total": total, "total_raw": total_raw,
                   "logr": rec.get("logr"), "baf": rec.get("baf"), "n_markers": rec.get("n_markers")}
            if seg["major"] is None and seg["total"] is None and seg["logr"] is None:
                seg_missing += 1
                continue
            segs_by[pid].append(seg)
            ensure(pid)["sources"].add("segments")
            track_pos(c, s, e)
    if seg_missing:
        warnings.append("%d segment rows skipped (missing patient, chromosome, start/end or any copy-number value)"
                        % seg_missing)
    lift_attrs(loaded["segments"])

    svs = []
    sv_missing = 0
    single_patient_tables = 0
    for ti, rows, spec in loaded["variants"]:
        for rec in rows:
            pid = pid_of(rec)
            if pid is None:
                if "patient" not in (spec.get("columns") or {}):
                    pid = "sample"
                    single_patient_tables += 1
                else:
                    sv_missing += 1
                    continue
            c = norm_chrom(rec.get("chrom"))
            pos = rec.get("pos")
            if c is None or pos is None:
                sv_missing += 1
                continue
            pos = int(round(pos))
            svtype = norm_svtype(rec.get("svtype"))
            if svtype is None and rec.get("alt"):
                alt = str(rec["alt"])
                svtype = norm_svtype(alt) if alt.startswith("<") else ("BND" if re.search(r"[\[\]]", alt) else None)
            if svtype is None and rec.get("id"):
                m = re.match(r"^(DEL|DUP|INV|BND|INS|CNV|TRA)", str(rec["id"]).upper())
                svtype = norm_svtype(m.group(1)) if m else None
            c2 = norm_chrom(rec.get("chrom2")) or c
            end = rec.get("end")
            pos2 = rec.get("pos2")
            end = None if end is None else int(round(end))
            pos2 = None if pos2 is None else int(round(pos2))
            if c2 != c:
                svtype = svtype if svtype in ("BND",) else ("BND" if svtype in (None, "INV", "DEL", "DUP") else svtype)
                if pos2 is None:
                    pos2 = end
                end = None
            elif end is None and pos2 is not None:
                end = pos2
            svlen = rec.get("svlen")
            if svlen is None and end is not None and svtype not in ("BND", "INS") and c2 == c:
                svlen = end - pos
            gt = rec.get("genotype")
            prec = rec.get("precise")
            svs.append({"patient": pid, "chrom": c, "pos": pos, "end": end, "chrom2": c2, "pos2": pos2,
                        "svtype": svtype or "OTHER", "svlen": None if svlen is None else abs(int(round(svlen))),
                        "qual": rec.get("qual"), "filter": None if rec.get("filter") is None else str(rec["filter"]),
                        "pe": rec.get("pe"), "sr": rec.get("sr"), "gt": None if gt is None else str(gt),
                        "id": None if rec.get("id") is None else str(rec["id"]),
                        "precise": None if prec is None else str(prec), "cn": rec.get("cn"),
                        "extra": rec.get("_extra")})
            ensure(pid)["sources"].add("variants")
            track_pos(c, pos, end)
            track_pos(c2, pos2)
    if sv_missing:
        warnings.append("%d variant rows skipped (missing patient, chromosome or position)" % sv_missing)
    if single_patient_tables:
        warnings.append("a variants table has no patient column; its rows were assigned to patient 'sample' "
                        "(map 'patient' with {\"value\": \"ID\"} or {\"from_filename\": \"regex\"})")
    lift_attrs(loaded["variants"])

    bins_by = defaultdict(list)
    bin_missing = 0
    for ti, rows, spec in loaded["bins"]:
        for rec in rows:
            pid = pid_of(rec) or "sample"
            c = norm_chrom(rec.get("chrom"))
            pos = rec.get("pos")
            if pos is None and rec.get("start") is not None:
                pos = (rec["start"] + (rec.get("end") or rec["start"])) / 2
            if c is None or pos is None:
                bin_missing += 1
                continue
            vals = (rec.get("logr"), rec.get("baf"), rec.get("cn"))
            if all(v is None for v in vals):
                bin_missing += 1
                continue
            bins_by[pid].append((c, int(round(pos)), vals))
            ensure(pid)["sources"].add("bins")
            track_pos(c, int(round(pos)))
    if bin_missing:
        warnings.append("%d bin rows skipped (missing chromosome/position or all values empty)" % bin_missing)
    lift_attrs(loaded["bins"])

    if not patients:
        raise MappingError("no patients found - check the 'patient' mapping")

    # ---------------- genome & chromosome indices ----------------
    genome_name, chroms = choose_genome(mapping.get("genome", "auto"), maxpos, warnings)
    cindex = {c: i for i, (c, _) in enumerate(chroms)}

    def drop_unknown(c):
        if c not in cindex:
            unknown_contigs[c] += 1
            return True
        return False

    for pid in list(segs_by):
        kept = [s for s in segs_by[pid] if not drop_unknown(s["chrom"])]
        kept.sort(key=lambda s: (cindex[s["chrom"]], s["start"]))
        segs_by[pid] = kept
    svs = [s for s in svs if not drop_unknown(s["chrom"])]
    for s in svs:
        if s["chrom2"] not in cindex:
            unknown_contigs[s["chrom2"]] += 1
            s["chrom2"] = None
    for pid in list(bins_by):
        kept = [b for b in bins_by[pid] if not drop_unknown(b[0])]
        kept.sort(key=lambda b: (cindex[b[0]], b[1]))
        bins_by[pid] = kept
    if unknown_contigs:
        top = ", ".join("%s (%d)" % kv for kv in unknown_contigs.most_common(6))
        warnings.append("rows on contigs outside the genome were left out of genome plots: %s" % top)

    # ---------------- purity / gof units ----------------
    pur = [p["purity"] for p in patients.values() if p["purity"] is not None]
    if pur and max(pur) > 1.5:
        for p in patients.values():
            if p["purity"] is not None:
                p["purity"] = p["purity"] / 100.0
        warnings.append("purity values look like percentages (max %.4g); divided by 100" % max(pur))
    gof = [p["goodness_of_fit"] for p in patients.values() if p["goodness_of_fit"] is not None]
    if gof and max(gof) <= 1.0:
        for p in patients.values():
            if p["goodness_of_fit"] is not None:
                p["goodness_of_fit"] = p["goodness_of_fit"] * 100.0
        warnings.append("goodness of fit given as a 0-1 fraction; shown as a percentage (ASCAT reports %)")
    bad = [p["id"] for p in patients.values() if p["purity"] is not None and not (0 <= p["purity"] <= 1.05)]
    if bad:
        warnings.append("purity outside 0-1 for %d patient(s): %s" % (len(bad), ", ".join(bad[:8])))

    # ---------------- attributes ----------------
    attr_cfg = mapping.get("attributes") or {}
    hidden = set(mapping.get("hide_attributes") or [])
    attributes = []
    for key in attr_order:
        if key in hidden:
            continue
        cfg = attr_cfg.get(key, {})
        vals = [p["attrs_raw"].get(key) for p in patients.values()]
        typ = cfg.get("type") or infer_attr_type(vals)
        if typ not in ("numeric", "category", "ordinal", "date", "text"):
            raise MappingError("attributes[%r].type must be numeric, category, ordinal, date or text" % key)
        order = cfg.get("order")
        if typ in ("category", "ordinal") and not order:
            present = [coerce_attr(v, typ) for v in vals if v is not None]
            present = [v for v in present if v is not None]
            if all(tableio.parse_number(v) is not None for v in present):
                order = sorted(set(present), key=lambda v: tableio.parse_number(v))
            elif typ == "ordinal":
                order = sorted(set(present), key=natural_key)
        attributes.append({"key": key, "label": cfg.get("label", key), "label_zh": cfg.get("label_zh"), "type": typ,
                           "order": order, "unit": cfg.get("unit"), "unit_zh": cfg.get("unit_zh"),
                           "description": cfg.get("description"), "description_zh": cfg.get("description_zh")})
    attr_types = {a["key"]: a["type"] for a in attributes}

    group_by = mapping.get("group_by")
    if group_by and group_by not in attr_types:
        warnings.append("group_by %r is not a patient attribute; ignored" % group_by)
        group_by = None
    if group_by is None:
        for a in attributes:
            if a["type"] in ("category", "ordinal"):
                n = len({p["attrs_raw"].get(a["key"]) for p in patients.values()} - {None})
                if 2 <= n <= 6:
                    group_by = a["key"]
                    break

    # ---------------- assemble output ----------------
    pids = list(patients)
    pindex = {pid: i for i, pid in enumerate(pids)}
    sv_by = defaultdict(list)
    for s in svs:
        sv_by[s["patient"]].append(s)

    y_limit = mapping.get("y_limit", 5)
    max_points = int(mapping.get("max_points_per_patient", 40000))
    out_patients, out_segments, out_bins = [], {}, {}
    downsampled = []
    for pid in pids:
        p = patients[pid]
        segs = segs_by.get(pid, [])
        metrics = ascat_metrics(segs) if segs else {}
        pts = sv_by.get(pid, [])
        if pts:
            metrics["n_sv"] = len(pts)
            metrics["n_sv_pass"] = sum(1 for s in pts if (s["filter"] or "").upper() == "PASS")
        ploidy, ploidy_src = p["ploidy"], "input"
        if ploidy is None and "ploidy_from_segments" in metrics:
            ploidy, ploidy_src = metrics["ploidy_from_segments"], "segments"
        elif ploidy is None:
            ploidy_src = None
        sex = p["sex"]
        sex_norm = SEX_SYNONYMS.get(str(sex).strip().upper(), SEX_SYNONYMS.get(str(sex).strip())) if sex else None
        out_patients.append({
            "id": pid,
            "attrs": {a["key"]: coerce_attr(p["attrs_raw"].get(a["key"]), a["type"]) for a in attributes},
            "purity": r4(p["purity"]), "ploidy": r4(ploidy), "ploidy_source": ploidy_src,
            "goodness_of_fit": r4(p["goodness_of_fit"]), "sex": sex_norm or (str(sex) if sex else None),
            "metrics": {k: r4(v) for k, v in metrics.items()},
            "has": {"segments": bool(segs), "variants": bool(pts), "bins": bool(bins_by.get(pid))},
        })
        if segs:
            out_segments[pid] = [[cindex[s["chrom"]], s["start"], s["end"], s["major"], s["minor"],
                                  r4(s["major_raw"]), r4(s["minor_raw"]), r4(s["total"]), r4(s["total_raw"]),
                                  r4(s["logr"]), r4(s["baf"])] for s in segs]
        b = bins_by.get(pid)
        if b:
            if len(b) > max_points:
                step = len(b) / float(max_points)
                b = [b[int(i * step)] for i in range(max_points)]
                downsampled.append(pid)
            out_bins[pid] = {"c": [cindex[x[0]] for x in b], "p": [x[1] for x in b],
                             "logr": [r4(x[2][0]) for x in b], "baf": [r4(x[2][1]) for x in b],
                             "cn": [r4(x[2][2]) for x in b]}
            for k in ("logr", "baf", "cn"):
                if all(v is None for v in out_bins[pid][k]):
                    del out_bins[pid][k]
    if downsampled:
        warnings.append("bins down-sampled to %d points for %d patient(s) to keep the HTML light"
                        % (max_points, len(downsampled)))

    sv_rows = []
    for s in svs:
        sv_rows.append([pindex[s["patient"]], cindex[s["chrom"]], s["pos"], s["end"],
                        None if s["chrom2"] is None else cindex[s["chrom2"]], s["pos2"], s["svtype"], s["svlen"],
                        r4(s["qual"]), s["filter"], r4(s["pe"]), r4(s["sr"]), s["gt"], s["id"], s["precise"],
                        r4(s["cn"])] + list(s["extra"] or []))

    missing_from_patient_table = [pid for pid in pids
                                  if loaded["patients"] and "patients" not in patients[pid]["sources"]]
    if missing_from_patient_table:
        warnings.append("%d patient(s) have genomic data but no row in the patients table: %s"
                        % (len(missing_from_patient_table), ", ".join(missing_from_patient_table[:10])))
    no_genomic = [pid for pid in pids if patients[pid]["sources"] == {"patients"}]
    if no_genomic and (segs_by or svs or bins_by):
        warnings.append("%d patient(s) in the patients table have no genomic rows: %s"
                        % (len(no_genomic), ", ".join(no_genomic[:10])))

    data = {
        "version": VERSION,
        "meta": {
            "title": mapping.get("title") or "Genomic report",
            "subtitle": mapping.get("subtitle") or "",
            "title_zh": mapping.get("title_zh"),
            "subtitle_zh": mapping.get("subtitle_zh"),
            "language": mapping.get("language", "en"),
            "generated": dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "genome": genome_name, "y_limit": y_limit, "group_by": group_by,
            "insights": mapping.get("insights") or [],
            "insights_zh": mapping.get("insights_zh") or [],
            "insights_label": mapping.get("insights_label"),
            "insights_label_zh": mapping.get("insights_label_zh"),
            "sources": sources, "warnings": warnings,
            "mapping": mapping,
        },
        "chroms": chroms,
        "attributes": attributes,
        "patients": out_patients,
        "segments": out_segments,
        "bins": out_bins,
        "sv_columns": ["patient", "chrom", "pos", "end", "chrom2", "pos2", "svtype", "svlen", "qual", "filter",
                       "pe", "sr", "gt", "id", "precise", "cn"],
        "sv_extra": sv_extra_names,
        "svs": sv_rows,
    }
    return data


# --------------------------------------------------------------------------- #
# cohort summary (printed for the agent; also --stats-json)
# --------------------------------------------------------------------------- #
def baseline_cn(p):
    """Total-CN baseline used for gain/loss calls: round(ploidy), 2 when unknown.

    Mirrors the report's JavaScript (cnState) so the printed summary and the
    cohort heatmap agree.
    """
    return max(1, int(math.floor(p["ploidy"] + 0.5))) if p.get("ploidy") else 2


def base_for_chrom(base, chrom, sex):
    # one X (and one Y) is the normal state in males
    if sex == "XY" and chrom in SEX_CHROMS:
        return max(1, int(math.floor(base / 2.0 + 0.5)))
    return base


def summarize(data):
    pts = data["patients"]
    chroms = data["chroms"]

    def med(xs):
        xs = sorted(x for x in xs if x is not None)
        if not xs:
            return None
        n = len(xs)
        return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2

    s = OrderedDict()
    s["patients"] = len(pts)
    s["with_segments"] = sum(p["has"]["segments"] for p in pts)
    s["with_variants"] = sum(p["has"]["variants"] for p in pts)
    s["with_bins"] = sum(p["has"]["bins"] for p in pts)
    s["genome"] = data["meta"]["genome"]
    s["purity_median"] = r4(med(p["purity"] for p in pts))
    s["ploidy_median"] = r4(med(p["ploidy"] for p in pts))
    wgd = [p["metrics"].get("WGD") for p in pts if "WGD" in p["metrics"]]
    if wgd:
        s["WGD_fraction"] = r4(sum(1 for w in wgd if w not in (0, None)) / len(wgd))
    for k in ("GI", "LOH", "FGA"):
        vals = [p["metrics"].get(k) for p in pts if p["metrics"].get(k) is not None]
        if vals:
            s[k + "_median"] = r4(med(vals))
    # chromosome-level recurrent events: >50% of a chromosome gained / lost vs round(ploidy)
    gain, loss, nseg = Counter(), Counter(), 0
    for p in pts:
        segs = data["segments"].get(p["id"])
        if not segs:
            continue
        nseg += 1
        base = baseline_cn(p)
        cov = defaultdict(lambda: [0, 0, 0])  # gain bp, loss bp, total bp
        for ci, st, en, mj, mn, mjr, mnr, tot, totr, lr, bf in segs:
            if tot is None:
                continue
            b = base_for_chrom(base, chroms[ci][0], p["sex"])
            L = en - st + 1
            cov[ci][2] += L
            if tot > b:
                cov[ci][0] += L
            elif tot < b:
                cov[ci][1] += L
        for ci, (g, l, t) in cov.items():
            if t and g / t > 0.5:
                gain[ci] += 1
            if t and l / t > 0.5:
                loss[ci] += 1
    if nseg:
        s["chromosome_gains"] = [["chr" + chroms[ci][0], r4(n / nseg)] for ci, n in gain.most_common(6)]
        s["chromosome_losses"] = [["chr" + chroms[ci][0], r4(n / nseg)] for ci, n in loss.most_common(6)]
    if data["svs"]:
        types = Counter(r[6] for r in data["svs"])
        s["sv_total"] = len(data["svs"])
        s["sv_pass"] = sum(1 for r in data["svs"] if (r[9] or "").upper() == "PASS")
        s["sv_types"] = dict(types.most_common())
        per = Counter(r[0] for r in data["svs"])
        s["sv_most"] = [[pts[i]["id"], n] for i, n in per.most_common(5)]
    s["attributes"] = [a["label"] + " (" + a["type"] + ")" for a in data["attributes"]]
    s["group_by"] = data["meta"]["group_by"]
    return s


def write_html(data, out_path):
    with open(TEMPLATE, encoding="utf-8") as fh:
        tpl = fh.read()
    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False, default=str)
    payload = payload.replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    title = data["meta"]["title"]
    title = title.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    html = tpl.replace("__REPORT_TITLE__", title)
    html = html.replace("__REPORT_DATA__", payload)
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    return len(html.encode("utf-8"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mapping")
    ap.add_argument("-o", "--out", help="output HTML (default: <mapping name>.html next to the mapping)")
    ap.add_argument("--check", action="store_true", help="validate and summarize only")
    ap.add_argument("--stats-json", help="write the cohort summary as JSON to this path")
    args = ap.parse_args(argv)

    try:
        with open(args.mapping, encoding="utf-8") as fh:
            mapping = json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        print("error: cannot read mapping %s: %s" % (args.mapping, exc), file=sys.stderr)
        return 2
    base = os.path.dirname(os.path.abspath(args.mapping))
    try:
        data = build(mapping, base)
    except MappingError as exc:
        print("mapping error: %s" % exc, file=sys.stderr)
        return 1
    except (OSError, ValueError) as exc:
        print("data error: %s" % exc, file=sys.stderr)
        return 1

    summary = summarize(data)
    print("summary:")
    for k, v in summary.items():
        print("  %s: %s" % (k, json.dumps(v, ensure_ascii=False) if not isinstance(v, str) else v))
    if data["meta"]["warnings"]:
        print("warnings:")
        for w in data["meta"]["warnings"]:
            print("  - " + w)
    if args.stats_json:
        with open(args.stats_json, "w", encoding="utf-8") as fh:
            json.dump(summary, fh, ensure_ascii=False, indent=1)
    if args.check:
        print("check ok (no HTML written)")
        return 0
    out = args.out or os.path.splitext(os.path.abspath(args.mapping))[0] + ".html"
    size = write_html(data, out)
    print("wrote %s (%.1f KB)" % (out, size / 1024.0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
