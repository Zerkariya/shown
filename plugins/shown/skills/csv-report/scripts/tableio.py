"""Read tabular files (CSV / TSV / TXT / XLSX / VCF, optionally gzipped) with the
Python standard library only, so the skill runs on any machine with python3.

Public API:
    read_table(path, sheet=None, header_row=None) -> Table
    list_sheets(path) -> [sheet names]            (xlsx only)
    parse_number(value) -> float | None
    MISSING                                         (set of NA-like tokens)
"""
from __future__ import annotations

import csv
import datetime as _dt
import gzip
import io
import os
import re
import zipfile
import xml.etree.ElementTree as ET

MISSING = {"", "na", "n/a", "nan", "null", "none", ".", "-", "--", "#n/a", "nd"}

_NUM_RE = re.compile(r"^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$")


class Table:
    """A parsed table: header names plus rows of strings (same length as header)."""

    def __init__(self, path, columns, rows, sheet=None, header_row=1, skipped_rows=0, fmt="csv", meta=None):
        self.path = path
        self.columns = columns
        self.rows = rows
        self.sheet = sheet
        self.header_row = header_row
        self.skipped_rows = skipped_rows
        self.format = fmt
        self.meta = meta or {}

    def column(self, name):
        i = self.columns.index(name)
        return [r[i] for r in self.rows]


def is_missing(value) -> bool:
    return value is None or str(value).strip().lower() in MISSING


def parse_number(value):
    """'1,234' -> 1234.0, '56%' -> 56.0, '3.2e5' -> 320000.0, NA-like -> None."""
    if value is None:
        return None
    s = str(value).strip()
    if s.lower() in MISSING:
        return None
    s = s.replace(",", "").replace(" ", "").replace(" ", "")
    if s.endswith("%"):
        s = s[:-1]
    if _NUM_RE.match(s):
        try:
            return float(s)
        except ValueError:
            return None
    return None


# --------------------------------------------------------------------------- #
# raw bytes / text helpers
# --------------------------------------------------------------------------- #
def _read_bytes(path):
    with open(path, "rb") as fh:
        data = fh.read()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data


def _decode(data: bytes) -> str:
    # utf-8 first (with BOM), then GB18030 for Chinese Excel exports, then latin-1.
    for enc in ("utf-8-sig", "gb18030"):
        try:
            return data.decode(enc)
        except UnicodeDecodeError:
            continue
    return data.decode("latin-1")


def _strip_gz(name: str) -> str:
    return name[:-3] if name.lower().endswith(".gz") else name


def detect_format(path: str) -> str:
    name = _strip_gz(os.path.basename(path)).lower()
    with open(path, "rb") as fh:
        head = fh.read(4)
    if head[:2] == b"PK" or name.endswith((".xlsx", ".xlsm")):
        return "xlsx"
    if name.endswith(".xls"):
        return "xls"
    if name.endswith(".bcf"):
        return "bcf"
    if name.endswith(".vcf"):
        return "vcf"
    return "delimited"


# --------------------------------------------------------------------------- #
# header detection
# --------------------------------------------------------------------------- #
def _nonempty(row):
    return [c for c in row if not is_missing(c)]


def guess_header_row(rows, scan=30) -> int:
    """Return the 0-based index of the most likely header row.

    Spreadsheets often carry title / note lines above the real header; the header
    is the first row that is about as wide as the widest row and mostly text.
    """
    sample = rows[:scan]
    if not sample:
        return 0
    widths = [len(_nonempty(r)) for r in sample]
    max_w = max(widths) or 1
    for i, row in enumerate(sample):
        cells = _nonempty(row)
        if len(cells) < max(1, 0.6 * max_w):
            continue
        texty = sum(1 for c in cells if parse_number(c) is None)
        if texty >= 0.6 * len(cells):
            return i
        # wide row of numbers: the file has no header at all
        return -1 if i == 0 else i
    return 0


def _finish(path, raw_rows, header_row, sheet, fmt, meta=None):
    raw_rows = [list(r) for r in raw_rows]
    # trim fully empty trailing rows
    while raw_rows and not _nonempty(raw_rows[-1]):
        raw_rows.pop()
    if header_row is None:
        h = guess_header_row(raw_rows)
    else:
        h = header_row - 1  # user-facing 1-based
    if h == -1:
        width = max((len(r) for r in raw_rows), default=0)
        header = ["column_%d" % (i + 1) for i in range(width)]
        body = raw_rows
        h_display, skipped = 0, 0
    else:
        header = raw_rows[h] if raw_rows else []
        body = raw_rows[h + 1:]
        h_display, skipped = h + 1, h
    width = max([len(header)] + [len(r) for r in body]) if (header or body) else 0
    names, seen = [], {}
    for i in range(width):
        name = str(header[i]).strip() if i < len(header) else ""
        name = re.sub(r"\s+", " ", name) or "column_%d" % (i + 1)
        if name in seen:
            seen[name] += 1
            name = "%s_%d" % (name, seen[name])
        else:
            seen[name] = 1
        names.append(name)
    rows = []
    for r in body:
        if not _nonempty(r):
            continue
        r = [("" if c is None else str(c).strip()) for c in r]
        r += [""] * (width - len(r))
        rows.append(r[:width])
    # drop columns that are entirely empty and unnamed
    keep = [i for i, n in enumerate(names)
            if not (n.startswith("column_") and all(is_missing(r[i]) for r in rows))]
    if len(keep) != len(names):
        names = [names[i] for i in keep]
        rows = [[r[i] for i in keep] for r in rows]
    # spreadsheet exports often end with a note line ("End of sheet", "注：...")
    meta = dict(meta or {})
    footer = []
    while rows and len(footer) < 5 and len(names) >= 3:
        cells = _nonempty(rows[-1])
        if len(cells) == 1 and _looks_like_note(cells[0]):
            footer.insert(0, cells[0])
            rows.pop()
        else:
            break
    if footer:
        meta["footer_rows_dropped"] = footer
    return Table(path, names, rows, sheet=sheet, header_row=h_display, skipped_rows=skipped, fmt=fmt, meta=meta)


_NOTE_START = re.compile(r"^(注|备注|说明|来源|数据来源|note|notes|source|end\b|total\b|\*|#)", re.I)


def _looks_like_note(cell: str) -> bool:
    cell = str(cell).strip()
    return bool(_NOTE_START.match(cell) or (re.search(r"\s", cell) and len(cell) >= 12) or len(cell) >= 40)


_DEC_COMMA = re.compile(r"^[+-]?\d+,\d+$")
_THOUSANDS = re.compile(r"^[+-]?\d{1,3}(,\d{3})+$")


def _fix_decimal_commas(t: "Table") -> None:
    """European exports write 0,55 for 0.55. Convert such columns in place.

    A column is converted when every value is a number, none uses '.', and at
    least one value cannot be a thousands separator (e.g. 0,0093 or 12,5).
    Values like 1,234 alone stay thousands-separated unless the same file
    already proved to use decimal commas and ',' is not the field delimiter.
    """
    cand, proven = [], False
    for j in range(len(t.columns)):
        vals = [r[j] for r in t.rows if not is_missing(r[j])]
        if not vals or any("." in v for v in vals):
            continue
        dec = [v for v in vals if _DEC_COMMA.match(v)]
        if not dec or not all(_DEC_COMMA.match(v) or _NUM_RE.match(v) for v in vals):
            continue
        unambiguous = any(not _THOUSANDS.match(v) for v in dec)
        proven = proven or unambiguous
        cand.append((j, unambiguous))
    converted = []
    for j, unambiguous in cand:
        if unambiguous or (proven and t.meta.get("delimiter") != ","):
            for r in t.rows:
                if _DEC_COMMA.match(r[j]):
                    r[j] = r[j].replace(",", ".")
            converted.append(t.columns[j])
    if converted:
        t.meta["decimal_comma_columns"] = converted


# --------------------------------------------------------------------------- #
# delimited text
# --------------------------------------------------------------------------- #
def _sniff_delimiter(text: str, path: str) -> str:
    lines = [l for l in text.splitlines()[:50] if l.strip() and not l.startswith("##")]
    sample = "\n".join(lines)
    counts = {d: sum(l.count(d) for l in lines) for d in ("\t", ",", ";", "|")}
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;|").delimiter
    except csv.Error:
        pass
    best = max(counts, key=counts.get)
    if counts[best] > 0:
        return best
    name = _strip_gz(path).lower()
    return "\t" if name.endswith((".tsv", ".txt", ".bed", ".seg")) else ","


def _read_delimited(path, header_row=None, delimiter=None):
    text = _decode(_read_bytes(path))
    delim = delimiter or _sniff_delimiter(text, path)
    # Skip '##' comment lines (VCF-like headers in tabular exports); keep '#CHROM'.
    lines = [l for l in text.splitlines() if not l.startswith("##")]
    if lines and lines[0].startswith("#") and delim in lines[0]:
        lines[0] = lines[0][1:]
    if delim == " ":
        rows = [re.split(r"\s+", l.strip()) for l in lines]
    else:
        rows = list(csv.reader(lines, delimiter=delim))
    t = _finish(path, rows, header_row, None, "delimited", {"delimiter": delim})
    _fix_decimal_commas(t)
    # ASCAT LogR/BAF files: first header cell is empty (row names column).
    if t.header_row and t.columns and t.columns[0] == "column_1" and t.rows and len(t.columns) > 2:
        t.columns[0] = "probe_id"
    return t


# --------------------------------------------------------------------------- #
# xlsx (Office Open XML) without openpyxl
# --------------------------------------------------------------------------- #
_NS = {
    "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "rel": "http://schemas.openxmlformats.org/package/2006/relationships",
}
_DATE_FMT_IDS = set(range(14, 23)) | {45, 46, 47}


def _col_index(ref: str) -> int:
    letters = re.match(r"[A-Z]+", ref.upper())
    n = 0
    for ch in letters.group(0) if letters else "A":
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def _xlsx_sheets(z):
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    targets = {r.get("Id"): r.get("Target") for r in rels.findall("rel:Relationship", _NS)}
    out = []
    for s in wb.findall("m:sheets/m:sheet", _NS):
        rid = s.get("{%s}id" % _NS["r"])
        target = targets.get(rid, "")
        target = target.lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        out.append((s.get("name"), target))
    return out


def list_sheets(path):
    if detect_format(path) != "xlsx":
        return []
    with zipfile.ZipFile(path) as z:
        return [name for name, _ in _xlsx_sheets(z)]


def _xlsx_date_styles(z):
    try:
        st = ET.fromstring(z.read("xl/styles.xml"))
    except KeyError:
        return set()
    custom = {}
    for nf in st.findall("m:numFmts/m:numFmt", _NS):
        code = re.sub(r'"[^"]*"|\[[^\]]*\]', "", nf.get("formatCode", "").lower())
        custom[int(nf.get("numFmtId"))] = bool(re.search(r"[dmy]", code)) and "general" not in code
    date_styles = set()
    xfs = st.find("m:cellXfs", _NS)
    if xfs is None:
        return date_styles
    for i, xf in enumerate(xfs.findall("m:xf", _NS)):
        fid = int(xf.get("numFmtId", "0"))
        if fid in _DATE_FMT_IDS or custom.get(fid):
            date_styles.add(i)
    return date_styles


def _excel_serial_to_date(v: float) -> str:
    base = _dt.datetime(1899, 12, 30)
    d = base + _dt.timedelta(days=v)
    if abs(v - round(v)) < 1e-9:
        return d.date().isoformat()
    return d.isoformat(sep=" ", timespec="minutes")


def _read_xlsx(path, sheet=None, header_row=None):
    with zipfile.ZipFile(path) as z:
        sheets = _xlsx_sheets(z)
        if not sheets:
            raise ValueError("%s: workbook has no sheets" % path)
        if sheet is None:
            name, target = sheets[0]
        else:
            match = [s for s in sheets if s[0] == str(sheet)]
            if not match and str(sheet).isdigit() and int(sheet) <= len(sheets):
                match = [sheets[int(sheet) - 1]]
            if not match:
                raise ValueError("%s: sheet %r not found; sheets are %s" % (path, sheet, [s[0] for s in sheets]))
            name, target = match[0]
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            sst = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in sst.findall("m:si", _NS):
                shared.append("".join(t.text or "" for t in si.iter("{%s}t" % _NS["m"])))
        date_styles = _xlsx_date_styles(z)
        root = ET.fromstring(z.read(target))
        rows = []
        for row in root.iter("{%s}row" % _NS["m"]):
            cells = {}
            for c in row.findall("m:c", _NS):
                ref = c.get("r")
                idx = _col_index(ref) if ref else len(cells)
                typ = c.get("t")
                v = c.find("m:v", _NS)
                text = v.text if v is not None else None
                if typ == "s" and text is not None:
                    val = shared[int(text)]
                elif typ == "inlineStr":
                    val = "".join(t.text or "" for t in c.iter("{%s}t" % _NS["m"]))
                elif typ == "b":
                    val = "TRUE" if text == "1" else "FALSE"
                else:
                    val = text or ""
                    if typ not in ("str", "e") and val and int(c.get("s", "0")) in date_styles:
                        try:
                            val = _excel_serial_to_date(float(val))
                        except ValueError:
                            pass
                    elif typ is None and val and re.match(r"^-?\d+\.0+$", val):
                        val = val.split(".")[0]
                cells[idx] = val
            rnum = int(row.get("r", len(rows) + 1))
            while len(rows) < rnum - 1:
                rows.append([])
            width = max(cells) + 1 if cells else 0
            rows.append([cells.get(i, "") for i in range(width)])
    return _finish(path, rows, header_row, name, "xlsx")


# --------------------------------------------------------------------------- #
# VCF (text) -> one row per record, INFO and FORMAT fields expanded
# --------------------------------------------------------------------------- #
def _read_vcf(path):
    text = _decode(_read_bytes(path))
    info_ids, fmt_ids, samples, header = [], [], [], None
    records = []
    for line in text.splitlines():
        if line.startswith("##INFO=<ID="):
            info_ids.append(line[11:].split(",", 1)[0])
        elif line.startswith("##FORMAT=<ID="):
            fmt_ids.append(line[13:].split(",", 1)[0])
        elif line.startswith("#CHROM"):
            header = line[1:].split("\t")
            samples = header[9:]
        elif line and not line.startswith("#"):
            records.append(line.split("\t"))
    if header is None:
        raise ValueError("%s: no #CHROM header line; is this a VCF?" % path)
    used_fmt, infos = [], []
    for rec in records:
        rec += [""] * (8 - len(rec))
        info = {}
        for item in rec[7].split(";"):
            if not item or item == ".":
                continue
            k, eq, v = item.partition("=")
            info[k] = v if eq else "TRUE"
            if k not in info_ids:
                info_ids.append(k)
        infos.append(info)
        if len(rec) > 8:
            for k in rec[8].split(":"):
                if k not in used_fmt:
                    used_fmt.append(k)
    used_info = [k for k in info_ids if any(k in d for d in infos)]
    fmt_cols = [k for k in fmt_ids if k in used_fmt] + [k for k in used_fmt if k not in fmt_ids]
    cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"] + ["INFO_" + k for k in used_info]
    for s in samples:
        cols += ["%s.%s" % (s, k) for k in fmt_cols]
    rows = []
    for rec, info in zip(records, infos):
        row = rec[:7] + [info.get(k, "") for k in used_info]
        keys = rec[8].split(":") if len(rec) > 8 else []
        for si in range(len(samples)):
            vals = rec[9 + si].split(":") if len(rec) > 9 + si else []
            d = dict(zip(keys, vals))
            row += [d.get(k, "") for k in fmt_cols]
        rows.append(row)
    width = len(cols)
    rows = [r + [""] * (width - len(r)) for r in rows]
    return Table(path, cols, rows, header_row=1, fmt="vcf", meta={"samples": samples})


# --------------------------------------------------------------------------- #
def read_table(path, sheet=None, header_row=None, delimiter=None) -> Table:
    fmt = detect_format(path)
    if fmt == "xlsx":
        return _read_xlsx(path, sheet=sheet, header_row=header_row)
    if fmt == "xls":
        raise ValueError("%s: legacy .xls is not supported; save it as .xlsx or .csv first" % path)
    if fmt == "bcf":
        raise ValueError("%s: BCF is binary; convert first: bcftools view %s > out.vcf" % (path, path))
    if fmt == "vcf":
        return _read_vcf(path)
    return _read_delimited(path, header_row=header_row, delimiter=delimiter)
