#!/usr/bin/env python3
"""Generate a simulated (fake) cohort in the shapes that Delly and ASCAT produce.

Outputs (in examples/fake/ by default):
  patients.xlsx    patient info sheet: a title row above the header, free-form column
                   names, dates, purity written as "58%" - deliberately messy
  segments.tsv     ASCAT .segments_raw.txt layout (sample chr startpos endpos nMajor nMinor nAraw nBraw)
  sv_calls.csv     Delly calls flattened with bcftools query (CHROM POS ID SVTYPE END CHR2 POS2 ...)
  snp_data.tsv.gz  per-SNP logR / BAF plus a Delly-style read-depth CN for 3 patients

All values are random; nothing here is real patient data.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import math
import os
import random
import zipfile
from xml.sax.saxutils import escape

HG38 = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
        133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
        58617616, 64444167, 46709983, 50818468, 156040895]
CHROMS = [str(i) for i in range(1, 23)] + ["X"]
# approximate hg38 centromere midpoints (Mb) - only used to place arm-level events
CENTRO = {"1": 123.4, "3": 90.9, "5": 48.8, "7": 60.1, "8": 45.2, "9": 43.0, "10": 39.8, "13": 17.7,
          "17": 25.1, "18": 18.5, "20": 28.1}
# recurrent arm events: (chrom, arm, kind, probability)
RECURRENT = [("8", "q", "gain", 0.55), ("1", "q", "gain", 0.45), ("17", "p", "loss", 0.5), ("3", "p", "loss", 0.35),
             ("13", "q", "loss", 0.3), ("20", "q", "gain", 0.35), ("9", "p", "loss", 0.3), ("5", "p", "gain", 0.25),
             ("18", "q", "loss", 0.3), ("7", "p", "gain", 0.3)]
CANCERS = ["Lung adenocarcinoma", "Colorectal", "Breast", "Gastric"]
REGIMENS = ["Chemotherapy", "Targeted therapy", "Immunotherapy", "Surgery + adjuvant chemo", "Radiotherapy", "Observation"]


def simulate_patient(rng, pid, sex, wgd):
    base = 2 if wgd else 1
    segs = []
    for ci, chrom in enumerate(CHROMS):
        length = HG38[ci]
        if chrom == "X" and sex == "XY":
            segs.append([chrom, 1, length, base, 0])
            continue
        # arm-level recurrent events
        state = [[1, length, base, base]]
        for c, arm, kind, p in RECURRENT:
            if c != chrom or rng.random() > p:
                continue
            cen = int(CENTRO.get(chrom, length / 2e6) * 1e6)
            a, b = (1, cen) if arm == "p" else (cen + 1, length)
            mj, mn = (base + rng.choice([1, 1, 2]), base) if kind == "gain" else rng.choice([(base, base - 1), (base, 0)])
            state = _paint(state, a, b, mj, mn)
        # focal noise: more for unstable genomes
        for _ in range(rng.randint(0, 3 + (3 if wgd else 0))):
            a = rng.randint(1, length - 2_000_000)
            size = int(10 ** rng.uniform(5.5, 7.6))
            b = min(length, a + size)
            r = rng.random()
            if r < 0.06:
                mj, mn = 0, 0                                  # homozygous deletion
                b = min(length, a + int(10 ** rng.uniform(5, 6.3)))
            elif r < 0.12:
                mj, mn = base + rng.randint(4, 8), base        # amplification
            elif r < 0.55:
                mj, mn = base + 1, base                        # gain
            elif r < 0.8:
                mj, mn = base, 0                               # LOH
            else:
                mj, mn = base, max(0, base - 1)                # loss
            state = _paint(state, a, b, mj, mn)
        for a, b, mj, mn in state:
            segs.append([chrom, a, b, mj, mn])
    return segs


def _paint(state, a, b, mj, mn):
    out = []
    for s, e, x, y in state:
        if e < a or s > b:
            out.append([s, e, x, y])
            continue
        if s < a:
            out.append([s, a - 1, x, y])
        out.append([max(s, a), min(e, b), mj, mn])
        if e > b:
            out.append([b + 1, e, x, y])
    merged = []
    for seg in out:
        if merged and merged[-1][2:] == seg[2:] and merged[-1][1] + 1 == seg[0]:
            merged[-1][1] = seg[1]
        else:
            merged.append(seg)
    return merged


def ploidy_of(segs):
    tot = sum((e - s + 1) * (mj + mn) for _, s, e, mj, mn in segs)
    size = sum(e - s + 1 for _, s, e, _, _ in segs)
    return tot / size


def simulate_svs(rng, pid, segs, n_random):
    rows = []
    counter = {"DEL": 0, "DUP": 0, "INV": 0, "BND": 0, "INS": 0}

    def add(svtype, chrom, pos, end, chrom2, pos2):
        counter[svtype] += 1
        pe = rng.randint(0, 40)
        sr = rng.randint(0, 25)
        precise = sr > 3
        qual = int(min(10000, 20 * (pe + sr) + rng.randint(0, 200)))
        filt = "PASS" if (pe + sr >= 6 and qual >= 100) else "LowQual"
        svlen = "" if svtype == "BND" else (rng.randint(50, 6000) if svtype == "INS" else end - pos)
        if svtype == "DEL":
            svlen = -(end - pos)
        gt = rng.choice(["0/1", "0/1", "0/1", "1/1"])
        rows.append([chrom, pos, "%s%08d" % (svtype, counter[svtype]), svtype, end, chrom2, pos2, svlen,
                     pe, sr, "PRECISE" if precise else "IMPRECISE", qual, filt, gt, pid])

    # breakpoints at CN changes
    for i in range(1, len(segs)):
        prev, cur = segs[i - 1], segs[i]
        if prev[0] != cur[0]:
            continue
        if cur[3] + cur[4] < prev[3] + prev[4]:
            end = cur[2]
            if end - cur[1] < 5e7:
                add("DEL", cur[0], cur[1], end, cur[0], end)
        elif cur[3] + cur[4] > prev[3] + prev[4]:
            end = cur[2]
            if end - cur[1] < 5e7:
                add("DUP", cur[0], cur[1], end, cur[0], end)
    for _ in range(n_random):
        ci = rng.randrange(len(CHROMS))
        chrom = CHROMS[ci]
        pos = rng.randint(1_000_000, HG38[ci] - 1_000_000)
        r = rng.random()
        if r < 0.35:
            end = min(HG38[ci] - 1, pos + int(10 ** rng.uniform(2.5, 6.5)))
            add("DEL", chrom, pos, end, chrom, end)
        elif r < 0.55:
            end = min(HG38[ci] - 1, pos + int(10 ** rng.uniform(3, 6.8)))
            add("INV", chrom, pos, end, chrom, end)
        elif r < 0.7:
            end = min(HG38[ci] - 1, pos + int(10 ** rng.uniform(3, 6.5)))
            add("DUP", chrom, pos, end, chrom, end)
        elif r < 0.85:
            cj = rng.randrange(len(CHROMS))
            pos2 = rng.randint(1_000_000, HG38[cj] - 1_000_000)
            add("BND", chrom, pos, pos2, CHROMS[cj], pos2)
        else:
            add("INS", chrom, pos, pos + 1, chrom, pos + 1)
    return rows


def simulate_snps(rng, pid, segs, purity, psi, step=150_000):
    rows = []
    for chrom, s, e, mj, mn in segs:
        pos = s + rng.randint(0, step)
        n = mj + mn
        while pos < e:
            denom = purity * psi + 2 * (1 - purity)
            logr = math.log2(max(1e-3, (purity * n + 2 * (1 - purity)) / denom)) + rng.gauss(0, 0.16)
            if rng.random() < 0.35:   # heterozygous SNP
                b = (1 - purity + purity * mn) / (2 - 2 * purity + purity * n) if n else 0.5
                baf = b if rng.random() < 0.5 else 1 - b
                baf = min(1, max(0, baf + rng.gauss(0, 0.035)))
            else:                      # homozygous
                baf = min(1, max(0, rng.choice([0, 1]) + rng.gauss(0, 0.012)))
                baf = round(baf, 4)
            cn = max(0, purity * n + 2 * (1 - purity) + rng.gauss(0, 0.22))
            rows.append([pid, "chr" + chrom, pos, round(logr, 4), round(baf, 4), round(cn, 3)])
            pos += step + rng.randint(-step // 3, step // 3)
    return rows


# --------------------------------------------------------------------------- #
# minimal xlsx writer (shared strings + a date style) to exercise the reader
# --------------------------------------------------------------------------- #
def _col(n):
    s = ""
    n += 1
    while n:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def write_xlsx(path, sheet_name, rows):
    strings, index = [], {}

    def sid(v):
        if v not in index:
            index[v] = len(strings)
            strings.append(v)
        return index[v]

    xml_rows = []
    for ri, row in enumerate(rows, start=1):
        cells = []
        for ci, v in enumerate(row):
            ref = "%s%d" % (_col(ci), ri)
            if v is None or v == "":
                continue
            if isinstance(v, dt.date):
                serial = (v - dt.date(1899, 12, 30)).days
                cells.append('<c r="%s" s="1"><v>%d</v></c>' % (ref, serial))
            elif isinstance(v, (int, float)) and not isinstance(v, bool):
                cells.append('<c r="%s"><v>%s</v></c>' % (ref, repr(v)))
            else:
                cells.append('<c r="%s" t="s"><v>%d</v></c>' % (ref, sid(str(v))))
        xml_rows.append('<row r="%d">%s</row>' % (ri, "".join(cells)))
    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    rel_ns = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
    files = {
        "[Content_Types].xml": (
            '<?xml version="1.0" encoding="UTF-8"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            '<Override PartName="/xl/sharedStrings.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml"/>'
            '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
            '</Types>'),
        "_rels/.rels": (
            '<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
            '</Relationships>'),
        "xl/workbook.xml": (
            '<?xml version="1.0" encoding="UTF-8"?><workbook %s %s><sheets><sheet name="%s" sheetId="1" r:id="rId1"/></sheets></workbook>'
            % (ns, rel_ns, escape(sheet_name))),
        "xl/_rels/workbook.xml.rels": (
            '<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>'
            '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/sharedStrings" Target="sharedStrings.xml"/>'
            '<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
            '</Relationships>'),
        "xl/styles.xml": (
            '<?xml version="1.0" encoding="UTF-8"?><styleSheet %s><fonts count="1"><font/></fonts><fills count="1"><fill/></fills>'
            '<borders count="1"><border/></borders><cellStyleXfs count="1"><xf/></cellStyleXfs>'
            '<cellXfs count="2"><xf numFmtId="0"/><xf numFmtId="14" applyNumberFormat="1"/></cellXfs></styleSheet>' % ns),
        "xl/sharedStrings.xml": (
            '<?xml version="1.0" encoding="UTF-8"?><sst %s count="%d" uniqueCount="%d">%s</sst>'
            % (ns, len(strings), len(strings), "".join("<si><t>%s</t></si>" % escape(s) for s in strings))),
        "xl/worksheets/sheet1.xml": (
            '<?xml version="1.0" encoding="UTF-8"?><worksheet %s><sheetData>%s</sheetData></worksheet>' % (ns, "".join(xml_rows))),
    }
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in files.items():
            z.writestr(name, data)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "fake"))
    ap.add_argument("--patients", type=int, default=40)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    os.makedirs(args.out, exist_ok=True)

    patient_rows = [["Cohort A - patient information (SIMULATED, demo only)"], [],
                    ["Patient No.", "Age", "Gender", "Cancer type", "Clinical stage", "Smoking history",
                     "Diagnosis date", "Tumor purity", "Ploidy", "GoF (%)", "Treatment", "Follow-up (months)", "Notes"]]
    seg_rows = [["sample", "chr", "startpos", "endpos", "nMajor", "nMinor", "nAraw", "nBraw"]]
    sv_rows = [["CHROM", "POS", "ID", "SVTYPE", "END", "CHR2", "POS2", "SVLEN", "PE", "SR", "PRECISION",
                "QUAL", "FILTER", "GT", "SAMPLE"]]
    snp_rows = [["Sample", "Chromosome", "Position", "LogR", "BAF", "CN"]]
    for i in range(1, args.patients + 1):
        pid = "P%03d" % i
        sex = rng.choice(["XX", "XY"])
        wgd = rng.random() < 0.35
        purity = round(rng.uniform(0.22, 0.95), 2)
        segs = simulate_patient(rng, pid, sex, wgd)
        psi = ploidy_of(segs)
        for chrom, s, e, mj, mn in segs:
            nar = max(0, mj + rng.gauss(0, 0.07))
            nbr = max(0, mn + rng.gauss(0, 0.07))
            seg_rows.append([pid, chrom, s, e, mj, mn, round(nar, 4), round(nbr, 4)])
        sv_rows += simulate_svs(rng, pid, segs, rng.randint(10, 90))
        if i <= 3:
            snp_rows += simulate_snps(rng, pid, segs, purity, psi)
        age = rng.randint(31, 84)
        stage = rng.choice(["I", "II", "II", "III", "III", "IV"])
        smoke = rng.choice(["Yes", "No", "No", "Unknown"])
        diag = dt.date(2019, 1, 1) + dt.timedelta(days=rng.randint(0, 5 * 365))
        note = rng.choice(["", "", "", "Relapse", "Multiple metastases", "Family history"])
        patient_rows.append([pid, age if rng.random() > 0.05 else "", "Male" if sex == "XY" else "Female",
                             rng.choice(CANCERS), stage, smoke, diag, "%d%%" % round(purity * 100),
                             round(psi, 2), round(rng.uniform(82, 99.5), 1), rng.choice(REGIMENS),
                             rng.randint(2, 60), note])

    write_xlsx(os.path.join(args.out, "patients.xlsx"), "Patients", patient_rows)
    with open(os.path.join(args.out, "segments.tsv"), "w") as fh:
        fh.write("\n".join("\t".join(map(str, r)) for r in seg_rows) + "\n")
    with open(os.path.join(args.out, "sv_calls.csv"), "w") as fh:
        fh.write("\n".join(",".join(map(str, r)) for r in sv_rows) + "\n")
    with gzip.open(os.path.join(args.out, "snp_data.tsv.gz"), "wt") as fh:
        fh.write("\n".join("\t".join(map(str, r)) for r in snp_rows) + "\n")
    print("wrote %d patients, %d segments, %d SVs, %d SNP rows to %s" % (
        args.patients, len(seg_rows) - 1, len(sv_rows) - 1, len(snp_rows) - 1, args.out))


if __name__ == "__main__":
    main()
