"""Regression tests for the csv-report skill (stdlib unittest, no dependencies).

    python3 -m unittest discover -s tests -v

Each test writes small synthetic inputs in different shapes (GBK CSV with title
rows, Delly VCF, Delly cov.gz + header-less seg.bed, "2+1" CN strings with Mb
positions, xlsx) and checks that the profiler and the builder handle them.
"""
from __future__ import annotations

import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True  # importing the scripts must not litter the skill folder

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "plugins", "shown", "skills", "csv-report", "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, os.path.join(ROOT, "examples"))

import build_report  # noqa: E402
import profile_table  # noqa: E402
import tableio  # noqa: E402

DELLY_VCF = """##fileformat=VCFv4.2
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="PE confidence interval around END">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for POS2">
##INFO=<ID=POS2,Number=1,Type=Integer,Description="Genomic position for CHR2">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Split-read support">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=RDCN,Number=1,Type=Integer,Description="Read-depth based copy-number estimate">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor1\tcontrol1
chr1\t1000000\tDEL00000001\tN\t<DEL>\t600\tPASS\tPRECISE;SVTYPE=DEL;END=1050000;PE=12;SR=8;SVLEN=-50000;CIEND=-5,5\tGT:GQ:RDCN\t0/1:99:1\t0/0:99:2
chr2\t5000000\tDUP00000001\tN\t<DUP>\t300\tPASS\tIMPRECISE;SVTYPE=DUP;END=5800000;PE=7;SVLEN=800000\tGT:GQ:RDCN\t0/1:80:3\t0/0:90:2
chr3\t2000000\tBND00000001\tN\tN[chr8:3000000[\t900\tPASS\tPRECISE;SVTYPE=BND;CHR2=chr8;POS2=3000000;PE=20;SR=11\tGT:GQ:RDCN\t0/1:99:2\t0/0:99:2
chr5\t700000\tINV00000001\tN\t<INV>\t40\tLowQual\tIMPRECISE;SVTYPE=INV;END=900000;PE=3;SVLEN=200000\tGT:GQ:RDCN\t0/1:20:2\t0/1:15:2
"""


def run(args):
    return subprocess.run([sys.executable] + args, capture_output=True, text=True)


class FixtureMixin:
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, name, text, encoding="utf-8", gz=False):
        path = os.path.join(self.d, name)
        data = text.encode(encoding)
        if gz:
            data = gzip.compress(data)
        with open(path, "wb") as fh:
            fh.write(data)
        return path

    def load(self, path):
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)

    def mapping(self, obj, name="mapping.json"):
        path = os.path.join(self.d, name)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(obj, fh, ensure_ascii=False)
        return path


class TableIOTests(FixtureMixin, unittest.TestCase):
    def test_gbk_csv_with_title_rows(self):
        text = "临床信息导出\n\n编号,姓名,年龄,性别,肿瘤纯度\nA01,张三,61,男,55%\nA02,李四,48,女,72%\n"
        p = self.write("c.csv", text, encoding="gb18030")
        t = tableio.read_table(p)
        self.assertEqual(t.columns, ["编号", "姓名", "年龄", "性别", "肿瘤纯度"])
        self.assertEqual(t.header_row, 3)
        self.assertEqual(t.rows[1], ["A02", "李四", "48", "女", "72%"])

    def test_headerless_bed(self):
        p = self.write("x.seg.bed", "chr1\t0\t1000000\tSEG1\t2.04\nchr1\t1000000\t2000000\tSEG2\t3.1\n")
        t = tableio.read_table(p)
        self.assertEqual(t.columns[:5], ["column_1", "column_2", "column_3", "column_4", "column_5"])
        self.assertEqual(len(t.rows), 2)

    def test_vcf_flattening(self):
        p = self.write("d.vcf", DELLY_VCF)
        t = tableio.read_table(p)
        self.assertIn("INFO_SVTYPE", t.columns)
        self.assertIn("tumor1.GT", t.columns)
        self.assertIn("control1.RDCN", t.columns)
        row = dict(zip(t.columns, t.rows[2]))
        self.assertEqual(row["INFO_CHR2"], "chr8")
        self.assertEqual(row["INFO_PRECISE"], "TRUE")
        self.assertEqual(row["INFO_SR"], "11")

    def test_xlsx_dates_and_title_row(self):
        import make_fake_data
        import datetime as dt
        p = os.path.join(self.d, "s.xlsx")
        make_fake_data.write_xlsx(p, "Sheet A", [["Title"], ["ID", "Dx", "Purity"], ["P1", dt.date(2021, 5, 3), 0.6]])
        self.assertEqual(tableio.list_sheets(p), ["Sheet A"])
        t = tableio.read_table(p)
        self.assertEqual(t.columns, ["ID", "Dx", "Purity"])
        self.assertEqual(t.rows[0], ["P1", "2021-05-03", "0.6"])


class ProfileTests(FixtureMixin, unittest.TestCase):
    def test_hints(self):
        p = self.write("s.tsv", "sample\tchr\tstartpos\tendpos\tnMajor\tnMinor\n"
                       + "".join("S%d\tchr%d\t%d\t%d\t2\t1\n" % (i % 3, i % 22 + 1, i * 1000000, i * 1000000 + 999999)
                                 for i in range(1, 60)))
        prof = profile_table.profile(tableio.read_table(p))
        hints = {c["name"]: c["hints"] for c in prof["column_profiles"]}
        self.assertIn("looks like chromosome", hints["chr"])
        self.assertIn("major_cn", hints["nMajor"])
        self.assertIn("patient", hints["sample"])
        self.assertEqual(prof["grouping"]["id_column"], "sample")


class BuildTests(FixtureMixin, unittest.TestCase):
    def test_mixed_inputs_end_to_end(self):
        self.write("clinical.csv",
                   "Patient,Name,Phone,Age,Sex,Tumour cellularity,Subtype\n"
                   "tumor1,Ann,555-1,61,F,0.55,Luminal A\n"
                   "T2,Bob,555-2,48,M,0.72,Basal\n"
                   "T3,Cy,555-3,55,F,0.4,Basal\n", encoding="gb18030")
        self.write("d.vcf", DELLY_VCF)
        cov = "chr\tstart\tend\ttumor1_uniqfrac\ttumor1_logR\ttumor1_CN\n" + "".join(
            "chr1\t%d\t%d\t0.98\t%s\t%s\n" % (i * 10000, (i + 1) * 10000, "0.01" if i % 7 else "NA", "2.02" if i % 7 else "NA")
            for i in range(300))
        self.write("t.cov.gz", cov, gz=True)
        self.write("t.seg.bed", "chr1\t0\t1500000\tSEG1\t2.04\nchr1\t1500000\t3000000\tSEG2\t3.1\n")
        # wide export: CN as "2+1", positions in Mb, patient attributes repeated on each row
        self.write("wide.csv", "pid,subtype,chrom,from_mb,to_mb,cn_state\n"
                   "T2,Basal,chr8,0.5,40.1,2+1\nT2,Basal,chr8,40.1,145.0,3+1\nT3,Basal,17,0,25,1+0\n")
        m = self.mapping({
            "title": "t", "language": "en", "genome": "hg38", "hide_attributes": ["Name", "Phone"],
            "tables": [
                {"file": "clinical.csv", "kind": "patients",
                 "columns": {"patient": "Patient", "purity": "Tumour cellularity", "sex": "Sex"}},
                {"file": "d.vcf", "kind": "variants",
                 "columns": {"patient": {"value": "tumor1"}, "chrom": "CHROM", "pos": "POS", "id": "ID",
                             "svtype": "INFO_SVTYPE", "end": "INFO_END", "chrom2": "INFO_CHR2", "pos2": "INFO_POS2",
                             "pe": "INFO_PE", "sr": "INFO_SR", "qual": "QUAL", "filter": "FILTER",
                             "genotype": "tumor1.GT", "cn": "tumor1.RDCN"},
                 "row_filter": [{"column": "tumor1.GT", "op": "not_in", "value": ["0/0", "./."]}]},
                {"file": "t.cov.gz", "kind": "bins",
                 "columns": {"patient": {"value": "tumor1"}, "chrom": "chr", "start": "start", "end": "end",
                             "cn": "tumor1_CN", "logr": "tumor1_logR"}},
                {"file": "t.seg.bed", "kind": "segments",
                 "columns": {"patient": {"value": "tumor1"}, "chrom": "column_1", "start": "column_2",
                             "end": "column_3", "total_raw": "column_5"}},
                {"file": "wide.csv", "kind": "segments", "attributes": ["subtype"],
                 "columns": {"patient": "pid", "chrom": "chrom", "start": "from_mb", "end": "to_mb",
                             "major_cn": "cn_state", "minor_cn": "cn_state"},
                 "transforms": {"start": {"op": "multiply", "value": 1e6}, "end": {"op": "multiply", "value": 1e6},
                                "major_cn": {"op": "regex", "pattern": "^(\\d+)\\D+(\\d+)$", "group": 1},
                                "minor_cn": {"op": "regex", "pattern": "^(\\d+)\\D+(\\d+)$", "group": 2}}},
            ]})
        data = build_report.build(self.load(m), self.d)
        ids = [p["id"] for p in data["patients"]]
        self.assertEqual(ids, ["tumor1", "T2", "T3"])
        keys = [a["key"] for a in data["attributes"]]
        self.assertNotIn("Name", keys)
        self.assertNotIn("Phone", keys)
        self.assertIn("Subtype", keys)
        self.assertIn("subtype", keys)
        # VCF: control-only INV is still kept for tumor1 (0/1), BND normalized with partner chrom
        svtypes = sorted(r[6] for r in data["svs"])
        self.assertEqual(svtypes, ["BND", "DEL", "DUP", "INV"])
        bnd = [r for r in data["svs"] if r[6] == "BND"][0]
        self.assertEqual(data["chroms"][bnd[4]][0], "8")
        # seg.bed total_raw becomes total CN; cov NA rows dropped
        segs = data["segments"]["tumor1"]
        self.assertEqual([s[7] for s in segs], [2.04, 3.1])
        self.assertEqual(len(data["bins"]["tumor1"]["p"]), 300 - len(range(0, 300, 7)))
        # "2+1" strings and Mb positions
        t2 = data["segments"]["T2"]
        self.assertEqual((t2[1][1], t2[1][3], t2[1][4]), (40100000, 3, 1))
        t2p = data["patients"][1]
        self.assertEqual(t2p["metrics"]["n_segs"], 2)
        self.assertEqual(t2p["sex"], "XY")
        # HTML is written and the payload cannot break out of its <script>
        out = os.environ.get("CSV_REPORT_KEEP_HTML") or os.path.join(self.d, "r.html")
        build_report.write_html(data, out)
        html = open(out, encoding="utf-8").read()
        self.assertNotIn("__REPORT_DATA__", html)
        self.assertEqual(html.count("</script>"), 2)

    def test_mapping_errors_suggest_columns(self):
        self.write("p.csv", "PatientID,Purity\nA,0.5\n")
        m = self.mapping({"tables": [{"file": "p.csv", "kind": "patients", "columns": {"patient": "Patient_ID"}}]})
        r = run([os.path.join(SCRIPTS, "build_report.py"), m, "--check"])
        self.assertEqual(r.returncode, 1)
        self.assertIn("'PatientID'", r.stderr)

    def test_percent_purity_and_genome_guess(self):
        self.write("p.csv", "id,purity\nA,55\nB,80\n")
        self.write("s.csv", "id,chr,s,e,cn\nA,1,1,248000000,2\nB,X,1,155000000,1\n")
        m = self.mapping({"tables": [
            {"file": "p.csv", "kind": "patients", "columns": {"patient": "id", "purity": "purity"}},
            {"file": "s.csv", "kind": "segments", "columns": {"patient": "id", "chrom": "chr", "start": "s", "end": "e", "total_cn": "cn"}}]})
        data = build_report.build(self.load(m), self.d)
        self.assertEqual(data["patients"][0]["purity"], 0.55)
        self.assertTrue(any("percent" in w for w in data["meta"]["warnings"]))
        self.assertEqual(data["meta"]["genome"], "hg38")

    def test_example_mapping_builds(self):
        fake = os.path.join(ROOT, "examples", "fake")
        if not os.path.exists(os.path.join(fake, "patients.xlsx")):
            self.skipTest("run examples/make_fake_data.py first")
        r = run([os.path.join(SCRIPTS, "build_report.py"), os.path.join(ROOT, "examples", "mapping.json"), "--check"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("patients: 40", r.stdout)


class ConsistencyTests(FixtureMixin, unittest.TestCase):
    """Warnings that point at likely mapping mistakes (found by a blind test with CNVkit-style data)."""

    def segs(self, rows):
        return self.write("s.csv", "id,chr,s,e,maj,min\n" + "".join("%s,%s,%d,%d,%d,%d\n" % r for r in rows))

    def build_with(self, patients_csv, cols, seg_rows):
        self.write("p.csv", patients_csv)
        self.segs(seg_rows)
        m = self.mapping({"genome": "hg19", "tables": [
            {"file": "p.csv", "kind": "patients", "columns": cols},
            {"file": "s.csv", "kind": "segments", "columns": {"patient": "id", "chrom": "chr", "start": "s", "end": "e",
                                                               "major_cn": "maj", "minor_cn": "min"}}]})
        return build_report.build(self.load(m), self.d)

    def test_ploidy_contradicting_segments(self):
        diploid = [("A", "1", 1, 249000000, 1, 1), ("A", "2", 1, 243000000, 1, 1)]
        data = self.build_with("id,ploidy\nA,3.9\n", {"patient": "id", "ploidy": "ploidy"}, diploid)
        self.assertTrue(any("mapped ploidy differs" in w for w in data["meta"]["warnings"]))
        data = self.build_with("id,ploidy\nA,2.0\n", {"patient": "id", "ploidy": "ploidy"}, diploid)
        self.assertFalse(any("mapped ploidy differs" in w for w in data["meta"]["warnings"]))

    def test_male_x_at_autosomal_level(self):
        rows = [("A", "1", 1, 249000000, 1, 1), ("A", "X", 1, 155000000, 1, 1)]
        data = self.build_with("id,sex\nA,M\n", {"patient": "id", "sex": "sex"}, rows)
        self.assertTrue(any("chrX is at the autosomal" in w for w in data["meta"]["warnings"]))
        rows = [("A", "1", 1, 249000000, 1, 1), ("A", "X", 1, 155000000, 1, 0)]
        data = self.build_with("id,sex\nA,M\n", {"patient": "id", "sex": "sex"}, rows)
        self.assertFalse(any("chrX is at the autosomal" in w for w in data["meta"]["warnings"]))

    def test_borderline_wgd(self):
        # 2+0 (WGD state) covers barely more than 1+1 -> borderline
        rows = [("A", "1", 1, 125000000, 2, 0), ("A", "1", 125000001, 249000000, 1, 1),
                ("A", "2", 1, 121000000, 1, 1), ("A", "2", 121000001, 243000000, 2, 0)]
        data = self.build_with("id\nA\n", {"patient": "id"}, rows)
        self.assertTrue(data["patients"][0]["metrics"].get("WGD_borderline"))
        self.assertTrue(any("borderline" in w for w in data["meta"]["warnings"]))

    def test_overlapping_sources(self):
        rows = [("A", "1", 1, 249000000, 1, 1), ("A", "1", 1000000, 50000000, 2, 1)]
        data = self.build_with("id\nA\n", {"patient": "id"}, rows)
        self.assertTrue(any("segments overlap" in w for w in data["meta"]["warnings"]))

    def test_facets_total_and_minor_only(self):
        self.write("f.tsv", "ID\tchrom\tstart\tend\ttcn.em\tlcn.em\nA\t23\t1\t1000000\t3\t1\nA\t1\t1\t2000000\t2\tNA\n")
        m = self.mapping({"genome": "hg19", "tables": [{"file": "f.tsv", "kind": "segments", "columns": {
            "patient": "ID", "chrom": "chrom", "start": "start", "end": "end", "total_cn": "tcn.em", "minor_cn": "lcn.em"}}]})
        data = build_report.build(self.load(m), self.d)
        segs = data["segments"]["A"]
        x = [s for s in segs if data["chroms"][s[0]][0] == "X"][0]
        self.assertEqual((x[3], x[4], x[7]), (2, 1, 3))
        one = [s for s in segs if data["chroms"][s[0]][0] == "1"][0]
        self.assertEqual((one[3], one[4], one[7]), (None, None, 2))


class ProfileHintTests(FixtureMixin, unittest.TestCase):
    def test_identifiers_and_arm_layout(self):
        arms = ["%d%s" % (c, a) for c in range(1, 8) for a in "pq"]
        p = self.write("arm.csv", "Sample,Patient name,MRN," + ",".join(arms) + "\n"
                       + "S1,Ann,M1," + ",".join("0" for _ in arms) + "\n"
                       + "S2,Bob,M2," + ",".join("-1" for _ in arms) + "\n"
                       + "S3,Cy,M3," + ",".join("1" for _ in arms) + "\n")
        prof = profile_table.profile(tableio.read_table(p))
        hints = {c["name"]: c["hints"] for c in prof["column_profiles"]}
        self.assertIn("patient", hints["Sample"])
        self.assertTrue(hints["Patient name"][0].startswith("DIRECT IDENTIFIER"))
        self.assertTrue(hints["MRN"][0].startswith("DIRECT IDENTIFIER"))
        self.assertTrue(any("relative call" in h for h in hints["1p"]))
        self.assertTrue(any("chromosome arms" in h for h in prof["layout_hints"]))

    def test_scripts_do_not_write_bytecode(self):
        cache = os.path.join(SCRIPTS, "__pycache__")
        shutil.rmtree(cache, ignore_errors=True)
        p = self.write("x.csv", "a,b\n1,2\n")
        r = run([os.path.join(SCRIPTS, "profile_table.py"), p])
        self.assertEqual(r.returncode, 0, r.stderr)
        m = self.mapping({"tables": [{"file": "x.csv", "kind": "patients", "columns": {"patient": "a"}}]})
        r = run([os.path.join(SCRIPTS, "build_report.py"), m, "--check"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("warnings: none", r.stdout)
        self.assertFalse(os.path.exists(cache))


if __name__ == "__main__":
    unittest.main()
