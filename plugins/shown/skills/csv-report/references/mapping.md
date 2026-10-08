# mapping.json reference

`build_report.py` reads one JSON file. Relative paths inside it resolve
against the folder that holds the mapping file.

```jsonc
{
  "title": "Cohort A - copy number & SV report",   // report heading (English)
  "title_zh": "队列 A 拷贝数与结构变异报告",          // optional Chinese variant for the 中文 toggle
  "subtitle": "40 patients · ASCAT + Delly",       // optional; also "subtitle_zh"
  "language": "en",                                 // default view: "en" (default) | "zh"; the page has a toggle
  "genome": "hg38",                                 // "hg19" | "hg38" | "chm13" | "auto" | {"name": "...", "chroms": [["1", 248956422], ...]}
  "group_by": "Cancer type",                        // attribute key (source column name) used for colours/sorting
  "y_limit": 5,                                     // ASCAT profile y-axis cap (ASCAT default 5)
  "max_points_per_patient": 40000,                  // bins are down-sampled above this
  "hide_attributes": ["Name"],                      // attribute keys never shown
  "insights": ["..."],                              // optional AI-written notes (from --check summary)
  "insights_zh": ["..."],                           // the same notes in Chinese
  "tables": [ /* one entry per input file (or glob) */ ],
  "attributes": { /* display config per attribute key */ }
}
```

## tables[]

| key | meaning |
|---|---|
| `file` | path, list of paths, or glob (`"ascat/*.segments.txt"`) |
| `sheet` | xlsx sheet name or 1-based index (default: first sheet) |
| `header_row` | 1-based header row (default: auto-detect, skipping title rows) |
| `delimiter` | force a delimiter (`","`, `"\t"`, `";"`, `"|"`, `" "`) |
| `kind` | `patients` / `segments` / `variants` / `bins` |
| `columns` | canonical field -> source column (see below) |
| `transforms` | canonical field -> transform op or list of ops |
| `attributes` | source columns to lift to patient attributes (non-`patients` tables); for `patients`, restricts the displayed attributes (default: every unmapped column) |
| `exclude` | source columns never used as attributes / extra columns |
| `row_filter` | list of `{"column", "op", "value"}`; ops `in`, `not_in`, `equals`, `not_equals`, `nonempty`, `gt`, `ge`, `lt`, `le` |

A `columns` value is either a column name, or an object:
- `{"column": "Sample"}`: same as the plain name
- `{"value": "P001"}`: a constant, e.g. a single-patient file
- `{"from_filename": "^(.+?)\\.segments"}`: regex group 1 of the file name; use
  it with a glob of per-sample files

### canonical fields per kind

| kind | required | optional |
|---|---|---|
| patients | `patient` | `purity` (0-1; values still above 1.5 after transforms are divided by 100, so a `divide` transform is optional and never applied twice), `ploidy`, `goodness_of_fit` (%), `sex` (XX/XY/男/女/M/F) |
| segments | `patient chrom start end` | `major_cn minor_cn` (integers, ASCAT nMajor/nMinor), `total_cn`, `major_raw minor_raw` (ASCAT nAraw/nBraw), `total_raw` (e.g. Delly float CN), `logr` (segment mean), `baf`, `n_markers` |
| variants | `chrom pos` | `patient end chrom2 pos2 svtype svlen qual filter pe sr genotype id precise alt cn` |
| bins | `chrom` and `pos` or `start` (+`end` -> midpoint) | `patient logr baf cn` |

Notes:
- `svtype` is normalized: DEL, DUP, INV, BND (TRA/CTX/translocation), INS,
  CNV. When it is missing, it is derived from `alt` (`<DEL>`, `N[chr2:...[`) or
  from the Delly ID prefix (`DEL00000012`). Different chromosomes for `chrom`
  and `chrom2` mean BND.
- `svlen` is stored as an absolute size; when missing it is `end - pos`.
- `filter`: `PASS` drives the "PASS only" toggle in the report.
- `variants` without a `patient` column are assigned to patient `sample`.
  Give `patient` a `value` instead.
- `patients` tables:
  - A mapped `sex` column is still shown as an attribute; purity, ploidy and
    goodness of fit are shown in their own places instead.
  - Several `patients` tables are merged by `patient`. For each field and
    attribute, the first non-empty value wins, in table order. A patient
    listed twice in the same table keeps its first row (with a warning).
- `ploidy` sets the gain/loss baseline (round(ploidy)). Leave it unmapped to
  derive it from the segments (length-weighted mean of the total CN).
  The builder warns when a mapped ploidy is more than 0.5 away from that.
- `sex` = XY halves the X/Y baseline. The builder warns when chrX in a male is
  at the autosomal level, which means the caller did not adjust X for sex.
  Patients without a known sex (XX/XY) have X/Y left out of gain/loss calls.
- `variants`: up to 12 unmapped columns are kept as extra columns in the SV
  table. Use `exclude` to drop the uninteresting ones.

## transforms

Each op is an object, applied in order to the raw cell (string) before
numeric parsing:

| op | example | effect |
|---|---|---|
| `regex` | `{"op": "regex", "pattern": "^(\\d+)\\+(\\d+)$", "group": 1}` | keep a capture group |
| `map` | `{"op": "map", "values": {"F": "XX", "M": "XY"}, "default": null}` | value lookup (case-insensitive fallback) |
| `replace` | `{"op": "replace", "old": ",", "new": ""}` | string replace |
| `strip_prefix` | `{"op": "strip_prefix", "value": "TCGA-"}` | |
| `multiply` / `divide` / `add` | `{"op": "multiply", "value": 1000000}` | numeric |
| `log2`, `abs`, `negate`, `round` | `{"op": "log2"}` / `{"op": "round", "digits": 0}` | numeric |

Splitting a `"2+1"` column into major / minor:

```json
"columns": {"major_cn": "CN state", "minor_cn": "CN state"},
"transforms": {
  "major_cn": {"op": "regex", "pattern": "^(\\d+)\\D+(\\d+)$", "group": 1},
  "minor_cn": {"op": "regex", "pattern": "^(\\d+)\\D+(\\d+)$", "group": 2}
}
```

## attributes{}

The keys are source column names; every key is optional. `label` (default:
the column name), `unit` and `description` are shown in the English view;
`label_zh`, `unit_zh` and `description_zh` are shown after clicking 中文.
Category values are shown as they appear in the data.

```json
"attributes": {
  "年龄":   {"label": "Age", "label_zh": "年龄", "type": "numeric", "unit": "y", "unit_zh": "岁"},
  "Stage": {"label_zh": "分期", "type": "ordinal", "order": ["I", "II", "III", "IV"]},
  "Dx date": {"label": "Diagnosis date", "label_zh": "诊断日期", "type": "date"},
  "Smoker": {"label_zh": "吸烟", "type": "category", "order": ["Yes", "No", "Unknown"]},
  "Notes": {"type": "text", "description": "free text from the clinic"}
}
```

When `type` is omitted it is inferred:
- numbers -> `numeric`, except small integer codes (at most 6 distinct values,
  all between 0 and 10, e.g. 0/1 flags or grades 1-4) -> `category`. Scores
  like HRD 32-42 stay numeric.
- ISO dates -> `date` (Excel date cells are converted)
- few distinct strings -> `category`
- anything else -> `text`

How each type is shown:
- numeric and date -> histogram
- category and ordinal -> bars
- text -> distinct count and top values

## Examples

### Clinical sheet only (any columns)

```json
{
  "title": "Patient overview", "language": "en",
  "tables": [{"file": "patients.xlsx", "kind": "patients",
              "columns": {"patient": "Patient ID", "purity": "Tumor purity (%)"},
              "exclude": ["Name", "Phone"]}]
}
```

### ASCAT output for many samples (per-sample files)

```json
{"file": "ascat/*.segments_raw.txt", "kind": "segments",
 "columns": {"patient": "sample", "chrom": "chr", "start": "startpos", "end": "endpos",
             "major_cn": "nMajor", "minor_cn": "nMinor", "major_raw": "nAraw", "minor_raw": "nBraw"}}
```

Purity and ploidy usually come from a separate `ascat.metrics` / summary table
(`kind: patients`, `purity`, `ploidy`, `goodness_of_fit`). Its remaining
columns (WGD, GI, LOH, ...) appear as attributes automatically.

### Delly VCF (text; convert BCF first: `bcftools view in.bcf > in.vcf`)

The reader flattens the VCF:
- INFO fields become `INFO_<ID>` columns
- per-sample FORMAT fields become `<sample>.<ID>` columns

```json
{"file": "delly.vcf", "kind": "variants",
 "columns": {"patient": {"value": "tumor1"}, "chrom": "CHROM", "pos": "POS", "id": "ID",
             "svtype": "INFO_SVTYPE", "end": "INFO_END", "chrom2": "INFO_CHR2", "pos2": "INFO_POS2",
             "svlen": "INFO_SVLEN", "pe": "INFO_PE", "sr": "INFO_SR", "qual": "QUAL", "filter": "FILTER",
             "precise": "INFO_PRECISE", "genotype": "tumor1.GT", "cn": "tumor1.RDCN"}}
```

For a multi-sample VCF, add one table entry per sample. Each entry uses the
same `file`, a different `patient` value, and that sample's `.GT` column. Add a
`row_filter` (e.g. `{"column": "tumor1.GT", "op": "not_in", "value": ["0/0", "./."]}`)
to keep only that sample's carriers.

### Delly read depth (`delly cnv -c out.cov.gz -u out.seg.bed`)

```json
{"file": "out.cov.gz", "kind": "bins",
 "columns": {"patient": {"value": "tumor1"}, "chrom": "chr", "start": "start", "end": "end",
             "cn": "tumor1_CN", "logr": "tumor1_logR"}},
{"file": "out.seg.bed", "kind": "segments",
 "columns": {"patient": {"value": "tumor1"}, "chrom": "column_1", "start": "column_2",
             "end": "column_3", "total_raw": "column_5"}}
```

`seg.bed` has no header, so its columns are named `column_1` … `column_5`.
