**English** | [简体中文](README.zh-CN.md)

# shown

Give an AI agent patient or cohort tables whose columns are **not fixed** (CSV, TSV, Excel, VCF). You get back **one HTML file** that opens with a double-click and works fully offline. The charts follow the output style of [Delly](https://github.com/dellytools/delly) (structural variants, read-depth copy number) and [ASCAT](https://github.com/VanLoo-lab/ascat) (allele-specific copy number, purity, ploidy).

This repo is a Claude Code / Codex **plugin** containing one skill, `csv-report`:

- **The AI (Claude or Codex)** reads the headers and values and works out what each column means: patient ID, chromosome, position, nMajor/nMinor, SV type, purity, clinical fields and so on. It then writes a `mapping.json`.
- **Fixed scripts** read the data through that mapping, compute ASCAT metrics, draw the charts, and write a single HTML file. They need only the `python3` standard library; there is nothing to install.

The report is in English, with a **中文** button that switches the whole page to Chinese. It also has a light/dark theme toggle.

A sample report is in [`examples/report.html`](examples/report.html); download it and double-click to open. All of its data is simulated.

## Install

Colleagues need access to this repository. If it is private, set up git credentials first.

### Claude Code

```text
/plugin marketplace add zerkariya/shown
/plugin install shown@shown
```

Or from a shell:

```bash
claude plugin marketplace add zerkariya/shown
claude plugin install shown@shown
```

Start a new session (or run `/reload-plugins`) after installing. To update later: `claude plugin marketplace update shown`.

### Codex

```bash
codex plugin marketplace add zerkariya/shown
codex plugin add shown@shown
```

Codex reads this repo's `.claude-plugin/marketplace.json` and `plugins/shown/.claude-plugin/plugin.json` directly, so both tools share the same files.

### Without the plugin system (either tool)

Copy the skill folder into your personal skills directory:

```bash
git clone https://github.com/zerkariya/shown.git
# Claude Code
cp -r shown/plugins/shown/skills/csv-report ~/.claude/skills/
# Codex
mkdir -p ~/.agents/skills && cp -r shown/plugins/shown/skills/csv-report ~/.agents/skills/
```

## Use

Describe what you want in Claude Code or Codex, for example:

> Turn `data/patients.xlsx` and `data/sv.csv` into a visual report

In Claude Code you can also invoke it explicitly: `/shown:csv-report data/patients.xlsx`. In Codex, mention the skill as `$csv-report`.

The agent then:

1. runs `profile_table.py` on every column
2. writes `report/mapping.json`
3. runs `build_report.py --check` and fixes any errors
4. builds `report/report.html`
5. tells you how each table was mapped, what it assumed, and which sensitive columns it hid

Keep `mapping.json`; when the data is updated you can rebuild with it directly.

### Writing the prompt

A one-liner works, and the agent asks when something is ambiguous. Each of
these lines saves a round trip, most useful first:

```text
Use the csv-report skill to turn these into a visual report:
- files: data/clinical.xlsx, data/ascat/*.segments.txt, data/delly_sv.csv
- what they are: clinical sheet; ASCAT segments (one file per sample); Delly SVs (bcftools query export)
- genome: hg38
- patient ID: "Sample" in the clinical sheet, the sample column in the ASCAT files, "SAMPLE" in the Delly table
- colour groups by "Histology"
- leave out: patient names, MRN, phone numbers
- title "XX cohort copy-number report", with a Chinese version
- write it to report/ and tell me how you mapped each file and what you assumed
```

### What input works

Column names, column order and the number of columns are free. What matters
is what one row stands for:

| to get | one row is | needs at least | nice to have |
|---|---|---|---|
| cohort overview | a patient | patient ID | any clinical / analysis columns (each gets a chart) |
| purity vs ploidy | a patient | patient ID, purity, ploidy | goodness of fit |
| ASCAT profile, CN heatmap | a CN segment | patient ID, chromosome, start, end, copy number | nMajor/nMinor, unrounded values, logR |
| SV charts | an SV | chromosome, position (+ patient ID) | end / second breakpoint, type, FILTER, read support |
| LogR / BAF / depth tracks | a bin or SNP | chromosome, position | logR, BAF, CN estimate |

Rules that do matter:
- **Header:** one header row. Title rows above it are fine; merged or
  multi-row headers are not.
- **IDs:** the same patient ID across files, or derivable from file names.
- **Coordinates:** all on one genome build.
- **Formats:** CSV / TSV / TXT / `.xlsx` / text VCF, optionally gzipped.
  - Convert `.xls` to `.xlsx` first.
  - Convert BCF with `bcftools view` first.
- **Other shapes** (one column per sample or per chromosome arm, gene-level
  tables) go through a small conversion script the agent writes.
- **Not drawn as dedicated charts yet:** gene-level results (expression,
  MAF mutations, enrichment) and survival curves.

### Run it by hand (no AI)

```bash
SKILL=plugins/shown/skills/csv-report
python3 $SKILL/scripts/profile_table.py data/*.csv data/*.xlsx      # column types and role hints
python3 $SKILL/scripts/build_report.py my_mapping.json --check      # validate the mapping
python3 $SKILL/scripts/build_report.py my_mapping.json -o report.html
```

The mapping format is documented in [`references/mapping.md`](plugins/shown/skills/csv-report/references/mapping.md). Ready-made mappings for native Delly / ASCAT output are in [`references/formats.md`](plugins/shown/skills/csv-report/references/formats.md).

## Supported data

| `kind` | one row is | examples |
|---|---|---|
| `patients` | a patient | Clinical sheets with any number of columns; every unmapped column gets a chart. ASCAT metrics / summary tables also fit here. |
| `segments` | a copy-number segment | ASCAT `*.segments.txt` / `segments_raw.txt`, Delly `seg.bed`, or any table with chromosome / start / end / copy number |
| `variants` | an SV / CNV call | Delly VCF (text), `bcftools query` exports, any SV list |
| `bins` | a bin / SNP | Delly `cov.gz`, ASCAT `Tumor_LogR.txt` / `Tumor_BAF.txt` |

- **Handled automatically:**
  - GBK and UTF-8 encodings
  - title rows above the header and note rows below the table in Excel
  - Excel dates
  - decimal commas (`0,55`)
  - chromosome spellings (`chr1`, `1`, `23`)
  - SV type synonyms (`deletion`, `TRA`, `<DEL>`, `缺失`)
  - purity written as a percentage
  - genome build inference (hg19 / hg38 / CHM13)
- **Handled by AI-written `transforms`:** positions in Mb, copy number packed into one cell as `"2+1"`, and similar irregular encodings.

## What the report shows

- **Overview**
  - stat tiles
  - AI notes (clearly labelled as AI-written)
  - ASCAT purity vs ploidy scatter
  - one distribution chart per patient field, however many fields there are
  - a sortable, searchable patient table
- **Copy number**
  - cohort gain/loss frequency across the genome
  - patients × genome heatmap
- **Structural variants**
  - SVs per patient, stacked by Delly SVTYPE
  - size and per-chromosome distributions
  - a call table you can filter by PASS and type
- **Patient**
  - ASCAT copy-number profile, on real chromosome coordinates:
    - rounded: red = nMajor, blue = nMinor
    - unrounded: purple = total CN, green = nMinor
  - Delly read-depth CN: black dots + green segment line
  - LogR / BAF: red points + blue segmented / fitted values
  - SV arcs
  - drag to zoom, double-click to reset
- **Data & mapping**
  - which column of which file fed each field
  - unused columns and warnings
  - metric definitions (WGD / GI / LOH as in `ascat.metrics`)
  - the raw `mapping.json`, so anyone can audit the AI's decisions

## Privacy

The HTML embeds the data it plots, so sharing a report means sharing that data. By default the agent leaves direct identifiers out of the report: names, national ID numbers, phone numbers, addresses, medical record numbers. You can double-check with `exclude` / `hide_attributes` in `mapping.json`.

## Repository layout

```
.claude-plugin/marketplace.json        # marketplace manifest (read by Claude Code and Codex)
plugins/shown/.claude-plugin/plugin.json
plugins/shown/skills/csv-report/
  SKILL.md                             # workflow instructions for the agent
  scripts/tableio.py                   # CSV/TSV/XLSX/VCF reader (stdlib only)
  scripts/profile_table.py             # column profiles + role hints
  scripts/build_report.py              # mapping -> data -> HTML
  assets/report_template.html          # report template (vanilla JS, no external dependencies)
  references/mapping.md                # mapping.json reference
  references/formats.md                # Delly / ASCAT output formats
examples/                              # fake-data generator, sample mapping, sample report
tests/                                 # python3 -m unittest discover -s tests
```

Regenerate the example:

```bash
python3 examples/make_fake_data.py
python3 plugins/shown/skills/csv-report/scripts/build_report.py examples/mapping.json -o examples/report.html
```
