---
name: csv-report
description: Turn patient or cohort tables whose columns are not known in advance (CSV, TSV, TXT, Excel .xlsx, VCF; gzipped too) into one self-contained offline HTML report styled after Delly (structural variants, read-depth copy number) and ASCAT (allele-specific copy-number profile, purity, ploidy). Use when someone wants to visualize, chart, summarize or build a report/dashboard from patient information, clinical tables, CNV / copy-number segments, SV calls, logR/BAF, or ASCAT / Delly results.
---

# csv-report: any table -> Delly / ASCAT style HTML report

The column layout of the input is **not fixed**. Your job is the part a script
cannot do: read the headers and values (any language, any naming), decide what
each column means, and write a `mapping.json`. Two bundled scripts do the rest
deterministically:

- `scripts/profile_table.py`: prints, for every column, its type, missing and
  unique counts, a value summary, and *hints* about the column's likely role.
- `scripts/build_report.py`: reads `mapping.json`, normalizes the data, computes
  ASCAT-style metrics, and writes one HTML file. The file has no network
  dependencies and opens by double-click.

The scripts live next to this file. In Claude Code this directory is
`${CLAUDE_SKILL_DIR}`. In other agents (Codex etc.) use the directory that
contains this SKILL.md. Below, `$SKILL` stands for that directory. Only
`python3` (3.8+) is needed; there are no packages to install.

## Workflow

1. **Find the inputs.** Use the files the user gave or pointed to. If none were
   given, look in the working directory for `*.csv *.tsv *.txt *.xlsx *.vcf(.gz)`.
   Ask only if that is ambiguous. Never modify the input files. Put everything
   you write (mapping, helper scripts, report) in an output folder such as
   `./report/`.

2. **Profile every input file:**
   ```bash
   python3 "$SKILL/scripts/profile_table.py" data/*.csv data/*.xlsx
   ```
   Every sheet of an .xlsx is profiled. The header row is auto-detected (title
   rows above it are skipped); override with `--header-row N`. Note rows at
   the end of a table ("End of sheet …", "注：…") are dropped, and decimal
   commas (`0,55`) are read as decimal points. The profile lists both as
   `note:` lines. Read the hints, but judge from the values: hints are
   heuristics.

3. **Decide what each table is and write `mapping.json`.** The full schema with
   examples is in `references/mapping.md`; read it the first time. Each input
   file becomes one entry in `tables` with a `kind`:

   | kind | one row is | required fields |
   |---|---|---|
   | `patients` | a patient / sample (clinical info, purity, ploidy, any extra columns) | `patient` |
   | `segments` | a copy-number segment (ASCAT `segments.txt`, Delly `seg.bed`, CNVkit, any CN table) | `patient chrom start end` + one of `major_cn`+`minor_cn` / `total_cn` / `logr` |
   | `variants` | an SV / CNV call (Delly VCF or `bcftools query` export, any SV table) | `chrom pos` |
   | `bins` | a genomic bin or SNP (Delly `cov.gz`, ASCAT LogR/BAF) | `chrom` + `pos` or `start` |

   How to decide:
   - Map *meaning*, not names. For example, `患者编号`, `Sample_ID`, `case` and
     `Tumor_Sample_Barcode` are all `patient`. `肿瘤纯度`, `cellularity` and
     `ACF` are all `purity`.
   - Every column of a `patients` table that you do not map becomes a displayed
     attribute automatically. Configure labels, types (`numeric`, `category`,
     `ordinal` with an `order`, `date`, `text`) and units in the top-level
     `attributes` object.
   - When patient-level columns repeat on every row of a genomic table (a
     "wide" export), list them in that table's `attributes`. They are lifted to
     the patient.
   - Fix units with `transforms`:
     - purity given in percent: `{"op": "divide", "value": 100}`. Optional:
       the builder also divides by 100 when purity is still above 1.5 after
       your transforms, so it never divides twice.
     - positions given in Mb: `multiply` by 1e6
     - a combined `"2+1"` CN cell or a `chr1:100-200` locus: `regex` with a
       `group` (the same source column can feed two fields)
     - ratios that should be logR: `log2`
   - Chromosome spelling (`chr1`, `1`, `23`, `chrX`) and SV type synonyms
     (`deletion`, `TRA`, `tandem_dup`, `<DEL>`, `缺失`) are normalized
     automatically.
   - **`ploidy` and `sex` change the copy-number calls.** Gain / loss / amp in
     the heatmap, frequency plot and summary are relative to round(ploidy);
     mapped `sex` = XY halves that baseline on X/Y. So:
     - map `ploidy` only from the same analysis as the segments (e.g. ASCAT's
       own ploidy). A ploidy from another pipeline that disagrees with the
       segments makes most of the genome look gained or lost. When in doubt,
       leave it unmapped and it is derived from the segments.
     - map `sex` only if the CN caller adjusted chrX for sex (ASCAT does).
       Otherwise males show a false X gain.
     - with `sex` unmapped, X and Y are not called at all; they show as no
       data instead of a false male "loss". The builder says so.
     - The builder warns about all of these; take those warnings seriously.
   - **One copy-number source per patient.** If two tables give segments for
     the same patients (e.g. CNVkit and ASCAT), map the one the user cares
     about, or ask. Overlapping segments trigger a warning.
   - Set `genome` (`hg19` / `hg38` / `chm13`) when you can tell; otherwise the
     builder infers it from coordinates and warns.
   - Set `group_by` to the clinical attribute that best splits the cohort (2-6
     groups, e.g. cancer type). It colours the purity/ploidy scatter and sorts
     the heatmap.
   - The report is English-first, with a one-click 中文 toggle in the page.
     - Write `title`, `subtitle`, attribute `label`s and `insights` in English.
     - Add `title_zh`, `subtitle_zh`, `label_zh` (and `unit_zh`) and
       `insights_zh` so the Chinese view is complete as well.
     - Set `"language": "zh"` only when the user wants Chinese as the default
       view.
   - **Privacy:** never display direct identifiers (patient names, national ID
     or phone numbers, addresses, MRNs) as attributes, and never use a name
     column as the patient ID. The profiler marks such columns
     `DIRECT IDENTIFIER?`. Put them in `exclude` / `hide_attributes`, unless
     the user explicitly asks for them. Tell the user which columns you left
     out.
   - Note every decision you are not sure about while you map: an ambiguous ID
     column, a unit, an opaque column name you interpreted, a preprocessing
     choice. Step 5 asks the user about them; do not settle them silently.

   `references/formats.md` describes the native Delly, ASCAT and CNVkit outputs
   and their ready-made mappings.

4. **Validate and iterate:**
   ```bash
   python3 "$SKILL/scripts/build_report.py" report/mapping.json --check
   ```
   Errors name the problem and suggest close column names. Fix the mapping and
   re-run until it passes. Then read the printed `warnings` and fix what is
   fixable. It always prints either a list or `warnings: none`. Warnings cover:
   - dropped rows and unknown contigs
   - percent purity and the genome guess
   - patients missing from a table
   - ploidy or sex contradicting the segments
   - borderline WGD calls
   - overlapping segments

   It then prints `questions for the user`, decisions the data checks cannot
   settle (ploidy vs segments, sex, overlapping sources, ID mismatches, genome
   build, …), each with a suggested default. Add `--stats-json
   report/summary.json` to keep the summary, warnings and questions as a file.

5. **Confirm with the user before building. This step is required.** Send
   ONE message, in the user's language, and wait for the answer:
   - **What you understood:** one line per file: what it is, what one row
     is, the key columns → roles, units, genome build.
   - **What you will leave out:** direct identifiers you will hide, columns
     you could not interpret, data that cannot be drawn.
   - **Numbered questions, each with your suggested default.** Include:
     - every item under `questions for the user` from step 4;
     - every uncertain decision you noted in step 3: ID columns that do not
       match across files, unclear units, opaque columns (`V3`, `F12`, pinyin
       abbreviations) whose meaning you inferred, how to handle a shape that
       needs preprocessing, and which attribute to `group_by` when several fit.
   - End with: reply with changes, or "OK" to accept all defaults.

   Rules:
   - Keep it short. Merge related points and ask at most about 7 questions.
     Skip what the user already told you or what the data settles beyond
     doubt.
   - Use your environment's structured question tool if it has one (e.g.
     `AskUserQuestion` in Claude Code, which takes up to 4 questions per
     call); otherwise ask in plain text. Do not build the report (step 7)
     until the user has answered.
   - Apply the answers to `mapping.json` and re-run `--check`. Ask again only
     about new questions the answers created.
   - Skip the confirmation only when the user said not to ask ("just do it",
     "不用问") or nobody can answer (a non-interactive run). Then use the
     defaults and list every one of them as an assumption in step 8.

6. **Optional: insights.** From the printed `summary` only, write 3-6 short,
   factual observations into `"insights"` (e.g. "chr8 gain in 43% of
   patients; WGD in 55%"), with the same points in Chinese in
   `"insights_zh"`. Do not speculate beyond the numbers and give no clinical
   advice; the report labels them as AI-written.
   - `gain_loss_rule` in the summary states exactly what `chromosome_gains` /
     `chromosome_losses` measure. Reuse its wording rather than saying "arm"
     or "focal".
   - Mention it when `WGD_borderline` is above 0.

7. **Build:**
   ```bash
   python3 "$SKILL/scripts/build_report.py" report/mapping.json -o report/report.html
   ```

8. **Report back** in the user's language with:
   - the output path
   - a one-line description of each table and how you mapped it
   - what the user confirmed, and any defaults you used without confirmation
   - the remaining warnings
   - how to open the report: double-click, works offline, and the file
     contains the data, so share it as carefully as the data itself

## When the shape does not fit

If a table cannot be expressed as one of the four kinds, write a small
preprocessing script in the output folder that produces a tidy CSV, then map
that CSV. Keep the script next to the report so the result can be reproduced.
Examples:
- CN encoded in a sheet layout
- several patients side by side in column blocks
- one column per chromosome arm (see below)

**Arm-level or other relative calls** (one row per patient, columns like
`1p 1q 2p …`, values −1/0/+1 or loss/neutral/gain; the profiler reports
`layout: wide layout …`):
- These are relative states, not copy numbers. Do not invent `total_cn`
  from them.
- Default: summarize them per patient (counts of gained / lost / altered arms,
  plus the lists of arms) in a tidy CSV, and map it as an extra `patients`
  table.
- Draw them on the genome as segments only when no real segment table exists
  and the user wants that. Arm boundaries need centromere positions for the
  right genome build (e.g. the UCSC cytoBand file). Ask the user for that
  file rather than guessing coordinates.
- If segments for the same patients also exist, check whether the two
  sources agree, and tell the user.

ASCAT's `Tumor_LogR.txt` / `Tumor_BAF.txt` (one column per sample) do **not**
need this: add one `bins` entry per sample with `"patient": {"value": "S1"}` and
`"logr": "S1"`.

## What the report contains

- **Overview (概览)**
  - stat tiles
  - AI notes
  - purity vs ploidy scatter
  - one chart per patient attribute, whatever the columns are
  - a sortable, searchable patient table
- **Copy number (拷贝数)**
  - cohort gain/loss frequency across the genome
  - a patients × genome heatmap (relative to round(ploidy))
- **Structural variants (结构变异)**
  - SVs per patient stacked by type
  - size and per-chromosome distributions
  - a filterable table of calls (VCF fields)
- **Patient (患者详情)**
  - ASCAT profile: rounded (red nMajor, blue nMinor) or unrounded (purple
    total, green minor), on real chromosome coordinates like
    `ascat.plotAdjustedAscatProfile`
  - Delly read-depth CN track: dots, plus a green segment line
  - LogR and BAF tracks: red points, blue segmented / fitted values (as in
    ASCAT's ASPCF plot)
  - SV arcs
  - zoom: drag on a track, or pick a chromosome
- **Data & mapping (数据与映射)**
  - which file and column fed which field
  - unused columns and warnings
  - metric definitions
  - the raw `mapping.json`, so anyone can audit what the AI decided

Every page has an English / 中文 toggle and a light / dark theme toggle.

Metrics follow `ascat.metrics`: WGD, GI, LOH, homdel_*, mode_majA/minA, and
FGA. They are computed from segments when allele-specific CN is available.
