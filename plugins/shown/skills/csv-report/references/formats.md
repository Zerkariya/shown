# Delly and ASCAT outputs: what they look like and how to map them

Use this to recognize native tool output among the user's files. Their own
CSVs are often re-exports of these with renamed or extra columns.

## Delly (github.com/dellytools/delly): SV and read-depth CNV caller

**SV calls** come as VCF/BCF, one record per SV.
- **Fixed columns:** `CHROM POS ID REF ALT QUAL FILTER`
  - ID: `DEL00000001`, `DUP…`, `INV…`, `BND…`, `INS…`, `CNV…`
  - ALT: `<DEL> <DUP> <INV> <INS> <BND> <CNV>`
  - FILTER: `PASS` or `LowQual`
- **INFO:**
  - type and position: `SVTYPE END SVLEN`, plus `CHR2 POS2` for inter-chromosomal BND
  - support: `PE SR` (paired-end / split-read support), `MAPQ SRMAPQ SRQ`
  - breakpoint: `CT` (connection type), `PRECISE`/`IMPRECISE`, `CIPOS CIEND`, `HOMLEN`
  - other: `CONSENSUS`, `SUBTYPE`
  - after `delly filter`: `SOMATIC RDRATIO AF …`
- **Per sample (FORMAT):**
  - `GT GQ FT`
  - `DR DV`: reference / variant read pairs
  - `RR RV`: reference / variant junction reads
  - `RC RCL RCR RDCN`
- Users usually flatten with `bcftools query`, e.g.
  `-f "%CHROM\t%POS\t%INFO/END\t%ID[\t%RDCN]\n"`. That gives columns
  `CHROM POS END ID <sample RDCN…>`.
- Map to `variants`: `chrom pos id svtype end chrom2 pos2 svlen pe sr qual
  filter genotype precise cn(RDCN)`.

**CNV mode** (`delly cnv -c out.cov.gz -u out.seg.bed -o cnv.bcf`):
- `out.cov.gz` (gzip, header):
  - columns: `chr start end <sample>_uniqfrac <sample>_logR <sample>_CN`
  - one row per window; `NA` in unmappable windows
  - map to `bins` (`cn`, `logr`)
- `out.seg.bed` (no header): `chr start end SEG<n> cn`, where `cn` is a float.
  Map to `segments` with `total_raw`.
- CNV calls in the BCF: `SVTYPE=CNV`, `END`, FORMAT `CN` (integer), `RDCN`
  (float), `CNL`, `RDSD`.

**Delly's own plots (R scripts):**
- `R/rd.R`
  - genome-wide copy-number scatter, one dot per window, black, alpha 0.2
  - green (`#31a354`) segment lines
  - y axis: 0–12 copies
  - one facet per chromosome
- `R/cnv.R`: per-CNV histograms of RDCN across samples
- `R/gcbias.R`: GC bias

The report's "Delly 读深拷贝数" track reproduces `rd.R`.

## ASCAT (github.com/VanLoo-lab/ascat): allele-specific copy number, purity, ploidy

**Text outputs** (tab-separated, with header):
- `<sample>.segments.txt`: `sample chr startpos endpos nMajor nMinor`
  (integers, rounded profile)
- `<sample>.segments_raw.txt`: adds `nAraw nBraw` (unrounded)
- `<sample>.LogR.PCFed.txt`, `<sample>.BAF.PCFed.txt`: segmented values per SNP
- inputs `Tumor_LogR.txt` / `Tumor_BAF.txt` (and `Germline_*`):
  - layout: `<probe id> Chr Position <sample1> <sample2> …`
  - the first header cell is empty; the reader names it `probe_id`
- per sample (`ascat.output`): `purity` (aberrant cell fraction), `ploidy`,
  `psi`, `goodnessOfFit` (%), non-aberrant / failed flags
- `ascat.metrics()` table, one row per sample:
  - QC: `sex tumour_mapd normal_mapd GC_correction_* RT_correction_*`
  - segmentation: `n_het_SNP n_segs_logR n_segs_BAF n_segs_logRBAF_diff frac_homo`
  - fit: `purity ploidy goodness_of_fit`
  - segment sizes: `size_intermediate_segments size_odd_segments n_segs segs_size n_segs_1kSNP`
  - homozygous deletions: `homdel_segs homdel_largest homdel_size homdel_fraction`
  - genome state: `LOH mode_minA mode_majA WGD GI`
  - the released TCGA summary also has `name patient cancer_type barcodeTumour barcodeNormal … QC`

**Mapping:**
- segments -> `segments` (`major_cn`=nMajor, `minor_cn`=nMinor,
  `major_raw`=nAraw, `minor_raw`=nBraw)
- metrics or summary -> `patients` (`purity`, `ploidy`, `goodness_of_fit`; the
  other columns become attributes)
- LogR/BAF -> `bins`, one table entry per sample column

**ASCAT plots, reproduced in the report:**
- `ASCATprofile.png`
  - rounded profile: red `#E03546` nMajor drawn at −0.1, blue `#3557E0`
    nMinor at +0.1
  - y axis 0–5; values above are capped and drawn paler
  - title `Ploidy: 2.33, purity: 58%, goodness of fit: 87.0%`
- `rawprofile.png`: unrounded, purple `#943CC3` total, green `#60AF36` minor
- `adjustedASCATprofile.png`: the same profiles on real chromosome lengths.
  The report always uses real coordinates.
- `tumour.png` / `ASPCF.png`
  - LogR (−2..2) and BAF (0..1) panels
  - red raw points, blue `#1b38ae` segmented values; BAF mirrored (b and 1−b)
- `sunrise.png`: purity × ploidy goodness-of-fit heatmap. It needs the
  distance matrix, which is not in tabular outputs. The report shows a
  purity-vs-ploidy cohort scatter instead.

**Metrics the builder recomputes from segments** (`ascat.metrics` definitions;
sex chromosomes excluded):
- `LOH`: fraction with nMinor = 0
- `mode_majA` / `mode_minA`: length-weighted mode, capped at 5
- `WGD`: mode_majA 1 -> 0, 2 -> 1, 3–5 -> "1+"
- `GI`: 1 − fraction at the WGD baseline (1+1 or 2+2)
- `homdel_*`: nMajor = nMinor = 0
- `FGA` (extra): fraction whose total CN differs from its mode
