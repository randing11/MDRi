# MDRi Explorer

MDRi is a single-cell-informed pan-cancer myeloid damage-response scoring
framework with a companion Shiny explorer.

It computes four scores from transcriptomic data:

- `MDR_injury`: inflammatory injury, inflammasome, hypoxia and remodeling.
- `MDR_resolution`: lipid/efferocytosis, resident repair and heme/iron handling.
- `MDR_apc_ifn`: antigen presentation and interferon activation.
- `MDR_axis`: `MDR_injury - MDR_resolution`.

## Quick Start

From R:

```r
shiny::runApp("F:/project/LYZY/single_pan_cancer/MDRi/inst/shiny")
```

Or on the server:

```r
shiny::runApp("/hwdata/home/longzq/project/Pancan_myeloid/MDRi/inst/shiny")
```

## Input Expression Format

Rows are genes and columns are samples.

```text
gene    Sample1    Sample2
NLRP3   5.2        3.1
INHBA   2.4        1.9
FOLR2   6.2        7.0
APOE    8.1        7.8
```

## Optional Clinical Format

```text
sample    OS.time    OS    PFI.time    PFI
Sample1   1200       1     800         1
Sample2   2300       0     1700        0
```

The app can merge clinical data either by `sample` or by TCGA-style `patient`
barcode using the first 12 characters of the sample ID.

## Included Reference

The MVP includes a TCGA pan-cancer MDRi reference generated from UCSC Xena Toil
expression and PanCanAtlas survival tables:

- TCGA MDRi scores and survival metadata.
- Pan-cancer Cox results.
- Median-split KM summary results.
- Example bulk expression and clinical files.

## Benchmark Module

The repository includes pan-cancer benchmark analyses comparing MDRi with
conventional immune infiltration estimates from TIMER, MCPcounter, xCell,
CIBERSORT, quanTIseq and EPIC.

Included benchmark files:

- `inst/extdata/benchmark/MDRi_infiltration_spearman_correlations.csv`
- `inst/extdata/benchmark/MDRi_infiltration_panTCGA_univariate_Cox.csv`
- `inst/extdata/benchmark/MDRi_multivariable_adjusted_for_benchmark_infiltration.csv`
- `inst/extdata/benchmark/benchmark_summary.txt`

The reproducible benchmark workflow is provided in:

- `inst/scripts/benchmark_tcga_infiltration.R`

Run it with:

```r
Rscript inst/scripts/benchmark_tcga_infiltration.R infiltration_estimation_for_tcga.csv.gz
```

The large merged TCGA infiltration table is intentionally not bundled to keep
the repository lightweight.

## Main Functions

The core functions live in `R/mdri_core.R`:

- `mdri_read_expression()`
- `mdri_score()`
- `mdri_gene_coverage()`
- `mdri_read_clinical()`
- `mdri_merge_clinical()`
- `mdri_cox()`
- `mdri_km_data()`

## Current Status

This is an MVP for manuscript packaging and interactive exploration. It includes
MDRi scoring, TCGA reference visualization, clinical survival analysis, and a
benchmark module against conventional immune infiltration estimates. A fuller
methods-paper version can further add external cohort benchmarking, time-dependent
AUC/C-index comparisons, and published TAM signature collections.
