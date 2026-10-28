# FAMasv

**Species-resolved analysis of household microbial sharing based on full-length 16S rRNA gene sequencing.**

FAMasv compares ASV sharing within household pairs with a species-specific
non-household background. It provides sharing rates, within-species Jaccard
similarity, permutation tests, multiple-testing correction and batch figure export.
Analyze **one tissue or sample type at a time** using three input tables.

[中文说明](README.zh-CN.md) · [Complete example](inst/examples/run_oral_example.R) · [Input data](inst/extdata/README.md) · [Data provenance](DATA_PROVENANCE.md)

![Household excess sharing in the de-identified oral example](man/figures/README-oral-landscape.png)

*Top 20 species ranked by excess sharing in the bundled oral dataset. Colour
indicates sharing FDR; point size indicates the number of double-carrier households.*

## Installation

Requires **R >= 4.1.0** and **ggplot2 >= 3.4.0**. The optional `phyloseq` package
is needed only for the phyloseq input adapter.

### Install from the GitHub source folder

On the repository page, choose **Code → Download ZIP**, extract the archive,
and open R or RStudio with the extracted package folder as the working directory.
This is the folder containing `DESCRIPTION`, `NAMESPACE` and `R/`.

```r
install.packages("remotes") # Once, if not already installed.
remotes::install_github("callAgene/FAMasv", dependencies = NA, upgrade = "never")
library(FAMasv)
```

### Install from an R source package

If you have separately downloaded `FAMasv_0.3.0.tar.gz`, run the following
in the directory containing that file:

```r
install.packages("ggplot2") # Skip if already installed.
install.packages("FAMasv_0.3.0.tar.gz", repos = NULL, type = "source")
library(FAMasv)
```

The GitHub **Download ZIP** archive and the R `.tar.gz` source package are
different formats; use the matching installation method above.

## Quick start

The package includes one **de-identified real oral dataset**: 50 PD-HC
households, 100 samples, 10,111 ASVs and 508 annotated taxa. These are
observational data, not simulated data or simulation ground truth.

```r
library(FAMasv)

input_dir <- system.file("extdata", package = "FAMasv")
input <- read_famasv(
  file.path(input_dir, "counts.csv"),
  file.path(input_dir, "metadata.csv"),
  file.path(input_dir, "taxonomy.csv")
)

# This example was already filtered using the original species-rank annotation.
parameters <- famasv_parameters(
  exclude_unclassified = FALSE,
  n_permutations = 10000,
  n_bootstrap = 10000
)
result <- run_famasv(input$counts, input$metadata, input$taxonomy, parameters)

# Complete tables and a small set of distribution plots for an initial preview.
selected <- select_famasv_species(result, significant_only = FALSE, top_n = 6)
export_famasv(result, "results_preview", species = selected$Species)
```

For the bundled example, use **`exclude_unclassified = FALSE`**. Species-rank
filtering has already been performed, and some full names retain an unresolved
upper-rank prefix such as `norank`. Filtering those names a second time would
remove taxa retained in the original analysis. For your own unfiltered taxonomy,
choose this setting according to its annotation scheme.

Use a fresh output directory for each export. `overwrite = TRUE` replaces
matching files but does not remove old files. The complete runnable workflow is
[run_oral_example.R](inst/examples/run_oral_example.R).

## Data Inputs

### Required files and dimensions

Dimensions below exclude the header row.

| File | Each row | Columns | Bundled example |
|---|---|---|---|
| [counts.csv](inst/extdata/counts.csv) | One ASV | `asv_id`, then one column per sample | 10,111 × 101 |
| [metadata.csv](inst/extdata/metadata.csv) | One sample | `sample_id`, `family_id`, `group` | 100 × 3 |
| [taxonomy.csv](inst/extdata/taxonomy.csv) | One ASV assignment | `asv_id`, `species` | 10,111 × 2 |

### 1. Counts: ASVs in rows, samples in columns

Selected rows and sample columns from the bundled count table:

| asv_id | ORAL_S0085 | ORAL_S0088 |
|---|---:|---:|
| ORAL_A00001 | 0 | 45 |
| ORAL_A00002 | 35 | 0 |
| ORAL_A00003 | 25 | 0 |

Use **non-negative integer reads**, not relative abundances or CLR values.
Zero indicates no reads. ASV IDs must identify the same exact sequence across
all samples. Upstream ASV inference, sequence filtering and taxonomic
classification are performed before FAMasv.

`read_famasv()` converts the first column to row names, giving a **10,111 × 100**
numeric matrix. Sample column names are retained unchanged.

### 2. Metadata: one sample per row

Two complete households from the same example:

| sample_id | family_id | group |
|---|---|---|
| ORAL_S0077 | ORAL_F050 | HC |
| ORAL_S0085 | ORAL_F050 | PD |
| ORAL_S0088 | ORAL_F022 | HC |
| ORAL_S0097 | ORAL_F022 | PD |

`sample_id` must match count-table columns exactly. Shared `family_id` values
identify the original household pair. Each household must have exactly one
member of each group. Incomplete pairs are recorded and excluded; duplicate
same-group samples cause an error. Extra metadata columns are allowed but are
not used as covariates. The default groups are PD and HC; change them with
`group_labels = c("case_label", "control_label")` in `run_famasv()`.

### 3. Taxonomy: one ASV assignment per row

| asv_id | species |
|---|---|
| ORAL_A00001 | Prevotella melaninogenica |
| ORAL_A00002 | Selenomonas sp. oral taxon 478 |
| ORAL_A00003 | Actinomyces graevenitzii |

Use complete taxon names. Several ASVs can map to one species; identical
labels are grouped together. Species-level labels can include named oral taxa
or unresolved upper-rank prefixes, as in the source annotation.

Counts and metadata must contain exactly the same sample IDs. Counts and
taxonomy must contain exactly the same ASV IDs. Row order can differ: matching
uses IDs rather than positions. Do not add spreadsheet row-number columns.

### Use your own data

```r
input <- read_famasv("data/counts.csv", "data/metadata.csv", "data/taxonomy.csv")
parameters <- famasv_parameters(exclude_unclassified = TRUE)
result <- run_famasv(input$counts, input$metadata, input$taxonomy,
                     parameters, group_labels = c("PD", "HC"))
```

For TSV inputs add `sep = "\t"`. For an existing single-tissue phyloseq object:

```r
input <- famasv_from_phyloseq(ps,
  family_col = "familyID", group_col = "CaseAConB",
  species_col = "species", genus_col = "genus")
result <- do.call(run_famasv, c(input, list(parameters = parameters)))
```

Use `genus_col = NULL` if species names are already complete and no genus
column is available. The adapter handles OTU-table orientation, taxonomy
prefixes and underscores. Inspect converted labels before analysis.

## Modification Guide

### Analysis settings

| Parameter | Default | Interpretation |
|---|---:|---|
| `detection_count` | 1 | Minimum reads to call an ASV present |
| `minimum_shared_asvs` | 1 | Minimum distinct identical ASVs to call a household sharing-positive |
| `minimum_codetected_families` | 5 | Minimum households with two carriers for formal testing; at least 3 |
| `n_permutations` | 10000 | Random non-household matching draws |
| `n_bootstrap` | 10000 | Household resamples for visualization; 0 disables |
| `random_seed` | 20260825 | Base permutation seed |
| `bootstrap_seed` | 20260908 | Base bootstrap seed |
| `fdr_cutoff` | 0.05 | Strict significance threshold: FDR < cutoff |
| `exclude_unclassified` | TRUE | Exclude labels containing unclassified, uncultured or norank; use FALSE for the bundled, already filtered example |

Blank taxonomy labels are always excluded. Two reads and two shared ASVs are
separate settings: `detection_count = 2` requires two reads for presence;
`minimum_shared_asvs = 2` requires two distinct identical ASVs for sharing.
The bundled input has no positive count below 2, so presence thresholds 1 and
2 give identical presence calls on this dataset.

### Batch plots and species selection

```r
# Every formally tested species: overviews, combined pages and individual plots.
export_famasv(result, "results_all", formats = c("pdf", "png"),
              species_per_page = 4, per_species = TRUE)

# Sharing-significant species only.
export_famasv(result, "results_significant", significant_only = TRUE)

# Filter by significance and effect size, retaining the original all-species FDR.
selected <- select_famasv_species(result, significant_only = TRUE,
                                  min_excess = 0.10, top_n = 10)
if (nrow(selected) > 0) {
  export_famasv(result, "results_filtered", species = selected$Species)
}

# Manually choose taxa and return a ggplot object.
p <- plot_famasv_species(result,
  species = c("Lachnospira pectinoschiza", "Prevotella jejuni"))
ggplot2::ggsave("selected_species.pdf", p, width = 10, height = 6)
```

![Observed-household and permuted-pair distributions](man/figures/README-oral-distributions.png)

*Three taxa with the largest excess sharing in the oral example. Grey:
non-household null; blue: household bootstrap. Solid black lines: observed
statistics; dashed lines: null means. FDR values are from the complete analysis.*

`min_excess = 0.10` means 10 percentage points. Use `metric = "jaccard"` for
Jaccard-based selection. To omit individual plots use `per_species = FALSE`.
PDF, PNG, TIFF and SVG are supported. Overview plots show the top 30 by default;
`plot_famasv_landscape(result, top_n = NULL)` shows all tested taxa.

The batch exporter covers all formally tested taxa unless a selection is
specified. Explicit names may include untested taxa, which are labelled as
such. Selecting taxa never recalculates P values or FDR. Reload the saved RDS
to change figures without rerunning permutations.

## Key Outputs

```text
results_preview/
├── tables/
│   ├── species_summary.csv
│   ├── selected_species.csv
│   ├── family_observed.csv
│   ├── pair_audit.csv
│   ├── permutation_audit.csv
│   └── input_audit.csv
├── figures/
│   ├── landscape_sharing.pdf / .png
│   ├── landscape_jaccard.pdf / .png
│   ├── distributions_page_001.pdf / .png
│   └── species_001_<species_name>.pdf / .png
├── famasv_result.rds
├── parameters.R
├── sessionInfo.txt
├── manifest.csv
└── README.txt
```

Full summary tables and overview plots are retained even for selected-species
exports. An empty significant selection still produces tables and overviews.

## Methodology

1. Identify complete two-group household pairs and call ASV presence using the
   read threshold.
2. For each species, retain **only households where both original members carry
   it**. If 7 of 46 households qualify, the null uses those 7 PD and 7 HC
   participants; carriers from the other households are excluded.
3. Calculate the fraction of eligible households sharing the required number
   of identical ASVs and mean ASV Jaccard similarity (intersection / union).
4. Reassign HC partners one-to-one, using each HC once per iteration and
   excluding every original household match. Matching schemes may repeat
   between iterations. Uniform derangements are sampled by rejection.
5. Define excess as observed minus the **mean** null statistic. The one-sided
   P value is `(1 + sum(null >= observed)) / (B + 1)`.
6. Apply BH correction across all formally tested species, separately for
   sharing and Jaccard. Bootstrap the original eligible households to display
   uncertainty of the observed statistics.

Permutation schedules are cached by carrier count with seed
`random_seed + carrier_count`. Bootstrap streams use sorted species index.
The caller's random-number state is restored. Null and bootstrap draws are
retained for reproducible plotting; memory and runtime increase with dataset
size and iteration count.

