# Oral example workflow

## Data Inputs

The sole bundled example: de-identified real oral data, 50 households, 100
samples, 10,111 ASVs and 508 annotated taxa. Three CSV tables contain only
replacement IDs, group labels, household membership, counts and taxon names.

## Methodology

Run 10,000 carrier-restricted non-household permutations and 10,000 household
bootstraps. The original species-rank filtering is already applied, so set
exclude_unclassified to FALSE to retain resolved taxa with upper-rank prefixes.
All taxa are analyzed before choosing any distribution plots.

## Key Outputs

Preview and selected-species directories contain full summaries, audit tables,
RDS results, parameter/version records and PDF/PNG figures. Optional complete
batch export adds a separate all_species directory.

## Modification Guide

Copy the R script to a writable working directory. Change parameters at the
top and run it after installing FAMasv. Set export_all_species to TRUE for all
formally tested taxa. Use a fresh output_root for a new run.

## Statistical Interpretation

Selected figures retain all-species FDR. Real data have no synthetic ground
truth, and ASV sharing alone does not establish transmission.

## Results-to-Manuscript Mapping

Use species_summary.csv for observed/expected/excess sharing and FDR. Cite
parameters.R and sessionInfo.txt for the actual settings and software versions.
