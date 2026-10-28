# CSV input workflow

## Data Inputs

Three CSV files: ASV-by-sample integer counts, metadata and taxonomy. Defaults
read the de-identified oral files bundled with the installed package. Edit `input_dir`
to use your own files after matching the schemas in the package README.

## Methodology

Single-tissue analyses use complete two-group household pairs, species-specific
double-carrier pools and one-to-one non-household permutations. Input counts
must already have undergone the intended upstream filtering.

## Key Outputs

Complete statistics, audit tables, resampling distributions, PDF/PNG figures,
parameters and session information are written to `output_dir`.

## Modification Guide

Edit the parameter section at the top of the R script. Counts require globally
consistent exact-sequence ASV IDs. Metadata needs `sample_id`, `family_id`,
`group`. Taxonomy needs `asv_id` and full `species` names.
The bundled data use exclude_unclassified = FALSE because species-rank filtering
was applied before export. Review this setting for your own annotations.
For subset plots, use `plot_famasv_species()` or `export_famasv(species = ...)`
after analysis, without rerunning or recalculating FDR.

## Statistical Interpretation

Read-count detection and the required number of shared ASVs are separate
thresholds. The analysis tests excess similarity conditional on both members
carrying a species; it does not demonstrate transmission.

## Results-to-Manuscript Mapping

Report carrier counts, observed and expected sharing, excess and original FDR
from `species_summary.csv`, with the recorded parameters and package versions.
