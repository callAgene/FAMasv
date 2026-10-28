# README figure generation

## Data Inputs

Saved result from the bundled real oral example, after running
inst/examples/run_oral_example.R from the repository root.

## Methodology

Select 20 taxa for the overview and three for distribution plots, ranked by
excess sharing. Keep the full-analysis FDR and stored resampling distributions.

## Key Outputs

Two PNG files in man/figures, linked by both READMEs.

## Modification Guide

Adjust result_file, output_dir and display counts at the top of the R script.
Run after installing the current package. Existing PNGs are replaced.

## Statistical Interpretation

Selection is for display only. Data are real observations, not simulation truth.

## Results-to-Manuscript Mapping

Use the complete result table for statistics and its parameter/session files
for reporting the run that produced the illustrations.
