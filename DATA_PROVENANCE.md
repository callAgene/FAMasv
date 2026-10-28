# Example-data provenance

## Dataset

The sole bundled example is de-identified observational oral microbiome data
from a PD-HC household cohort supplied by the maintainer. It comprises 50
complete households, 100 samples, 10,111 observed ASVs and 508 annotated taxa.
It is not simulated data and has no synthetic ground-truth labels.

## Preparation

1. Select oral PD and HC samples and retain households with exactly one member
   of each group.
2. Apply the source analysis's species-rank annotation filter: exclude blank,
   unclassified, uncultured and norank species-rank assignments.
3. Retain ASVs with at least one nonzero count in the selected samples.
4. Independently replace sample, household and ASV IDs with randomly assigned
   example identifiers and randomize table order. Preserve paired membership,
   group labels and every retained count.
5. Export only counts.csv, metadata.csv and taxonomy.csv. No correspondence
   table, original phyloseq object or sequence file is included.

Some full taxon labels preserve unresolved upper-rank prefixes while their
species-rank annotation is retained by the source analysis. Therefore the
bundled workflow uses `exclude_unclassified = FALSE` to avoid filtering those
full names again. This yields 86 formally tested taxa with default numerical
settings. Presence thresholds 1 and 2 are equivalent for these files because
the smallest nonzero count is 2.

## Contents and omissions

Included: replacement IDs, PD/HC labels, household membership, taxon names and
ASV read counts. Excluded: original identifiers, nucleotide sequences, dates,
demographic and other clinical variables, and identifier correspondence tables.
De-identification does not imply complete anonymity: disease labels and
microbial profiles remain part of the example.

## Reproducibility

The distributed CSVs are fixed inputs. Their SHA-256 checksums are recorded in
`inst/extdata/SHA256SUMS`. Statistical random seeds are specified by the package
parameters. Independently changing household numbering/order can change finite
Monte Carlo draws even when observed statistics and the null-model definition
remain identical.

The example preparation was checked against the source analysis: eligible
household counts, observed sharing and observed Jaccard were unchanged. The
formal-test set and sharing-significance classification were also preserved
in the validated 10,000-permutation run. This is a software validation result,
not an additional clinical finding.

## Licence and contact

Code and bundled example data are distributed under the MIT licence. See
`LICENSE.md` and `inst/extdata/LICENSE`. Maintainer:
<bingzhang@ucas.ac.cn>.
