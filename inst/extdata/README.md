# De-identified real oral example

This is the only example dataset distributed with FAMasv. It contains real
observational data, not a simulation.

| File | Rows, excluding header | Columns |
|---|---:|---:|
| counts.csv | 10,111 ASVs | 101: asv_id + 100 samples |
| metadata.csv | 100 samples | 3: sample_id, family_id, group |
| taxonomy.csv | 10,111 ASVs | 2: asv_id, species |

There are 50 complete PD-HC households and 508 annotated taxa. Retained counts
are unchanged from the source dataset; the smallest nonzero count is 2.

Use `famasv_parameters(exclude_unclassified = FALSE)` for this example.
Species-rank filtering has already been applied; full labels can retain an
unresolved upper-rank prefix despite a retained species-level assignment.

Sample, household and ASV IDs were replaced independently. Original identifiers,
sequences, additional clinical metadata and ID correspondence tables are not
included. Disease groups, paired relationships and microbial profiles remain;
de-identification is not a guarantee of complete anonymity.

`SHA256SUMS` identifies the fixed CSV version. Data are distributed under the
MIT licence; see LICENSE in this directory. Maintainer:
bingzhang@ucas.ac.cn.
