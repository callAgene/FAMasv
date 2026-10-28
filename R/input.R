#' Read the three input tables
#' @param counts_file CSV with asv_id as its first column, followed by samples.
#' @param metadata_file CSV with sample_id, family_id and group (one tissue).
#' @param taxonomy_file CSV with asv_id and species (full species labels).
#' @param sep Field separator, default comma; use '\\t' for TSV.
#' @return A list with counts (ASVs in rows), metadata and taxonomy.
#' @export
#' @examples
#' path <- system.file("extdata", package = "FAMasv")
#' input <- read_famasv(file.path(path, "counts.csv"),
#'   file.path(path, "metadata.csv"), file.path(path, "taxonomy.csv"))
read_famasv <- function(counts_file, metadata_file, taxonomy_file, sep = ",") {
  counts <- utils::read.table(counts_file, header = TRUE, sep = sep,
    check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
  if (names(counts)[1] != "asv_id") stop("The first counts column must be asv_id.")
  if (anyDuplicated(counts$asv_id)) stop("Duplicate ASV IDs in counts.")
  ids <- counts$asv_id
  counts <- as.matrix(counts[, -1, drop = FALSE])
  rownames(counts) <- ids
  list(counts = counts,
    metadata = utils::read.table(metadata_file, header = TRUE, sep = sep,
      check.names = FALSE, stringsAsFactors = FALSE, comment.char = ""),
    taxonomy = utils::read.table(taxonomy_file, header = TRUE, sep = sep,
      check.names = FALSE, stringsAsFactors = FALSE, comment.char = ""))
}

#' Convert an existing phyloseq object to the FAMasv input format
#' @param object A phyloseq object containing counts, taxonomy and metadata.
#' @param family_col,group_col Metadata column names for household and group.
#' @param species_col,genus_col Taxonomy column names. genus_col may be NULL if
#'   species_col already contains full names.
#' @return A list accepted by run_famasv via do.call.
#' @export
famasv_from_phyloseq <- function(object, family_col = "familyID",
    group_col = "CaseAConB",
    species_col = "species", genus_col = "genus") {
  if (!requireNamespace("phyloseq", quietly = TRUE)) stop("Install phyloseq to use this adapter.")
  counts <- methods::as(phyloseq::otu_table(object), "matrix")
  if (!phyloseq::taxa_are_rows(object)) counts <- t(counts)
  metadata <- data.frame(phyloseq::sample_data(object), check.names = FALSE)
  taxonomy <- data.frame(phyloseq::tax_table(object), check.names = FALSE)
  if (!all(c(family_col, group_col) %in% names(metadata)))
    stop("Requested metadata columns are missing.")
  if (!all(c(species_col, genus_col) %in% names(taxonomy)))
    stop("Requested taxonomy columns are missing.")
  species <- gsub("_", " ", sub("^[a-z]__", "", as.character(taxonomy[[species_col]])))
  if (!is.null(genus_col)) {
    genus <- gsub("_", " ", sub("^[a-z]__", "", as.character(taxonomy[[genus_col]])))
    valid <- !is.na(species) & nzchar(trimws(species)) &
      !is.na(genus) & nzchar(trimws(genus))
    add_genus <- valid & !startsWith(species, paste0(genus, " "))
    species[add_genus] <- paste(genus[add_genus], species[add_genus])
  }
  list(counts = counts, metadata = data.frame(sample_id = rownames(metadata),
    family_id = metadata[[family_col]], group = metadata[[group_col]]),
    taxonomy = data.frame(asv_id = rownames(taxonomy), species = species))
}

famasv_validate_input <- function(counts, metadata, taxonomy) {
  counts <- as.matrix(counts)
  if (!is.numeric(counts) || length(counts) == 0 || any(!is.finite(counts)) ||
      any(counts < 0) || any(abs(counts - round(counts)) > 1e-8))
    stop("counts must contain finite, non-negative integer read counts.")
  for (ids in list(rownames(counts), colnames(counts))) {
    if (is.null(ids) || anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids))
      stop("counts must have unique, non-empty ASV row names and sample column names.")
  }
  if (!is.data.frame(metadata) ||
      !all(c("sample_id", "family_id", "group") %in% names(metadata)))
    stop("metadata requires sample_id, family_id and group columns.")
  if (!is.data.frame(taxonomy) || !all(c("asv_id", "species") %in% names(taxonomy)))
    stop("taxonomy requires asv_id and species columns.")
  for (field in c("sample_id", "family_id", "group")) {
    metadata[[field]] <- as.character(metadata[[field]])
    if (anyNA(metadata[[field]]) || any(!nzchar(trimws(metadata[[field]]))))
      stop("metadata ", field, " contains missing or empty values.")
  }
  taxonomy$asv_id <- as.character(taxonomy$asv_id)
  taxonomy$species <- trimws(as.character(taxonomy$species))
  if (anyNA(taxonomy$asv_id) || anyDuplicated(taxonomy$asv_id) ||
      anyDuplicated(metadata$sample_id)) stop("Duplicate or missing input IDs.")
  if (!setequal(colnames(counts), metadata$sample_id))
    stop("Counts and metadata must contain exactly the same sample IDs.")
  if (!setequal(rownames(counts), taxonomy$asv_id))
    stop("Counts and taxonomy must contain exactly the same ASV IDs.")
  list(counts = counts, metadata = metadata[match(colnames(counts), metadata$sample_id), ],
    taxonomy = taxonomy[match(rownames(counts), taxonomy$asv_id), ])
}
