#' Analyze species sharing within paired households
#'
#' For each species, the null uses only households in which BOTH
#' original members carry the species. Each shuffled assignment uses every
#' group-2 participant once and excludes every original household match.
#' Expected sharing is the null MEAN. One-sided P values use the plus-one
#' correction. BH correction is separate for each metric, across
#' all species meeting the carrier-household threshold. Bootstrap distributions
#' describe sampling uncertainty of observed households; they do not determine
#' P values. The result is conditional household similarity, not transmission.
#'
#' @param counts Numeric matrix of integer counts, ASVs in rows, samples in columns.
#' @param metadata Data frame with sample_id, family_id and group for one tissue.
#' @param taxonomy Data frame with asv_id and full species labels in species.
#'   Identical species labels are treated as the same species. ASV IDs must
#'   identify the same exact sequence across all samples.
#' @param parameters A list returned by famasv_parameters.
#' @param group_labels Two group values in metadata, in focal/partner order.
#' @return A plain list containing summary, family_data, distributions,
#'   pair_audit, input_audit, permutation_audit, parameters and group_labels.
#'   Rates and Jaccard values are fractions from zero to one. Excess values are
#'   observed minus null mean. Untested species have missing P values and FDR.
#' @export
#' @examples
#' path <- system.file("extdata", package = "FAMasv")
#' input <- read_famasv(file.path(path, "counts.csv"),
#'   file.path(path, "metadata.csv"), file.path(path, "taxonomy.csv"))
#' result <- run_famasv(input$counts, input$metadata, input$taxonomy,
#'   famasv_parameters(n_permutations = 99, n_bootstrap = 0,
#'     exclude_unclassified = FALSE))
#' result$summary[, c("Species", "Sharing_FDR")]
run_famasv <- function(counts, metadata, taxonomy,
    parameters = famasv_parameters(), group_labels = c("PD", "HC")) {
  parameters <- do.call(famasv_parameters, parameters)
  if (!is.character(group_labels) || length(group_labels) != 2L ||
      anyNA(group_labels) || any(!nzchar(group_labels)) || anyDuplicated(group_labels))
    stop("group_labels must contain two distinct non-empty group names.")
  input <- famasv_validate_input(counts, metadata, taxonomy)
  counts <- input$counts
  metadata <- input$metadata
  taxonomy <- input$taxonomy
  valid <- !is.na(taxonomy$species) & nzchar(taxonomy$species)
  if (parameters$exclude_unclassified) valid <- valid & !grepl(
    "unclassified|uncultured|norank", taxonomy$species, ignore.case = TRUE)
  if (!any(valid)) stop("No ASVs with eligible species labels remain.")
  input_audit <- data.frame(Item = c("Input samples", "Input ASVs",
    "Retained ASVs", "Excluded ASVs", "Samples outside selected groups"),
    N = c(ncol(counts), nrow(counts), sum(valid), sum(!valid),
      sum(!metadata$group %in% group_labels)))
  counts <- counts[valid, , drop = FALSE]
  taxonomy <- taxonomy[valid, , drop = FALSE]
  species_names <- sort(unique(taxonomy$species))
  if (!any(metadata$group %in% group_labels)) stop("No samples have the requested group labels.")

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(famasv_restore_seed(had_seed, old_seed), add = TRUE)
  permutation_cache <- list()
  audit_rows <- list()
  result_rows <- list()
  family_rows <- list()
  distributions <- list()
  result_index <- 0L

  current_metadata <- metadata[metadata$group %in% group_labels, , drop = FALSE]
  family_names <- sort(unique(current_metadata$family_id))
  family_count <- table(factor(current_metadata$family_id, levels = family_names),
    factor(current_metadata$group, levels = group_labels))
  if (any(family_count > 1L)) stop("Multiple samples from the same group in a household. ",
    "Use one tissue and resolve technical replicates before analysis.")
  complete <- rowSums(family_count == 1L) == 2L
  pair_audit <- data.frame(Family_ID = family_names,
    N_group1 = as.integer(family_count[, 1]), N_group2 = as.integer(family_count[, 2]),
    Included = complete, Reason = ifelse(complete, "Complete pair", "Missing partner"))
  family_names <- family_names[complete]
  if (!length(family_names)) stop("No complete households remain.")
  group1 <- current_metadata[current_metadata$group == group_labels[1], , drop = FALSE]
  group2 <- current_metadata[current_metadata$group == group_labels[2], , drop = FALSE]
  sample1 <- group1$sample_id[match(family_names, group1$family_id)]
  sample2 <- group2$sample_id[match(family_names, group2$family_id)]
  presence1 <- t(counts[, sample1, drop = FALSE]) >= parameters$detection_count
  presence2 <- t(counts[, sample2, drop = FALSE]) >= parameters$detection_count

  for (species_index in seq_along(species_names)) {
    species <- species_names[species_index]
    asv_index <- which(taxonomy$species == species)
    group1_current <- presence1[, asv_index, drop = FALSE]
    group2_current <- presence2[, asv_index, drop = FALSE]
    carrier1 <- rowSums(group1_current) > 0
    carrier2 <- rowSums(group2_current) > 0
    eligible <- which(carrier1 & carrier2)
    n <- length(eligible)
    formal <- n >= parameters$minimum_codetected_families
    observed <- expected <- excess <- p_value <- z_value <- c(NA_real_, NA_real_)
    null <- bootstrap <- data.frame(Sharing = numeric(), Jaccard = numeric())
    if (n > 0L) {
      a <- group1_current[eligible, , drop = FALSE]
      b <- group2_current[eligible, , drop = FALSE]
      intersection <- a %*% t(b)
      union <- outer(rowSums(a), rowSums(b), "+") - intersection
      binary <- intersection >= parameters$minimum_shared_asvs
      jaccard <- intersection / union
      observed_binary <- diag(binary)
      observed_jaccard <- diag(jaccard)
      observed <- c(mean(observed_binary), mean(observed_jaccard))
      family_rows[[length(family_rows) + 1L]] <- data.frame(Species = species, Family_ID = family_names[eligible],
        Group1_sample = sample1[eligible], Group2_sample = sample2[eligible],
        Group1_ASVs = rowSums(a), Group2_ASVs = rowSums(b),
        N_shared_ASVs = diag(intersection), Shared_exact_ASV = as.integer(observed_binary),
        Within_species_ASV_Jaccard = observed_jaccard)
      if (formal) {
        key <- as.character(n)
        if (is.null(permutation_cache[[key]])) {
          permutations <- famasv_derangements(n, parameters$n_permutations,
            parameters$random_seed + n)
          permutation_cache[[key]] <- permutations
          audit_rows[[key]] <- data.frame(Carrier_N = n,
            N_permutations = nrow(permutations),
            Unique_pairings = nrow(unique(permutations)),
            No_true_household_pairing = all(sweep(permutations, 2, seq_len(n), "!=") ),
            Every_partner_used_once = all(apply(permutations, 1,
              function(x) identical(sort(x), seq_len(n)))))
        }
        permutations <- permutation_cache[[key]]
        row_index <- matrix(seq_len(n), nrow = nrow(permutations), ncol = n, byrow = TRUE)
        null_index <- cbind(as.vector(row_index), as.vector(permutations))
        null <- data.frame(
          Sharing = rowMeans(matrix(binary[null_index], nrow(permutations), n)),
          Jaccard = rowMeans(matrix(jaccard[null_index], nrow(permutations), n)))
        expected <- colMeans(null)
        excess <- observed - expected
        null_sd <- vapply(null, stats::sd, numeric(1))
        z_value[null_sd > 1e-12] <- excess[null_sd > 1e-12] / null_sd[null_sd > 1e-12]
        p_value <- c((1 + sum(null$Sharing >= observed[1])) / (nrow(null) + 1),
          (1 + sum(null$Jaccard >= observed[2])) / (nrow(null) + 1))
        if (parameters$n_bootstrap > 0L) {
          set.seed(parameters$bootstrap_seed + species_index)
          bootstrap_index <- matrix(sample.int(n, parameters$n_bootstrap * n,
            replace = TRUE), nrow = parameters$n_bootstrap, ncol = n)
          bootstrap <- data.frame(
            Sharing = rowMeans(matrix(observed_binary[bootstrap_index], parameters$n_bootstrap, n)),
            Jaccard = rowMeans(matrix(observed_jaccard[bootstrap_index], parameters$n_bootstrap, n)))
        }
      }
    }
    result_index <- result_index + 1L
    distributions[[result_index]] <- list(Species = species,
      null = null, bootstrap = bootstrap)
    result_rows[[result_index]] <- data.frame(Species = species,
      N_total_ASVs = length(asv_index), N_complete_families = length(family_names),
      N_group1_carriers = sum(carrier1), N_group2_carriers = sum(carrier2),
      Group1_prevalence = mean(carrier1), Group2_prevalence = mean(carrier2),
      N_codetected_real_families = n,
      Observed_exact_ASV_sharing_rate = observed[1],
      Expected_exact_ASV_sharing_rate = unname(expected[1]),
      Excess_exact_ASV_sharing = unname(excess[1]), Sharing_Z = z_value[1],
      Sharing_permutation_P = p_value[1],
      Observed_mean_ASV_Jaccard = observed[2], Expected_mean_ASV_Jaccard = unname(expected[2]),
      Excess_ASV_Jaccard = unname(excess[2]), Jaccard_Z = z_value[2],
      Jaccard_permutation_P = p_value[2], Formal_analysis = formal)
  }
  summary <- do.call(rbind, result_rows)
  summary$Sharing_FDR <- summary$Jaccard_FDR <- NA_real_
  eligible <- summary$Formal_analysis
  summary$Sharing_FDR[eligible] <- stats::p.adjust(summary$Sharing_permutation_P[eligible], "BH")
  summary$Jaccard_FDR[eligible] <- stats::p.adjust(summary$Jaccard_permutation_P[eligible], "BH")
  summary$Sharing_significant <- !is.na(summary$Sharing_FDR) & summary$Sharing_FDR < parameters$fdr_cutoff
  summary$Jaccard_significant <- !is.na(summary$Jaccard_FDR) & summary$Jaccard_FDR < parameters$fdr_cutoff
  empty_family <- data.frame(Species = character(), Family_ID = character(),
    Group1_sample = character(), Group2_sample = character(), Group1_ASVs = integer(),
    Group2_ASVs = integer(), N_shared_ASVs = integer(), Shared_exact_ASV = integer(),
    Within_species_ASV_Jaccard = numeric())
  empty_audit <- data.frame(Carrier_N = integer(), N_permutations = integer(),
    Unique_pairings = integer(), No_true_household_pairing = logical(), Every_partner_used_once = logical())
  list(summary = summary,
    family_data = if (length(family_rows)) do.call(rbind, family_rows) else empty_family,
    distributions = distributions, pair_audit = pair_audit,
    permutation_audit = if (length(audit_rows)) do.call(rbind, audit_rows) else empty_audit,
    input_audit = input_audit, parameters = parameters, group_labels = group_labels,
    package_version = "0.3.0")
}
