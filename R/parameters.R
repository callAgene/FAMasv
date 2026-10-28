#' Analysis parameters
#'
#' The defaults reproduce the primary species-sharing definitions in the
#' originating research notebooks. Detection is a count threshold, not a
#' relative-abundance threshold. Presence does not imply strain identity.
#' @param detection_count Minimum reads per sample for an ASV to be present.
#' @param minimum_codetected_families Minimum households in which both members
#'   carry the species; must be at least three. The default is five.
#' @param n_permutations Number of sampled derangements per carrier count.
#' @param n_bootstrap Household resamples for visualization; zero disables them.
#' @param random_seed Seed for permutation generation; carrier count is added.
#' @param bootstrap_seed Seed for household bootstrap distributions.
#' @param minimum_shared_asvs Minimum distinct identical ASVs defining sharing.
#'   The default is one. Changing this does not change the Jaccard definition.
#' @param fdr_cutoff Threshold used to label significant results.
#' @param exclude_unclassified Exclude blank, unclassified, uncultured and norank
#'   species labels. Blank species labels are always excluded.
#' @return A named list of validated parameters.
#' @export
#' @examples
#' parameters <- famasv_parameters(n_permutations = 99, n_bootstrap = 99)
famasv_parameters <- function(detection_count = 1L,
    minimum_codetected_families = 5L, n_permutations = 10000L,
    n_bootstrap = 10000L, random_seed = 20260825L,
    bootstrap_seed = 20260908L, minimum_shared_asvs = 1L,
    fdr_cutoff = 0.05, exclude_unclassified = TRUE) {
  values <- list(detection_count = detection_count,
    minimum_codetected_families = minimum_codetected_families,
    n_permutations = n_permutations, n_bootstrap = n_bootstrap,
    random_seed = random_seed, bootstrap_seed = bootstrap_seed,
    minimum_shared_asvs = minimum_shared_asvs)
  lower <- c(1, 3, 2, 0, 0, 0, 1)
  for (i in seq_along(values)) {
    x <- values[[i]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x != floor(x) || x < lower[i] || x > .Machine$integer.max - 100000L) {
      stop(names(values)[i], " must be a finite integer >= ", lower[i], ".")
    }
  }
  if (!is.numeric(fdr_cutoff) || length(fdr_cutoff) != 1L ||
      !is.finite(fdr_cutoff) || fdr_cutoff <= 0 || fdr_cutoff > 1) {
    stop("fdr_cutoff must lie in (0, 1].")
  }
  if (!is.logical(exclude_unclassified) || length(exclude_unclassified) != 1L ||
      is.na(exclude_unclassified)) stop("exclude_unclassified must be TRUE or FALSE.")
  c(values, list(fdr_cutoff = fdr_cutoff,
    exclude_unclassified = exclude_unclassified))
}

# Keep reproducible analyses from changing the caller's random-number stream.
famasv_restore_seed <- function(had_seed, old_seed) {
  if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    rm(".Random.seed", envir = .GlobalEnv)
}

famasv_derangements <- function(n, n_permutations, seed) {
  set.seed(seed)
  target <- seq_len(n)
  permutations <- matrix(NA_integer_, n_permutations, n)
  accepted <- 0L
  while (accepted < n_permutations) {
    candidate <- sample.int(n)
    if (all(candidate != target)) {
      accepted <- accepted + 1L
      permutations[accepted, ] <- candidate
    }
  }
  permutations
}
