# Small controlled fixtures for statistical unit tests only.
make_test_counts <- function(n_families = 20L, seed = 101L) {
  if (length(n_families) != 1L || !is.numeric(n_families) || !is.finite(n_families) ||
      n_families < 6 || n_families != floor(n_families)) stop("n_families must be an integer >= 6.")
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed > .Machine$integer.max || seed != floor(seed)) stop("seed must be a non-negative integer.")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  species <- c("Simulibacter familiaris", "Simulibacter moderatus",
    "Simulibacter randomus", "Simulibacter communis", "Simulibacter separatus",
    "Simulibacter unilateralis", "Simulibacter rarus", "unclassified bacterium")
  scenario <- c("family_specific", "moderate_sharing", "independent_draws",
    "ubiquitous_identical", "group_specific_ASVs", "one_sided_carriers",
    "fewer_than_five_co_carriers", "excluded_annotation")
  taxonomy <- data.frame(asv_id = sprintf("ASV%04d", seq_len(8L * n_families)),
    species = rep(species, each = n_families))
  metadata <- expand.grid(group = c("PD", "HC"),
    family_id = sprintf("F%03d", seq_len(n_families)),
    stringsAsFactors = FALSE)
  metadata$sample_id <- sprintf("Sample%04d", seq_len(nrow(metadata)))
  metadata <- metadata[, c("sample_id", "family_id", "group")]
  counts <- matrix(0L, nrow(taxonomy), nrow(metadata),
    dimnames = list(taxonomy$asv_id, metadata$sample_id))
  for (f in seq_len(n_families)) {
    columns <- which(metadata$family_id == sprintf("F%03d", f))
    for (s in seq_along(species)) {
      chosen <- switch(s,
        list(f, f),
        list(f, if (stats::runif(1) < 0.65) f else sample.int(n_families, 1)),
        list(sample.int(n_families, 2), sample.int(n_families, 2)),
        list(1L, 1L),
        list(1L, 2L),
        if (f %% 2L == 0) list(f, integer()) else list(integer(), f),
        if (f <= 3L) list(f, f) else list(integer(), integer()),
        list(1L, 1L))
      for (g in 1:2) {
        rows <- (s - 1L) * n_families + chosen[[g]]
        if (length(rows)) counts[rows, columns[g]] <- sample(2:100, length(rows), replace = TRUE)
      }
    }
  }
  list(counts = counts, metadata = metadata, taxonomy = taxonomy,
    truth = data.frame(species = species, scenario = scenario))
}
