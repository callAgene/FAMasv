#' Export analysis tables and batch or selected-species plots
#'
#' species = NULL exports all formally tested species by default. Explicit
#' names can include untested species. Significant-only export uses the original
#' all-species FDR. The complete summary is always exported, even for selected
#' plots. An existing nonempty output folder requires overwrite = TRUE.
#' @inheritParams select_famasv_species
#' @param output_dir Output folder.
#' @param formats Any of pdf, png, tiff, svg. PDF and PNG are the defaults.
#' @param per_species Export a separate two-metric figure for every selection.
#' @param species_per_page Maximum species panels per combined page.
#' @param dpi Raster resolution.
#' @param overwrite Allow writing into a nonempty output folder.
#' @return Invisibly, a manifest of output file paths and types.
#' @export
export_famasv <- function(result, output_dir, species = NULL,
    significant_only = FALSE, metric = c("sharing", "jaccard"),
    formats = c("pdf", "png"), per_species = TRUE, species_per_page = 4L,
    dpi = 300, overwrite = FALSE) {
  metric <- match.arg(metric)
  formats <- unique(match.arg(formats, c("pdf", "png", "tiff", "svg"), several.ok = TRUE))
  if (length(species_per_page) != 1L || !is.numeric(species_per_page) ||
      !is.finite(species_per_page) || species_per_page < 1 ||
      species_per_page != floor(species_per_page)) stop("species_per_page must be a positive integer.")
  if (length(dpi) != 1L || !is.numeric(dpi) || !is.finite(dpi) || dpi < 72) stop("dpi must be >= 72.")
  selected <- select_famasv_species(result, species, metric = metric,
    significant_only = significant_only)
  if (is.null(species)) selected <- selected[selected$Formal_analysis, , drop = FALSE]
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) && !overwrite)
    stop("output_dir is not empty. Choose another folder or set overwrite = TRUE.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  table_dir <- file.path(output_dir, "tables")
  figure_dir <- file.path(output_dir, "figures")
  dir.create(table_dir, showWarnings = FALSE)
  dir.create(figure_dir, showWarnings = FALSE)
  tables <- list(species_summary = result$summary, selected_species = selected,
    family_observed = result$family_data, pair_audit = result$pair_audit,
    permutation_audit = result$permutation_audit, input_audit = result$input_audit)
  manifest <- data.frame(file = character(), type = character())
  for (name in names(tables)) {
    path <- file.path(table_dir, paste0(name, ".csv"))
    utils::write.csv(tables[[name]], path, row.names = FALSE, na = "")
    manifest <- rbind(manifest, data.frame(file = path, type = "table"))
  }
  saveRDS(result, file.path(output_dir, "famasv_result.rds"))
  dput(result$parameters, file = file.path(output_dir, "parameters.R"))
  writeLines(c("FAMasv output", "",
    "Rates and Jaccard values in tables are fractions (0-1); plot axes use percentages.",
    "Excess = observed minus mean permutation expectation.",
    "FDR is computed over all tested species per metric before selection.",
    "Non-household null: one-to-one HC reassignment within true double-carrier households.",
    "True household matches are excluded. Assignments can repeat across permutations.",
    "Bootstrap distributions visualize uncertainty; permutation distributions determine P values.",
    "Exact-ASV similarity does not establish transmission or whole-genome strain identity.",
    paste0("Selected species: ", nrow(selected))), file.path(output_dir, "README.txt"))

  # This helper is reused for every figure and format.
  save_plot <- function(plot, stem, height) {
    paths <- character()
    for (format in formats) {
      path <- file.path(figure_dir, paste0(stem, ".", format))
      device <- switch(format, pdf = grDevices::pdf, png = grDevices::png,
        tiff = grDevices::tiff, svg = grDevices::svg)
      ggplot2::ggsave(path, plot = plot, device = device, width = 10,
        height = height, units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
      paths <- c(paths, path)
    }
    paths
  }
  for (m in c("sharing", "jaccard")) {
    paths <- save_plot(plot_famasv_landscape(result, metric = m),
      paste0("landscape_", m), 8)
    manifest <- rbind(manifest, data.frame(file = paths, type = "figure"))
  }
  if (nrow(selected)) {
    names <- selected$Species
    pages <- split(names, ceiling(seq_along(names) / species_per_page))
    for (page in seq_along(pages)) {
      p <- plot_famasv_species(result, pages[[page]])
      paths <- save_plot(p, sprintf("distributions_page_%03d", page),
        max(3.5, 2.2 * length(pages[[page]]) + 1.2))
      manifest <- rbind(manifest, data.frame(file = paths, type = "figure"))
    }
    if (per_species) for (s in seq_along(names)) {
      label <- gsub("[^A-Za-z0-9]+", "_", names[s])
      paths <- save_plot(plot_famasv_species(result, names[s]),
        sprintf("species_%03d_%s", s, substr(label, 1, 80)), 3.6)
      manifest <- rbind(manifest, data.frame(file = paths, type = "figure"))
    }
  }
  # Record versions after plotting has loaded ggplot2 and its dependencies.
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  manifest <- rbind(manifest, data.frame(file = file.path(output_dir,
    c("famasv_result.rds", "parameters.R", "sessionInfo.txt", "README.txt")), type = "record"))
  utils::write.csv(manifest, file.path(output_dir, "manifest.csv"), row.names = FALSE)
  invisible(manifest)
}
