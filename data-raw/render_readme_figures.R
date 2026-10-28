# Run from the repository root after installing FAMasv.
# ============================================================
# Parameters
# ============================================================
result_file <- file.path("FAMasv_oral_example", "preview", "famasv_result.rds")
output_dir <- file.path("man", "figures")
overview_species <- 20L
distribution_species <- 3L

library(FAMasv)
if (!file.exists(result_file)) {
  stop("Run inst/examples/run_oral_example.R first or change result_file.")
}
result <- readRDS(result_file)
selected <- select_famasv_species(result, significant_only = FALSE,
                                  top_n = distribution_species)
overview <- plot_famasv_landscape(result, top_n = overview_species)
distributions <- plot_famasv_species(result, selected$Species)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(output_dir, "README-oral-landscape.png"), overview,
                width = 10, height = 7, dpi = 160, bg = "white")
ggplot2::ggsave(file.path(output_dir, "README-oral-distributions.png"), distributions,
                width = 10, height = 8, dpi = 160, bg = "white")
