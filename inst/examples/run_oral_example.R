# ============================================================
# 0. Parameters
# ============================================================
input_dir <- system.file("extdata", package = "FAMasv")
output_root <- "FAMasv_oral_example"
n_permutations <- 10000L
n_bootstrap <- 10000L
detection_count <- 1L
minimum_shared_asvs <- 1L
minimum_codetected_families <- 5L
fdr_cutoff <- 0.05
exclude_unclassified <- FALSE # Bundled taxa were filtered on species-rank annotation.
preview_species_n <- 6L
selected_species <- c("Lachnospira pectinoschiza", "Prevotella jejuni")
export_all_species <- FALSE # TRUE additionally exports every formally tested taxon.
figure_formats <- c("pdf", "png")

# ============================================================
# 1. Read the single de-identified oral example
# ============================================================
library(FAMasv)
input <- read_famasv(file.path(input_dir, "counts.csv"),
                    file.path(input_dir, "metadata.csv"),
                    file.path(input_dir, "taxonomy.csv"))

# ============================================================
# 2. Analyze all taxa before selecting figures
# ============================================================
parameters <- famasv_parameters(detection_count = detection_count,
  minimum_shared_asvs = minimum_shared_asvs,
  minimum_codetected_families = minimum_codetected_families,
  n_permutations = n_permutations, n_bootstrap = n_bootstrap,
  fdr_cutoff = fdr_cutoff, exclude_unclassified = exclude_unclassified)
result <- run_famasv(input$counts, input$metadata, input$taxonomy, parameters)

# ============================================================
# 3. Preview, manual selection and optional complete batch export
# ============================================================
preview <- select_famasv_species(result, significant_only = FALSE,
                                 top_n = preview_species_n)
export_famasv(result, file.path(output_root, "preview"),
  species = preview$Species, formats = figure_formats, species_per_page = 3L)
export_famasv(result, file.path(output_root, "selected"),
  species = selected_species, formats = figure_formats)
if (export_all_species) {
  export_famasv(result, file.path(output_root, "all_species"), formats = figure_formats)
}
print(result$summary[result$summary$Formal_analysis,
  c("Species", "N_codetected_real_families", "Excess_exact_ASV_sharing", "Sharing_FDR")])
