# ============================================================
# 0. Parameters: replace these paths with your own input files
# ============================================================
input_dir <- system.file("extdata", package = "FAMasv")
output_dir <- "FAMasv_csv_output"
group_labels <- c("PD", "HC")
detection_count <- 1L
minimum_shared_asvs <- 1L
minimum_codetected_families <- 5L
n_permutations <- 10000L
n_bootstrap <- 10000L
fdr_cutoff <- 0.05
exclude_unclassified <- FALSE # Bundled input is already filtered; review for your own data.

# ============================================================
# 1. Load package and read data
# ============================================================
library(FAMasv)
input <- read_famasv(file.path(input_dir, "counts.csv"),
                    file.path(input_dir, "metadata.csv"),
                    file.path(input_dir, "taxonomy.csv"))

# ============================================================
# 2. Run analysis
# ============================================================
parameters <- famasv_parameters(detection_count = detection_count,
  minimum_shared_asvs = minimum_shared_asvs,
  minimum_codetected_families = minimum_codetected_families,
  n_permutations = n_permutations, n_bootstrap = n_bootstrap,
  fdr_cutoff = fdr_cutoff, exclude_unclassified = exclude_unclassified)
result <- run_famasv(input$counts, input$metadata, input$taxonomy,
                     parameters, group_labels = group_labels)

# ============================================================
# 3. Export all tested species
# ============================================================
export_famasv(result, output_dir)
selected_species <- select_famasv_species(result, metric = "sharing",
  significant_only = TRUE, min_excess = 0, top_n = 10)
print(selected_species)
