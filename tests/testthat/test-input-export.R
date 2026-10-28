test_that("real oral example has only de-identified analysis fields", {
  path <- system.file("extdata", package = "FAMasv")
  input <- read_famasv(file.path(path, "counts.csv"), file.path(path, "metadata.csv"),
    file.path(path, "taxonomy.csv"))
  expect_equal(dim(input$counts), c(10111L, 100L))
  expect_identical(names(input$metadata), c("sample_id", "family_id", "group"))
  expect_identical(names(input$taxonomy), c("asv_id", "species"))
  expect_true(all(grepl("^ORAL_S[0-9]{4}$", input$metadata$sample_id)))
  expect_true(all(grepl("^ORAL_F[0-9]{3}$", input$metadata$family_id)))
  expect_true(all(grepl("^ORAL_A[0-9]{5}$", input$taxonomy$asv_id)))
  expect_true(all(table(input$metadata$family_id, input$metadata$group) == 1))
  expect_equal(length(unique(input$taxonomy$species)), 508)
  expect_identical(colnames(input$counts), input$metadata$sample_id)
  expect_identical(rownames(input$counts), input$taxonomy$asv_id)
  expect_true(all(input$counts >= 0 & input$counts == floor(input$counts)))
  expect_false("simulate_famasv" %in% getNamespaceExports("FAMasv"))
})

test_that("invalid inputs fail informatively", {
  input <- make_test_counts(6)
  expect_equal(dim(input$counts), c(48L, 12L))
  expect_identical(names(input$metadata), c("sample_id", "family_id", "group"))
  parameters <- famasv_parameters(n_permutations = 19, n_bootstrap = 0)
  expect_error(run_famasv(input$counts / 7, input$metadata, input$taxonomy, parameters), "integer")
  metadata <- input$metadata
  metadata$family_id[3:4] <- metadata$family_id[1]
  expect_error(run_famasv(input$counts, metadata, input$taxonomy, parameters), "Multiple samples")
  metadata <- input$metadata[-1, ]
  result <- run_famasv(input$counts[, -1], metadata, input$taxonomy, parameters)
  expect_equal(sum(!result$pair_audit$Included), 1)
  expect_error(run_famasv(input$counts, metadata, input$taxonomy, parameters), "same sample IDs")
  expect_error(famasv_parameters(minimum_codetected_families = 2), "integer")
})

test_that("degenerate and untested plots build without warnings", {
  input <- make_test_counts(6)
  result <- run_famasv(input$counts, input$metadata, input$taxonomy,
    famasv_parameters(n_permutations = 29, n_bootstrap = 20))
  for (species in c("Simulibacter familiaris", "Simulibacter rarus",
                    "Simulibacter unilateralis", "Simulibacter randomus")) {
    expect_warning(ggplot2::ggplot_build(plot_famasv_species(result, species)), NA)
  }
  expect_error(plot_famasv_species(result, "Unknown species"), "not found")
})

test_that("exports support all, selected and empty significant subsets", {
  input <- make_test_counts(6)
  parameters <- famasv_parameters(n_permutations = 29, n_bootstrap = 20)
  result <- run_famasv(input$counts, input$metadata, input$taxonomy, parameters)
  directory <- tempfile("famasv-test-")
  on.exit(unlink(directory, recursive = TRUE))
  manifest <- export_famasv(result, directory,
    species = c("Simulibacter familiaris", "Simulibacter rarus"), formats = "pdf")
  expect_true(all(file.exists(manifest$file)))
  selected <- read.csv(file.path(directory, "tables", "selected_species.csv"))
  expect_setequal(selected$Species, c("Simulibacter familiaris", "Simulibacter rarus"))
  expect_equal(nrow(read.csv(file.path(directory, "tables", "species_summary.csv"))), 7)
  expect_false("Habitat" %in% names(result$summary))
  expect_false(any(grepl("habitat", basename(manifest$file))))
  expect_error(export_famasv(result, directory), "not empty")
  result$parameters$fdr_cutoff <- 1e-10
  empty_dir <- file.path(directory, "empty")
  expect_warning(export_famasv(result, empty_dir, significant_only = TRUE, formats = "pdf"), NA)
  expect_equal(nrow(read.csv(file.path(empty_dir, "tables", "selected_species.csv"))), 0)
  expect_equal(length(list.files(file.path(empty_dir, "figures"))), 2)
})
