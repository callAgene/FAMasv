test_that("hand-calculated sharing controls distinguish family specificity", {
  input <- make_test_counts(6)
  result <- run_famasv(input$counts, input$metadata, input$taxonomy,
    famasv_parameters(n_permutations = 99, n_bootstrap = 50))
  tab <- result$summary
  family <- tab[tab$Species == "Simulibacter familiaris", ]
  common <- tab[tab$Species == "Simulibacter communis", ]
  separate <- tab[tab$Species == "Simulibacter separatus", ]
  expect_equal(family$Observed_exact_ASV_sharing_rate, 1)
  expect_equal(family$Expected_exact_ASV_sharing_rate, 0)
  expect_equal(family$Observed_mean_ASV_Jaccard, 1)
  expect_equal(family$Expected_mean_ASV_Jaccard, 0)
  expect_equal(family$Sharing_permutation_P, 0.01)
  expect_equal(family$Jaccard_permutation_P, 0.01)
  expect_true(is.na(family$Sharing_Z))
  expect_equal(common$Observed_exact_ASV_sharing_rate, 1)
  expect_equal(common$Expected_exact_ASV_sharing_rate, 1)
  expect_equal(common$Excess_exact_ASV_sharing, 0)
  expect_equal(common$Sharing_permutation_P, 1)
  expect_equal(separate$Observed_exact_ASV_sharing_rate, 0)
  expect_equal(separate$Expected_exact_ASV_sharing_rate, 0)
  expect_equal(separate$Sharing_permutation_P, 1)
  expect_equal(tab$N_codetected_real_families[tab$Species == "Simulibacter unilateralis"], 0)
  expect_true(is.na(tab$Sharing_FDR[tab$Species == "Simulibacter rarus"]))
  expect_false("unclassified bacterium" %in% tab$Species)
  expect_true(all(result$permutation_audit$No_true_household_pairing))
  expect_true(all(result$permutation_audit$Every_partner_used_once))
})

test_that("non-household pool excludes carriers whose original partner is absent", {
  metadata <- data.frame(sample_id = paste0("s", 1:8),
    family_id = rep(paste0("F", 1:4), each = 2), group = rep(c("PD", "HC"), 4))
  counts <- matrix(0L, 3, 8, dimnames = list(paste0("a", 1:3), metadata$sample_id))
  counts[1, 1:2] <- 2L
  counts[2, 3:4] <- 2L
  counts[3, 5:6] <- 2L
  counts[1, 8] <- 2L # Extra HC carrier cannot enter the background: partner is absent.
  taxonomy <- data.frame(asv_id = rownames(counts), species = "Testus species")
  parameters <- famasv_parameters(minimum_codetected_families = 3,
    n_permutations = 99, n_bootstrap = 0)
  result <- run_famasv(counts, metadata, taxonomy, parameters)
  expect_equal(result$summary$N_group2_carriers, 4)
  expect_equal(result$summary$N_codetected_real_families, 3)
  expect_equal(result$summary$Expected_exact_ASV_sharing_rate, 0)
  expect_equal(result$permutation_audit$Unique_pairings, 2)
  expect_false("F4" %in% result$family_data$Family_ID)
})

test_that("Jaccard averages within-species sets and respects both possible derangements", {
  metadata <- data.frame(sample_id = paste0("s", 1:6),
    family_id = rep(paste0("F", 1:3), each = 2), group = rep(c("PD", "HC"), 3))
  # Family sets: {a,b}/{a}; {b}/{b,c}; {c}/{a,c}.
  counts <- cbind(c(2,2,0), c(2,0,0), c(0,2,0), c(0,2,2), c(0,0,2), c(2,0,2))
  dimnames(counts) <- list(c("a", "b", "c"), metadata$sample_id)
  taxonomy <- data.frame(asv_id = rownames(counts), species = "Testus species")
  result <- run_famasv(counts, metadata, taxonomy,
    famasv_parameters(minimum_codetected_families = 3, n_permutations = 999, n_bootstrap = 20))
  expect_equal(result$summary$Observed_mean_ASV_Jaccard, 1/2)
  null <- result$distributions[[1]]$null
  # The two legal cycles have means (1/3 + 0 + 0)/3 and (1/3 + 0 + 1/2)/3.
  expect_equal(sort(unique(round(null$Jaccard, 10))), round(c(1/9, 5/18), 10))
  expect_equal(sort(unique(null$Sharing)), c(1/3, 2/3))
  expect_equal(result$summary$Jaccard_permutation_P, 1/1000)
})

test_that("read threshold and number of distinct shared ASVs are separate", {
  input <- make_test_counts(6)
  low <- famasv_parameters(n_permutations = 29, n_bootstrap = 0)
  result1 <- run_famasv(input$counts, input$metadata, input$taxonomy, low)
  low$minimum_shared_asvs <- 2
  result2 <- run_famasv(input$counts, input$metadata, input$taxonomy, low)
  name <- "Simulibacter familiaris"
  expect_equal(result2$summary$Observed_exact_ASV_sharing_rate[result2$summary$Species == name], 0)
  expect_equal(result1$summary$Observed_mean_ASV_Jaccard, result2$summary$Observed_mean_ASV_Jaccard)
  low$minimum_shared_asvs <- 1
  low$detection_count <- 101
  result3 <- run_famasv(input$counts, input$metadata, input$taxonomy, low)
  expect_true(all(result3$summary$N_codetected_real_families == 0))
  expect_true(all(is.na(result3$summary$Sharing_FDR)))
})

test_that("BH is separate by metric and retained after selection", {
  input <- make_test_counts(8)
  result <- run_famasv(input$counts, input$metadata, input$taxonomy,
    famasv_parameters(n_permutations = 99, n_bootstrap = 0))
  tab <- subset(result$summary, Formal_analysis)
  expect_equal(tab$Sharing_FDR, p.adjust(tab$Sharing_permutation_P, "BH"))
  expect_equal(tab$Jaccard_FDR, p.adjust(tab$Jaccard_permutation_P, "BH"))
  selected <- select_famasv_species(result, species = "Simulibacter familiaris")
  expect_equal(selected$Sharing_FDR, tab$Sharing_FDR[tab$Species == selected$Species])
})

test_that("reproducibility is independent of input row order and caller RNG", {
  input <- make_test_counts(6)
  parameters <- famasv_parameters(n_permutations = 49, n_bootstrap = 25)
  set.seed(500)
  old_seed <- .Random.seed
  result1 <- run_famasv(input$counts, input$metadata, input$taxonomy, parameters)
  expect_identical(.Random.seed, old_seed)
  result2 <- run_famasv(input$counts[6:1, 12:1], input$metadata[12:1, ],
    input$taxonomy[6:1, ], parameters)
  expected <- result1$summary[result1$summary$Species == "Simulibacter familiaris", ]
  rownames(expected) <- NULL
  # Restricting analysis to a single species changes the multiple-testing family.
  expect_equal(result2$summary[, 1:19], expected[, 1:19])
  result3 <- run_famasv(input$counts[nrow(input$counts):1, 12:1], input$metadata[12:1, ],
    input$taxonomy[nrow(input$taxonomy):1, ], parameters)
  expect_identical(result1$summary, result3$summary)
  expect_identical(result1$distributions, result3$distributions)
})
