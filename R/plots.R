#' Select species after the complete analysis
#'
#' Filtering never recalculates P values or FDR. Explicit species names are
#' matched exactly. Set significant_only = FALSE to keep named nonsignificant
#' or untested species. Ranking and top_n are applied after filtering.
#' @param result Output of run_famasv.
#' @param species Optional vector of full species names.
#' @param metric Either sharing or jaccard.
#' @param significant_only Require FDR below fdr_cutoff.
#' @param fdr_cutoff FDR cutoff; NULL uses the analysis parameter.
#' @param min_excess Optional minimum excess (a fraction, not percentage points).
#' @param top_n Optional maximum number of rows after filtering.
#' @return A filtered result table retaining the original FDR values.
#' @export
select_famasv_species <- function(result, species = NULL,
    metric = c("sharing", "jaccard"), significant_only = TRUE,
    fdr_cutoff = NULL, min_excess = NULL, top_n = NULL) {
  metric <- match.arg(metric)
  if (!is.list(result) || !is.data.frame(result$summary)) stop("Expected a run_famasv result.")
  dat <- result$summary
  if (!is.null(species)) {
    missing <- setdiff(species, dat$Species)
    if (length(missing)) stop("Species not found: ", paste(missing, collapse = "; "))
    dat <- dat[dat$Species %in% species, , drop = FALSE]
  }
  if (length(significant_only) != 1L || !is.logical(significant_only) || is.na(significant_only))
    stop("significant_only must be TRUE or FALSE.")
  if (is.null(fdr_cutoff)) fdr_cutoff <- result$parameters$fdr_cutoff
  if (length(fdr_cutoff) != 1L || !is.numeric(fdr_cutoff) ||
      !is.finite(fdr_cutoff) || fdr_cutoff <= 0 || fdr_cutoff > 1) stop("Invalid fdr_cutoff.")
  fdr_column <- if (metric == "sharing") "Sharing_FDR" else "Jaccard_FDR"
  excess_column <- if (metric == "sharing") "Excess_exact_ASV_sharing" else "Excess_ASV_Jaccard"
  if (significant_only) dat <- dat[!is.na(dat[[fdr_column]]) & dat[[fdr_column]] < fdr_cutoff, , drop = FALSE]
  if (!is.null(min_excess)) {
    if (!is.numeric(min_excess) || length(min_excess) != 1L || !is.finite(min_excess)) stop("Invalid min_excess.")
    dat <- dat[!is.na(dat[[excess_column]]) & dat[[excess_column]] >= min_excess, , drop = FALSE]
  }
  dat <- dat[order(-dat[[excess_column]], dat$Species, na.last = TRUE), , drop = FALSE]
  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) ||
        top_n < 1 || top_n != floor(top_n)) stop("top_n must be a positive integer.")
    dat <- utils::head(dat, top_n)
  }
  rownames(dat) <- NULL
  dat
}

#' Plot the species sharing landscape
#' @inheritParams select_famasv_species
#' @return A ggplot object. Excess is shown in percentage points.
#' @export
plot_famasv_landscape <- function(result, metric = c("sharing", "jaccard"),
    top_n = 30L, significant_only = FALSE) {
  metric <- match.arg(metric)
  dat <- select_famasv_species(result, metric = metric,
    significant_only = significant_only, top_n = top_n)
  column <- if (metric == "sharing") "Excess_exact_ASV_sharing" else "Excess_ASV_Jaccard"
  fdr <- if (metric == "sharing") "Sharing_FDR" else "Jaccard_FDR"
  dat <- dat[is.finite(dat[[column]]), , drop = FALSE]
  if (!nrow(dat)) return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
    label = "No tested species meet the selection criteria.") + ggplot2::theme_void())
  dat$Excess_percent <- dat[[column]] * 100
  dat$Significance <- ifelse(dat[[fdr]] < result$parameters$fdr_cutoff,
    "FDR below cutoff", "FDR at or above cutoff")
  dat$Species <- factor(dat$Species, levels = rev(unique(dat$Species)))
  ggplot2::ggplot(dat, ggplot2::aes(x = Excess_percent, y = Species)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = Excess_percent,
      yend = Species), colour = "grey80", linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(colour = Significance,
      size = N_codetected_real_families)) +
    ggplot2::scale_colour_manual(values = c("FDR below cutoff" = "#B84446",
      "FDR at or above cutoff" = "#8995A1")) +
    ggplot2::labs(x = paste0("Excess ", metric, " (percentage points)"), y = NULL,
      size = "Co-carrier households", colour = NULL,
      title = "Household-associated ASV sharing",
      subtitle = paste0("Observed minus mean non-household expectation; FDR cutoff = ",
        result$parameters$fdr_cutoff)) + ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom")
}

#' Plot permutation and household-bootstrap distributions for selected species
#'
#' Both sharing and Jaccard are shown. A constant distribution is represented
#' by a vertical line instead of an artificial density. Display smoothing does
#' not change statistics. Species below the testing threshold remain explicitly
#' labelled untested. No P values or FDR are recalculated for selected plots.
#' @inheritParams select_famasv_species
#' @param sharing_bandwidth,jaccard_bandwidth Density bandwidths in percentage
#'   points. The defaults follow the research plotting notebooks.
#' @return A ggplot object with one row per species and two metric columns.
#' @export
plot_famasv_species <- function(result, species,
    sharing_bandwidth = 8, jaccard_bandwidth = 6) {
  if (!length(species)) stop("Select at least one species.")
  for (bw in list(sharing_bandwidth, jaccard_bandwidth)) {
    if (!is.numeric(bw) || length(bw) != 1L || !is.finite(bw) || bw <= 0) stop("Bandwidths must be positive.")
  }
  rows <- select_famasv_species(result, species, significant_only = FALSE)
  if (!nrow(rows)) stop("No species remain.")
  rows <- rows[order(match(rows$Species, species)), , drop = FALSE]
  density_rows <- list()
  point_rows <- list()
  line_rows <- list()
  label_rows <- list()
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, ]
    panel <- paste0(row$Species, "\nn = ", row$N_codetected_real_families)
    index <- which(vapply(result$distributions, function(x)
      identical(x$Species, row$Species), logical(1)))
    distribution <- result$distributions[[index]]
    for (m in c("Sharing", "Jaccard")) {
      metric <- if (m == "Sharing") "ASV sharing rate" else "Mean ASV Jaccard"
      observed <- if (m == "Sharing") row$Observed_exact_ASV_sharing_rate else row$Observed_mean_ASV_Jaccard
      expected <- if (m == "Sharing") row$Expected_exact_ASV_sharing_rate else row$Expected_mean_ASV_Jaccard
      p <- if (m == "Sharing") row$Sharing_permutation_P else row$Jaccard_permutation_P
      fdr <- if (m == "Sharing") row$Sharing_FDR else row$Jaccard_FDR
      label <- if (row$Formal_analysis) sprintf("P = %.3g; FDR = %.3g", p, fdr) else
        paste0("Not tested: fewer than ", result$parameters$minimum_codetected_families, " co-carrier households")
      label_rows[[length(label_rows) + 1L]] <- data.frame(Panel = panel, Metric = metric, Label = label)
      for (kind in c("null", "bootstrap")) {
        values <- distribution[[kind]][[m]] * 100
        if (!length(values)) next
        name <- if (kind == "null") "Non-household null" else "Household bootstrap"
        if (length(unique(values)) == 1L) {
          point_rows[[length(point_rows) + 1L]] <- data.frame(Panel = panel, Metric = metric,
            Value = values[1], Distribution = name)
        } else {
          d <- stats::density(values, bw = if (m == "Sharing") sharing_bandwidth else jaccard_bandwidth,
            from = 0, to = 100, n = 256)
          density_rows[[length(density_rows) + 1L]] <- data.frame(Panel = panel, Metric = metric,
            Value = d$x, Density = d$y, Distribution = name)
        }
      }
      for (statistic in c("Observed", "Null mean")) {
        value <- if (statistic == "Observed") observed else expected
        if (is.finite(value)) line_rows[[length(line_rows) + 1L]] <- data.frame(
          Panel = panel, Metric = metric, Value = value * 100, Statistic = statistic)
      }
    }
  }
  labels <- do.call(rbind, label_rows)
  panel_levels <- unique(labels$Panel)
  labels$Panel <- factor(labels$Panel, levels = panel_levels)
  p <- ggplot2::ggplot() +
    ggplot2::geom_blank(data = labels, ggplot2::aes(x = 0, y = 0))
  if (length(density_rows)) {
    densities <- do.call(rbind, density_rows)
    densities$Panel <- factor(densities$Panel, levels = panel_levels)
    p <- p + ggplot2::geom_area(data = densities,
      ggplot2::aes(x = Value, y = Density, fill = Distribution, colour = Distribution),
      position = "identity", alpha = 0.28, linewidth = 0.5)
  }
  if (length(point_rows)) {
    points <- do.call(rbind, point_rows)
    points$Panel <- factor(points$Panel, levels = panel_levels)
    p <- p + ggplot2::geom_vline(data = points,
      ggplot2::aes(xintercept = Value, colour = Distribution), linewidth = 2, alpha = 0.45)
  }
  if (length(line_rows)) {
    lines <- do.call(rbind, line_rows)
    lines$Panel <- factor(lines$Panel, levels = panel_levels)
    p <- p + ggplot2::geom_vline(data = lines,
      ggplot2::aes(xintercept = Value, linetype = Statistic), linewidth = 0.45)
  }
  colours <- c("Non-household null" = "#7F8891", "Household bootstrap" = "#315EAD")
  if (length(density_rows)) p <- p + ggplot2::scale_fill_manual(values = colours)
  if (length(density_rows) || length(point_rows))
    p <- p + ggplot2::scale_colour_manual(values = colours)
  if (length(line_rows)) p <- p + ggplot2::scale_linetype_manual(
    values = c("Observed" = "solid", "Null mean" = "dashed"))
  p + ggplot2::geom_text(data = labels, ggplot2::aes(x = 2, y = Inf, label = Label),
      hjust = 0, vjust = 1.3, size = 2.7) +
    ggplot2::facet_grid(Panel ~ Metric, scales = "free_y") +
    ggplot2::scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.28))) +
    ggplot2::labs(x = "Statistic (%)", y = "Density", colour = NULL, fill = NULL, linetype = NULL,
      subtitle = "One-sided permutation P; bootstrap shows uncertainty of observed households",
      caption = "Constant distributions are vertical lines. Display smoothing does not determine significance.") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom",
      strip.text.y = ggplot2::element_text(angle = 0, size = 8))
}

utils::globalVariables(c("Excess_percent", "Species", "Significance",
  "N_codetected_real_families", "Value", "Density", "Distribution", "Statistic", "Label"))
