###
# proteomics_data_processing.R
#
# - Visualize data
# - Filter high-missingness proteins
# - Filter low-variance proteins
# - Generate bioimpedance deltas (delta_PBF and delta_PSMM)
# - Export processed data and patient metadata for linear regression script
###

library(dplyr)
library(tidyr)
library(ggplot2)
library(matrixStats)
source(here::here("common_functions.R"))


###
# To run code as a function:
run_proteomics_data_processing <- function(subset_HCT_type, low_variance_cutoff) {
###
# subset_HCT_type <- "ALL"
# low_variance_cutoff <- proteomics_low_variance_cutoff

# Validate input
if (low_variance_cutoff < 0 || low_variance_cutoff > 1) {
  stop("low_variance_cutoff must be between 0 and 1 (e.g. 0.1 for 10%)")
}


###
# Functions
###
filter_low_variance_proteins <- function(data_df, ptid_visit_colnames, low_variance_cutoff, report_file, output_dir, filename_suffix, font_size=5, id_colname="OlinkID") {
  protein_vars <- apply(data_df[, ptid_visit_colnames], 1, var, na.rm = TRUE)
  names(protein_vars) <- data_df[[id_colname]]

  if (any(is.na(protein_vars))) {
    warning("Some proteins have NA variance and will be removed.")
    protein_vars <- protein_vars[!is.na(protein_vars)]
  }

  # Plot variance distribution with candidate cutoffs
  cut_probs <- c(0.01, 0.05, 0.10, 0.15, 0.25)
  cuts <- quantile(protein_vars, probs = cut_probs)
  line_colors <- c("#919191", "#686868", "#4d4d4d", "#363636", "#000000")
  cutoff_match_idx <- which(cut_probs == low_variance_cutoff)
  if (length(cutoff_match_idx) > 0) {
    line_colors[cutoff_match_idx] <- "#f86161"
  }
  cut_df <- data.frame(
    label = factor(paste0(names(cuts), ": ", round(cuts, 4)),
                   levels = paste0(names(cuts), ": ", round(cuts, 4))),
    value = as.numeric(cuts),
    color = line_colors,
    stringsAsFactors = FALSE
  )
  p_var <- ggplot(data.frame(variance = protein_vars), aes(x = .data$variance)) +
    geom_histogram(bins = 300, fill = "lightgray", color = "black", linewidth = 0.2) +
    geom_vline(data = cut_df, aes(xintercept = .data$value, color = .data$label),
               linetype = "dashed", linewidth = 0.8) +
    scale_color_manual(values = setNames(line_colors, levels(cut_df$label)), name = NULL) +
    coord_cartesian(xlim = c(0, max(cuts["25%"] * 2, 3))) +
    labs(x = "Variance", y = "Count") +
    theme_bw(base_size = font_size) +
    theme(
      legend.position = "right",
      panel.grid = element_blank()
    ) +
    guides(color = guide_legend(byrow = TRUE, override.aes = list(linetype = "solid", linewidth = 1)))
  ggsave(
    filename = file.path(output_dir, paste0("protein_variance_distribution", filename_suffix, ".png")),
    plot = p_var,
    width = 8, height = 6, dpi = 150
  )

  variance_threshold <- quantile(protein_vars, probs = low_variance_cutoff)
  proteins_to_keep <- names(protein_vars)[protein_vars >= variance_threshold]
  proteins_removed <- nrow(data_df) - length(proteins_to_keep)
  data_df <- data_df %>% dplyr::filter(!!sym(id_colname) %in% proteins_to_keep)
  message(paste0("Filtered out ", proteins_removed, " low-variance proteins (bottom ", low_variance_cutoff * 100,
                "%). Kept ", length(proteins_to_keep), " proteins."))
  append_report_line(report_file, "")
  append_report_line(report_file, "--- Low-Variance Protein Filtering ---")
  append_report_line(report_file, paste0("Variance threshold (", low_variance_cutoff * 100, "th percentile): ", round(variance_threshold, 6)))
  append_report_line(report_file, paste0("Proteins removed: ", proteins_removed))
  append_report_line(report_file, paste0("Proteins kept: ", length(proteins_to_keep)))

  return(data_df)
}


plot_outcome_var_boxplot <- function(patient_meta, output_dir, filename_suffix, font_size, timepoint_order = c("Baseline", "6m", "12m"), outcome_var = "VO2peak") {
  outcome_long <- patient_meta %>%
    dplyr::select(PTID, dplyr::starts_with(paste0(outcome_var, "_"))) %>%
    tidyr::pivot_longer(cols = dplyr::starts_with(paste0(outcome_var, "_")), names_to = "visit", values_to = outcome_var) %>%
    mutate(visit = sub(paste0(outcome_var, "_"), "", visit)) %>%
    mutate(visit = factor(visit, levels = timepoint_order))

  outcome_var_label <- gsub("VO2", "VO\u2082", outcome_var)
  p <- ggplot(outcome_long, aes(x = visit, y = !!sym(outcome_var))) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.4, fill = "#CFE3F2") +
    geom_jitter(width = 0.15, height = 0, alpha = 0.8, size = 1.6, color = "#2E5C8A") +
    labs(x = "Timepoint", y = paste0(outcome_var_label, " (ml/kg/min)")) +
    theme_bw(base_size = font_size) +
    theme(plot.title = element_blank(), panel.grid = element_blank())

  ggsave(file.path(output_dir, paste0(outcome_var, "_distribution_by_timepoint", filename_suffix, ".png")), p, width = 9, height = 7, dpi = 600)
}

plot_outcome_var_trendlines <- function(patient_meta, output_dir, filename_suffix, font_size, delta_cutoff = 2, y_range=c(5,48), outcome_var = "VO2peak") {
  # Determine trend direction per patient using delta (Baseline -> 12m only)
  # Up: delta > delta_cutoff, Static: -delta_cutoff <= delta <= delta_cutoff, Down: delta < -delta_cutoff
  trend_df <- patient_meta %>%
    dplyr::select(PTID, dplyr::starts_with(paste0(outcome_var, "_Baseline")), dplyr::starts_with(paste0(outcome_var, "_12m"))) %>%
    dplyr::mutate(
      delta = dplyr::case_when(
        !is.na(!!sym(paste0(outcome_var, "_Baseline"))) & !is.na(!!sym(paste0(outcome_var, "_12m"))) ~ !!sym(paste0(outcome_var, "_12m")) - !!sym(paste0(outcome_var, "_Baseline")),
        TRUE ~ NA_real_
      ),
      trend_group = dplyr::case_when(
        !is.na(delta) & delta < -delta_cutoff ~ "down",
        !is.na(delta) & delta >= -delta_cutoff & delta <= delta_cutoff ~ "static",
        !is.na(delta) & delta >  delta_cutoff ~ "up",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(PTID, trend_group)

  # Only show Baseline and 12m (exclude 6m)
  plot_visits <- c("Baseline", "12m")

  outcome_var_long_all <- patient_meta %>%
    dplyr::select(PTID, dplyr::starts_with(paste0(outcome_var, "_Baseline")), dplyr::starts_with(paste0(outcome_var, "_12m"))) %>%
    tidyr::pivot_longer(cols = dplyr::starts_with(paste0(outcome_var, "_")), names_to = "visit", values_to = outcome_var) %>%
    mutate(visit = sub(paste0(outcome_var, "_"), "", visit)) %>%
    dplyr::filter(!is.na(!!sym(outcome_var)), visit %in% plot_visits) %>%
    mutate(visit = factor(visit, levels = plot_visits)) %>%
    dplyr::left_join(trend_df, by = "PTID")

  outcome_var_label <- gsub("VO2", "VO\u2082", outcome_var)
  group_configs <- list(
    list(group = "down",   color = "#2166AC", panel_title = paste0("\u0394", outcome_var_label, " < ", -delta_cutoff)),
    list(group = "static", color = "gray",      panel_title = paste0("|\u0394", outcome_var_label, "| \u2264 ", delta_cutoff)),
    list(group = "up",     color = "indianred", panel_title = paste0("\u0394", outcome_var_label, " > ", delta_cutoff))
  )

  plot_list <- list()
  for (i in seq_along(group_configs)) {
    cfg <- group_configs[[i]]
    outcome_var_long_subset <- outcome_var_long_all %>% dplyr::filter(trend_group == cfg$group)
    n_patients <- sum(trend_df$trend_group == cfg$group, na.rm = TRUE)

    p <- ggplot(outcome_var_long_all, aes(x = visit, y = !!sym(outcome_var))) +
      geom_line(data = outcome_var_long_subset, aes(group = PTID), color = cfg$color, alpha = 0.5, linewidth = 1.2) +
      geom_point(data = outcome_var_long_subset, alpha = 0.5, size = 3, color = "black") +
      coord_cartesian(ylim = y_range) +
      annotate("text", x = -Inf, y = Inf, label = paste0("n=", n_patients),
               hjust = -0.1, vjust = 1.5, size = font_size / ggplot2::.pt, color = "black") +
      labs(x="",
      y = if (i == 1) paste0(outcome_var_label, " (ml/kg/min)") else NULL,
           title = cfg$panel_title) +
      theme_bw(base_size = font_size) +
      theme(
               plot.title    = element_text(size = font_size, hjust = 0.5, color = "black"),
        panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid    = element_blank(),
        axis.title    = element_text(size = font_size, color = "black"),
        axis.text     = element_text(size = font_size, color = "black"),
        axis.text.y        = if (i == 1) element_text(size = font_size, color = "black") else element_blank(),
        axis.ticks         = element_line(color = "black", linewidth = 0.5),
        axis.ticks.y       = if (i == 1) element_line(color = "black", linewidth = 0.5) else element_blank(),
        axis.ticks.length  = unit(5, "pt"),
        plot.margin        = if (i == 1) margin(5, 0, 5, 5) else if (i == length(group_configs)) margin(5, 5, 5, 0) else margin(5, 0, 5, 0)
      )

    plot_list[[i]] <- p
  }

  combined_plot <- patchwork::wrap_plots(plot_list, nrow = 1)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    file.path(output_dir, paste0(outcome_var, "_trendlines", filename_suffix, ".png")),
    combined_plot, width = 9, height = 6, dpi = 600
  )
}

plot_delta_outcome_var_histogram <- function(patient_metadata_df, output_dir, filename_suffix, font_size, plot_title = "", cutoff_labels = data.frame(cutoff = c(-2,2)), outcome_var = "VO2peak") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  delta_configs <- list(
    list(col = paste0("delta_", outcome_var, "_6m_Baseline"),  label = "6m - Baseline",  file_tag = "_delta_6m", border_color = "#000000"),
    list(col = paste0("delta_", outcome_var, "_12m_Baseline"), label = "12m - Baseline", file_tag = "_delta_12m", border_color = "#000000")
  )

  for (cfg in delta_configs) {
    vals <- patient_metadata_df[[cfg$col]]
    if (is.null(vals) || all(is.na(vals))) next

    outcome_var_label <- gsub("VO2", "VO\u2082", outcome_var)
    plot_distribution_histogram(
      data         = vals,
      title        = plot_title,
      xlab         = paste0("\u0394", outcome_var_label, " (", cfg$label, ")"),
      filename     = file.path(output_dir, paste0("delta_", tolower(outcome_var), cfg$file_tag, filename_suffix, ".png")),
      text_size    = font_size,
      bin_num = 20,
      vertical_cutoff_labels = cutoff_labels,
      region_colors = c("below" = "#2166AC", "between" = "gray", "above" = "indianred"),
      width = 6,
      height = 6,
      y_label = "Number of Patients"
    )
  }
}

plot_outcome_var_violin_plot <- function(patient_metadata_df, output_dir, filename_suffix, font_size, timepoint_order = c("Baseline", "6m", "12m"), colors = timepoint_palette, plot_title="", outcome_var = "VO2peak") {
  outcome_var_long <- patient_metadata_df %>%
    dplyr::select(PTID, paste0(outcome_var, "_Baseline"), paste0(outcome_var, "_6m"), paste0(outcome_var, "_12m")) %>%
    tidyr::pivot_longer(cols = starts_with(paste0(outcome_var, "_")), names_to = "visit", values_to = outcome_var) %>%
    mutate(visit = sub(paste0(outcome_var, "_"), "", visit)) %>%
    mutate(visit = factor(visit, levels = timepoint_order)) %>%
    dplyr::filter(!is.na(!!sym(outcome_var)))

  # Apply colors
  p <- ggplot(outcome_var_long, aes(x = visit, y = !!sym(outcome_var), fill = visit, color = visit)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.6, size = 1.4) +
    geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.8, fill = "white", color = "black", alpha = 0.3) +
    scale_fill_manual(values = setNames(colors, levels(outcome_var_long$visit))) +
    scale_color_manual(values = setNames(colors, levels(outcome_var_long$visit))) +
    coord_cartesian(ylim = c(0, 50)) +
    labs(x = "Timepoint", y = gsub("VO2", "VO\u2082", outcome_var)) +
    theme_bw(base_size = font_size) +
    theme(
      plot.title   = element_blank(),
      axis.title   = element_text(size = font_size),
      axis.text    = element_text(size = font_size),
      legend.text  = element_text(size = font_size),
      legend.title = element_text(size = font_size),
      panel.grid   = element_blank(),
      legend.position = "none"
    )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(output_dir, paste0(tolower(outcome_var), "_violin_by_timepoint", filename_suffix, ".png")), p, width = 9, height = 7, dpi = 600)
}

plot_outcome_var_violin_allo_vs_auto <- function(patient_metadata_df, output_dir, filename_suffix, font_size, timepoint_order = c("Baseline", "6m", "12m"), colors = c("ALLO" = "#5fa2b6", "AUTO" = "#0a3468"), outcome_var = "VO2peak") {
  outcome_var_long <- patient_metadata_df %>%
    dplyr::filter(!is.na(TRANSPLANT_TYPE) & TRANSPLANT_TYPE %in% c("ALLO", "AUTO")) %>%
    dplyr::select(PTID, TRANSPLANT_TYPE, paste0(outcome_var, "_Baseline"), paste0(outcome_var, "_6m"), paste0(outcome_var, "_12m")) %>%
    tidyr::pivot_longer(cols = starts_with(paste0(outcome_var, "_")), names_to = "visit", values_to = outcome_var) %>%
    mutate(visit = sub(paste0(outcome_var, "_"), "", visit)) %>%
    mutate(visit = factor(visit, levels = timepoint_order)) %>%
    dplyr::filter(!is.na(!!sym(outcome_var))) %>%
    dplyr::mutate(TRANSPLANT_TYPE = factor(TRANSPLANT_TYPE, levels = c("ALLO", "AUTO")))

  # One-sided Wilcoxon rank-sum test per timepoint: H_a: ALLO < AUTO
  test_results <- dplyr::bind_rows(lapply(timepoint_order, function(tp) {
    tp_data <- outcome_var_long %>% dplyr::filter(visit == tp)
    allo_vals <- tp_data %>% dplyr::filter(TRANSPLANT_TYPE == "ALLO") %>% dplyr::pull(!!sym(outcome_var))
    auto_vals <- tp_data %>% dplyr::filter(TRANSPLANT_TYPE == "AUTO") %>% dplyr::pull(!!sym(outcome_var))
    if (length(allo_vals) >= 3 && length(auto_vals) >= 3) {
      wt <- wilcox.test(allo_vals, auto_vals, alternative = "less", exact = FALSE)
      tibble::tibble(visit = tp, n_allo = length(allo_vals), n_auto = length(auto_vals), p_value = wt$p.value)
    } else {
      tibble::tibble(visit = tp, n_allo = length(allo_vals), n_auto = length(auto_vals), p_value = NA_real_)
    }
  }))
  test_results <- test_results %>%
    dplyr::mutate(
      p_adj = stats::p.adjust(p_value, method = "BH"),
      sig_label = dplyr::case_when(
        is.na(p_adj)  ~ "NA",
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ "ns"
      ),
      visit = factor(visit, levels = timepoint_order)
    )

  p <- ggplot(outcome_var_long, aes(x = visit, y = !!sym(outcome_var), fill = TRANSPLANT_TYPE, color = TRANSPLANT_TYPE)) +
    geom_violin(trim = FALSE, alpha = 0.4, position = position_dodge(0.8)) +
    geom_jitter(aes(group = TRANSPLANT_TYPE), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), alpha = 0.6, size = 1.4) +
    geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.8, fill = "white", color = "black", alpha = 0.3, position = position_dodge(0.8)) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    coord_cartesian(ylim = c(0, 50)) +
    labs(x = "Timepoint", y = gsub("VO2", "VO\u2082", outcome_var), fill = "HCT Type", color = "HCT Type") +
    theme_bw(base_size = font_size) +
    theme(
      plot.title   = element_blank(),
      axis.title   = element_text(size = font_size),
      axis.text    = element_text(size = font_size),
      legend.text  = element_text(size = font_size),
      legend.title = element_text(size = font_size),
      panel.grid   = element_blank()
    )

  # Add significance annotations per timepoint
  y_max <- 50
  y_bracket <- y_max * 0.88
  bracket_height <- y_max * 0.02
  for (i in seq_len(nrow(test_results))) {
    x_pos <- as.numeric(test_results$visit[i])
    label_text <- paste0(test_results$sig_label[i], "\np=", signif(test_results$p_adj[i], 2))
    p <- p +
      annotate("segment", x = x_pos - 0.2, xend = x_pos + 0.2, y = y_bracket, yend = y_bracket, linewidth = 0.4) +
      annotate("segment", x = x_pos - 0.2, xend = x_pos - 0.2, y = y_bracket, yend = y_bracket - bracket_height, linewidth = 0.4) +
      annotate("segment", x = x_pos + 0.2, xend = x_pos + 0.2, y = y_bracket, yend = y_bracket - bracket_height, linewidth = 0.4) +
      annotate("text", x = x_pos, y = y_bracket + bracket_height * 1.5, label = label_text, size = font_size / ggplot2::.pt, vjust = 0, lineheight = 0.9)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(output_dir, paste0(outcome_var, "_violin_allo_vs_auto", filename_suffix, ".png")), p, width = 10, height = 7, dpi = 600)

  # Save test results table
  write.csv(test_results, file.path(output_dir, paste0(outcome_var, "_allo_vs_auto_wilcoxon", filename_suffix, ".csv")), row.names = FALSE)
  invisible(test_results)
}

plot_outcome_var_violin_plot_one_timepoint <- function(patient_metadata_df, output_dir, filename_suffix, font_size,
                                                   timepoint = "Baseline",
                                                   color_all  = "#494c4e",
                                                   color_auto = "#4c488d",
                                                   color_allo = "#5fa2b6",
                                                   outcome_var = "VO2peak") {
  outcome_var_col <- paste0(outcome_var, "_", timepoint)
  if (!outcome_var_col %in% colnames(patient_metadata_df)) {
    message(paste0("Warning: Column '", outcome_var_col, "' not found in patient metadata. Skipping single-timepoint violin plot."))
    return(invisible(NULL))
  }

  all_df <- patient_metadata_df %>%
    dplyr::select(PTID, !!sym(outcome_var) := all_of(outcome_var_col)) %>%
    dplyr::filter(!is.na(!!sym(outcome_var))) %>%
    dplyr::mutate(group = "All Patients")

  auto_df <- patient_metadata_df %>%
    dplyr::filter(TRANSPLANT_TYPE == "AUTO") %>%
    dplyr::select(PTID, !!sym(outcome_var) := all_of(outcome_var_col)) %>%
    dplyr::filter(!is.na(!!sym(outcome_var))) %>%
    dplyr::mutate(group = "AUTO")

  allo_df <- patient_metadata_df %>%
    dplyr::filter(TRANSPLANT_TYPE == "ALLO") %>%
    dplyr::select(PTID, !!sym(outcome_var) := all_of(outcome_var_col)) %>%
    dplyr::filter(!is.na(!!sym(outcome_var))) %>%
    dplyr::mutate(group = "ALLO")

  y_lim <- c(0, 50)

  # --- Subplot 1: All Patients ---
  p1 <- ggplot(all_df, aes(x = group, y = !!sym(outcome_var), fill = group, color = group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.6, size = 1.4) +
    geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.8, fill = "white", color = "black", alpha = 0.3) +
    scale_fill_manual(values = c("All Patients" = color_all)) +
    scale_color_manual(values = c("All Patients" = color_all)) +
    coord_cartesian(ylim = y_lim) +
    labs(x = "", y = paste0(timepoint, " ", gsub("VO2", "VO\u2082", outcome_var))) +
    theme_bw(base_size = font_size) +
    theme(
      plot.title      = element_blank(),
      axis.title      = element_text(size = font_size, color = "black"),
      axis.text       = element_text(size = font_size, color = "black"),
      legend.position    = "none",
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid         = element_blank(),
      axis.ticks         = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length  = unit(5, "pt"),
      plot.margin        = margin(5, 0, 5, 5)
    )

  # --- One-sided Wilcoxon rank-sum test: H_a: ALLO < AUTO ---
  allo_vals <- allo_df[[outcome_var]]
  auto_vals <- auto_df[[outcome_var]]
  wt_p <- if (length(allo_vals) >= 3 && length(auto_vals) >= 3) {
    wilcox.test(allo_vals, auto_vals, alternative = "less", exact = FALSE)$p.value
  } else {
    NA_real_
  }
  sig_label <- dplyr::case_when(
    is.na(wt_p)  ~ "NA",
    wt_p < 0.001 ~ "***",
    wt_p < 0.01  ~ "**",
    wt_p < 0.05  ~ "*",
    TRUE         ~ "ns"
  )
  label_text   <- paste0(sig_label, "\np=", signif(wt_p, 2))
  y_bracket      <- y_lim[2] * 0.88
  bracket_height <- y_lim[2] * 0.02

  # --- Subplot 2: AUTO vs ALLO with significance bracket ---
  type_df <- dplyr::bind_rows(auto_df, allo_df) %>%
    dplyr::mutate(group = factor(group, levels = c("AUTO", "ALLO")))

  p2 <- ggplot(type_df, aes(x = group, y = !!sym(outcome_var), fill = group, color = group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.6, size = 1.4) +
    geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.8, fill = "white", color = "black", alpha = 0.3) +
    scale_fill_manual(values = c("AUTO" = color_auto, "ALLO" = color_allo)) +
    scale_color_manual(values = c("AUTO" = color_auto, "ALLO" = color_allo)) +
    coord_cartesian(ylim = y_lim) +
    annotate("segment", x = 1, xend = 2, y = y_bracket, yend = y_bracket, linewidth = 0.4) +
    annotate("segment", x = 1, xend = 1, y = y_bracket, yend = y_bracket - bracket_height, linewidth = 0.4) +
    annotate("segment", x = 2, xend = 2, y = y_bracket, yend = y_bracket - bracket_height, linewidth = 0.4) +
    annotate("text", x = 1.5, y = y_bracket + bracket_height * 1.5,
             label = label_text, size = font_size / ggplot2::.pt, vjust = 0, lineheight = 0.9) +
    labs(x = "", y = NULL) +
    theme_bw(base_size = font_size) +
    theme(
      plot.title      = element_blank(),
      axis.title      = element_text(size = font_size, color = "black"),
      axis.text.x          = element_text(size = font_size, color = "black"),
      axis.text.y          = element_blank(),
      panel.border         = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.ticks.x         = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length.x  = unit(5, "pt"),
      axis.ticks.y         = element_blank(),
      axis.ticks.length.y  = unit(0, "pt"),
      axis.line.y          = element_blank(),
      legend.position      = "none",
      panel.grid           = element_blank(),
      plot.margin          = margin(5, 5, 5, 0)
    )

  combined_plot <- patchwork::wrap_plots(p1, p2, nrow = 1, widths = c(1, 2))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    file.path(output_dir, paste0(outcome_var, "_violin_", tolower(timepoint), "_all_auto_allo", filename_suffix, ".png")),
      combined_plot, width = 6, height = 6, dpi = 600
  )
}


###
# Values
###
input_folder <- get_input_folder()
output_folder <- get_output_folder()
output_script_folder <- "proteomics_data_processing"

# Protein data input; rows of OlinkIDs, columns of PTID_visit, values of NPX_mean:
batch_correction_subfolder <- "batch_correction_output" # in output_folder
proteomics_input_filename <- "proteomics_wide_corrected.csv"
# Patient metadata input: columns include PTID, covariates:
# e.g. AgeatHCT, BMI, gender, Diabetes, Hypertension, TRANSPLANT_TYPE
# - numerical covariates imputed with median below
# - categorical covariates imputed with mode below
patient_metadata_filename <- "Patient_data.csv"
# VO2 data; columns of PTID, visit, cpet_vo2_adjusted_num
vo2_data_filename <- "CRESTpfizer_VO2peak_100125_filtered.csv"
# Bioimpedance data; columns of PTID, BFM, LBM, SMM, BMI, PBF
bioimpedance_filename <- "CLEAN_Baseline_BIA_03.10.2025.xlsx"
# Protein annotations for correlation plot labels
olink_annotations_filename <- "olink_mapped.csv"

feature_id_colname <- prot_id_colname

# Set options (font_size comes from config.R)
pca_var_top_pct <- proteomics_pca_var_top_pct # keep top n percent in PCA for selected proteins


###
# Take into account HCT type subsets
###
if (subset_HCT_type == "ALLO" || subset_HCT_type == "AUTO") {
  output_script_folder <- paste0(subset_HCT_type, "_", output_script_folder)
  filename_suffix <- paste0("_", subset_HCT_type)
} else if (subset_HCT_type == "ALL") {
  filename_suffix <- ""
} else {
  stop("Invalid subset_HCT_type specified. Please choose 'ALLO', 'AUTO', or 'ALL'.")
}


###
# Organize Directories
###
output_dir <- file.path(output_folder, output_script_folder)
protein_distribution_dir <- file.path(output_dir, "protein_distribution")
sample_distribution_dir <- file.path(output_dir, "sample_protein_distributions")
vo2_distribution_output_dir <- file.path(output_dir, "vo2_distributions")
pca_output_dir <- file.path(output_dir, "pca_plots")
top_pct_output_dir <- file.path(output_dir, paste0("pca_plots_top", pca_var_top_pct, "pct_variance"))
bioimp_corr_output_dir <- file.path(output_dir, "bioimpedance_correlations")

ensure_dir(top_pct_output_dir)
ensure_dir(output_dir)
ensure_dir(protein_distribution_dir)
ensure_dir(sample_distribution_dir)
ensure_dir(vo2_distribution_output_dir)
ensure_dir(pca_output_dir)
ensure_dir(bioimp_corr_output_dir)

# Initialize report file
report_file <- file.path(output_dir, paste0("prot_data_processing_report", filename_suffix, ".txt"))
if (file.exists(report_file)) file.remove(report_file)
append_report_line(report_file, paste0("Data Processing Report: subset_HCT_type = ", subset_HCT_type))
append_report_line(report_file, paste0("Run date: ", Sys.time()))
append_report_line(report_file, paste0("Low-variance cutoff: bottom ", low_variance_cutoff * 100, "%"))
append_report_line(report_file, "==================================================")


###
# Import Data
###
# Proteomics Data
data_df <- read.csv(file.path(output_folder, batch_correction_subfolder, proteomics_input_filename))
ptid_visit_colnames <- setdiff(colnames(data_df), feature_id_colname)
append_report_line(report_file, "")
append_report_line(report_file, paste0("Input data: ", nrow(data_df), " proteins x ", length(ptid_visit_colnames), " samples"))

# Patient Metadata, VO2 data
patient_metadata_df <- read.csv(file.path(input_folder, patient_metadata_filename))
vo2_data_df <- read.csv(file.path(output_folder, batch_correction_subfolder, vo2_data_filename))

# Patient Bioimpedance Data
bioimp_data_df <- read.xlsx(file.path(input_folder, bioimpedance_filename))

# Protein Annotation
protein_annotations_df <- read.csv(file.path(input_folder, olink_annotations_filename))


###
# Filter out patients without TRANSPLANT_TYPE info
###
transplant_filter_result <- remove_ptids_without_transplant_type(patient_metadata_df, data_df, ptid_visit_colnames, report_file)
patient_metadata_df <- transplant_filter_result$patient_metadata_df
data_df <- transplant_filter_result$data_df
ptid_visit_colnames <- transplant_filter_result$ptid_visit_colnames


###
# Format patient metadata - impute covariates
###
# Impute numerical and categorical covariates with median and mode, respectively
patient_metadata_df <- impute_covariates_patient_metadata(patient_metadata_df, report_file, covars_num, covars_cat)


###
# Format Bioimpedance Data and Join with Patient Metadata
###
patient_metadata_df <- join_patient_metadata_with_bioimpedance(patient_metadata_df, bioimp_data_df)


###
# Filter Low-Variance Proteins (full cohort, before any patient subsetting)
###
data_df <- filter_low_variance_proteins(data_df, ptid_visit_colnames, low_variance_cutoff, report_file, protein_distribution_dir, filename_suffix)


###
# Apply HCT_type Subsetting if Specified
###
if (subset_HCT_type == "ALLO" || subset_HCT_type == "AUTO") {
  # Get PTIDs matching the specified transplant type from patient metadata
  subset_ptids <- patient_metadata_df %>%
    dplyr::filter(TRANSPLANT_TYPE == subset_HCT_type) %>%
    dplyr::pull(PTID)

  # Filter columns of data_df to keep only matching PTID_visit samples + OlinkID
  keep_cols <- c(feature_id_colname, ptid_visit_colnames[
    sub("^X", "", sub("_.*$", "", ptid_visit_colnames)) %in% as.character(subset_ptids)
  ])
  data_df <- data_df[, keep_cols, drop = FALSE]

  message(paste0("Subsetted proteomics data to ", length(keep_cols) - 1,
                 " samples for HCT type: ", subset_HCT_type))
  append_report_line(report_file, "")
  append_report_line(report_file, paste0("--- HCT Type Subsetting: ", subset_HCT_type, " ---"))
  append_report_line(report_file, paste0("Samples retained: ", length(keep_cols) - 1))
  ptid_visit_colnames <- setdiff(colnames(data_df), feature_id_colname)
}
append_report_line(report_file, "")
append_report_line(report_file, paste0("Final dimensions for analysis: ", nrow(data_df), " proteins x ", ncol(data_df) - 1, " samples"))


###
# Plot NPX Distributions
###
avg_prot_abundance <- rowMeans(data_df[ , -which(names(data_df) == feature_id_colname)], na.rm = TRUE) # NAs not imputed yet, ignore them for now.
plot_distribution_histogram(
  data     = avg_prot_abundance,
  title    = "Protein Distribution, Post-Filtering and Batch Correction",
  xlab     = "Average NPX per Protein",
  filename = file.path(protein_distribution_dir, paste0("protein_post_filtering_distribution", filename_suffix, ".png")),
  text_size = font_size,
  fill_color = "#4A90E2",
  border_color = "#2E5C8A"
)

generate_sample_wise_feature_distribution_plots(
  data_df = data_df,
  output_dir = sample_distribution_dir,
  filename_stem = "protein_distribution_samples_all_timepoints",
  filename_suffix = filename_suffix,
  id_colname = feature_id_colname,
  y_axis_label = "NPX_mean",
  single_plot = TRUE
)

generate_sample_median_variability_plots(
  data_df = data_df,
  output_dir = sample_distribution_dir,
  filename_stem = "sample_median_variability_by_timepoint",
  filename_suffix = filename_suffix,
  id_colname = feature_id_colname,
  font_size = font_size,
  add_stats = TRUE
)


###
# Generate VO2 deltas
###
vo2_wide <- data_processing_generate_outcome_var_deltas(vo2_data_df)


###
# Generate Percent_Change_VO2
###
vo2_wide <- vo2_wide %>%
  mutate(
    pct_change_VO2peak_6m_Baseline = if_else(!is.na(VO2peak_Baseline) & !is.na(VO2peak_6m) & VO2peak_Baseline != 0, (VO2peak_6m - VO2peak_Baseline) / abs(VO2peak_Baseline) * 100, NA_real_),
    pct_change_VO2peak_12m_Baseline = if_else(!is.na(VO2peak_Baseline) & !is.na(VO2peak_12m) & VO2peak_Baseline != 0, (VO2peak_12m - VO2peak_Baseline) / abs(VO2peak_Baseline) * 100, NA_real_)
  )


###
# Create Processed Data
###
# Export processed proteomics data and patient metadata
# 1. Proteomics measurements (long format, no VO2 data; columns PTID, visit, OlinkID, NPX_mean)
prot_measurements <- data_df %>%
  tidyr::pivot_longer(
    cols = -!!sym(feature_id_colname),
    names_to = "sample_id",
    values_to = "NPX_mean"
  ) %>%
  tidyr::separate(sample_id, into = c("PTID", "visit"), sep = "_", remove = FALSE, extra = "merge", fill = "right") %>%
  mutate(PTID = sub("^X", "", PTID)) %>%  # Remove leading X if present
  mutate(PTID = as.integer(PTID)) %>%
  dplyr::select(PTID, visit, !!sym(feature_id_colname), NPX_mean)

write.csv(prot_measurements, 
        file.path(output_dir, paste0("proteomics_data_processed", filename_suffix, ".csv")), 
        row.names = FALSE)

# 2. Proteomics measurements (wide format, columns PTID, rows OlinkID, values NPX_mean)
data_df_wide_out <- data_df %>%
  dplyr::rename_with(~ sub("^X(\\d)", "\\1", .x), .cols = -!!sym(feature_id_colname))
write.csv(data_df_wide_out, 
        file.path(output_dir, paste0("proteomics_data_processed_wide", filename_suffix, ".csv")), 
        row.names = FALSE)

# 3. Patient metadata (one row per patient, includes VO2 timepoints)
# Create this BEFORE PCA plotting so it's available for sample_meta
# Append vo2_wide to patient_metadata_df
patient_metadata_df_final <- patient_metadata_df %>%
  dplyr::left_join(vo2_wide, by = "PTID") %>%
  dplyr::select(PTID, TRANSPLANT_TYPE, everything())

write.csv(patient_metadata_df_final,
          file.path(output_dir, paste0("proteomics_patient_metadata", filename_suffix, ".csv")),
          row.names = FALSE)


###
# Plot VO2peak Distributions Per Timepoint
###
plot_outcome_var_boxplot(patient_metadata_df_final, vo2_distribution_output_dir, filename_suffix, font_size, timepoint_order = timepoint_order, outcome_var = "VO2peak")
plot_outcome_var_trendlines(patient_metadata_df_final, vo2_distribution_output_dir, filename_suffix, font_size, outcome_var = "VO2peak")


###
# Plot VO2peak Violin Plots and Delta VO2peakHistograms (All, and ALLO- and AUTO-separated)
###
if (subset_HCT_type == "ALL") {
  # All patients together
  plot_outcome_var_violin_plot(patient_metadata_df_final, vo2_distribution_output_dir, filename_suffix, font_size, timepoint_order = timepoint_order, outcome_var = "VO2peak")
  plot_delta_outcome_var_histogram(patient_metadata_df_final, vo2_distribution_output_dir, filename_suffix, font_size, outcome_var = "VO2peak")

  # Patients with TRANSPLANT_TYPE == ALLO
  plot_outcome_var_violin_plot(
    patient_metadata_df_final %>% dplyr::filter(TRANSPLANT_TYPE == "ALLO"),
    vo2_distribution_output_dir, paste0(filename_suffix, "_ALLO"), font_size, timepoint_order = timepoint_order, plot_title = "Allogeneic", outcome_var = "VO2peak"
  )
  plot_delta_outcome_var_histogram(
    patient_metadata_df_final %>% dplyr::filter(TRANSPLANT_TYPE == "ALLO"),
    vo2_distribution_output_dir, paste0(filename_suffix, "_ALLO"), font_size, plot_title = "Allogeneic", outcome_var = "VO2peak"
  )

  # Patients with TRANSPLANT_TYPE == AUTO
  plot_outcome_var_violin_plot(
    patient_metadata_df_final %>% dplyr::filter(TRANSPLANT_TYPE == "AUTO"),
    vo2_distribution_output_dir, paste0(filename_suffix, "_AUTO"), font_size, timepoint_order = timepoint_order, plot_title = "Autologous", outcome_var = "VO2peak"
  )
  plot_delta_outcome_var_histogram(
    patient_metadata_df_final %>% dplyr::filter(TRANSPLANT_TYPE == "AUTO"),
    vo2_distribution_output_dir, paste0(filename_suffix, "_AUTO"), font_size, plot_title = "Autologous", outcome_var = "VO2peak"
  )

}
# Skip if violin plotting if subset_HCT_type != "ALL"


###
# Compare ALLO vs AUTO VO2peak distributions per timepoint
###
if (subset_HCT_type == "ALL") {
  plot_outcome_var_violin_allo_vs_auto(patient_metadata_df_final, vo2_distribution_output_dir, filename_suffix, font_size, timepoint_order = timepoint_order, outcome_var = "VO2peak")
  plot_outcome_var_violin_plot_one_timepoint(patient_metadata_df_final, vo2_distribution_output_dir, filename_suffix, font_size, timepoint = "Baseline", outcome_var = "VO2peak")
}


###
# Create a PCA plot of the proteomics data (post-filtering) to Check for Possible Outliers / Batch Effects
###
sample_meta <- build_data_for_pca(data_df, patient_metadata_df_final, id_colname=feature_id_colname)

# Run PCA plot generation:
# (1) Run PCA plot generation using all proteins:
run_all_pca_plot_generation(
  data_df = data_df,
  sample_meta = sample_meta,
  output_dir = pca_output_dir,
  filename_suffix = filename_suffix,
  subset_HCT_type = subset_HCT_type,
  font_size = font_size
)

# (2) Run PCA plot generation using only the top pca_var_top_pct% highest-variance proteins:
variance_by_protein <- data_df %>%
  dplyr::select(-!!sym(feature_id_colname)) %>%
  as.matrix() %>%
  rowVars(na.rm = TRUE)
variance_threshold <- quantile(variance_by_protein, probs = (100-pca_var_top_pct)/100, na.rm = TRUE)
data_df_high_variance <- data_df[variance_by_protein > variance_threshold, ]
run_all_pca_plot_generation(
  data_df = data_df_high_variance,
  sample_meta = sample_meta,
  output_dir = top_pct_output_dir,
  filename_suffix = paste0(filename_suffix, "_highVariance"),
  subset_HCT_type = subset_HCT_type,
  font_size = font_size
)


###
# Compare protein variance across timepoints
###
variance_by_timepoint <- data_df %>%
  pivot_longer(cols = -!!sym(feature_id_colname), names_to = "sample", values_to = "NPX") %>%
  separate(sample, into = c("PTID", "visit"), sep = "_", extra = "merge") %>%
  group_by(!!sym(feature_id_colname), visit) %>%
  summarise(variance = var(NPX, na.rm = TRUE), .groups = "drop") %>%
  group_by(visit) %>%
  summarise(
    mean_variance = mean(variance, na.rm = TRUE),
    median_variance = median(variance, na.rm = TRUE),
    sd_variance = sd(variance, na.rm = TRUE)
  )

write.csv(variance_by_timepoint,
          file.path(output_dir, paste0("Variance_By_Timepoint", filename_suffix, ".csv")),
          row.names = FALSE)

message("\n✓ Proteomic Data Processing complete\n")
append_report_line(report_file, "")
append_report_line(report_file, "==================================================")
append_report_line(report_file, "Data processing complete.")
append_report_line(report_file, "==================================================")


###
# Calculate Correlation of PBF and PSMM with VO2peak at Baseline
###
corr_vars <- list(
  list(var = "PBF",  label = "% Body Fat"),
  list(var = "PSMM", label = "% Skeletal Muscle Mass")
)

for (cv in corr_vars) {
  d_corr <- patient_metadata_df_final %>%
    dplyr::filter(!is.na(.data[[cv$var]]), !is.na(VO2peak_Baseline))

  if (nrow(d_corr) < 3) {
    message(cv$label, " vs VO2peak_Baseline: insufficient data (n=", nrow(d_corr), ")")
    next
  }

  ct <- cor.test(d_corr[[cv$var]], d_corr$VO2peak_Baseline, method = "pearson")
  r_val <- round(ct$estimate, 3)
  p_val <- formatC(ct$p.value, format = "e", digits = 2)
  annot_label <- paste0("r = ", r_val, "\np = ", p_val, "\nn = ", nrow(d_corr))

  p_corr <- ggplot(d_corr, aes(x = .data[[cv$var]], y = VO2peak_Baseline)) +
    geom_point(size = 2.5, alpha = 0.7, color = "#4A90E2") +
    geom_smooth(method = "lm", se = TRUE, color = "#2E5C8A", fill = "#4A90E2", alpha = 0.15) +
    annotate("text", x = Inf, y = Inf, label = annot_label,
             hjust = 1.1, vjust = 1.5, size = font_size * 0.2, color = "black") +
    labs(
      x = cv$label,
      y = expression("VO"[2]*"peak (mL/kg/min)")
    ) +
    theme_bw(base_size = font_size*0.7)

  out_file <- file.path(bioimp_corr_output_dir,
                        paste0("corr_", cv$var, "_vs_VO2peak_Baseline", filename_suffix, ".png"))
  ggsave(out_file, plot = p_corr, width = 6, height = 5, dpi = 600)
  message("Saved: ", basename(out_file))
}

}


# Only run when executed directly (not when sourced from another script)
if (sys.nframe() == 0L) {
  run_proteomics_data_processing(subset_HCT_type = "ALL", low_variance_cutoff = proteomics_low_variance_cutoff)
}