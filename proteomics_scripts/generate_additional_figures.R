run_generate_additional_figures <- function() {

######################################################################
# Protein-specific scatter plots for VO2peak vs NPX_mean
######################################################################

generate_scatter_plots <- function() {

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
source("common_functions.R")


###
# Values
###
proteins_of_interest <- c(
"IL6",
"LEP",
"CDHR2",
"FABP4",
"FGF21",
"ASGR1",
"ADAM12",
"ITGA11",
"OXT",
"CTSD",
"GDF15",
"DUSP13A",
"CDH2",
"IGSF3",
"KIT",
"CGREF1"
)


input_folder               <- get_input_folder()
output_folder              <- get_output_folder()
output_script_folder       <- "additional_figures"

data_processing_folder <- "proteomics_data_processing" # in output_folder
proteomics_input_filename  <- "proteomics_data_processed.csv"
proteomics_metadata_filename <- "proteomics_patient_metadata.csv"
olink_annotations_filename <- "olink_mapped.csv"

subset_HCT_type <- "ALL"


###
# Take into account HCT type subsets
###
if (subset_HCT_type %in% c("ALLO", "AUTO")) {
  data_processing_folder  <- paste0(subset_HCT_type, "_", data_processing_folder)
  proteomics_input_filename      <- gsub(".csv", paste0("_", subset_HCT_type, ".csv"), proteomics_input_filename)
  proteomics_metadata_filename   <- gsub(".csv", paste0("_", subset_HCT_type, ".csv"), proteomics_metadata_filename)
  output_script_folder           <- paste0(subset_HCT_type, "_", output_script_folder)
  filename_suffix                <- paste0("_", subset_HCT_type)
} else if (subset_HCT_type == "ALL") {
  filename_suffix <- ""
} else {
  stop("Invalid subset_HCT_type. Choose 'ALLO', 'AUTO', or 'ALL'.")
}


###
# Organize Directories
###
output_dir         <- file.path(output_folder, output_script_folder)
scatter_all_dir    <- file.path(output_dir, "Scatter_All_Timepoints")
scatter_tp_dir     <- file.path(output_dir, "Scatter_Per_Timepoint")

ensure_dir(output_dir)
ensure_dir(scatter_all_dir)
ensure_dir(scatter_tp_dir)


###
# Import Data
###
prot_measurements <- read.csv(
  file.path(output_folder, data_processing_folder, proteomics_input_filename)
)
patient_meta <- read.csv(
  file.path(output_folder, data_processing_folder, proteomics_metadata_filename)
)
olink_annotations <- read.csv(
  file.path(input_folder, olink_annotations_filename)
)


###
# Build visit-matched VO2 long table
# Maps each PTID × visit to the corresponding VO2peak value
###
vo2_long <- patient_meta %>%
  dplyr::select(PTID, VO2peak_Baseline, VO2peak_6m, VO2peak_12m) %>%
  tidyr::pivot_longer(
    cols      = c(VO2peak_Baseline, VO2peak_6m, VO2peak_12m),
    names_to  = "visit",
    values_to = "VO2peak"
  ) %>%
  dplyr::mutate(
    visit = dplyr::recode(visit,
      "VO2peak_Baseline" = "Baseline",
      "VO2peak_6m"       = "6m",
      "VO2peak_12m"      = "12m"
    )
  )


###
# Map gene symbols → OlinkIDs using olink_annotations
###
olink_id_map <- olink_annotations %>%
  dplyr::filter(SYMBOL %in% proteins_of_interest) %>%
  dplyr::select(OlinkID, SYMBOL) %>%
  dplyr::distinct()

missing_symbols <- setdiff(proteins_of_interest, olink_id_map$SYMBOL)
if (length(missing_symbols) > 0) {
  message("Warning: No OlinkID found for symbol(s): ", paste(missing_symbols, collapse = ", "))
}


###
# Merge NPX measurements with visit-matched VO2
###
plot_df_full <- prot_measurements %>%
  dplyr::filter(OlinkID %in% olink_id_map$OlinkID) %>%
  dplyr::left_join(olink_id_map, by = "OlinkID") %>%
  # Average NPX_mean across OlinkIDs for the same SYMBOL × PTID × visit
  # (handles symbols with multiple OlinkIDs, e.g. technical replicates)
  dplyr::group_by(SYMBOL, PTID, visit) %>%
  dplyr::summarise(NPX_mean = mean(NPX_mean, na.rm = TRUE), .groups = "drop") %>%
  dplyr::left_join(vo2_long, by = c("PTID", "visit")) %>%
  dplyr::filter(!is.na(VO2peak), !is.na(NPX_mean)) %>%
  dplyr::mutate(visit = factor(visit, levels = timepoint_order))


###
# Generate scatter plots
###
generate_feature_scatter_plot <- function(df_prot, prot, timepoint_order,
                                          timepoint_palette, scatter_all_dir,
                                          scatter_tp_dir, filename_suffix) {

  if (nrow(df_prot) == 0) {
    message("No data found for protein: ", prot, " — skipping.")
    return(invisible(NULL))
  }

  # ---------------------------------------------------------------
  # 1. Combined plot — all timepoints, points colored by timepoint,
  #    with a per-timepoint linear trend line and correlation labels
  # ---------------------------------------------------------------

  # Compute per-timepoint Pearson r and n for annotation
  cor_labels <- df_prot %>%
    dplyr::group_by(visit) %>%
    dplyr::summarise(
      r = cor(NPX_mean, VO2peak, use = "complete.obs", method = "pearson"),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(label = sprintf("%s: r = %.2f (n = %d)", visit, r, n))

  p_all <- ggplot(df_prot, aes(x = NPX_mean, y = VO2peak, color = visit)) +
    geom_point(size = 2.5, alpha = 0.75) +
    geom_smooth(method  = "lm", se = TRUE,
                aes(group = visit, fill = visit),
                linewidth = 0.8, alpha = 0.15) +
    geom_text(
      data    = cor_labels,
      mapping = aes(label = label, color = visit),
      x       = -Inf, y = Inf,
      hjust   = -0.05,
      vjust   = seq(1.5, by = 1.6, length.out = nrow(cor_labels)),
      size    = 3.5,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = timepoint_palette, name = "Timepoint") +
    scale_fill_manual(values  = timepoint_palette, guide = "none") +
    labs(
      title = paste0(prot, " \u2014 NPX vs VO<sub>2</sub>max (All Timepoints)"),
      x     = paste0(prot, " NPX (normalized protein abundance)"),
      y     = "VO<sub>2</sub>max (mL/kg/min)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title      = element_markdown(hjust = 0.5),
      axis.title.y    = element_markdown(),
      legend.position = "right"
    )

  ggsave(
    filename = file.path(scatter_all_dir,
                         paste0(prot, "_NPX_vs_VO2peak_AllTimepoints",
                                filename_suffix, ".png")),
    plot   = p_all,
    width  = 7, height = 5, dpi = 600
  )

  # ---------------------------------------------------------------
  # 2. Per-timepoint plots — one plot per timepoint
  # ---------------------------------------------------------------
  for (tp in timepoint_order) {
    df_tp <- df_prot %>% dplyr::filter(visit == tp)

    if (nrow(df_tp) < 2) {
      message("  Insufficient data for ", prot, " at ", tp, " — skipping.")
      next
    }

    r_val <- cor(df_tp$NPX_mean, df_tp$VO2peak, use = "complete.obs",
                 method = "pearson")
    n_val <- nrow(df_tp)

    tp_color <- timepoint_palette[[tp]]

    p_tp <- ggplot(df_tp, aes(x = NPX_mean, y = VO2peak)) +
      geom_point(size = 2.5, alpha = 0.8, color = tp_color) +
      geom_smooth(method  = "lm", se = TRUE,
                  color   = tp_color, fill = tp_color,
                  linewidth = 0.8, alpha = 0.15) +
      annotate("text",
               x = -Inf, y = Inf,
               label  = sprintf("r = %.2f\nn = %d", r_val, n_val),
               hjust  = -0.15, vjust = 1.4,
               size   = 4, color = "grey30") +
      labs(
        title = paste0(prot, " \u2014 NPX vs VO<sub>2</sub>max (", tp, ")"),
        x     = paste0(prot, " NPX (normalized)"),
        y     = "VO<sub>2</sub>max (mL/kg/min)"
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title   = element_markdown(hjust = 0.5),
        axis.title.y = element_markdown()
      )

    tp_subdir <- file.path(scatter_tp_dir, tp)
    ensure_dir(tp_subdir)

    ggsave(
      filename = file.path(tp_subdir,
                           paste0(prot, "_NPX_vs_VO2peak_", tp,
                                  filename_suffix, ".png")),
      plot   = p_tp,
      width  = 5.5, height = 5, dpi = 600
    )
  }

  message("Saved scatter plots for: ", prot)
}


###
# Loop over proteins and generate plots
###
for (prot in proteins_of_interest) {
  df_prot <- plot_df_full %>% dplyr::filter(SYMBOL == prot)
  generate_feature_scatter_plot(
    df_prot         = df_prot,
    prot            = prot,
    timepoint_order = timepoint_order,
    timepoint_palette = timepoint_palette,
    scatter_all_dir  = scatter_all_dir,
    scatter_tp_dir   = scatter_tp_dir,
    filename_suffix  = filename_suffix
  )
}

message("All scatter plots saved to: ", output_dir)
}

generate_scatter_plots()


######################################################################
# Generate Top VO2peak-correlating Proteins List
######################################################################
calculate_top_correlating_proteins <- function() {

  library(dplyr)
  library(tidyr)
  source("common_functions.R")

  ###
  # Values
  ###
  input_folder               <- get_input_folder()
  output_folder              <- get_output_folder()
  output_script_folder       <- "additional_figures"

  data_processing_folder       <- "proteomics_data_processing" # in output_folder
  prot_df_filename           <- "proteomics_data_processed.csv"
  patient_metadata_filename  <- "proteomics_patient_metadata.csv"
  prot_feature_metadata_filename <- "olink_mapped.csv"

  n_top_labels     <- 30
  font_size        <- 14


  ###
  # Organize Directories
  ###
  output_dir <- file.path(output_folder, output_script_folder, "Top_Correlating_Proteins")
  ensure_dir(output_dir)


  ###
  # Import Data
  ###
  prot_long    <- read.csv(file.path(output_folder, data_processing_folder, prot_df_filename))
  patient_meta <- read.csv(file.path(output_folder, data_processing_folder, patient_metadata_filename))
  prot_feature_meta <- read.csv(file.path(input_folder, prot_feature_metadata_filename))


  ###
  # Build visit-matched VO2 long table
  ###
  vo2_long <- patient_meta %>%
    dplyr::select(PTID, VO2peak_Baseline, VO2peak_6m, VO2peak_12m) %>%
    tidyr::pivot_longer(
      cols      = c(VO2peak_Baseline, VO2peak_6m, VO2peak_12m),
      names_to  = "visit",
      values_to = "VO2peak"
    ) %>%
    dplyr::mutate(
      visit = dplyr::recode(visit,
        "VO2peak_Baseline" = "Baseline",
        "VO2peak_6m"       = "6m",
        "VO2peak_12m"      = "12m"
      )
    )

  # Build sample_meta: sample_id, visit, VO2peak
  sample_meta <- vo2_long %>%
    dplyr::mutate(sample_id = paste0(PTID, "_", visit)) %>%
    dplyr::select(sample_id, visit, VO2peak) %>%
    dplyr::filter(!is.na(VO2peak))


  ###
  # Pivot proteomics to wide format: OlinkID x samples (PTID_visit)
  ###
  prot_wide <- prot_long %>%
    dplyr::mutate(sample_id = paste0(PTID, "_", visit)) %>%
    dplyr::select(OlinkID, sample_id, NPX_mean) %>%
    tidyr::pivot_wider(
      names_from  = sample_id,
      values_from = NPX_mean,
      values_fn   = mean
    )


  ###
  # Cross-Sectional: protein abundance vs VO2peak (all timepoints + per-timepoint)
  # Always generated regardless of jobs_to_run_linear_regression
  ###

  # All timepoints combined
  generate_correlation_plot(
    df           = prot_wide,
    sample_meta  = sample_meta,
    y_var        = "VO2peak",
    filename     = file.path(output_dir, "protein_cor_VO2peak_AllTimepoints.png"),
    id_name      = "OlinkID",
    feature_meta = prot_feature_meta,
    label_col    = "SYMBOL",
    visit_filter = NULL,
    n_top_labels = n_top_labels,
    text_size    = font_size,
    y_label      = "Pearson r,\nProtein Abundance vs. VO\u2082peak"
  )

  # Per-timepoint
  for (tp in timepoint_order) {
    generate_correlation_plot(
      df           = prot_wide,
      sample_meta  = sample_meta,
      y_var        = "VO2peak",
      filename     = file.path(output_dir, paste0("protein_cor_VO2peak_", tp, ".png")),
      id_name      = "OlinkID",
      feature_meta = prot_feature_meta,
      label_col    = "SYMBOL",
      visit_filter = tp,
      n_top_labels = n_top_labels,
      text_size    = font_size,
      y_label      = paste0("Pearson r,\n", tp, " Protein Abundance vs. VO\u2082peak")
    )
  }


  ###
  # Build VO2 delta table for Delta, Lagged_Association, and Double_Delta jobs
  ###
  vo2_delta <- patient_meta %>%
    dplyr::select(PTID, delta_VO2peak_6m_Baseline, delta_VO2peak_12m_Baseline) %>%
    dplyr::filter(!is.na(delta_VO2peak_6m_Baseline) | !is.na(delta_VO2peak_12m_Baseline))

  # Helper: pivot NPX deltas to wide format for a given timepoint delta
  # Returns OlinkID x PTID columns with delta NPX values
  compute_npx_deltas_wide <- function(prot_long, baseline_visit = "Baseline", followup_visit) {
    prot_bl <- prot_long %>%
      dplyr::filter(visit == baseline_visit) %>%
      dplyr::select(PTID, OlinkID, NPX_Baseline = NPX_mean)
    prot_fu <- prot_long %>%
      dplyr::filter(visit == followup_visit) %>%
      dplyr::select(PTID, OlinkID, NPX_followup = NPX_mean)
    delta_long <- dplyr::inner_join(prot_bl, prot_fu, by = c("PTID", "OlinkID")) %>%
      dplyr::mutate(delta_NPX = NPX_followup - NPX_Baseline) %>%
      dplyr::mutate(sample_id = as.character(PTID)) %>%
      dplyr::select(OlinkID, sample_id, delta_NPX)
    delta_wide <- delta_long %>%
      tidyr::pivot_wider(names_from = sample_id, values_from = delta_NPX, values_fn = mean)
    return(delta_wide)
  }


  ###
  # Double_Delta: delta protein abundance vs delta VO2 at 6m and 12m
  ###
  if ("Double_Delta" %in% proteomics_jobs_to_run_linear_regression) {
    dd_dir <- file.path(output_dir, "Double_Delta")
    ensure_dir(dd_dir)

    dd_configs <- list(
      list(followup = "6m",  vo2_col = "delta_VO2peak_6m_Baseline"),
      list(followup = "12m", vo2_col = "delta_VO2peak_12m_Baseline")
    )

    for (cfg in dd_configs) {
      delta_prot_wide <- compute_npx_deltas_wide(prot_long, followup_visit = cfg$followup)
      dd_sample_meta <- vo2_delta %>%
        dplyr::mutate(sample_id = as.character(PTID)) %>%
        dplyr::select(sample_id, !!rlang::sym(cfg$vo2_col)) %>%
        dplyr::filter(!is.na(!!rlang::sym(cfg$vo2_col)))

      generate_correlation_plot(
        df           = delta_prot_wide,
        sample_meta  = dd_sample_meta,
        y_var        = cfg$vo2_col,
        filename     = file.path(dd_dir, paste0("protein_cor_DoubleDelta_", cfg$followup, ".png")),
        id_name      = "OlinkID",
        feature_meta = prot_feature_meta,
        label_col    = "SYMBOL",
        visit_filter = NULL,
        n_top_labels = n_top_labels,
        text_size    = font_size,
        y_label      = paste0("Pearson r, Delta Protein Abundance\nvs. Delta VO\u2082peak (", cfg$followup, "-Baseline)")
      )
    }
    message("Double_Delta correlation plots saved.")
  }


  ###
  # Delta: protein abundance at follow-up vs delta VO2 at 6m and 12m
  ###
  if ("Delta" %in% proteomics_jobs_to_run_linear_regression) {
    delta_dir <- file.path(output_dir, "Delta")
    ensure_dir(delta_dir)

    delta_configs <- list(
      list(visit_filter = "6m",  vo2_col = "delta_VO2peak_6m_Baseline"),
      list(visit_filter = "12m", vo2_col = "delta_VO2peak_12m_Baseline")
    )

    for (cfg in delta_configs) {
      delta_sample_meta <- vo2_delta %>%
        dplyr::mutate(
          sample_id = paste0(PTID, "_", cfg$visit_filter),
          visit = cfg$visit_filter
        ) %>%
        dplyr::select(sample_id, visit, !!rlang::sym(cfg$vo2_col)) %>%
        dplyr::filter(!is.na(!!rlang::sym(cfg$vo2_col)))

      generate_correlation_plot(
        df           = prot_wide,
        sample_meta  = delta_sample_meta,
        y_var        = cfg$vo2_col,
        filename     = file.path(delta_dir, paste0("protein_cor_Delta_", cfg$visit_filter, ".png")),
        id_name      = "OlinkID",
        feature_meta = prot_feature_meta,
        label_col    = "SYMBOL",
        visit_filter = cfg$visit_filter,
        n_top_labels = n_top_labels,
        text_size    = font_size,
        y_label      = paste0("Pearson r, ", cfg$visit_filter,
                              " Protein Abundance\nvs. Delta VO\u2082peak (", cfg$visit_filter, "-Baseline)")
      )
    }
    message("Delta correlation plots saved.")
  }


  ###
  # Lagged_Association: protein at earlier timepoint vs delta VO2 at later timepoint
  #   Baseline protein -> delta VO2 at 6m
  #   Baseline protein -> delta VO2 at 12m
  #   6m protein       -> delta VO2 at 12m
  ###
  if ("Lagged_Association" %in% proteomics_jobs_to_run_linear_regression) {
    lagged_dir <- file.path(output_dir, "Lagged_Association")
    ensure_dir(lagged_dir)

    lagged_configs <- list(
      list(prot_visit = "Baseline", vo2_col = "delta_VO2peak_6m_Baseline",  label = "Baseline_to_delta6m"),
      list(prot_visit = "Baseline", vo2_col = "delta_VO2peak_12m_Baseline", label = "Baseline_to_delta12m"),
      list(prot_visit = "6m",       vo2_col = "delta_VO2peak_12m_Baseline", label = "6m_to_delta12m")
    )

    for (cfg in lagged_configs) {
      lagged_sample_meta <- vo2_delta %>%
        dplyr::mutate(
          sample_id = paste0(PTID, "_", cfg$prot_visit),
          visit = cfg$prot_visit
        ) %>%
        dplyr::select(sample_id, visit, !!rlang::sym(cfg$vo2_col)) %>%
        dplyr::filter(!is.na(!!rlang::sym(cfg$vo2_col)))

      # Build descriptive y-axis label for lagged association
      vo2_delta_label <- if (cfg$vo2_col == "delta_VO2peak_6m_Baseline") {
        "Delta VO\u2082peak (6m-Baseline)"
      } else {
        "Delta VO\u2082peak (12m-Baseline)"
      }

      generate_correlation_plot(
        df           = prot_wide,
        sample_meta  = lagged_sample_meta,
        y_var        = cfg$vo2_col,
        filename     = file.path(lagged_dir, paste0("protein_cor_Lagged_", cfg$label, ".png")),
        id_name      = "OlinkID",
        feature_meta = prot_feature_meta,
        label_col    = "SYMBOL",
        visit_filter = cfg$prot_visit,
        n_top_labels = n_top_labels,
        text_size    = font_size,
        y_label      = paste0("Pearson r, ", cfg$prot_visit,
                              " Protein Abundance\nvs. ", vo2_delta_label)
      )
    }
    message("Lagged_Association correlation plots saved.")
  }

  message("All top-correlating protein plots saved to: ", output_dir)
}

calculate_top_correlating_proteins()
}


# Only run when executed directly (not when sourced from another script)
if (sys.nframe() == 0L) {
  run_generate_additional_figures()
}


######################################################################
# Ad hoc figure generation: for HCT-adjusted run, generate heatmap for the 16 concordant proteins
######################################################################
# Run this to generate a heatmap for the 16 concordant proteins from non-HCT-adjusted analysis, using the HCT-adjusted analysis results
# - requires having run the HCT-adjusted analysis first, to access those files
generate_ad_hoc_topn_heatmap <- function(job_folders, result_folder, excel_result_filename_prefix, top_prots, output_folder, output_subfolder) {
  library(openxlsx)
  source("common_functions.R")

  # Build results_list: one data frame per timepoint, read from Excel
  results_list <- setNames(
    lapply(job_folders, function(folder) {
      read.xlsx(
        file.path(result_folder, folder,
                  paste0(excel_result_filename_prefix, folder, ".xlsx")),
        sheet = "all"
      )
    }),
    job_folders
  )

  out_dir <- file.path(output_folder, output_subfolder)
  ensure_dir(out_dir)

  plot_protein_time_series_heatmap(
    results_list = results_list,
    proteins     = top_prots,
    timepoints   = job_folders,
    outputdir    = out_dir,
    filename     = "Protein_Time_Series_Heatmap_TopN_AdHoc.png",
    pval_col     = "FDR",
    metric_col   = "beta",
    sig_cutoff   = 0.05,
    sort_type    = "input"
  )
}

run_ad_hoc_figure_gen <- FALSE
if (run_ad_hoc_figure_gen) {
  output_folder <- "output"
  output_subfolder <- "additional_figures" # place heatmap figure here

  # plot in this order on the y-axis
  top_prots <- c("ITGA11",
                 "KIT",
                 "DUSP13A",
                 "OXT",
                 "FGF21",
                 "LEP",
                 "IL6",
                 "GDF15",
                 "IGSF3",
                 "CGREF1",
                 "CTSD",
                 "CDHR2",
                 "ADAM12",
                 "FABP4",
                 "CDH2",
                 "ASGR1")
  result_folder <- file.path(output_folder, "proteomics_linear_regression_enrichment")
  excel_result_filename_prefix <- "Linear_Protein_Results_Annotated_" # ie "Linear_Protein_Results_Annotated_Baseline.xlsx"

  generate_ad_hoc_topn_heatmap(timepoint_order, result_folder, excel_result_filename_prefix, top_prots, output_folder, output_subfolder)
}