###
# Proteomics Extreme Subsets Differential Expression Plots
###
library(dplyr)
library(openxlsx)
library(ggplot2)
source(here::here("common_functions.R"))
# Seed set in common_functions.R via config.R

###
# To run code as a function:
run_proteomics_top_vs_bottom_diff_exp_plots <- function(subset_HCT_type = "ALL") {
###


###
# Values
###
output_folder <- get_output_folder()
# sub-folder of output folder
enrichment_output_folder <- "proteomics_top_vs_bottom_enrichment"
# sub-folders in enrichment output folder
contrasts <- timepoint_order
filename_suffix <- ""


output_script_folder <- "proteomics_top_vs_bottom_plots"


###
# Take into account HCT type subsets
###
if (subset_HCT_type == "ALLO" || subset_HCT_type == "AUTO") {
  enrichment_output_folder <- paste0(subset_HCT_type, "_", enrichment_output_folder)
  output_script_folder <- paste0(subset_HCT_type, "_", output_script_folder)
  filename_suffix <- paste0("_", subset_HCT_type)
} else if (subset_HCT_type == "ALL") {
  filename_suffix <- ""
} else {
  stop("Invalid subset_HCT_type specified. Please choose 'ALLO', 'AUTO', or 'ALL'.")
}

timeseries_subfolder <- "Time_Series_Heatmaps"
volcano_subfolder <- "Volcano_Plots"
barplot_subfolder <- "Bar_Plots"

approach_contrast_combos <- list(
    "Top_vs_Bottom" = timepoint_order
)


###
# Organize Directories
###
output_dir <- file.path(output_folder, output_script_folder)
ensure_dir(output_dir)
ensure_dir(file.path(output_dir, timeseries_subfolder))
ensure_dir(file.path(output_dir, volcano_subfolder))
ensure_dir(file.path(output_dir, barplot_subfolder))

# Time series sub-type subfolders (one per data source / generation method)
ts_proteins_subfolder      <- file.path(timeseries_subfolder, "Proteins")
ts_gsea_gobp_subfolder     <- file.path(timeseries_subfolder, "GSEA_GOBP")
ts_ora_gobp_subfolder      <- file.path(timeseries_subfolder, "ORA_GOBP")
ts_gsea_reactome_subfolder <- file.path(timeseries_subfolder, "GSEA_Reactome")
ts_ora_reactome_subfolder  <- file.path(timeseries_subfolder, "ORA_Reactome")
ensure_dir(file.path(output_dir, ts_proteins_subfolder))
ensure_dir(file.path(output_dir, ts_gsea_gobp_subfolder))
ensure_dir(file.path(output_dir, ts_ora_gobp_subfolder))
ensure_dir(file.path(output_dir, ts_gsea_reactome_subfolder))
ensure_dir(file.path(output_dir, ts_ora_reactome_subfolder))


###
# Import Enrichment Results
###
# Store all linear regression results
lin_regr_results_list <- list()
# Store all enrichment results
enrichment_results_list <- list()
# Store additional enrichment types extracted from Summary_Pathways
enrichment_gsea_gobp_results_list <- list()
enrichment_ora_gobp_results_list <- list()
enrichment_gsea_reactome_results_list <- list()
enrichment_ora_reactome_results_list <- list()

for (contrast in contrasts) {
  contrast_folder <- file.path(output_folder, enrichment_output_folder, contrast)
  # All protein results
  lin_regr_results_list[[contrast]] <- read.xlsx(file.path(contrast_folder, paste0("Linear_Protein_Results_Annotated_", contrast, filename_suffix, ".xlsx")), sheet = "all")
  # All enrichment results
  enrichment_results_list[[contrast]] <- list(
    Enrich = read.xlsx(file.path(contrast_folder, paste0("Summary_Pathways_", contrast, filename_suffix, ".xlsx")), sheet = "all")
  )
}


###
# Clean up pathway names in enrichment results
###
# Normalize all pathway names and propagate cleaned names to both
# `pathway` and `Pathway` columns so downstream functions are consistent
for (contrast in contrasts) {
  enrich_df <- enrichment_results_list[[contrast]]$Enrich
  enrich_df <- enrich_df %>%
    mutate(
      Pathway = mapply(normalize_pathway_name, DB, Pathway),
      pathway = Pathway
    )
  enrichment_results_list[[contrast]]$Enrich <- enrich_df
}


###
# Extract additional enrichment types from Summary_Pathways
###
for (contrast in contrasts) {
  enrich_all <- enrichment_results_list[[contrast]]$Enrich
  enrichment_gsea_gobp_results_list[[contrast]]     <- list(GSEA_GOBP     = extract_gsea_enrich(enrich_all, "GOBP"))
  enrichment_gsea_reactome_results_list[[contrast]] <- list(GSEA_Reactome = extract_gsea_enrich(enrich_all, "Reactome"))
  enrichment_ora_gobp_results_list[[contrast]]      <- list(ORA_GOBP      = extract_ora_enrich(enrich_all, "GOBP"))
  enrichment_ora_reactome_results_list[[contrast]]  <- list(ORA_Reactome  = extract_ora_enrich(enrich_all, "Reactome"))
}


###
# GSEA Bar Plots
###
for (contrast in contrasts) {
  plot_gsea_combined(enrichment_gsea_gobp_results_list[[contrast]]$GSEA_GOBP,
                     file.path(output_dir, barplot_subfolder),
                     contrast,
                     filename_suffix = filename_suffix,
                     font_size = font_size*1.1)
}


###
# Volcano Plots with Top Points Labeled
###
for (contrast in contrasts) {
  generate_volcano_plot(
    data = lin_regr_results_list[[contrast]],
    pval_col = "FDR",
    beta_col = "logFC",
    contrast = contrast,
    p_thresh = 0.05,
    beta_thresh = 0.0,
    output_dir = file.path(output_dir, volcano_subfolder),
    feature_type = "Protein",
    label_col = "SYMBOL",
    filename_suffix = paste0(filename_suffix, "_top_vs_bottom"),
    n_label_limit = 30,
    font_size = font_size
  )
}

# Create versions with -log10(raw p-value) on y-axis
for (contrast in contrasts) {
  generate_volcano_plot(
    data = lin_regr_results_list[[contrast]],
    pval_col = "p",
    beta_col = "logFC",
    contrast = contrast,
    p_thresh = 0.1,
    beta_thresh = 0.0,
    output_dir = file.path(output_dir, volcano_subfolder),
    feature_type = "Protein",
    label_col = "SYMBOL",
    filename_suffix = paste0(filename_suffix, "_top_vs_bottom_raw_pval"),
    n_label_limit = 30,
    font_size = font_size
  )
}


###
# Time Series Heatmap, Individual Proteins
###
top_n_list <- c(20, 50)

for (approach in names(approach_contrast_combos)) {
  tps <- approach_contrast_combos[[approach]]
  run_protein_time_series_heatmap_generation(
    results_list = lin_regr_results_list,
    tps_for_focus = tps,
    top_n_list = top_n_list,
    outputdir = file.path(output_dir, ts_proteins_subfolder),
    filename_suffix = filename_suffix,
    timepoints_label = approach,
    metric_col = "logFC",
    pval_col = "FDR",
    pval_cutoff = 0.05
  )
}


###
# Annotation Time Series Heatmaps — Top n by timpoint in a contrast
###
annotation_sets <- list(
  GSEA_GOBP      = enrichment_gsea_gobp_results_list,
  ORA_GOBP       = enrichment_ora_gobp_results_list,
  GSEA_Reactome  = enrichment_gsea_reactome_results_list,
  ORA_Reactome   = enrichment_ora_reactome_results_list
)

annotation_subdirs <- list(
  GSEA_GOBP      = ts_gsea_gobp_subfolder,
  ORA_GOBP       = ts_ora_gobp_subfolder,
  GSEA_Reactome  = ts_gsea_reactome_subfolder,
  ORA_Reactome   = ts_ora_reactome_subfolder
)

# For each focus timepoint, select the top 20 lowest-padj pathways from that
# timepoint's enrichment results and plot them across all main timepoints.
# This mirrors the tp_focus approach used for individual protein heatmaps.
n_top_annotations <- 20

run_tp_focus_annotation_time_series_heatmap_generation(
  annotation_sets = annotation_sets,
  annotation_subdirs = annotation_subdirs,
  tps_for_focus = approach_contrast_combos$Top_vs_Bottom,
  timepoints = approach_contrast_combos$Top_vs_Bottom,
  n_top_annotations = n_top_annotations,
  output_dir = output_dir,
  filename_suffix = filename_suffix
)


###
# Annotation Time Series Heatmaps, Top 20 Overall in a Contrast
###
for (enrich_name in names(annotation_sets)) {
  enrich_list <- annotation_sets[[enrich_name]]

  for (approach in names(approach_contrast_combos)) {
    plot_annotation_time_series_heatmap(enrich_list, approach_contrast_combos[[approach]],
                                        file.path(output_dir, annotation_subdirs[[enrich_name]],
                                        paste0(approach, "_Time_Series_Heatmap_", enrich_name, filename_suffix, ".png")),
                                        enrichment_type = enrich_name,
                                        font_size = 12)
  }
}

}


# Only run when executed directly (not when sourced from another script)
if (sys.nframe() == 0L) {
  run_proteomics_top_vs_bottom_diff_exp_plots(subset_HCT_type = "ALL")
}