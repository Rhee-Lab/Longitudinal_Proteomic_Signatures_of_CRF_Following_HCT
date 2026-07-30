############################################################
## PROTEOMICS → VO2 LINEAR REGRESSION ANALYSIS
## - Linear regression per protein
## - Plot p-value histograms and beta coefficient distributions
## - Save results as Excel files
############################################################

library(dplyr)
library(tidyr)
library(openxlsx)
source(here::here("common_functions.R"))
source(here::here("config.R"))
# Seed set in common_functions.R via config.R


###
# To run code as a function:
run_proteomics_linear_regression <- function(subset_HCT_type, covars, covars_to_remove_at_baseline, jobs_to_run) {
###


###
# Functions
###
calculate_npx_deltas <- function(prot_measurements) {
  # Pivot NPX data to wide format by visit
  prot_wide <- prot_measurements %>%
    pivot_wider(
      id_cols = c(PTID, OlinkID),
      names_from = visit,
      values_from = NPX_mean,
      names_prefix = "NPX_"
    )
  
  # Calculate NPX deltas from Baseline
  prot_deltas <- prot_wide %>%
    mutate(
      delta_NPX_6m_Baseline = NPX_6m - NPX_Baseline,
      delta_NPX_12m_Baseline = NPX_12m - NPX_Baseline
    ) %>%
    dplyr::select(PTID, OlinkID, delta_NPX_6m_Baseline, delta_NPX_12m_Baseline)
  
  return(prot_deltas)
}


###
# Values
###
input_folder <- get_input_folder()
output_folder <- get_output_folder()
output_script_folder <- "proteomics_linear_regression"

data_processing_folder <- "proteomics_data_processing" # in output_folder
proteomics_input_filename <- "proteomics_data_processed.csv"
patient_metadata_filename <- "proteomics_patient_metadata.csv"


###
# Take into account HCT type subsets
###
if (subset_HCT_type == "ALLO" || subset_HCT_type == "AUTO") {
  # Adjust all relevant file paths and names for ALLO or AUTO subset
  data_processing_folder <- paste0(subset_HCT_type, "_", data_processing_folder)
  proteomics_input_filename <- gsub(".csv", paste0("_", subset_HCT_type, ".csv"), proteomics_input_filename)
  patient_metadata_filename <- gsub(".csv", paste0("_", subset_HCT_type, ".csv"), patient_metadata_filename)
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
ensure_dir(output_dir)


###
# Import Data
###
# Processed protein NPX measurements and patient metadata in proteomics dataset
prot_measurements <- read.csv(file.path(output_folder, data_processing_folder, proteomics_input_filename))
patient_meta <- read.csv(file.path(output_folder, data_processing_folder, patient_metadata_filename))


###
# Calculate NPX Deltas
###
# Calculate change in NPX from Baseline to follow-up timepoints
npx_deltas <- calculate_npx_deltas(prot_measurements)


###
# Fit linear regression models per contrast
###
# Cross-sectional contrasts, "Baseline", "6m", "12m":
if ("Cross-Sectional" %in% jobs_to_run) {
  covars_for_baseline <- setdiff(covars, covars_to_remove_at_baseline)
  res_Baseline <- run_regression_contrast(prot_measurements, patient_meta, "Baseline", "VO2peak_Baseline", covars_for_baseline,
                                          id_col = "OlinkID", measurement_col = "NPX_mean",
                                          visit_filter = "Baseline", impute_covars = FALSE, return_statistic = TRUE)
  res_6m  <- run_regression_contrast(prot_measurements, patient_meta, "6m", "VO2peak_6m", covars,
                                     id_col = "OlinkID", measurement_col = "NPX_mean",
                                     visit_filter = "6m", impute_covars = FALSE, return_statistic = TRUE)
  res_12m <- run_regression_contrast(prot_measurements, patient_meta, "12m", "VO2peak_12m", covars,
                                     id_col = "OlinkID", measurement_col = "NPX_mean",
                                     visit_filter = "12m", impute_covars = FALSE, return_statistic = TRUE)
}

# Delta contrasts, "Delta_6m", "Delta_12m"
if ("Delta" %in% jobs_to_run) {
  res_d6  <- run_regression_contrast(prot_measurements, patient_meta, "Delta_6m", "delta_VO2peak_6m_Baseline", covars,
                                     id_col = "OlinkID", measurement_col = "NPX_mean",
                                     visit_filter = "6m", impute_covars = FALSE, return_statistic = TRUE)
  res_d12 <- run_regression_contrast(prot_measurements, patient_meta, "Delta_12m", "delta_VO2peak_12m_Baseline", covars,
                                     id_col = "OlinkID", measurement_col = "NPX_mean",
                                     visit_filter = "12m", impute_covars = FALSE, return_statistic = TRUE)
}

# Pct contrasts (like Delta but using percent change in VO2peak), "Pct_Change_6m", "Pct_Change_12m"
if ("Pct_Change" %in% jobs_to_run) {
  res_pct_d6  <- run_regression_contrast(prot_measurements, patient_meta, "Pct_Change_6m", "pct_change_VO2peak_6m_Baseline", covars,
                                         id_col = "OlinkID", measurement_col = "NPX_mean",
                                         visit_filter = "6m", impute_covars = FALSE, return_statistic = TRUE)
  res_pct_d12 <- run_regression_contrast(prot_measurements, patient_meta, "Pct_Change_12m", "pct_change_VO2peak_12m_Baseline", covars,
                                         id_col = "OlinkID", measurement_col = "NPX_mean",
                                         visit_filter = "12m", impute_covars = FALSE, return_statistic = TRUE)
}

# Lagged Association contrasts, "Baseline_Delta_6m", "Baseline_Delta_12m", "6m_Delta_12m"
if ("Lagged_Association" %in% jobs_to_run) {
  res_bl_d6  <- run_regression_contrast(prot_measurements, patient_meta, "Baseline_Delta_6m", "delta_VO2peak_6m_Baseline", covars,
                                        id_col = "OlinkID", measurement_col = "NPX_mean",
                                        visit_filter = "Baseline", impute_covars = FALSE, return_statistic = TRUE)
  res_bl_d12 <- run_regression_contrast(prot_measurements, patient_meta, "Baseline_Delta_12m", "delta_VO2peak_12m_Baseline", covars,
                                        id_col = "OlinkID", measurement_col = "NPX_mean",
                                        visit_filter = "Baseline", impute_covars = FALSE, return_statistic = TRUE)
  res_6m_d12 <- run_regression_contrast(prot_measurements, patient_meta, "6m_Delta_12m", "delta_VO2peak_12m_Baseline", covars,
                                        id_col = "OlinkID", measurement_col = "NPX_mean",
                                        visit_filter = "6m", impute_covars = FALSE, return_statistic = TRUE)
}

# Double Delta contrasts, "Double_Delta_6m", "Double_Delta_12m"
if ("Double_Delta" %in% jobs_to_run) {
  res_npx_d6  <- run_regression_contrast_delta(npx_deltas, patient_meta, "Double_Delta_6m", "delta_VO2peak_6m_Baseline", covars,
                                               id_col = "OlinkID", delta_col = "delta_NPX_6m_Baseline",
                                               return_statistic = TRUE)
  res_npx_d12 <- run_regression_contrast_delta(npx_deltas, patient_meta, "Double_Delta_12m", "delta_VO2peak_12m_Baseline", covars,
                                               id_col = "OlinkID", delta_col = "delta_NPX_12m_Baseline",
                                               return_statistic = TRUE)
}


###
# Save results per contrast as excel file
###
# Cross-sectional analysis results
if ("Cross-Sectional" %in% jobs_to_run) {
  write.xlsx(res_Baseline, file.path(output_dir, paste0("Linear_Regression_Results_Baseline", filename_suffix, ".xlsx")), sheetName = "Baseline", rowNames = FALSE)
  write.xlsx(res_6m, file.path(output_dir, paste0("Linear_Regression_Results_6m", filename_suffix, ".xlsx")), sheetName = "6m",        rowNames = FALSE)
  write.xlsx(res_12m, file.path(output_dir, paste0("Linear_Regression_Results_12m", filename_suffix, ".xlsx")), sheetName = "12m",       rowNames = FALSE)
}

# Delta analysis results
if ("Delta" %in% jobs_to_run) {
  write.xlsx(res_d6, file.path(output_dir, paste0("Linear_Regression_Results_Delta_6m", filename_suffix, ".xlsx")), sheetName = "Delta_6m",  rowNames = FALSE)
  write.xlsx(res_d12, file.path(output_dir, paste0("Linear_Regression_Results_Delta_12m", filename_suffix, ".xlsx")), sheetName = "Delta_12m", rowNames = FALSE)
}

# Pct change analysis results
if ("Pct_Change" %in% jobs_to_run) {
  write.xlsx(res_pct_d6, file.path(output_dir, paste0("Linear_Regression_Results_Pct_Change_6m", filename_suffix, ".xlsx")), sheetName = "Pct_Change_6m", rowNames = FALSE)
  write.xlsx(res_pct_d12, file.path(output_dir, paste0("Linear_Regression_Results_Pct_Change_12m", filename_suffix, ".xlsx")), sheetName = "Pct_Change_12m", rowNames = FALSE)
}

# Lagged association analysis results
if ("Lagged_Association" %in% jobs_to_run) {
  write.xlsx(res_bl_d6, file.path(output_dir, paste0("Linear_Regression_Results_Baseline_Delta_6m", filename_suffix, ".xlsx")), sheetName = "Baseline_Delta_6m", rowNames = FALSE)
  write.xlsx(res_bl_d12, file.path(output_dir, paste0("Linear_Regression_Results_Baseline_Delta_12m", filename_suffix, ".xlsx")), sheetName = "Baseline_Delta_12m", rowNames = FALSE)
  write.xlsx(res_6m_d12, file.path(output_dir, paste0("Linear_Regression_Results_6m_Delta_12m", filename_suffix, ".xlsx")), sheetName = "6m_Delta_12m", rowNames = FALSE)
}

# Double Delta analysis results
if ("Double_Delta" %in% jobs_to_run) {
  write.xlsx(res_npx_d6, file.path(output_dir, paste0("Linear_Regression_Results_Double_Delta_6m", filename_suffix, ".xlsx")), sheetName = "Delta_NPX_6m", rowNames = FALSE)
  write.xlsx(res_npx_d12, file.path(output_dir, paste0("Linear_Regression_Results_Double_Delta_12m", filename_suffix, ".xlsx")), sheetName = "Delta_NPX_12m", rowNames = FALSE)
}


###
# Generate p-value Histograms
###
# Look for histograms with high values at low p-values
# Flat histograms indicates likely no signal
# Red line = signal.cutoff
# Blue line = signal.cutoff adjsuted for multiple testing (Bonferroni)
subfolder_pvalue_histograms <- file.path(output_dir, "PValue_Histograms")
ensure_dir(subfolder_pvalue_histograms)

if ("Cross-Sectional" %in% jobs_to_run) {
  pvalue_histogram(res_Baseline$p, subfolder_pvalue_histograms, "Baseline", filename_suffix=filename_suffix)
  pvalue_histogram(res_6m$p, subfolder_pvalue_histograms, "6m", filename_suffix=filename_suffix)
  pvalue_histogram(res_12m$p, subfolder_pvalue_histograms, "12m", filename_suffix=filename_suffix)
}

if ("Delta" %in% jobs_to_run) {
  pvalue_histogram(res_d6$p, subfolder_pvalue_histograms, "Delta_6m", filename_suffix=filename_suffix)
  pvalue_histogram(res_d12$p, subfolder_pvalue_histograms, "Delta_12m", filename_suffix=filename_suffix)
}

if ("Pct_Change" %in% jobs_to_run) {
  pvalue_histogram(res_pct_d6$p, subfolder_pvalue_histograms, "Pct_Change_6m", filename_suffix=filename_suffix)
  pvalue_histogram(res_pct_d12$p, subfolder_pvalue_histograms, "Pct_Change_12m", filename_suffix=filename_suffix)
}

if ("Baseline_Delta" %in% jobs_to_run) {
  pvalue_histogram(res_bl_d6$p, subfolder_pvalue_histograms, "Baseline_Delta_6m", filename_suffix=filename_suffix)
  pvalue_histogram(res_bl_d12$p, subfolder_pvalue_histograms, "Baseline_Delta_12m", filename_suffix=filename_suffix)
}

if ("6m_Delta_12m" %in% jobs_to_run) {
  pvalue_histogram(res_6m_d12$p, subfolder_pvalue_histograms, "6m_Delta_12m", filename_suffix=filename_suffix)
}

if ("Double_Delta" %in% jobs_to_run) {
  pvalue_histogram(res_npx_d6$p, subfolder_pvalue_histograms, "Double_Delta_6m", filename_suffix=filename_suffix)
  pvalue_histogram(res_npx_d12$p, subfolder_pvalue_histograms, "Double_Delta_12m", filename_suffix=filename_suffix)
}


###
# Plot beta coefficient distributions
###
subfolder_beta_distributions <- file.path(output_dir, "Beta_Distributions")
ensure_dir(subfolder_beta_distributions)

if ("Cross-Sectional" %in% jobs_to_run) {
  plot_beta_distribution(res_Baseline, subfolder_beta_distributions, "Baseline", filename_suffix=filename_suffix)
  plot_beta_distribution(res_6m, subfolder_beta_distributions, "6m", filename_suffix=filename_suffix)
  plot_beta_distribution(res_12m, subfolder_beta_distributions, "12m", filename_suffix=filename_suffix)
}

if ("Delta" %in% jobs_to_run) {
  plot_beta_distribution(res_d6, subfolder_beta_distributions, "Delta_6m", filename_suffix=filename_suffix)
  plot_beta_distribution(res_d12, subfolder_beta_distributions, "Delta_12m", filename_suffix=filename_suffix)
}

if ("Pct_Change" %in% jobs_to_run) {
  plot_beta_distribution(res_pct_d6, subfolder_beta_distributions, "Pct_Change_6m", filename_suffix=filename_suffix)
  plot_beta_distribution(res_pct_d12, subfolder_beta_distributions, "Pct_Change_12m", filename_suffix=filename_suffix)
}

if ("Baseline_Delta" %in% jobs_to_run) {
  plot_beta_distribution(res_bl_d6, subfolder_beta_distributions, "Baseline_Delta_6m", filename_suffix=filename_suffix)
  plot_beta_distribution(res_bl_d12, subfolder_beta_distributions, "Baseline_Delta_12m", filename_suffix=filename_suffix)
}

if ("6m_Delta_12m" %in% jobs_to_run) {
  plot_beta_distribution(res_6m_d12, subfolder_beta_distributions, "6m_Delta_12m", filename_suffix=filename_suffix)
}

if ("Double_Delta" %in% jobs_to_run) {
  plot_beta_distribution(res_npx_d6, subfolder_beta_distributions, "Double_Delta_6m", filename_suffix=filename_suffix)
  plot_beta_distribution(res_npx_d12, subfolder_beta_distributions, "Double_Delta_12m", filename_suffix=filename_suffix)
}


###
# Generate Venn Diagrams of Sig. Proteins Between Contrasts
###
# Default: use FDR<0.05 to define sig proteins
output_dir_venn <- file.path(output_dir, "Venn_Diagrams")
ensure_dir(output_dir_venn)

# Load Olink protein annotations for Excel summaries
olink_mapped_venn <- read.csv(file.path(input_folder, "olink_mapped.csv")) %>%
  mutate(ENTREZID = as.character(ENTREZID), OlinkID = as.character(OlinkID))

# Cross-sectional contrasts, FDR<0.05
if ("Cross-Sectional" %in% jobs_to_run) {
  cross_sectional_res_list <- list(
  "Baseline" = res_Baseline,
  "6m"       = res_6m,
  "12m"      = res_12m
  )

  # Cross-sectional view, any beta direction, FDR<0.05
  for (dir in c("any", "positive", "negative")) {
    generate_sig_venn_diagram(
      res_list = cross_sectional_res_list,
      title = NULL,
      output_dir = output_dir_venn,
      filename_suffix = filename_suffix,
      beta_dir = dir,
      protein_mapping = olink_mapped_venn,
      font_size = font_size
    )
    }
  
  # Generate a table image with the proteins significant across the 3 contrasts,
  # with a column of the protein descriptions,
  # and columns for adjusted p-value and beta coefficient (at baseline)
  get_sig_symbols_cs <- function(res) {
    res %>%
      filter(FDR < 0.05) %>%
      left_join(olink_mapped_venn %>% dplyr::select(OlinkID, SYMBOL) %>% distinct(), by = "OlinkID") %>%
      filter(!is.na(SYMBOL), SYMBOL != "") %>%
      group_by(SYMBOL) %>%
      arrange(p, .by_group = TRUE) %>%
      dplyr::filter(dplyr::row_number() == 1L) %>%
      ungroup()
  }

  sig_bl_sym  <- get_sig_symbols_cs(res_Baseline)
  sig_6m_sym  <- get_sig_symbols_cs(res_6m)
  sig_12m_sym <- get_sig_symbols_cs(res_12m)

  common_symbols_cs <- Reduce(intersect, list(sig_bl_sym$SYMBOL, sig_6m_sym$SYMBOL, sig_12m_sym$SYMBOL))

  if (length(common_symbols_cs) > 0) {
    table_df_cs <- sig_bl_sym %>%
      filter(SYMBOL %in% common_symbols_cs) %>%
      left_join(olink_mapped_venn %>% dplyr::select(OlinkID, ENTREZID) %>% distinct(), by = "OlinkID") %>%
      mutate(ENTREZID = as.character(ENTREZID))

    descs_cs <- get_gene_descriptions(table_df_cs$ENTREZID)
    beta_col_nm <- "\u03b2 (Baseline)"
    table_df_cs <- table_df_cs %>%
      mutate(Description = ifelse(!is.na(ENTREZID) & ENTREZID %in% names(descs_cs), descs_cs[ENTREZID], "")) %>%
      mutate(Description = ifelse(is.na(Description), "", Description)) %>%
      transmute(
        Protein              = SYMBOL,
        Description          = Description,
        !!beta_col_nm        := round(beta, 1),
        `p-adj (Baseline)`   = formatC(FDR, format = "e", digits = 1)
      )

    # Order rows to match heatmap that sorts by decreasing baseline beta
    table_df_cs <- table_df_cs %>%
      arrange(desc(!!sym(beta_col_nm)))

    # Build per-cell text-color matrix: β column colored red/blue by sign
    beta_col_idx <- which(colnames(table_df_cs) == beta_col_nm)
    fg_matrix <- matrix("black", nrow = nrow(table_df_cs), ncol = ncol(table_df_cs))
    fg_matrix[, beta_col_idx] <- ifelse(table_df_cs[[beta_col_nm]] > 0, "#B2182B", "#2166AC")

    # Build alternating row background matrix (white / light gray)
    row_bg <- rep(c("white", "gray94"), length.out = nrow(table_df_cs))
    fill_matrix <- matrix(row_bg, nrow = nrow(table_df_cs), ncol = ncol(table_df_cs))

    tbl_theme_cs <- gridExtra::ttheme_minimal(
      base_size = 11,
      core    = list(
        bg_params = list(fill = fill_matrix, col = NA),
        fg_params = list(col = fg_matrix, hjust = 0, x = 0.05)
      ),
      colhead = list(fg_params = list(hjust = 0, x = 0.05))
    )

    tbl_grob_cs <- gridExtra::tableGrob(table_df_cs, rows = NULL, theme = tbl_theme_cs)
    nat_w <- as.numeric(grid::convertWidth(tbl_grob_cs$widths, "mm", valueOnly = TRUE))
    tbl_grob_cs$widths <- unit(nat_w, "null")
    table_img_file_cs <- file.path(output_dir_venn, paste0("Table_Sig_All_CrossSectional_Contrasts", filename_suffix, ".png"))
    h_cs <- max(2, nrow(table_df_cs) * 0.35 + 1)
    grDevices::png(table_img_file_cs, width = 8, height = h_cs, units = "in", res = 300, bg = "white")
    grid::grid.draw(tbl_grob_cs)
    grDevices::dev.off()
    message("Table image saved: ", basename(table_img_file_cs))
  } else {
    message("No proteins significant (FDR < 0.05) in all 3 cross-sectional contrasts.")
  }
}

# Delta contrasts, FDR<0.05
if ("Delta" %in% jobs_to_run) {
  delta_res_list <- list(
    "Delta_6m" = res_d6,
    "Delta_12m" = res_d12
  )
  for (dir in c("any", "positive", "negative")) {
    generate_sig_venn_diagram(
      res_list = delta_res_list,
      title = NULL,
      output_dir = output_dir_venn,
      filename_suffix = filename_suffix,
      beta_dir = dir,
      protein_mapping = olink_mapped_venn,
      font_size = font_size
    )
  }
}

# Pct change contrasts, FDR<0.05
if ("Pct_Change" %in% jobs_to_run) {
  pct_res_list <- list(
    "Pct_Change_6m" = res_pct_d6,
    "Pct_Change_12m" = res_pct_d12
  )
  for (dir in c("any", "positive", "negative")) {
    generate_sig_venn_diagram(
      res_list = pct_res_list,
      title = NULL,
      output_dir = output_dir_venn,
      filename_suffix = filename_suffix,
      beta_dir = dir,
      protein_mapping = olink_mapped_venn,
      font_size = font_size
    )
  }
}

# Lagged association contrasts, FDR<0.05
if ("Lagged_Association" %in% jobs_to_run) {
  lagged_res_list <- list(
    "Baseline_Delta_6m" = res_bl_d6,
    "Baseline_Delta_12m" = res_bl_d12,
    "6m_Delta_12m" = res_6m_d12
  )
  for (dir in c("any", "positive", "negative")) {
    generate_sig_venn_diagram(
      res_list = lagged_res_list,
      title = NULL,
      output_dir = output_dir_venn,
      filename_suffix = filename_suffix,
      beta_dir = dir,
      protein_mapping = olink_mapped_venn,
      font_size = font_size
    )
  }
}

# Double delta contrasts, choose raw p-value<0.05
# Justified with p-value histogram left-enrichment
if ("Double_Delta" %in% jobs_to_run) {
  double_delta_res_list <- list(
    "Double_Delta_6m" = res_npx_d6,
    "Double_Delta_12m" = res_npx_d12
  )
  for (dir in c("any", "positive", "negative")) {
    generate_sig_venn_diagram(
      res_list = double_delta_res_list,
      title = NULL,
      output_dir = output_dir_venn,
      filename_suffix = filename_suffix,
      beta_dir = dir,
      protein_mapping = olink_mapped_venn,
      pval_col = "p",
      pval_cutoff = 0.05,
      font_size = font_size
    )
  }
}


}


# Only run when executed directly (not when sourced from another script)
if (sys.nframe() == 0L) {
  run_proteomics_linear_regression(subset_HCT_type = "ALL", covars = c("AgeatHCT", "gender"), covars_to_remove_at_baseline = NULL, jobs_to_run = proteomics_jobs_to_run_linear_regression)
}