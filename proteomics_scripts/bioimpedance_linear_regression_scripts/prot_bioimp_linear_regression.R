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
# Seed set in common_functions.R via config.R


###
# To run code as a function:
run_biomp_proteomics_linear_regression <- function(subset_HCT_type, covars, covars_to_remove_at_baseline, outcome_variable_list) {
###

###
# Functions
###
run_contrast <- function(df, patient_meta, label, outcome_var, covars, visit_filter = NULL) {
  
  message("Running regression for: ", label)
  
  # Filter by visit if specified
  if (!is.null(visit_filter)) {
    d <- df %>% filter(visit == visit_filter)
  } else {
    d <- df
  }
  
  if (nrow(d) == 0) return(NULL)
  
  # Join patient metadata and prepare data before grouping by OlinkID
  d <- d %>% 
    left_join(patient_meta, by = "PTID") %>%
    mutate(
      gender = factor(gender),
      Diabetes = factor(Diabetes),
      Hypertension = factor(Hypertension)
    ) %>%
    filter(!is.na(.data[[outcome_var]]))  # Filter out missing outcomes
  
  # Convert any remaining character covariates to factors
  for (cv in covars) {
    if (cv %in% colnames(d) && is.character(d[[cv]])) d[[cv]] <- factor(d[[cv]])
  }
  
  if (nrow(d) == 0) return(NULL)
  
  # Pre-filter covariates: remove those with >50% missing across all data
  valid_covars <- covars[sapply(d[covars], function(x) mean(is.na(x)) <= 0.5)]
  
  # Impute covariate values
  for (cv in valid_covars) {
    if (is.factor(d[[cv]])) {
      d[[cv]] <- impute_factor_mode(d[[cv]])
    } else {
      d[[cv]] <- impute_numeric_mean(d[[cv]])
    }
  }
  
  # Group by OlinkID and fit models
  result <- d %>%
    group_by(OlinkID) %>%
    group_modify(~fit_one_protein(.x, outcome_var, valid_covars, return_statistic = TRUE)) %>%
    ungroup() %>%
    mutate(Protein = OlinkID, contrast = label)
  
  # Add FDR correction (Benjamini-Hochberg)
  if (nrow(result) > 0) {
    result <- result %>%
      mutate(FDR = p.adjust(p, method = "BH"))
  }

  # Filter by ascending p
  result <- result %>% arrange(p)
  
  return(result)
}


###
# Values
###
input_folder <- get_input_folder()
output_folder <- get_output_folder()
output_script_folder <- "prot_bioimp_linear_regression"

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
subfolder_pvalue_histograms <- file.path(output_dir, "PValue_Histograms")
ensure_dir(subfolder_pvalue_histograms)
subfolder_beta_distributions <- file.path(output_dir, "Beta_Distributions")
ensure_dir(subfolder_beta_distributions)
output_dir_residuals <- file.path(output_dir, "Residual_Diagnostics")
ensure_dir(output_dir_residuals)


###
# Import Data
###
# Processed protein NPX measurements and patient metadata in proteomics dataset
prot_measurements <- read.csv(file.path(output_folder, data_processing_folder, proteomics_input_filename))
patient_meta <- read.csv(file.path(output_folder, data_processing_folder, patient_metadata_filename))


###
# Initialise residual analysis report
###
# QQ plots per contrast + summary statistics appended to report.txt.
# See common_functions.R for the underlying functions.
residual_report_file <- file.path(output_dir, paste0("report", filename_suffix, ".txt"))
init_residual_report(residual_report_file,
                     script_label    = "prot_bioimp_linear_regression.R (protein NPX -> bioimpedance)",
                     covars          = covars,
                     subset_HCT_type = subset_HCT_type)

residual_summaries <- list()


###
# Fit linear regression models per contrast
###
# Cross-sectional contrast, "Baseline" only
# For each bioimpedance variable in outcome_variable_list (config.R)
for (var in outcome_variable_list) {
  # Fit model
  covars_for_baseline <- setdiff(covars, covars_to_remove_at_baseline)
  res_Baseline <- run_contrast(prot_measurements, patient_meta, "Baseline",  var, covars_for_baseline, visit_filter = "Baseline")

  # Save results in excel
  write.xlsx(res_Baseline, file.path(output_dir, paste0("Linear_Regression_Results_Baseline_", var, filename_suffix, ".xlsx")), sheetName = "Baseline", rowNames = FALSE)

  # Generate p-value histogram
  pvalue_histogram(res_Baseline$p, subfolder_pvalue_histograms, paste0("Baseline_", var), filename_suffix=filename_suffix)

  # Plot beta coefficient distribution
  plot_beta_distribution(res_Baseline, subfolder_beta_distributions, paste0("Baseline_", var), filename_suffix=filename_suffix)

  # Residual analysis: QQ plots + summary statistics in report.txt.
  # impute_covars = TRUE mirrors run_contrast() above, which imputes covariates
  # and fits via fit_one_protein() (min-imputation of missing NPX).
  residual_summaries[[paste0("Baseline_", var)]] <- run_contrast_residual_diagnostics(
    prot_measurements, patient_meta, paste0("Baseline_", var), var, covars_for_baseline,
    output_dir = output_dir_residuals, report_file = residual_report_file,
    id_col = "OlinkID", measurement_col = "NPX_mean", visit_filter = "Baseline",
    impute_covars = TRUE, filename_suffix = filename_suffix)
}


###
# Combined "Summary - All Contrasts" table
###
append_residual_summary_table(residual_report_file, residual_summaries)
message("Residual analysis report written: ", residual_report_file)


}


# Only run when executed directly (not when sourced from another script)
if (sys.nframe() == 0L) {
  run_biomp_proteomics_linear_regression(subset_HCT_type = "ALL", covars = c("AgeatHCT", "gender"), covars_to_remove_at_baseline = NULL, outcome_variable_list = outcome_variable_list)
}