###
# Run all proteomics scripts
###
# Batch Correction and Data Processing
# source(here::here("proteomics_scripts", "proteomics_batch_correction_chain_bridging.R"))
source(here::here("proteomics_scripts", "proteomics_batch_correction_ComBat.R"))
source(here::here("proteomics_scripts", "proteomics_data_processing.R"))
# Linear Regression
source(here::here("proteomics_scripts", "VO2peak_linear_regression_scripts", "proteomics_linear_regression.R"))
source(here::here("proteomics_scripts", "VO2peak_linear_regression_scripts", "proteomics_linear_regression_enrichment.R"))
source(here::here("proteomics_scripts", "VO2peak_linear_regression_scripts", "proteomics_linear_regression_plots.R"))
# Extreme Subsets, Linear Regression
source(here::here("proteomics_scripts", "VO2peak_linear_regression_scripts", "proteomics_top_vs_bottom_diff_exp.R"))
source(here::here("proteomics_scripts", "VO2peak_linear_regression_scripts", "proteomics_top_vs_bottom_diff_exp_enrichment.R"))
source(here::here("proteomics_scripts", "VO2peak_linear_regression_scripts", "proteomics_top_vs_bottom_diff_exp_plots.R"))

# Additional Figures (Proteomics, Metabolomics)
source(here::here("proteomics_scripts", "generate_additional_figures.R"))
# Bioimpedance analyses
source(here::here("proteomics_scripts", "bioimpedance_linear_regression_scripts", "prot_bioimp_linear_regression.R"))
source(here::here("proteomics_scripts", "bioimpedance_linear_regression_scripts", "prot_bioimp_linear_regression_enrichment.R"))
source(here::here("proteomics_scripts", "bioimpedance_linear_regression_scripts", "prot_bioimp_linear_regression_plots.R"))


#####################################################################
###
# Save Job Info
###
# Save a job text file in output folder with values used
job_info <- list(
  choices_HCT_types = choices_HCT_types,
  covars = covars,
  covars_to_remove_at_baseline = covars_to_remove_at_baseline,
  proteomics_lod_filtering_cutoff = proteomics_lod_filtering_cutoff,
  proteomics_low_variance_cutoff = proteomics_low_variance_cutoff
)
job_info_text <- paste(
  sapply(names(job_info), function(nm) {
    paste0(nm, ": ", paste(job_info[[nm]], collapse = ", "))
  }),
  collapse = "\n"
)
writeLines(job_info_text, con = here::here("output", "proteomics_analysis_job_info.txt"))


###
# Run batch correction
###
cat("Begin running proteomics batch correction script...\n")
cat("==================================================\n")
print(Sys.time())
run_proteomics_batch_correction_ComBat(lod_filtering_cutoff = proteomics_lod_filtering_cutoff)
cat("Finished batch correction.\n")
print(Sys.time())


###
# Run data processing
###
for (job_type in choices_HCT_types) {
  cat(paste0("Running proteomics data processing script for HCT type: ", job_type, "\n"))

  run_proteomics_data_processing(subset_HCT_type = job_type, low_variance_cutoff = proteomics_low_variance_cutoff)

  cat(paste0("Finished data processing for HCT type: ", job_type, "\n"))
  print(Sys.time())
}


###
# Run linear regression
###
for (job_type in choices_HCT_types) {
  cat(paste0("Running proteomics linear regression script for HCT type: ", job_type, "\n"))

  if (job_type != "ALL") {
    covars_to_remove_at_baseline_temp <- NULL
  } else {
    covars_to_remove_at_baseline_temp <- covars_to_remove_at_baseline
  }

  run_proteomics_linear_regression(subset_HCT_type = job_type, covars = covars, covars_to_remove_at_baseline = covars_to_remove_at_baseline_temp, jobs_to_run = proteomics_jobs_to_run_linear_regression)

  cat(paste0("Finished linear regression for HCT type: ", job_type, "\n"))
  print(Sys.time())
}


###
# Run linear regression enrichment
###
for (job_type in choices_HCT_types) {
  cat(paste0("Running proteomics linear regression enrichment script for HCT type: ", job_type, "\n"))

  run_proteomics_linear_regression_enrichment(subset_HCT_type = job_type, jobs_to_run = proteomics_jobs_to_run_linear_regression)

  cat(paste0("Finished proteomics linear regression enrichment script for HCT type: ", job_type, "\n"))
  print(Sys.time())
}


###
# Run linear regression plots
###
for (job_type in choices_HCT_types) {
  cat(paste0("Running proteomics linear regression plots script for HCT type: ", job_type, "\n"))

  run_proteomics_linear_regression_plots(subset_HCT_type = job_type, jobs_to_run = proteomics_jobs_to_run_linear_regression)

  cat(paste0("Finished proteomics linear regression plots script for HCT type: ", job_type, "\n"))
  print(Sys.time())
}


###
# Run extreme subsets differential expression
###
for (job_type in choices_HCT_types) {
  cat(paste0("Running extreme subset differential expression for HCT type: ", job_type, "\n"))

  if (job_type != "ALL") {
    covars_to_remove_at_baseline_temp <- NULL
  } else {
    covars_to_remove_at_baseline_temp <- covars_to_remove_at_baseline
  }

  run_top_vs_bottom_diff_exp(covars, covars_to_remove_at_baseline_temp, subset_HCT_type = job_type)
  run_proteomics_top_vs_bottom_diff_exp_enrichment(subset_HCT_type = job_type)
  run_proteomics_top_vs_bottom_diff_exp_plots(subset_HCT_type = job_type)

  cat(paste0("Finished extreme subset analysis for HCT type: ", job_type, "\n"))
  print(Sys.time())
}


###
# Run bioimpedance linear regression analyses
###
for (job_type in choices_HCT_types) {
    cat(paste0("Running bioimpedance proteomics linear regression script for HCT type: ", job_type, "\n"))

    if (job_type != "ALL") {
      covars_to_remove_at_baseline_temp <- NULL
    } else {
      covars_to_remove_at_baseline_temp <- covars_to_remove_at_baseline
    }

    run_biomp_proteomics_linear_regression(subset_HCT_type = job_type, covars = covars, covars_to_remove_at_baseline = covars_to_remove_at_baseline_temp, outcome_variable_list = outcome_variable_list)

    cat(paste0("Finished bioimpedance proteomics linear regression for HCT type: ", job_type, "\n"))
    print(Sys.time())
  }

# Bioimpedance, linear regression enrichment
for (job_type in choices_HCT_types) {
  cat(paste0("Running bioimpedance proteomics linear regression enrichment script for HCT type: ", job_type, "\n"))

  run_bioimp_proteomics_linear_regression_enrichment(subset_HCT_type = job_type, outcome_variable_list = outcome_variable_list)

  cat(paste0("Finished bioimpedance proteomics linear regression enrichment script for HCT type: ", job_type, "\n"))
  print(Sys.time())
}

# Bioimpedance, linear regression plotting
for (job_type in choices_HCT_types) {
  cat(paste0("Running bioimpedance proteomics linear regression plots script for HCT type: ", job_type, "\n"))

  run_bioimp_proteomics_linear_regression_plots(subset_HCT_type = job_type, outcome_variable_list = outcome_variable_list)

  cat(paste0("Finished proteomics linear regression plots script for HCT type: ", job_type, "\n"))
  print(Sys.time())
}


###
# Run additional figures script
###
if ("ALL" %in% choices_HCT_types) {
  run_generate_additional_figures()
}