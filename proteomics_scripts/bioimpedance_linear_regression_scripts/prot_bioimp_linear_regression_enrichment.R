###
# Proteomics linear regression pathway enrichment
# - Continue from linear regression
# - Perform enrichment analysis only (no plotting)
###

library(dplyr)
library(openxlsx)
library(AnnotationDbi)
library(org.Hs.eg.db)
source(here::here("common_functions.R"))
# Seed set in common_functions.R via config.R

###
# To run code as a function:
run_bioimp_proteomics_linear_regression_enrichment <- function(subset_HCT_type, outcome_variable_list, use_t_statistic = FALSE) {
###


###
# Function
###
run_all_linear_regression_pathways_bioimp <- function(base_dir, filename_suffix, use_t_statistic = FALSE) {
  # Baseline PBF:
  if ("PBF" %in% outcome_variable_list) {
    run_pathway_enrichment(res_Baseline_PBF, "Baseline_PBF",  olink_mapped, base_dir, filename_suffix, use_t_statistic = use_t_statistic)
  }
  # Baseline PSMM:
  if ("PSMM" %in% outcome_variable_list) {
    run_pathway_enrichment(res_Baseline_PSMM, "Baseline_PSMM",  olink_mapped, base_dir, filename_suffix, use_t_statistic = use_t_statistic)
  }

  build_feature_summary(reg_list, file.path(base_dir, paste0("Prot_Bioimp_Baseline_Linear_Protein_Summary", filename_suffix, ".csv")))
}


###
# Values
###
input_folder <- get_input_folder()
olink_annotations_filename <- "olink_mapped.csv"
msigdb_cache_folder <- file.path(input_folder, "msigdb_annotations")

output_folder <- get_output_folder()
output_reg_results_folder <- "prot_bioimp_linear_regression"
output_script_folder <- "prot_bioimp_linear_regression_enrichment"

summary_list <<- list() # pathway summaries per contrast
reg_list <<- list() # regression tables per contrast


###
# Take into account HCT type subsets
###
if (subset_HCT_type == "ALLO" || subset_HCT_type == "AUTO") {
  # Adjust all relevant file paths and names for ALLO or AUTO subset
  output_reg_results_folder <- paste0(subset_HCT_type, "_", output_reg_results_folder)
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
ensure_dir(msigdb_cache_folder)


###
# Import Regression Results
###
# Read excel files for linear regression results
if ("PBF" %in% outcome_variable_list) {
  res_Baseline_PBF <- read.xlsx(file.path(output_folder, output_reg_results_folder, paste0("Linear_Regression_Results_Baseline_PBF", filename_suffix, ".xlsx")))
}
if ("PSMM" %in% outcome_variable_list) {
  res_Baseline_PSMM <- read.xlsx(file.path(output_folder, output_reg_results_folder, paste0("Linear_Regression_Results_Baseline_PSMM", filename_suffix, ".xlsx")))
}


###
# Import Olink-Protein Mapping
###
olink_mapped <- read.csv(file.path(input_folder, olink_annotations_filename))

olink_mapped <- olink_mapped %>%
  mutate(
    ENTREZID = as.character(ENTREZID),
    OlinkID  = as.character(OlinkID)
  )


###
# Load or download MSigDB sets
###
msig_h_raw    <- get_msigdb_set("H", NULL, msigdb_cache_folder)
msig_re_raw   <- get_msigdb_set("C2", "CP:REACTOME", msigdb_cache_folder)
msig_bp_raw   <- get_msigdb_set("C5", "BP", msigdb_cache_folder)
msig_kegg_raw <- get_msigdb_set("C2", "CP:KEGG_MEDICUS", msigdb_cache_folder)

hallmark_sets <<- msig_to_list(msig_h_raw)
react_sets    <<- msig_to_list(msig_re_raw)
gobp_sets     <<- msig_to_list(msig_bp_raw)
kegg_sets     <<- msig_to_list(msig_kegg_raw)

t2g_h    <<- make_t2g(msig_h_raw)
t2g_re   <<- make_t2g(msig_re_raw)
t2g_bp   <<- make_t2g(msig_bp_raw)
t2g_kegg <<- make_t2g(msig_kegg_raw)


###
# Run All Pathway Enrichments for Linear Regression Analysis
###
run_all_linear_regression_pathways_bioimp(output_dir, filename_suffix, use_t_statistic = use_t_statistic)

}


# Only run when executed directly (not when sourced from another script)
if (sys.nframe() == 0L) {
  run_bioimp_proteomics_linear_regression_enrichment("ALL", outcome_variable_list = outcome_variable_list)
}
