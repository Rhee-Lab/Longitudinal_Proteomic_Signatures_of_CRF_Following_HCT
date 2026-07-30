###############################################################################
# config.R — Central configuration for CREST HCT VO2 Study analysis
#
# All user-configurable parameters are defined here. Individual analysis
# scripts read from this file so that changing a value here propagates
# everywhere.
#
# Sourced automatically by common_functions.R; no need to source separately.
###############################################################################


###
# Reproducibility
###
random_seed <- 1234


###
# Directory layout (relative names only; resolved by helpers in common_functions.R)
###
input_dir_name  <- "input"
output_dir_name <- "output"


###
# HCT type subsets to analyse
# Valid values: "ALL", "ALLO", "AUTO" (or any combination)
###
choices_HCT_types <- c("ALL")


###
# Timepoints to consider and their order
###
timepoint_order <- c("Baseline", "6m", "12m")

###
# Covariates used in regression models
###
# Option 1:
covars <- c("AgeatHCT", "gender")
covars_num <- c("AgeatHCT")
covars_cat <- c("gender")
# Covariates to drop at Baseline for certain contrasts (leave empty if none):
covars_to_remove_at_baseline <- c()

# # Option 2:
# # Consider "AgeatHCT", "gender", "TRANSPLANT_TYPE" as covariates in linear regression models
# covars <- c("AgeatHCT", "gender", "TRANSPLANT_TYPE")
# covars_num <- c("AgeatHCT")
# covars_cat <- c("gender", "TRANSPLANT_TYPE")
# # Covariates to drop at Baseline for certain contrasts (leave empty if none):
# covars_to_remove_at_baseline <- c("TRANSPLANT_TYPE")


###
# Proteomics — Batch Correction & Filtering
###
prot_id_colname <- "OlinkID"
# LOD filtering: % of samples below LOD to trigger removal (0.50 = stringent, 0.75 = lenient)
proteomics_lod_filtering_cutoff <- 0.75

# Low-variance filtering: remove proteins below this variance percentile
proteomics_low_variance_cutoff <- 0.10

# PCA: keep top n% highest-variance proteins for a secondary PCA view
proteomics_pca_var_top_pct <- 10


###
# Proteomics - Linear Regression
###
proteomics_jobs_to_run_linear_regression <- c("Cross-Sectional", "Double_Delta")


###
# Metabolomics — Filtering
###
metab_id_colname <- "feature_label"
# Maximum fraction of missing values allowed per metabolite
metabolomics_missingness_cutoff <- 0.20

# Minimum max peak intensity across samples to keep a metabolite
metabolomics_low_max_peak_intensity_cutoff <- 500

# Low-variance filtering: remove metabolites below this variance percentile
metabolomics_low_variance_cutoff <- 0.05

# PCA: keep top n% highest-variance metabolites for a secondary PCA view
metabolomics_pca_var_top_pct <- 10


###
# Metabolomics - Linear Regression
###
metabolomics_jobs_to_run_linear_regression <- c("Cross-Sectional", "Double_Delta")
# Available job types: "Cross-Sectional", "Delta", "Pct_Change", "Lagged_Association", "Double_Delta"


###
# Plotting
###
font_size <- 25
font_style <- "Roboto"
run_for_metab <- FALSE

# GO:BP category colors for enrichment plots
gobp_palette <- c(
	"Immune / Inflammation"       = "#f47a00",
	"Metabolism"                  = "#007191",
	"Cell Cycle / Proliferation"  = "#009E73",
	"Cell Death / Stress"         = "#9c5177",
	"ECM / Adhesion"              = "#999999",
	"Morphogenesis / Remodeling"  = "#E69F00",
	"Signaling"                   = "#CC79A7",
	"Muscle"                      = "#62c8d3",
	"Other"                       = "#414141"
)

# Timepoint colors
timepoint_palette <- c("Baseline" = "#0072B2", "6m" = "#E69F00", "12m" = "#CC79A7")
blue_plot_color <- "#2166AC"
red_plot_color <- "#B2182B"


###
# Outcome Variable
###
# Proteomics Bioimpedance
outcome_variable_list <- c("PBF", "PSMM")