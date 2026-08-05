############################################################
## Protoemics --> Extreme Patient Subsets Differential Expression Analysis
## - split patient population into quarters based on VO2peak, per timepoint
## - differential expression analysis between top and bottom 25% of patients
## - use limma
## - Plot p-value histograms and log2 fold-change distributions
## - Save results as Excel files
############################################################

library(dplyr)
library(tidyr)
library(openxlsx)
library(limma)
library(ggplot2)
source(here::here("common_functions.R"))
# Seed set in common_functions.R via config.R


###
# To run code as a function:
run_top_vs_bottom_diff_exp <- function(covars, covars_to_remove_at_baseline, subset_HCT_type = "ALL") {
###
# covars <- c("AgeatHCT", "gender", "TRANSPLANT_TYPE")
# # ^ For this contrast, TRANSPLANT_TYPE should be adjusted,
# # since the distribution of transplant types may differ
# # between the top and bottom patient subsets.
# # However, at baseline, TRANSPLANT_TYPE should not be adjusted.
# covars_to_remove_at_baseline <- c("TRANSPLANT_TYPE")


###
# Functions
###
plot_VO2peak_distributions <- function(patient_meta, top_patients_list,
                                      bottom_patients_list, tp, output_dir,
                                      VO2peak_col_name, patient_ID_col_name,
                                      font_size = 18, y_range=c(0,48),
                                      colors=c("steelblue", "indianred")) {
    VO2peak_data <- patient_meta %>%
        filter(.data[[patient_ID_col_name]] %in% c(top_patients_list, bottom_patients_list)) %>%
        dplyr::select(all_of(patient_ID_col_name), all_of(VO2peak_col_name)) %>%
        mutate(Group = ifelse(.data[[patient_ID_col_name]] %in% top_patients_list, "Top 25%", "Bottom 25%"))

    n_labels <- VO2peak_data %>%
        group_by(Group) %>%
        summarise(n = n(), y = max(.data[[VO2peak_col_name]]), .groups = "drop") %>%
        mutate(label = paste0("n=", n))

    ggplot(VO2peak_data, aes(x = Group, y = .data[[VO2peak_col_name]], fill = Group)) +
        # Outlier points are hidden because geom_jitter already draws every sample.
        geom_boxplot(outlier.shape = NA, fill = colors) +
        geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
        geom_text(data = n_labels, aes(x = Group, y = y, label = label), vjust = -0.5, size = font_size / 3, inherit.aes = FALSE) +
        labs(title = "",
            x = "Patient Group",
                y = vo2_label(paste0("VO2peak at ", tp))) +
        coord_cartesian(ylim = y_range) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size = font_size),
        axis.text.y = element_text(color = "black", size = font_size),
        axis.title.x = element_text(color = "black", size = font_size),
        axis.title.y = element_text(color = "black", size = font_size),
        legend.position = "none"
        )
    ggsave(filename = file.path(output_dir, paste0("VO2peak_Distributions_", tp, ".png")), width = 6, height = 6, dpi=600)
}

subset_patients <- function(patient_meta, tp,
                            patient_ID_col_name, outcome_var_prefix,
                            num_top_outliers_to_remove,
                            num_bottom_outliers_to_remove, top_grp_name, bottom_grp_name) {
    VO2peak_col_name <- paste0(outcome_var_prefix, tp)

    # Prep patients list ranked by VO2peak at tp.
    # Distinct on patient ID avoids duplicate IDs reappearing after outlier trimming.
    ranked_patient_meta <- patient_meta %>%
        filter(!is.na(.data[[VO2peak_col_name]])) %>%
        arrange(.data[[VO2peak_col_name]]) %>%
        distinct(.data[[patient_ID_col_name]], .keep_all = TRUE)

    # Before dividing, remove outlier patients (need to run with 0s first to determine outliers based on distribution)
    num_top_outliers <- num_top_outliers_to_remove[[tp]]
    num_bottom_outliers <- num_bottom_outliers_to_remove[[tp]]
    if (is.null(num_top_outliers) || is.na(num_top_outliers)) num_top_outliers <- 0L
    if (is.null(num_bottom_outliers) || is.na(num_bottom_outliers)) num_bottom_outliers <- 0L
    num_top_outliers <- as.integer(max(0L, num_top_outliers))
    num_bottom_outliers <- as.integer(max(0L, num_bottom_outliers))

    n_ranked <- nrow(ranked_patient_meta)
    start_idx <- num_bottom_outliers + 1L
    end_idx <- n_ranked - num_top_outliers

    if (start_idx > end_idx) {
        stop(paste0("Outlier removal removed all patients for ", tp,
                    ". Check num_top_outliers_to_remove / num_bottom_outliers_to_remove."))
    }

    patients_ranked <- ranked_patient_meta[[patient_ID_col_name]][start_idx:end_idx]

    # Determine top and bottom patient lists
    n_quarter <- floor(length(patients_ranked) * 0.25)
    bottom_patients_list <- head(patients_ranked, n_quarter)
    top_patients_list <- tail(patients_ranked, n_quarter)

    # Group into reference metadata for limma.
    # Add VO2peak_group column.
    patient_metadata_grouped <- patient_meta %>%
        filter(.data[[patient_ID_col_name]] %in% c(top_patients_list, bottom_patients_list)) %>%
        mutate(VO2peak_group = ifelse(.data[[patient_ID_col_name]] %in% top_patients_list, top_grp_name, bottom_grp_name)) %>%
        mutate(VO2peak_group = factor(VO2peak_group, levels = c(bottom_grp_name, top_grp_name)))


    ###
    # Sanity Check: Plot VO2peak distributions in top and bottom patient subsets, for the timepoint
    ###
    plot_VO2peak_distributions(patient_meta = patient_meta, top_patients_list = top_patients_list, bottom_patients_list = bottom_patients_list, tp = tp, output_dir = VO2peak_distrib_output_dir, VO2peak_col_name = VO2peak_col_name, patient_ID_col_name = patient_ID_col_name)

    return(patient_metadata_grouped)
}

append_results_with_metadata <- function(results_df, metadata_df, metadata_id_col = "OlinkID") {
  # Add annotations to limma results
  results_df[[metadata_id_col]] <- rownames(results_df)
  
  # Create unique metadata by metadata_id_col to avoid duplicates
  unique_metadata <- metadata_df %>%
    distinct(!!sym(metadata_id_col), .keep_all = TRUE)
  
  combined_df <- merge(results_df, unique_metadata, by = metadata_id_col, all.x = TRUE)
  
  # Deduplicate by SYMBOL: keep the row with the lowest P.Value per SYMBOL
  # (consistent with how common_functions.R handles multiple OlinkIDs mapping to the same gene)
  if ("SYMBOL" %in% colnames(combined_df)) {
    combined_df <- combined_df %>%
      group_by(SYMBOL) %>%
      arrange(P.Value, .by_group = TRUE) %>%
      filter(row_number() == 1L) %>%
      ungroup()
  }
  
  # Rearrange columns. Additional column names come directly from limma.
  desired_order <- c(metadata_id_col, "P.Value", "adj.P.Val", "logFC")
  desired_order <- desired_order[desired_order %in% colnames(combined_df)]
  remaining_cols <- setdiff(colnames(combined_df), desired_order)
  combined_df <- combined_df[, c(desired_order, remaining_cols)]
  
  return(combined_df)
}

save_prot_results <- function(results_df, prot_metadata, output_filepath, protein_identifier = "OlinkID", pval_cutoff = 0.05) { 
  all_prots <- append_results_with_metadata(
    results_df   = results_df,
    metadata_df  = prot_metadata,
    metadata_id_col = protein_identifier
  )
  
  # Sort rows by p-value
  all_prots <- all_prots %>%
    arrange(P.Value)
  
  # Save all results as excel
  write.xlsx(all_prots, file = output_filepath, rowNames = FALSE)
  message("All results saved: ", output_filepath, " (", nrow(all_prots), " total proteins)")
}


###
# Values
###
input_folder <- get_input_folder()
output_folder <- get_output_folder()
output_script_folder <- "proteomics_top_vs_bottom_diff_exp_limma"

data_processing_folder <- "proteomics_data_processing" # in output_folder
proteomics_input_filename <- "proteomics_data_processed.csv"
proteomics_metadata_filename <- "proteomics_patient_metadata.csv"
protein_metadata_filename <- "olink_mapped.csv" # in input_folder

outcome_var_prefix <- "VO2peak_"
patient_ID_col_name <- "PTID"
prot_id_col_name <- "OlinkID"
tp_col_name <- "visit"
abund_col_name <- "NPX_mean"

top_grp_name <- "High VO2peak"
bottom_grp_name <- "Low VO2peak"

# Modify based on VO2peak distributions
num_top_outliers_to_remove <- list(
    "Baseline" = 0,
    "6m" = 0,
    "12m" = 0
)
num_bottom_outliers_to_remove <- list(
    "Baseline" = 0,
    "6m" = 0,
    "12m" = 0
)


###
# Take into account HCT type subsets
###
if (subset_HCT_type == "ALLO" || subset_HCT_type == "AUTO") {
  data_processing_folder <- paste0(subset_HCT_type, "_", data_processing_folder)
  proteomics_input_filename <- gsub(".csv", paste0("_", subset_HCT_type, ".csv"), proteomics_input_filename)
  proteomics_metadata_filename <- gsub(".csv", paste0("_", subset_HCT_type, ".csv"), proteomics_metadata_filename)
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

VO2peak_distrib_output_dir <- file.path(output_dir, "VO2peak_Distributions")
ensure_dir(VO2peak_distrib_output_dir)

p_value_histogram_output_dir <- file.path(output_dir, "PValue_Histograms")
ensure_dir(p_value_histogram_output_dir)

residual_diagnostics_output_dir <- file.path(output_dir, "Residual_Diagnostics")
ensure_dir(residual_diagnostics_output_dir)

# Residual analysis report (QQ plots + summary statistics).
# See common_functions.R for the underlying functions.
residual_report_file <- file.path(output_dir, paste0("report", filename_suffix, ".txt"))
init_residual_report(residual_report_file,
                     script_label    = "proteomics_top_vs_bottom_diff_exp.R (limma, top vs bottom VO2peak)",
                     covars          = covars,
                     subset_HCT_type = subset_HCT_type)
residual_summaries <- list()


res_limma <- list() # to store limma result dataframes per timepoint


###
# Import Data
###
# Processed protein NPX measurements and patient metadata in proteomics dataset
prot_measurements <- read.csv(file.path(output_folder, data_processing_folder, proteomics_input_filename))
patient_meta <- read.csv(file.path(output_folder, data_processing_folder, proteomics_metadata_filename))
prot_metadata <- read.csv(file.path(input_folder, protein_metadata_filename))

# Convert covariates to factors as appropriate for limma design matrix
if ("TRANSPLANT_TYPE" %in% covars) {
    patient_meta$TRANSPLANT_TYPE <- factor(patient_meta$TRANSPLANT_TYPE, levels = c("ALLO", "AUTO"))
}
if ("gender" %in% covars) {
    patient_meta$gender <- factor(patient_meta$gender, levels = c("F", "M"))
}

# Make patientIDs character type
patient_meta[[patient_ID_col_name]] <- as.character(patient_meta[[patient_ID_col_name]])
prot_measurements[[patient_ID_col_name]] <- as.character(prot_measurements[[patient_ID_col_name]])


###
# Perform Differential Expression Analysis per Timepoint
###
for (tp in timepoint_order) {
    if (!is.null(covars_to_remove_at_baseline) && tp == "Baseline") {
        covars_for_tp <- setdiff(covars, covars_to_remove_at_baseline)
    } else {
        covars_for_tp <- covars
    }

    # Generate limma design matrix, VO2peak_group
    # will be the metadata col indicating top vs bottom patient subset
    if (length(covars_for_tp) > 0) {
        covar_terms <- paste(covars_for_tp, collapse = " + ")
        design_formula <- as.formula(paste("~ VO2peak_group + ", covar_terms))
    } else {
        design_formula <- as.formula("~ VO2peak_group")
    }


    ###
    # Generate Top and Bottom Patient Subsets, for the Timepoint
    ###
    patient_metadata_grouped <- subset_patients(patient_meta = patient_meta,
                                                tp = tp,
                                                patient_ID_col_name = patient_ID_col_name,
                                                outcome_var_prefix = outcome_var_prefix,
                                                num_top_outliers_to_remove = num_top_outliers_to_remove,
                                                num_bottom_outliers_to_remove = num_bottom_outliers_to_remove,
                                                top_grp_name = top_grp_name,
                                                bottom_grp_name = bottom_grp_name)

    ###
    # Save Lists of Patients in Subsets
    ###
    patient_metadata_grouped_filename_out <- paste0("patient_metadata_grouped_", tp, ".xlsx")
    write.xlsx(patient_metadata_grouped, file.path(output_dir, patient_metadata_grouped_filename_out))


    ###
    # Prepare Data for Differential Expression Analysis
    ###
    # Convert data to wide format for limma analysis (rows = OlinkID, cols = PTID)
    prot_wide_at_tp <- prot_measurements %>%
        filter(.data[[tp_col_name]] == tp) %>%
        pivot_wider(
            id_cols = all_of(prot_id_col_name),
            names_from = all_of(patient_ID_col_name),
            values_from = all_of(abund_col_name)
        )
    # Some patients in the top/bottom subsets may lack proteomics at this timepoint;
    # restrict patient_metadata_grouped to patients that are actually present in prot_wide_at_tp.
    available_ptids <- setdiff(colnames(prot_wide_at_tp), prot_id_col_name)
    missing_ptids <- setdiff(patient_metadata_grouped[[patient_ID_col_name]], available_ptids)
    if (length(missing_ptids) > 0) {
        message("Note: ", length(missing_ptids), " patient(s) in top/bottom subsets have no proteomics data at ",
                tp, " and will be excluded: ", paste(missing_ptids, collapse = ", "))
        patient_metadata_grouped <- patient_metadata_grouped %>%
            filter(.data[[patient_ID_col_name]] %in% available_ptids)
    }
    # Ensure that the columns in prot_wide_at_tp match the order of patients in patient_metadata_grouped
    prot_wide_at_tp <- prot_wide_at_tp %>%
        dplyr::select(all_of(prot_id_col_name), all_of(patient_metadata_grouped[[patient_ID_col_name]]))


    ###
    # Perform Differential Expression Analysis using limma
    ###
    patient_metadata_grouped$VO2peak_group <- factor(patient_metadata_grouped$VO2peak_group, levels = c(bottom_grp_name, top_grp_name))
    prot_mat <- as.matrix(prot_wide_at_tp[,-1]) # Remove OlinkID column for matrix
    rownames(prot_mat) <- prot_wide_at_tp[[prot_id_col_name]]

    print(paste0("Running limma differential expression analysis for ", tp))

    # Build design matrix
    design_tp <- model.matrix(design_formula, data = patient_metadata_grouped)

    # Fit and extract the results
    fit_tp <- lmFit(prot_mat, design_tp) %>%
        eBayes()
    results_tp <- topTable(
        fit_tp,
        coef = "VO2peak_groupHigh VO2peak",
        adjust.method="fdr",
        number = Inf,
        sort.by = "none"
    )

    # Save results
    res_limma[[tp]] <- results_tp
    save_prot_results(
        results_df = results_tp,
        prot_metadata = prot_metadata,
        output_filepath = file.path(output_dir, paste0("Limma_Results_Top_vs_Bottom_VO2peak_", tp, ".xlsx")),
        protein_identifier = prot_id_col_name
    )

    # P-value histogram
    pvalue_histogram(
        pvalues = results_tp$P.Value,
        outputdir = p_value_histogram_output_dir,
        label = paste0("Top_vs_Bottom_", tp)
    )

    # Residual analysis: QQ plot + summary statistics in report.txt
    residual_summaries[[tp]] <- run_limma_residual_diagnostics(
        prot_mat        = prot_mat,
        design          = design_tp,
        label           = paste0("Top_vs_Bottom_", tp),
        output_dir      = residual_diagnostics_output_dir,
        report_file     = residual_report_file,
        filename_suffix = filename_suffix,
        font_size       = font_size
    )
}


###
# Combined "Summary - All Contrasts" table
###
append_residual_summary_table(residual_report_file, residual_summaries)
message("Residual analysis report written: ", residual_report_file)


}


# Only run when executed directly (not when sourced from another script)
if (sys.nframe() == 0L) {
  run_top_vs_bottom_diff_exp(covars = c("AgeatHCT", "gender"), covars_to_remove_at_baseline = c(), subset_HCT_type = "ALL")
}