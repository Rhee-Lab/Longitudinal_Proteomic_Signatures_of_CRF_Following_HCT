# ==========================================
# Full Proteomics Batch Correction Workflow
# ==========================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(sva)
library(openxlsx)
source(here::here("common_functions.R"))


###
# To run code as a function:
run_proteomics_batch_correction_ComBat <- function(lod_filtering_cutoff) {
###

# Validate input
if (lod_filtering_cutoff < 0 || lod_filtering_cutoff > 1) {
  stop("lod_filtering_cutoff must be between 0 and 1 (e.g. 0.75 for 75%)")
}


###
# Functions
###
read_prot <- function(f) {
  read_delim(f, delim = ";", col_types = cols(.default = "c"), trim_ws = TRUE)
}

parse_PTID_visit_col <- function(proteomics_df) {
  proteomics_df <- proteomics_df %>%
    mutate(
      # Extract PTID from position [2] and visit from position [3]
      .split = strsplit(SampleID, "_"),
      PTID = sapply(.split, function(x) if (length(x) >= 2) x[2] else NA_character_),
      visit = sapply(.split, function(x) if (length(x) >= 3) x[3] else NA_character_)
    ) %>%
    dplyr::select(-.split) %>%
    mutate(visit = trimws(visit)) %>%
    mutate(visit = case_when(
      visit %in% c("Year 1", "Year1") ~ "12m",
      visit %in% c("Month 6", "Month6", "month 6") ~ "6m",
      TRUE ~ visit
    )) %>%
    mutate(PTID_visit = paste(PTID, visit, sep="_"))
  
  return(proteomics_df)
}

keep_only_in_common_assays <- function(proteomics_clean, report_file) {
  assays_by_plate <- table(proteomics_clean$PlateID, proteomics_clean$Assay, proteomics_clean$Assay_Warning)
  # In-common assays (assays missing in 0 plates)
  # Assay ~ specific protein measurement
  # Assay_Warning ~ if EXCLUDED in any plate, exclude that assay in all plates
  excluded_matrix <- assays_by_plate[, , "EXCLUDED"]
  assays_to_exclude <- colnames(excluded_matrix)[colSums(excluded_matrix) > 0]
  total_assays_prefiltering <- length(unique(proteomics_clean$OlinkID))
  assays_remaining <- total_assays_prefiltering - length(assays_to_exclude)

  append_report_line(report_file, "")
  append_report_line(report_file, paste("Total unique OlinkIDs before filtering:", total_assays_prefiltering))
  append_report_line(report_file, paste0("Number of assays removed due to EXCLUDED warning: ", length(assays_to_exclude)))
  append_report_line(report_file, paste0("Remaining: ", assays_remaining, " (", round((length(assays_to_exclude) / total_assays_prefiltering) * 100, 2), "% removed)"))

  # Keep only in-common assays
  common_assays <- setdiff(unique(proteomics_clean$Assay), assays_to_exclude)
  proteomics_clean <- proteomics_clean %>%
    filter(Assay %in% common_assays)
  return(proteomics_clean)
}

apply_lod_filter <- function(proteomics_df, cutoff, report_file, output_dir, font_size) {
  total_assays_prefiltering <- length(unique(proteomics_df$OlinkID))

  # Remove assays with NA LOD
  append_report_line(report_file, "")
  num_assays_with_na_lod <- proteomics_df %>%
    group_by(Assay) %>%
    summarise(has_na_lod = any(is.na(LOD)), .groups = "drop") %>%
    filter(has_na_lod) %>%
    nrow()
  append_report_line(report_file, paste0("Number of unique assays with NA LOD values removed: ", num_assays_with_na_lod))
  proteomics_df <- proteomics_df %>%
    filter(!is.na(LOD))

  # Calculate number and percent of NPX values below LOD
  num_entries_total <- nrow(proteomics_df)
  num_entries_below_lod <- sum(proteomics_df$NPX < proteomics_df$LOD, na.rm = TRUE)
  percent_below_lod <- (num_entries_below_lod / num_entries_total) * 100
  append_report_line(report_file, paste0("Number of NPX values below LOD: ", num_entries_below_lod, " out of ", num_entries_total, " total entries (", round(percent_below_lod, 2), "%)"))
  append_report_line(report_file, "Note: NPX values below LOD were not removed or imputed in this script, but this information is important for understanding data quality and may be relevant for downstream analyses.")

  # Identify assays with >= cutoff% of values below LOD
  proteomics_df$below_lod_flag <- proteomics_df$NPX < proteomics_df$LOD

  # Plot distribution of per-assay % below LOD with candidate cutoff lines
  all_assay_lod_pct <- proteomics_df %>%
    group_by(Assay) %>%
    summarise(percent_below_lod = mean(below_lod_flag, na.rm = TRUE) * 100, .groups = "drop")

  candidate_cutoffs <- c(50, 75)
  cutoff_colors <- c("#919191", "#363636")
  n_assays_total <- nrow(all_assay_lod_pct)
  cutoff_labels <- sapply(seq_along(candidate_cutoffs), function(i) {
    pct_removed <- mean(all_assay_lod_pct$percent_below_lod >= candidate_cutoffs[i]) * 100
    paste0(candidate_cutoffs[i], "% LOD cutoff: ", round(pct_removed, 1), "% of proteins removed")
  })
  # Highlight the active cutoff in red
  active_idx <- which(candidate_cutoffs == cutoff * 100)
  if (length(active_idx) > 0) cutoff_colors[active_idx] <- "#f86161"

  cutoff_df <- data.frame(
    xintercept = candidate_cutoffs,
    label = cutoff_labels,
    color = cutoff_colors,
    stringsAsFactors = FALSE
  )

  p_lod <- ggplot(all_assay_lod_pct, aes(x = percent_below_lod)) +
    geom_histogram(bins = 40, fill = "lightgray", color = "black") +
    geom_vline(data = cutoff_df, aes(xintercept = xintercept, color = color),
               linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
    scale_color_identity() +
    annotate("text", x = max(all_assay_lod_pct$percent_below_lod) * 0.95,
             y = Inf, vjust = seq(1.5, by = 1.5, length.out = nrow(cutoff_df)),
             label = cutoff_df$label, color = cutoff_df$color,
             hjust = 1, size = font_size * 0.28) +
    labs(
      title = "Distribution of % Values Below LOD per Assay",
      x = "% of values below LOD",
      y = "Count"
    ) +
    theme_bw() +
    theme(text = element_text(size = font_size))

  ggsave(file.path(output_dir, "assay_percent_below_lod_distribution.png"),
         p_lod, width = 10, height = 7, dpi = 600)

  assays_below_lod <- proteomics_df %>%
    group_by(Assay) %>%
    summarise(
      percent_below_lod = mean(below_lod_flag, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    filter(percent_below_lod >= cutoff * 100)
  append_report_line(report_file, "")
  append_report_line(report_file, paste("Assays with >=", cutoff * 100, "% of values below LOD:", nrow(assays_below_lod)))

  # Export LOD-filtered assay details to Excel
  if (nrow(assays_below_lod) > 0) {
    lod_filtered_summary <- proteomics_df %>%
      filter(Assay %in% assays_below_lod$Assay) %>%
      group_by(Assay, OlinkID) %>%
      summarise(
        n_total        = n(),
        n_below_lod    = sum(below_lod_flag, na.rm = TRUE),
        percent_below_lod = mean(below_lod_flag, na.rm = TRUE) * 100,
        .groups = "drop"
      ) %>%
      arrange(desc(percent_below_lod))

    lod_filtered_per_plate <- proteomics_df %>%
      filter(Assay %in% assays_below_lod$Assay) %>%
      group_by(Assay, OlinkID, PlateID) %>%
      summarise(
        n_total        = n(),
        n_below_lod    = sum(below_lod_flag, na.rm = TRUE),
        percent_below_lod = mean(below_lod_flag, na.rm = TRUE) * 100,
        .groups = "drop"
      ) %>%
      arrange(Assay, PlateID)

    xlsx_file <- file.path(output_dir, "LOD_filter_summary.xlsx")
    write.xlsx(
      list(
        "Summary"           = lod_filtered_summary,
        "Per_Plate_Breakdown" = lod_filtered_per_plate
      ), xlsx_file)
  }

  # Remove these assays from proteomics_df
  n_rows_before_lod_filter <- nrow(proteomics_df)
  proteomics_df <- proteomics_df %>%
    anti_join(assays_below_lod %>% dplyr::select(Assay), by = "Assay") %>%
    dplyr::select(-below_lod_flag)
  assays_remaining_after_lod_filter <- length(unique(proteomics_df$OlinkID))
  append_report_line(report_file, paste0("OlinkIDs remaining after filtering: ", assays_remaining_after_lod_filter, " (", round((nrow(assays_below_lod) / total_assays_prefiltering) * 100, 2), "% removed)"))

  return(proteomics_df)
}

generate_vo2_filtered <- function(vo2_df, proteomics_clean, report_file, timepoint_order, output_dir) {
  vo2_df_total_prefilter <- nrow(vo2_df)
  # Remove rows with no value in cpet_vo2_adjusted_num or visit columns
  vo2_df_filtered <- vo2_df[!is.na(vo2_df$cpet_vo2_adjusted_num) & !is.na(vo2_df$visit), ]
  append_report_line(report_file, "")
  append_report_line(report_file, "For VO2peak coverage assessment, removed rows with missing cpet_vo2_adjusted_num or visit values.")
  append_report_line(report_file, paste("Remaining rows:", nrow(vo2_df_filtered)))
  append_report_line(report_file, paste("Removed rows:", vo2_df_total_prefilter - nrow(vo2_df_filtered)))

  # Convert numeric visit codes (0, 1, ...) to visit labels using timepoint_order vector
  visit_map <- setNames(timepoint_order, seq_along(timepoint_order) - 1)
  vo2_df_filtered <- vo2_df_filtered %>%
    mutate(visit = {
      mapped <- visit_map[as.character(visit)]
      ifelse(!is.na(mapped), mapped, as.character(visit))
    })

  # Export filtered version to output directory
  write.csv(vo2_df_filtered, file.path(output_dir, "CRESTpfizer_VO2peak_100125_filtered.csv"), row.names = FALSE)

  # Filter out PTID-visit combos with no VO2 data (before batch correction)
  vo2_ptid_visit_list <- paste(vo2_df_filtered$PTID, vo2_df_filtered$visit, sep="_")

  # Count unique PTID_visits per visit BEFORE VO2 filtering
  pre_filter_counts <- proteomics_clean %>%
    group_by(visit) %>%
    summarise(count_before = n_distinct(PTID_visit), .groups = "drop")

  pre_filter_rows <- n_distinct(proteomics_clean$PTID_visit)
  proteomics_clean <- proteomics_clean %>%
    filter(PTID_visit %in% vo2_ptid_visit_list)
  post_filter_rows <- n_distinct(proteomics_clean$PTID_visit)

  # Count unique PTID_visits per visit AFTER VO2 filtering
  post_filter_counts <- proteomics_clean %>%
    group_by(visit) %>%
    summarise(count_after = n_distinct(PTID_visit), .groups = "drop")

  # Merge and calculate removed counts per visit
  filter_comparison <- pre_filter_counts %>%
    full_join(post_filter_counts, by = "visit") %>%
    mutate(count_before = ifelse(is.na(count_before), 0, count_before),
          count_after = ifelse(is.na(count_after), 0, count_after),
          removed = count_before - count_after)

  append_report_line(report_file, "")
  append_report_line(report_file, paste0("Removed ", pre_filter_rows - post_filter_rows, " PTID-visit combos from proteomics_clean that had no matching PTID-visit in VO2 data (pre-batch correction).\nRemaining unique PTID_visits: ", post_filter_rows))
  append_report_line(report_file, "")
  append_report_line(report_file, "PTID_visit counts removed per visit due to VO2 filtering (pre-batch correction):")
  for (i in seq_len(nrow(filter_comparison))) {
    append_report_line(report_file, paste0("  Visit: ", filter_comparison$visit[i], ", Before: ", filter_comparison$count_before[i], ", After: ", filter_comparison$count_after[i], ", Removed: ", filter_comparison$removed[i]))
  }
  return(vo2_df_filtered)
}

generate_pca_plot_batch_correction <- function(df, batch_vector, output_path, title_prefix = "PCA") {
  # Transpose matrix so samples are rows (required for PCA)
  df_t <- t(df)
  
  # Perform PCA
  pca_result <- prcomp(df_t, center = TRUE, scale. = TRUE)
  
  # Extract variance explained
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100
  
  # Create data frame for plotting
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Sample = rownames(pca_result$x),
    Batch = batch_vector
  )
  
  # Create PCA plot
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(
      title = paste(title_prefix, "- Colored by Batch"),
      x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 2), "%)")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    ) +
    scale_color_brewer(palette = "Paired")
  
  # Save plot
  ggsave(output_path, p, width = 10, height = 7, dpi = 600)
}


###
# Values
###
input_folder <- get_input_folder()
input_subfolder <- "Raw_Olink_NPX_data"
output_folder <- get_output_folder()
output_subfolder <- "batch_correction_output"

# Input filenames
proteomics_files <- c(
  "CardioTox_CoH_1-4_EXTENDED_NPX_2025-03-17.csv",
  "CardioTox_CoH_PM5_EXTENDED_NPX_2025-03-17.csv",
  "CardioTox_CoH_PM7_8_EXTENDED_NPX_2025-03-17.csv",
  "CardioToxPM9_EXTENDED_NPX_2025-04-01.csv",
  "CardioToxPM10_6rpt_EXTENDED_NPX_2025-04-02.csv"
)

vo2_filename <- "CRESTpfizer_VO2_100125.csv"
patient_data_filename <- "Patient_data.csv"

# Set options (font_size comes from config.R)
proteomics_clean_cols_keep <- c("PTID", "visit", "PTID_visit", "PlateID", "SampleID", "tubeID", "OlinkID", "NPX")
combat_model_covariates <- c("gender", "BMI", "AgeatHCT", "TRANSPLANT_TYPE", "visit", "cpet_vo2_adjusted_num")


###
# Organize Directories
###
output_dir <- file.path(output_folder, output_subfolder)
ensure_dir(output_dir)

# Initialize report file
report_file <- file.path(output_dir, "batch_correction_report.txt")
if (file.exists(report_file)) file.remove(report_file)


###
# Import and combine proteomics data
###
append_report_line(report_file, "==================================================")
append_report_line(report_file, "Filtering Steps")
append_report_line(report_file, "==================================================")
proteomics <- proteomics_files %>%
  lapply(function(f) read_prot(file.path(input_folder, input_subfolder, f))) %>%
  bind_rows()
append_report_line(report_file, "")
append_report_line(report_file, paste("Total rows:", nrow(proteomics)))

# Clean up column types and whitespace
proteomics <- proteomics %>%
  mutate(across(c(SampleID, PlateID, Assay, OlinkID), ~ trimws(as.character(.))),
         NPX = as.numeric(gsub(",", ".", NPX)))


###
# Remove control and bridging samples
###
# Remove rows that include "CTRL", "CONTROL", "Pool", or "Bridging" in SampleID (case-insensitive)
proteomics_clean <- proteomics %>%
  filter(!grepl("CTRL|CONTROL|Pool|Bridging", SampleID, ignore.case = TRUE)) %>%
  filter(!grepl("Control", Sample_Type, ignore.case = TRUE)) %>%
  filter(!Assay %in% c("Amplification", "Extension", "Incubation"))

append_report_line(report_file, paste("Rows removed from proteomics (controls + bridging):", nrow(proteomics) - nrow(proteomics_clean)))


###
# Parse ID columns
###
# Parse PTID and visit
proteomics_clean <- parse_PTID_visit_col(proteomics_clean)
# Parse tubeID
proteomics_clean <- proteomics_clean %>%
  mutate(tubeID = sapply(strsplit(as.character(SampleID), "_"), function(x) x[1]))


###
# Identify common assays across all plates
###
proteomics_clean <- keep_only_in_common_assays(proteomics_clean, report_file)


###
# Generate NPX distribution with input data (pre-filtering/batch correction)
###
# Use plot_distribution_histogram
avg_prot_abundance_prefilter <- proteomics_clean %>%
  group_by(OlinkID) %>%
  summarise(mean_NPX = mean(NPX, na.rm = TRUE), .groups = "drop")
plot_distribution_histogram(
  data     = avg_prot_abundance_prefilter$mean_NPX,
  title    = "Protein Distribution, Pre-Filtering and Batch Correction",
  xlab     = "Average NPX per Protein",
  filename = file.path(output_dir, "average_protein_abundance_pre_filtering_and_batch_correction.png"),
  text_size = font_size
)


###
# Remove assays with >= lod_filtering_cutoff % of values below LOD
###
proteomics_clean <- apply_lod_filter(proteomics_clean, lod_filtering_cutoff, report_file = report_file, output_dir = output_dir, font_size = font_size)


###
# Generate NPX distribution pre-batch correction
###
# Use plot_distribution_histogram
avg_prot_abundance <- proteomics_clean %>%
  group_by(OlinkID) %>%
  summarise(mean_NPX = mean(NPX, na.rm = TRUE), .groups = "drop")
plot_distribution_histogram(
  data     = avg_prot_abundance$mean_NPX,
  title    = "Protein Distribution, Pre-Batch Correction",
  xlab     = "Average NPX per Protein",
  filename = file.path(output_dir, "average_protein_abundance_pre_batch_correction.png"),
  text_size = font_size
)


###
# Prepare data for batch correction
###
# Keep only necessary columns (retain PTID and visit for covariate model matrix)
# Note: column selection is done here (after LOD filtering) because apply_lod_filter
# requires columns like Assay and LOD that are not needed downstream
proteomics_clean <- proteomics_clean %>%
  dplyr::select(all_of(proteomics_clean_cols_keep))

# Check how many NPX values are NA
num_npx_na <- sum(is.na(proteomics_clean$NPX))
append_report_line(report_file, "")
append_report_line(report_file, paste("Number of NPX values that are NA:", num_npx_na, "out of", nrow(proteomics_clean), "rows.\nIf no NAs, then NAs were handled previously. NAs due to missing assays was handled in this script."))


###
# Import VO2peak data and filter proteomics to samples with VO2 data
# (ComBat model matrix cannot handle NAs, so filter before batch correction)
###
# Generate vo2_df_filtered, with rows with missing VO2 values removed
# Use vo2_df_filtered to filter proteomics_clean to samples with VO2 data (pre-batch correction)
vo2_df <- read.csv(file.path(input_folder, vo2_filename))
vo2_df_filtered <- generate_vo2_filtered(vo2_df, proteomics_clean, report_file = report_file, timepoint_order = timepoint_order, output_dir = output_dir)

# Create PTID_visit_PlateID column to keep plate-specific measurements separate for ComBat
proteomics_clean <- proteomics_clean %>%
  mutate(PTID_visit_PlateID = paste(PTID_visit, PlateID, sep = "_PLATE_"))

# Create mapping of PTID_visit_PlateID to PlateID for batch information
sample_plate_mapping <- proteomics_clean %>%
  dplyr::select(PTID_visit_PlateID, PlateID) %>%
  distinct()

# Create wide matrix for ComBat: proteins as rows, samples as columns
proteomics_wide_for_combat <- proteomics_clean %>%
  dplyr::select(PTID_visit_PlateID, OlinkID, NPX) %>%
  pivot_wider(names_from = PTID_visit_PlateID, values_from = NPX)

# Extract protein names and create numeric matrix, combat_matrix
protein_names <- proteomics_wide_for_combat$OlinkID
combat_matrix <- proteomics_wide_for_combat %>%
  dplyr::select(-OlinkID) %>%
  as.matrix()
rownames(combat_matrix) <- protein_names

# Ensure sample_plate_mapping is in same order as matrix columns
sample_plate_mapping <- sample_plate_mapping %>%
  filter(PTID_visit_PlateID %in% colnames(combat_matrix)) %>%
  arrange(match(PTID_visit_PlateID, colnames(combat_matrix)))

# Create batch vector (PlateID for each sample)
batch <- sample_plate_mapping$PlateID

append_report_line(report_file, "==================================================")
append_report_line(report_file, "ComBat Batch Correction")
append_report_line(report_file, "==================================================")
append_report_line(report_file, "")
append_report_line(report_file, paste0("ComBat matrix dimensions: ", nrow(combat_matrix), " proteins x ", ncol(combat_matrix), " samples"))
append_report_line(report_file, paste0("Batch levels: ", paste(unique(batch), collapse=", ")))
append_report_line(report_file, "Samples per batch:")
batch_table <- table(batch)
for (i in seq_along(batch_table)) {
  append_report_line(report_file, paste0("  ", names(batch_table)[i], ": ", batch_table[i]))
}

# Check for minimum batch size before ComBat
batch_counts <- table(batch)
small_batches <- names(batch_counts)[batch_counts < 2]
if (length(small_batches) > 0) {
  append_report_line(report_file, paste0("WARNING: Batch(es) with <2 samples: ", paste(small_batches, collapse = ", "), ". ComBat may produce unstable estimates."))
}

# Remove proteins with zero or near-zero variance
protein_vars <- apply(combat_matrix, 1, var, na.rm = TRUE)
zero_var_proteins <- which(protein_vars < 1e-10 | is.na(protein_vars))
if (length(zero_var_proteins) > 0) {
  append_report_line(report_file, paste("Removing", length(zero_var_proteins), "proteins with zero/near-zero variance"))
  combat_matrix <- combat_matrix[-zero_var_proteins, ]
}

append_report_line(report_file, paste0("Final matrix dimensions for ComBat: ", nrow(combat_matrix), " proteins x ", ncol(combat_matrix), " samples"))

# Generate PCA plot before batch correction
pca_before_path <- file.path(output_dir, "pca_before_batch_correction.png")
generate_pca_plot_batch_correction(combat_matrix, batch, pca_before_path, "PCA Before Batch Correction")


###
# Build model matrix for ComBat
# Includes biological covariates and variables of interest
# Ensures ComBat preserves variation associated with these variables
# while only removing batch (plate) effects.
###
append_report_line(report_file, "")
append_report_line(report_file, "Building model matrix (vars of interest + covariates) for ComBat...")

# Import patient info data
patient_data <- read.csv(file.path(input_folder, patient_data_filename))
patient_data$PTID <- as.character(patient_data$PTID)

# Impute numerical covariates with median value if NA
covars_num <- c("BMI", "Height_CM", "Weight_KG", "AgeatHCT")
covars_num <- covars_num[covars_num %in% colnames(patient_data)]
for (v in covars_num) {
  med <- median(patient_data[[v]], na.rm = TRUE)
  patient_data[[v]][is.na(patient_data[[v]])] <- med
}

# Use VO2 data already imported and filtered above (vo2_df_filtered)
vo2_data_for_model <- vo2_df_filtered %>%
  mutate(
    PTID = as.character(PTID),
    PTID_visit = paste(PTID, visit, sep = "_")
  )

# Build sample metadata aligned with combat_matrix columns
sample_metadata <- proteomics_clean %>%
  dplyr::select(PTID_visit_PlateID, PTID, visit, PTID_visit) %>%
  distinct() %>%
  filter(PTID_visit_PlateID %in% colnames(combat_matrix)) %>%
  arrange(match(PTID_visit_PlateID, colnames(combat_matrix)))

# Join with patient info (on PTID) — deduplicate to prevent row expansion
sample_metadata <- sample_metadata %>%
  left_join(
    patient_data %>%
      dplyr::select(PTID, gender, BMI, AgeatHCT, TRANSPLANT_TYPE) %>%
      distinct(PTID, .keep_all = TRUE),
    by = "PTID"
  )

# Join with VO2 data (on PTID_visit) — deduplicate to prevent row expansion
sample_metadata <- sample_metadata %>%
  left_join(
    vo2_data_for_model %>%
      dplyr::select(PTID_visit, cpet_vo2_adjusted_num) %>%
      distinct(PTID_visit, .keep_all = TRUE),
    by = "PTID_visit"
  )

# Report covariate coverage (check ALL model covariates for NAs)
n_total <- nrow(sample_metadata)

append_report_line(report_file, paste0("  Total samples in combat matrix: ", n_total))
for (cov in combat_model_covariates) {
  n_miss <- sum(is.na(sample_metadata[[cov]]))
  append_report_line(report_file, paste0("  Samples missing ", cov, ": ", n_miss))
}

# Identify samples with complete covariate data across all model variables
complete_rows <- complete.cases(sample_metadata[, combat_model_covariates, drop = FALSE])
n_incomplete <- sum(!complete_rows)

if (n_incomplete > 0) {
  incomplete_ids <- sample_metadata$PTID_visit_PlateID[!complete_rows]
  append_report_line(report_file, paste0("  WARNING: Removing ", n_incomplete, " samples with missing covariate data from combat_matrix and model matrix."))
  for (id in incomplete_ids) {
    append_report_line(report_file, paste0("    Removed: ", id))
  }

  # Filter sample_metadata, combat_matrix, and batch vector to keep only complete cases
  sample_metadata <- sample_metadata[complete_rows, ]
  keep_cols <- sample_metadata$PTID_visit_PlateID
  combat_matrix <- combat_matrix[, keep_cols, drop = FALSE]
  batch <- sample_plate_mapping$PlateID[
    match(sample_metadata$PTID_visit_PlateID, sample_plate_mapping$PTID_visit_PlateID)
  ]

  append_report_line(report_file, paste0("  Updated ComBat matrix dimensions: ", nrow(combat_matrix), " proteins x ", ncol(combat_matrix), " samples"))
}

# Ensure categorical variables are factors
sample_metadata <- sample_metadata %>%
  mutate(
    gender = as.factor(gender),
    TRANSPLANT_TYPE = as.factor(TRANSPLANT_TYPE),
    visit = factor(visit, levels = timepoint_order)
  )

# Build model matrix with biological covariates
mod <- model.matrix(~ gender + BMI + AgeatHCT + TRANSPLANT_TYPE + visit + cpet_vo2_adjusted_num,
                    data = sample_metadata)
model_formula_str <- "~ gender + BMI + AgeatHCT + TRANSPLANT_TYPE + visit + cpet_vo2_adjusted_num"

append_report_line(report_file, paste0("  Model formula: ", model_formula_str))
append_report_line(report_file, paste0("  Model matrix dimensions: ", nrow(mod), " samples x ", ncol(mod), " covariates"))
append_report_line(report_file, paste0("  Model matrix columns: ", paste(colnames(mod), collapse = ", ")))

# Verify alignment with combat_matrix
stopifnot("Model matrix rows must equal combat_matrix columns" = nrow(mod) == ncol(combat_matrix))


###
# Apply ComBat batch correction
###
append_report_line(report_file, "")
append_report_line(report_file, "Applying ComBat batch correction with biological covariate model...")
combat_corrected <- ComBat(dat = combat_matrix, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)
append_report_line(report_file, "ComBat correction complete.")

# Generate PCA plot after batch correction
pca_after_path <- file.path(output_dir, "pca_after_batch_correction.png")
generate_pca_plot_batch_correction(combat_corrected, batch, pca_after_path, "PCA After Batch Correction")

# Convert back to long format
proteomics_corrected <- combat_corrected %>%
  as.data.frame() %>%
  tibble::rownames_to_column("OlinkID") %>%
  pivot_longer(cols = -OlinkID, names_to = "PTID_visit_PlateID", values_to = "NPX_corrected")

# Add back PlateID, PTID_visit, and SampleID metadata
proteomics_corrected <- proteomics_corrected %>%
  left_join(
    proteomics_clean %>% 
      dplyr::select(PTID_visit_PlateID, PlateID, PTID_visit, SampleID, tubeID) %>% 
      distinct(PTID_visit_PlateID, .keep_all = TRUE),
    by = "PTID_visit_PlateID"
  )

append_report_line(report_file, "")
append_report_line(report_file, paste0("Proteomics corrected data dimensions: ", nrow(proteomics_corrected), " rows"))

# Save the full data table to compare directly with other script run
write.csv(proteomics_corrected, file.path(output_dir, "proteomics_post_correction.csv"), row.names = FALSE)


###
# Aggregate duplicates (average) across plates and pivot wide
###
# Count unique samples before aggregation (plate-specific samples)
num_samples_before_aggregating_duplicates <- n_distinct(proteomics_corrected$PTID_visit_PlateID)

# Average across plates for the same PTID_visit after batch correction
# This is the correct time to do this averaging - after ComBat has corrected for batch effects
proteomics_wide_corrected <- proteomics_corrected %>%
  group_by(PTID_visit, OlinkID) %>%
  summarise(NPX_corrected = mean(NPX_corrected, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = OlinkID, values_from = NPX_corrected)

# Count unique samples after aggregation (final unique PTID_visits)
num_samples_after_aggregating_duplicates <- nrow(proteomics_wide_corrected)

append_report_line(report_file, "")
append_report_line(report_file, paste("Number of unique samples with plate designation before aggregation:", num_samples_before_aggregating_duplicates))
append_report_line(report_file, paste("Number of unique PTID_visit samples after aggregating across plates:", num_samples_after_aggregating_duplicates))
append_report_line(report_file, paste("Number of sample entries collapsed:", num_samples_before_aggregating_duplicates - num_samples_after_aggregating_duplicates))


###
# Check PTID_visit coverage after batch correction
###
append_report_line(report_file, "")
append_report_line(report_file, "PTID_visit coverage after filtering:")
coverage_summary <- proteomics_wide_corrected %>%
  mutate(visit = sapply(strsplit(PTID_visit, "_"), function(x) x[2])) %>%
  group_by(visit) %>%
  summarise(Unique_PTIDs = n_distinct(PTID_visit)) # print number of unique PTID_visits per visit after filtering
for (i in seq_len(nrow(coverage_summary))) {
  append_report_line(report_file, paste0("  Visit: ", coverage_summary$visit[i], ", Unique PTID_visits: ", coverage_summary$Unique_PTIDs[i]))
}


###
# Export batch-corrected proteomics data
###
# Print any PTID_visit values that are duplicated (should be none)
duplicated_ptid_visits <- proteomics_wide_corrected$PTID_visit[duplicated(proteomics_wide_corrected$PTID_visit)]
append_report_line(report_file, "")
if (length(duplicated_ptid_visits) > 0) {
  append_report_line(report_file, "Warning: Found duplicated PTID_visit values after processing:")
  for (ptid in duplicated_ptid_visits) {
    append_report_line(report_file, paste0("  ", ptid))
  }
} else {
  append_report_line(report_file, "No duplicated PTID_visit values found after processing.")
}

# Transpose: proteins as rows, samples as columns
proteomics_export <- proteomics_wide_corrected %>%
  pivot_longer(cols = -PTID_visit, names_to = "OlinkID", values_to = "NPX_mean") %>%
  pivot_wider(names_from = PTID_visit, values_from = NPX_mean)

write.csv(proteomics_export, file.path(output_dir, "proteomics_wide_corrected.csv"), row.names = FALSE)

append_report_line(report_file, "")
append_report_line(report_file, "==================================================")
append_report_line(report_file, "Batch correction completed successfully.")
append_report_line(report_file, "==================================================")

}


# Only run when executed directly (not when sourced from another script)
if (sys.nframe() == 0L) {
  run_proteomics_batch_correction_ComBat(lod_filtering_cutoff = proteomics_lod_filtering_cutoff)
}