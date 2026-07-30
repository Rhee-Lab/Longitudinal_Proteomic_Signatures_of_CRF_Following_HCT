###
# Proteomics linear regression Plots
###
library(dplyr)
library(openxlsx)
library(ggplot2)
source(here::here("common_functions.R"))
# Seed set in common_functions.R via config.R

# For network generation:
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(igraph)
library(ggraph)


###
# To run code as a function:
run_proteomics_linear_regression_plots <- function(subset_HCT_type, jobs_to_run) {
###


###
# Functions
###
plot_pathway_focus_network <- function(gsea_res, prot_assoc_df, n_top = 10, pathway_pval_col = "padj",
                                       pathway_pval_cutoff = 0.1, gene_pval_col = "FDR",
                                       gene_pval_cutoff = 0.05, title = "", font_size = 16,
                                       category_color = "black", category_fontface = "bold",
                                       item_color = "grey10", select_pathways = NULL) {
  ###
  # Plot networks with a focus on pathway results, showing associated protein connectedness.
  # Inputs:
  # - gsea_res: Dataframe of GSEA results with columns for pathway name, p-value, NES, and leading-edge genes
  # - prot_assoc_df: Dataframe of association results with columns for gene symbol, beta, and p-value
  # - n_top: Number of top pathways to include based on pathway_pval_col
  # - pathway_pval_col: Column name in gsea_res to use for pathway significance filtering
  # - pathway_pval_cutoff: P-value cutoff for including pathways
  # - gene_pval_col: Column name in prot_assoc_df to use for gene significance filtering
  # - gene_pval_cutoff: P-value cutoff for labeling genes
  # - title: Plot title
  # - font_size: Base font size for plot
  # - category_color: Color for pathway text
  # - category_fontface: Font face for pathway text
  # - item_color: Color for gene node text
  # - select_pathways: Optional vector of specific pathway names to include (overrides n_top and pathway_pval_cutoff)
  ###
  if (!is.null(select_pathways)) {
    top_paths <- gsea_res %>%
      filter(pathway %in% select_pathways)
  } else {
    top_paths <- gsea_res %>%
      filter(!is.na(.data[[pathway_pval_col]]), !is.na(Genes), Genes != "") %>%
      filter(.data[[pathway_pval_col]] <= pathway_pval_cutoff) %>%
      arrange(.data[[pathway_pval_col]]) %>%
      slice_head(n = n_top)
  }

  if (nrow(top_paths) == 0) return(NULL)

  # Build pathway -> gene edge list from leading-edge gene strings
  edges <- top_paths %>%
    rowwise() %>%
    mutate(gene = list(trimws(strsplit(Genes, ",\\s*")[[1]]))) %>%
    unnest(gene) %>%
    filter(nchar(gene) > 0) %>%
    ungroup() %>%
    dplyr::select(from = pathway, to = gene)

  if (nrow(edges) == 0) return(NULL)

  # Pathway nodes: coloured by NES, sized by -log10(padj)
  # Gene nodes: coloured by regression beta, sized by -log10(FDR); label those meeting p-adj cutoff
  pathway_nodes <- top_paths %>%
    transmute(name = pathway, node_type = "Pathway", value = NES,
              node_size = -log10(pmax(.data[[pathway_pval_col]], 1e-300)),
              sig_label = FALSE)

  gene_nodes <- tibble(name = unique(edges$to), node_type = "Protein") %>%
    left_join(
      prot_assoc_df %>% dplyr::select(all_of(c("SYMBOL", "beta", gene_pval_col))) %>% distinct(),
      by = c("name" = "SYMBOL")
    ) %>%
    dplyr::rename(value = beta) %>%
    mutate(
      node_size = ifelse(is.na(.data[[gene_pval_col]]), 1, -log10(pmax(.data[[gene_pval_col]], 1e-300))),
      sig_label = !is.na(.data[[gene_pval_col]]) & .data[[gene_pval_col]] <= gene_pval_cutoff
    )

  nodes <- bind_rows(pathway_nodes, gene_nodes)

  nodes <- nodes %>%
    mutate(label_wrapped = ifelse(node_type == "Pathway",
           stringr::str_wrap(normalize_pathway_name_publish(NULL, name), width = 25), name))

  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

  set.seed(random_seed)
  layout_coords <- create_layout(g, layout = "fr", niter = 1500,
                                 start.temp = sqrt(igraph::vcount(g)) * 3)
  layout_coords$x <- layout_coords$x * 1.7
  layout_coords$y <- layout_coords$y * 1.7

  # Combine all labels into one ggrepel layer so pathway and gene labels repel
  # each other — separate geom_node_text(repel=TRUE) layers are unaware of each
  # other and will freely overlap.
  # Include ALL nodes with empty labels for non-labeled ones so that ggrepel
  # repels pathway/protein labels away from every node position, not just
  # the labeled ones. This pushes text into surrounding whitespace.
  label_data <- as.data.frame(layout_coords) %>%
    dplyr::mutate(
      has_label      = node_type == "Pathway" | (node_type == "Protein" & sig_label),
      display_label  = dplyr::case_when(
        node_type == "Pathway"             ~ label_wrapped,
        node_type == "Protein" & sig_label ~ name,
        TRUE                               ~ ""
      ),
      label_fontface = ifelse(node_type == "Pathway", category_fontface, "plain"),
      label_size     = ifelse(has_label, font_size * 0.3, font_size * 0.3)
    )

  p <- suppressWarnings(
    ggraph(layout_coords) +
      geom_edge_link(alpha = 0.2, color = "#383838", width = 0.6) +
      geom_node_point(aes(filter = node_type == "Pathway", shape = node_type, size = node_size), fill = "darkgoldenrod1", color = "black") +
      geom_node_point(aes(filter = node_type == "Protein", fill = value, shape = node_type, size = node_size), color = "black") +
      ggrepel::geom_text_repel(
        data               = label_data,
        aes(x = x, y = y, label = display_label, fontface = label_fontface),
        size               = I(label_data$label_size),
        color              = scales::alpha("black", 0.7),
        bg.color           = scales::alpha("white", 0.5),
        bg.r               = 0.15,
        force              = 15,   force_pull        = 0.5,
        box.padding        = 1, point.padding     = 0.4,
        max.overlaps       = Inf,
        max.iter           = 20000,
        segment.color      = "#2b2b2b", segment.size = 0.3,
        min.segment.length = 0.2,
        seed               = 42,
        lineheight         = 0.85,
        show.legend        = FALSE
      ) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, na.value = "grey70", name = expression(beta)) +
      scale_shape_manual(values = c("Pathway" = 23, "Protein" = 21), name = "Node Type", guide = guide_legend(override.aes = list(size = 5, color = "black", fill = c("darkgoldenrod1", "white")))) +
      scale_size_continuous(name = "-log10(p-adj)", range = c(3, 10)) +
      theme_graph(base_family = font_style) +
      theme(
        legend.text      = element_text(size = font_size),
        legend.title     = element_text(size = font_size),
        legend.key.size  = unit(1.2, "lines"),
        legend.spacing.y = unit(0.3, "lines")
      )
  )
  return(p)
}


###
# Values
###
output_folder <- get_output_folder()
# sub-folder of output folder
enrichment_output_folder <- "proteomics_linear_regression_enrichment"

output_script_folder <- "proteomics_linear_regression_plots"
timeseries_subfolder <- "Time_Series_Heatmaps"
volcano_subfolder <- "Volcano_Plots"
barplot_subfolder <- "Bar_Plots"
pathway_focus_network_subfolder  <- "Pathway_Focus_Network_Plots"

approach_contrast_combos <- list("Cross_sectional" = timepoint_order,
                      "Delta" = c("Delta_6m", "Delta_12m"),
                      "Pct_Change" = c("Pct_Change_6m", "Pct_Change_12m"),
                      "Lagged_Association" = c("Baseline_Delta_6m", "Baseline_Delta_12m", "6m_Delta_12m"),
                      "Double_Delta_Comparison" = c("Double_Delta_6m", "Double_Delta_12m"))

# Map jobs_to_run names to approach_contrast_combos keys
jobs_to_approach_map <- c(
  "Cross-Sectional"    = "Cross_sectional",
  "Delta"              = "Delta",
  "Pct_Change"         = "Pct_Change",
  "Lagged_Association" = "Lagged_Association",
  "Double_Delta"       = "Double_Delta_Comparison"
)

# Filter approach_contrast_combos and contrasts based on jobs_to_run
active_approaches <- unname(jobs_to_approach_map[intersect(names(jobs_to_approach_map), jobs_to_run)])
approach_contrast_combos <- approach_contrast_combos[intersect(names(approach_contrast_combos), active_approaches)]
contrasts <- unname(unlist(approach_contrast_combos))


###
# Take into account HCT type subsets
###
if (subset_HCT_type == "ALLO" || subset_HCT_type == "AUTO") {
  # Adjust all relevant file paths and names for ALLO or AUTO subset
  enrichment_output_folder <- paste0(subset_HCT_type, "_", enrichment_output_folder)
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
ensure_dir(file.path(output_dir, timeseries_subfolder))
ensure_dir(file.path(output_dir, volcano_subfolder))
ensure_dir(file.path(output_dir, barplot_subfolder))
ensure_dir(file.path(output_dir, pathway_focus_network_subfolder))

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
# Time Series Heatmap, Individual Proteins
###
top_n_list <- c(30, 50)

for (approach in names(approach_contrast_combos)) {
  tps <- approach_contrast_combos[[approach]]
  run_protein_time_series_heatmap_generation(
    results_list = lin_regr_results_list,
    tps_for_focus = tps,
    top_n_list = top_n_list,
    outputdir = file.path(output_dir, ts_proteins_subfolder),
    filename_suffix = filename_suffix,
    timepoints_label = approach,
    pval_col = "FDR",
    pval_cutoff = 0.05
  )
}


###
# Time Series Heatmap — Proteins significant (FDR < 0.05) across ALL cross-sectional contrasts
###
if ("Cross-Sectional" %in% jobs_to_run) {
  cs_tps <- approach_contrast_combos$Cross_sectional  # c("Baseline", "6m", "12m")
  sig_in_all_cs <- Reduce(
    intersect,
    lapply(cs_tps, function(tp) {
      lin_regr_results_list[[tp]] %>%
        filter(FDR < 0.05) %>%
        pull(SYMBOL)
    })
  )
  if (length(sig_in_all_cs) > 0) {
    # Order proteins by Baseline beta (descending) so positive associations appear first
    sig_in_all_cs_sorted <- lin_regr_results_list[["Baseline"]] %>%
      filter(SYMBOL %in% sig_in_all_cs) %>%
      arrange(desc(beta)) %>%
      pull(SYMBOL)
    plot_protein_time_series_heatmap(
      results_list   = lin_regr_results_list,
      proteins       = sig_in_all_cs_sorted,
      timepoints     = cs_tps,
      outputdir      = file.path(output_dir, ts_proteins_subfolder),
      filename       = paste0("Protein_Time_Series_Heatmap_Sig_All_CS", filename_suffix, ".png"),
      filename_suffix = filename_suffix,
      sig_cutoff     = 0.05,
      sort_type = "hierarchical clustering"
    )
    # Also plot option without the sig. dots
    plot_protein_time_series_heatmap(
      results_list   = lin_regr_results_list,
      proteins       = sig_in_all_cs_sorted,
      timepoints     = cs_tps,
      outputdir      = file.path(output_dir, ts_proteins_subfolder),
      filename       = paste0("Protein_Time_Series_Heatmap_Sig_All_CS_No_Sig_Dots", filename_suffix, ".png"),
      filename_suffix = filename_suffix,
      sig_cutoff     = NULL,
      sort_type = "hierarchical clustering"
    )
    # Also plot with proteins ordered by descending baseline beta
    plot_protein_time_series_heatmap(
      results_list   = lin_regr_results_list,
      proteins       = sig_in_all_cs_sorted,
      timepoints     = cs_tps,
      outputdir      = file.path(output_dir, ts_proteins_subfolder),
      filename       = paste0("Protein_Time_Series_Heatmap_Sig_All_CS_Ordered_Baseline_Beta", filename_suffix, ".png"),
      filename_suffix = filename_suffix,
      sig_cutoff     = 0.05,
      sort_type = "Baseline"
    )
  } else {
    message("No proteins are significant (FDR < 0.05) across all cross-sectional contrasts.")
  }
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

if ("Cross-Sectional" %in% jobs_to_run) {
  run_tp_focus_annotation_time_series_heatmap_generation(
    annotation_sets = annotation_sets,
    annotation_subdirs = annotation_subdirs,
    tps_for_focus = approach_contrast_combos$Cross_sectional,
    timepoints = approach_contrast_combos$Cross_sectional,
    n_top_annotations = n_top_annotations,
    output_dir = output_dir,
    filename_suffix = filename_suffix
  )
}


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
                                        font_size = 12,
                                        sig_cutoff = 0.05)
  }
}


###
# Volcano Plots with Top Points Labeled
###
volcano_plot_subfolders <- c("padj_0.05", "raw_pval", "padj_0.1", "padj_0.2")
for (subfolder in volcano_plot_subfolders) {
  ensure_dir(file.path(output_dir, volcano_subfolder, subfolder))
}
for (contrast in contrasts) {
  temp <- lin_regr_results_list[[contrast]]
  # Rename FDR column to p-adj for volcano plot function
  colnames(temp)[colnames(temp) == "FDR"] <- "p-adj"
  generate_volcano_plot(
    data = temp,
    pval_col = "p-adj",
    beta_col = "beta",
    contrast = contrast,
    p_thresh = 0.05,
    beta_thresh = 0.0,
    output_dir = file.path(output_dir, volcano_subfolder, "padj_0.05"),
    feature_type = "Protein",
    label_col = "SYMBOL",
    filename_suffix = filename_suffix,
    n_label_limit = 30,
    font_size = font_size
  )
}

# Create versions with -log10(raw p-value) on y-axis
for (contrast in contrasts) {
  temp <- lin_regr_results_list[[contrast]]
  generate_volcano_plot(
    data = temp,
    pval_col = "p",
    beta_col = "beta",
    contrast = contrast,
    p_thresh = 0.05,
    beta_thresh = 0.0,
    output_dir = file.path(output_dir, volcano_subfolder, "raw_pval"),
    feature_type = "Protein",
    label_col = "SYMBOL",
    filename_suffix = paste0(filename_suffix, "_raw_pval"),
    font_size = font_size,
    n_label_limit = 30
  )
}

# Create versions with -log10(p-adj) on y-axis and cutoff of 0.1
for (contrast in contrasts) {
  temp <- lin_regr_results_list[[contrast]]
  colnames(temp)[colnames(temp) == "FDR"] <- "p-adj"
  generate_volcano_plot(
    data = temp,
    pval_col = "p-adj",
    beta_col = "beta",
    contrast = contrast,
    p_thresh = 0.1,
    beta_thresh = 0.0,
    output_dir = file.path(output_dir, volcano_subfolder, "padj_0.1"),
    feature_type = "Protein",
    label_col = "SYMBOL",
    filename_suffix = paste0(filename_suffix, "_padj_0.1"),
    n_label_limit = 20,
    y_top_padding = 0.4,
    font_size = font_size
  )
}


# Create versions with -log10(p-adj) on y-axis and cutoff of 0.2
for (contrast in contrasts) {
  temp <- lin_regr_results_list[[contrast]]
  colnames(temp)[colnames(temp) == "FDR"] <- "p-adj"
  generate_volcano_plot(
    data = temp,
    pval_col = "p-adj",
    beta_col = "beta",
    contrast = contrast,
    p_thresh = 0.2,
    beta_thresh = 0.0,
    output_dir = file.path(output_dir, volcano_subfolder, "padj_0.2"),
    feature_type = "Protein",
    label_col = "SYMBOL",
    filename_suffix = paste0(filename_suffix, "_padj_0.2"),
    font_size = font_size
  )
}


###
# Pathway-focused Network
###
# Generate plot for specific pathways of interest for baseline cross-sectional contrast
if ("Cross-Sectional" %in% jobs_to_run) {
  pathways_of_interest_baseline_network <- c(
    # Manually curated, select pathways with p-adj<0.1 using main set of job parameters
    # Down
    "INFLAMMATORY RESPONSE",
    "LYMPHOCYTE DIFFERENTIATION",
    "REGULATION OF IMMUNE SYSTEM PROCESS",
    "TISSUE REMODELING",
    "FATTY ACID CATABOLIC PROCESS",
    "LIPID LOCALIZATION",
    # Up
    "AXON DEVELOPMENT",
    "PHASIC SMOOTH MUSCLE CONTRACTION",
    "GLYCEROLIPID METABOLIC PROCESS",
    "DEVELOPMENTAL CELL GROWTH",
    "MUSCLE CELL MIGRATION",
    "TRIGLYCERIDE METABOLIC PROCESS"
  )
  gsea_res_baseline <- enrichment_gsea_gobp_results_list[["Baseline"]]$GSEA_GOBP
  prot_assoc_df_baseline <- lin_regr_results_list[["Baseline"]]
  p_baseline <- plot_pathway_focus_network(gsea_res_baseline,
                                           prot_assoc_df_baseline,
                                           title = "Baseline Cross-Sectional GSEA GO:BP, Network of Selected Sig. Pathways (p-adj<0.1)",
                                           pathway_pval_cutoff = 0.1,
                                           select_pathways = pathways_of_interest_baseline_network)
  if (!is.null(p_baseline)) {
    ggsave(filename = file.path(output_dir, pathway_focus_network_subfolder,
                                paste0("Pathway_Focus_Network_GOBP_Baseline_Selected_Pathways", filename_suffix, ".png")),
          plot = p_baseline, width = 14, height = 10, dpi=600)
  }

  # Also generate a GSEA GO:BP Bar Plot that corresponds to pathways_of_interest_baseline_network
  plot_gsea_combined(gsea_res_baseline,
                     file.path(output_dir, barplot_subfolder),
                     "Baseline",
                     filename_suffix = paste0(filename_suffix, "_Selected_Pathways"),
                     select_pathways = pathways_of_interest_baseline_network,
                     font_size = font_size*1.1)

  # Generate plot for specific pathways of interest for 12m cross-sectional contrast
  pathways_of_interest_12m_network <- c(
    # Down
    "TISSUE REMODELING",
    "CELL ADHESION",
    "RESPONSE TO BACTERIUM",
    # Up
    "SARCOMERE ORGANIZATION",
    "MYOFIBRIL ASSEMBLY"
  )
  gsea_res_12m <- enrichment_gsea_gobp_results_list[["12m"]]$GSEA_GOBP
  prot_assoc_df_12m <- lin_regr_results_list[["12m"]]
  p_12m <- plot_pathway_focus_network(gsea_res_12m,
                                      prot_assoc_df_12m,
                                      title = "12m Cross-Sectional GSEA GO:BP, Network of Selected Sig. Pathways (p-adj<0.05)",
                                      pathway_pval_cutoff = 0.05,
                                      select_pathways = pathways_of_interest_12m_network)
  if (!is.null(p_12m)) {
    ggsave(filename = file.path(output_dir, pathway_focus_network_subfolder,
                                paste0("Pathway_Focus_Network_GOBP_12m_Selected_Pathways", filename_suffix, ".png")),
          plot = p_12m, width = 14, height = 10, dpi=600)
  }

}


}


# Only run when executed directly (not when sourced from another script)
if (sys.nframe() == 0L) {
  run_proteomics_linear_regression_plots(subset_HCT_type = "ALL", jobs_to_run = c("Cross-Sectional", "Double_Delta"))
}