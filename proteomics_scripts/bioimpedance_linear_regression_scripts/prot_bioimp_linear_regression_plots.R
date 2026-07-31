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
run_bioimp_proteomics_linear_regression_plots <- function(subset_HCT_type, outcome_variable_list) {
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
           stringr::str_wrap(name, width = 25), name))

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
      ggrepel::geom_label_repel(
        data               = label_data,
        aes(x = x, y = y, label = display_label, fontface = label_fontface),
        size               = I(label_data$label_size),
        color              = "black",
        fill               = alpha("white", 0.5),
        label.size         = NA,
        force              = 25,  force_pull        = 0.06,
        box.padding        = 1.2, point.padding     = 0.6,
        max.overlaps       = Inf,
        max.iter           = 20000,
        segment.color      = "#2b2b2b", segment.size = 0.3,
        min.segment.length = 0.2,
        seed               = 42,
        lineheight         = 0.85,
        show.legend        = FALSE
      ) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                           midpoint = 0, na.value = "grey70", name = expression(beta)) +
      scale_shape_manual(values = c("Pathway" = 23, "Protein" = 21), name = "Node Type", guide = guide_legend(override.aes = list(size = 5, color = "black", fill = c("darkgoldenrod1", "white")))) +
      scale_size_continuous(name = expression(-log[10](p[adj])), range = c(3, 10)) +
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
enrichment_output_folder <- "prot_bioimp_linear_regression_enrichment"

output_script_folder <- "prot_bioimp_linear_regression_plots"
volcano_subfolder <- "Volcano_Plots"
barplot_subfolder <- "Bar_Plots"

contrasts <- list()
for (var in outcome_variable_list) {
  contrasts[[paste0("Cross_sectional_", var)]] <- paste0("Baseline_", var)
}


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
ensure_dir(file.path(output_dir, volcano_subfolder))
ensure_dir(file.path(output_dir, barplot_subfolder))


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


}



# Only run when executed directly (not when sourced from another script)
if (sys.nframe() == 0L) {
  run_bioimp_proteomics_linear_regression_plots(subset_HCT_type = "ALL", outcome_variable_list = outcome_variable_list)
}