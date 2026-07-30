###
# Proteomics Diff Exp for Extreme Subsets, pathway enrichment
# - Continue from extreme subset differential expression script
# - Perform enrichment analysis only (no plotting)
###

library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(fgsea)
library(clusterProfiler)
library(openxlsx)
source(here::here("common_functions.R"))
# Seed set in common_functions.R via config.R


###
# To run code as a function:
run_proteomics_top_vs_bottom_diff_exp_enrichment <- function(subset_HCT_type = "ALL") {
###


###
# Functions
###
run_pathway_enrichment_diff_exp <- function(res_df, label, mapping,
                                   base_dir, filename_suffix="", regr_type="Linear", lfc_colname = "logFC", pval_colname = "P.Value", padj_colname = "adj.P.Val", ora_p_col = "p", ora_p_cutoff = 0.05) {
  
  if (is.null(res_df)) return(NULL)
  outdir <- file.path(base_dir, label) # create subfolder per contrast
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  message("\n=== PATHWAY ANALYSIS: ", label, " ===")
  
  # Use OlinkID if present; fallback to Protein for backward compatibility
  id_col <- if ("OlinkID" %in% names(res_df)) "OlinkID" else "Protein"
  
  # join annotation columns that are not already present in res_df
  cols_to_add <- setdiff(colnames(mapping), c("OlinkID", colnames(res_df)))
  if (length(cols_to_add) > 0) {
    mapping_subset <- mapping %>% dplyr::select(OlinkID, dplyr::all_of(cols_to_add))
    res_df <- res_df %>% left_join(mapping_subset, by = setNames("OlinkID", id_col))
  }
  
  reg <- res_df %>%
    filter(!is.na(ENTREZID)) %>%
    mutate(
      ENTREZID = as.character(ENTREZID),
      logFC     = as.numeric(.data[[lfc_colname]]),
      p        = as.numeric(.data[[pval_colname]]),
      FDR      = as.numeric(.data[[padj_colname]])
    ) %>%
    group_by(ENTREZID) %>%
    summarise(
      logFC    = mean(logFC, na.rm = TRUE),
      p       = suppressWarnings(min(p, na.rm = TRUE)),
      FDR     = suppressWarnings(min(FDR, na.rm = TRUE)),
      SYMBOL  = dplyr::first(SYMBOL),
      UniProt = dplyr::first(UniProt),
      OlinkID = dplyr::first(OlinkID),
      .groups = "drop"
    ) %>%
    filter(is.finite(logFC), is.finite(p))
  
  if (nrow(reg) < 10) {
    message("  Too few mapped genes – skipping enrichment.")
    return(NULL)
  }
  
  # Add gene descriptions
  message("  Fetching gene descriptions...")
  gene_desc <- get_gene_descriptions(reg$ENTREZID)
  reg <- reg %>%
    mutate(
      Description = gene_desc[ENTREZID],
      Description = ifelse(is.na(Description), "", Description)
    )
  
  reg_raw <- reg  # preserve pre-sort version for return

  # Sort reg by ascending p-value and reorder columns for better readability
  # Note: after summarise(), the p-value column has been renamed to "p"
  reg <- reg %>% 
    arrange(p) %>%
    dplyr::select(SYMBOL, Description, ENTREZID, UniProt, OlinkID, logFC, p, FDR)
  
  # Write Excel file with two sheets: all results and significant results (FDR < 0.1)
  reg_sig <- reg %>% 
    filter(FDR < 0.1) %>%
    mutate(
      logFC = signif(logFC, 2),
      p = signif(p, 2),
      FDR = signif(FDR, 2)
    ) %>%
    dplyr::select(SYMBOL, Description, logFC, p, FDR)
  
  write.xlsx(
    list(
      all = reg,
      sig = reg_sig
    ),
    file = file.path(outdir, paste0(regr_type, "_Protein_Results_Annotated_", label, filename_suffix, ".xlsx"))
  )
  
  # Create ranked gene list for GSEA: sign(logFC) * -log10(p)
  # Positive ranks = upregulated genes with low p-values
  # Negative ranks = downregulated genes with low p-values
  ranks <- sign(reg$logFC) * -log10(pmax(reg$p, .Machine$double.xmin))
  names(ranks) <- reg$ENTREZID
  
  # Verify ranks are valid (no NA, NaN, or Inf values)
  if (any(!is.finite(ranks))) {
    message("  Warning: Some ranks are not finite. Filtering...")
    valid_ranks <- is.finite(ranks)
    ranks <- ranks[valid_ranks]
  }
  
  # Filter out genes with no rank
  if (length(ranks) == 0) {
    message("  No genes with valid ranks – skipping GSEA.")
    return(NULL)
  }
  
  # Create ENTREZID to SYMBOL mapping for leadingEdge conversion
  entrez_to_symbol <- setNames(reg$SYMBOL, reg$ENTREZID)
  
  genes_up   <- reg %>% filter(!!sym(ora_p_col) < ora_p_cutoff, logFC > 0) %>% pull(ENTREZID)
  genes_down <- reg %>% filter(!!sym(ora_p_col) < ora_p_cutoff, logFC < 0) %>% pull(ENTREZID)
  universe   <- reg$ENTREZID
  
  gsea_db <- list(
    Hallmark = hallmark_sets,
    Reactome = react_sets,
    GOBP     = gobp_sets,
    KEGG     = kegg_sets
  )
  
  t2g_db <- list(
    Hallmark = t2g_h,
    Reactome = t2g_re,
    GOBP     = t2g_bp,
    KEGG     = t2g_kegg
  )
  
  all_summaries <- list()
  
  for (db in names(gsea_db)) {
    message("  → ", db)
    
    path_sets <- gsea_db[[db]]
    t2g       <- t2g_db[[db]]
    
    if (length(path_sets) == 0) {
      message("    no gene sets for ", db)
      next
    }
    
    ## ---------- GSEA ----------
    valid_sets <- lapply(path_sets, function(g) intersect(g, names(ranks)))
    valid_sets <- valid_sets[lengths(valid_sets) > 1]
    
    # use fgsea
    if (length(valid_sets) > 0) {
      set.seed(random_seed)
      gsea_res <- try(
        suppressWarnings(
          fgsea::fgseaMultilevel(
            pathways    = valid_sets,
            stats       = ranks,
            minSize     = 5,
            maxSize     = 500,
            nPermSimple = 5000,
            eps         = 1e-10   # Set minimum p-value threshold to prevent underestimation
          )
        ),
        silent = TRUE
      )
      
      if (!inherits(gsea_res, "try-error")) {
        gsea_res <- as.data.frame(gsea_res)
        
        # Filter out pathways with NA log2err (unreliable p-values)
        if ("log2err" %in% colnames(gsea_res)) {
          n_total <- nrow(gsea_res)
          gsea_res <- gsea_res %>% filter(!is.na(log2err))
          n_filtered <- n_total - nrow(gsea_res)
          if (n_filtered > 0) {
            message("    Filtered ", n_filtered, " pathways with unreliable p-values")
          }
        }
        
        # Convert leadingEdge ENTREZ IDs to gene symbols
        if ("leadingEdge" %in% colnames(gsea_res)) {
          gsea_res$leadingEdge <- vapply(
            gsea_res$leadingEdge,
            function(x) {
              symbols <- entrez_to_symbol[x]
              symbols <- symbols[!is.na(symbols)]
              paste(symbols, collapse = ", ")
            },
            FUN.VALUE = character(1)
          )
        }
        
        # Save GSEA result per annotation type (db)
        # Sort by ascending p-value
        gsea_res <- gsea_res %>% arrange(pval)
        # Only save GOBP results to reduce clutter, since GOBP is used in GSEA GOBP bar plotting
        if (db == "GOBP") {
          write.csv(gsea_res,
                    file.path(outdir, paste0("GSEA_", db, "_", label, filename_suffix, ".csv")),
                    row.names = FALSE)
        }
        
        sig <- gsea_res %>%
          mutate(
            NES  = as.numeric(NES),
            pval = as.numeric(pval)
          )
        
        gsea_summary <- sig %>%
          mutate(
            Pathway   = pathway,
            Direction = ifelse(NES > 0, "UP", "DOWN")
          ) %>%
          transmute(
            Source      = "GSEA",
            DB          = db,
            Pathway_raw = pathway,
            Pathway,
            Direction,
            PValue      = pval,
            FDR         = padj,
            Stat        = NES,
            Genes       = leadingEdge
          )
        
        if (db == "GOBP" && nrow(gsea_summary) > 0) {
          gsea_summary <- gsea_summary %>%
            mutate(
              Pathway    = normalize_pathway_name(db, Pathway),
              GOBP_Group = assign_gobp_group(Pathway)
            )
        } else if (nrow(gsea_summary) > 0) {
          gsea_summary <- gsea_summary %>%
            mutate(Pathway = normalize_pathway_name(db, Pathway))
        }
        
        if (nrow(gsea_summary) > 0)
          all_summaries <- append(all_summaries, list(gsea_summary))
      }
    }
    
    ## ---------- ORA ----------
    do_ora <- function(glist, direction) {
      if (length(glist) < 3) return(NULL)
      if (is.null(t2g) || nrow(t2g) == 0) return(NULL)
      
      t2g_use <- t2g %>% filter(gene %in% universe)
      if (nrow(t2g_use) == 0) return(NULL)
      
      # Use enricher from clusterProfiler for ORA
      enr <- try(
        clusterProfiler::enricher(
          gene          = glist,
          TERM2GENE     = t2g_use,
          universe      = universe,
          pAdjustMethod = "BH",
          minGSSize     = 3
        ),
        silent = TRUE
      )
      
      if (inherits(enr, "try-error") || is.null(enr)) return(NULL)
      
      df <- as.data.frame(enr)
      if (nrow(df) == 0) return(NULL)
      df2 <- df  # No p-value cutoff; save all enrichment results
      
      # # Uncomment to save ORA result per annotation type (db) and direction
      # # Sort by ascending p-value
      # df2_sorted <- df2 %>% arrange(pvalue)
      # write.csv(df2_sorted,
      #           file.path(outdir, paste0("ORA_", db, "_", direction, "_", label, filename_suffix, ".csv")),
      #           row.names = FALSE)
      
      out <- df2 %>%
        mutate(
          Pathway = Description,
          Genes = vapply(
            strsplit(geneID, "/"),
            function(x) {
              symbols <- entrez_to_symbol[x]
              symbols <- symbols[!is.na(symbols)]
              paste(symbols, collapse = ", ")
            },
            FUN.VALUE = character(1)
          )
        ) %>%
        transmute(
          Source      = "ORA",
          DB          = db,
          Direction   = direction,
          Pathway_raw = Description,
          Pathway     = normalize_pathway_name(db, Pathway),
          PValue      = pvalue,
          FDR         = p.adjust,
          Stat        = NA_real_,
          Genes       = Genes
        )
      
      if (db == "GOBP" && nrow(out) > 0) {
        out <- out %>%
          mutate(GOBP_Group = assign_gobp_group(Pathway))
      }
      
      out
    }
    
    all_summaries <- append(all_summaries, list(do_ora(genes_up,   "UP")))
    all_summaries <- append(all_summaries, list(do_ora(genes_down, "DOWN")))
  }
  
  ## ---------- Save summaries ----------
  all_summaries <- all_summaries[!sapply(all_summaries, is.null)]
  
  final_raw <- NULL
  if (length(all_summaries) > 0) {
    final_raw <- bind_rows(all_summaries)

    # Sort by ascending PValue
    final <- final_raw %>% arrange(PValue)
    
    # Write Excel file with two sheets: all results and significant results (FDR < 0.1)
    # Filter sig sheet to GSEA GOBP only, rename and format columns
    final_sig <- final %>% 
      filter(FDR < 0.1, Source == "GSEA", DB == "GOBP") %>%
      mutate(
        NES = signif(Stat, 2),
        PValue = signif(PValue, 2),
        FDR = signif(FDR, 2)
      ) %>%
      dplyr::select(Pathway, NES, PValue, FDR)
      
    write.xlsx(
      list(
        all = final,
        sig = final_sig
      ),
      file = file.path(outdir, paste0("Summary_Pathways_", label, filename_suffix, ".xlsx"))
    )
  }
  
  message("✓ Completed contrast: ", label)
  invisible(list(reg = reg_raw, summary = final_raw))
}


###
# Values
###
input_folder <- get_input_folder()
olink_annotations_filename <- "olink_mapped.csv"
msigdb_cache_folder <- file.path(input_folder, "msigdb_annotations")

output_folder <- get_output_folder()
output_reg_results_folder <- "proteomics_top_vs_bottom_diff_exp_limma"
output_script_folder <- "proteomics_top_vs_bottom_enrichment"

summary_list <- list() # pathway summaries per contrast
reg_list     <- list() # regression tables per contrast


###
# Take into account HCT type subsets
###
if (subset_HCT_type == "ALLO" || subset_HCT_type == "AUTO") {
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
# Import Differential Expression Results
###
res_Baseline <- read.xlsx(file.path(output_folder, output_reg_results_folder, "Limma_Results_Top_vs_Bottom_VO2peak_Baseline.xlsx"))
res_6m       <- read.xlsx(file.path(output_folder, output_reg_results_folder, "Limma_Results_Top_vs_Bottom_VO2peak_6m.xlsx"))
res_12m      <- read.xlsx(file.path(output_folder, output_reg_results_folder, "Limma_Results_Top_vs_Bottom_VO2peak_12m.xlsx"))


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
# Run All Pathway Enrichments
###
res_Baseline_out <- run_pathway_enrichment_diff_exp(res_Baseline, "Baseline", olink_mapped, output_dir, filename_suffix = filename_suffix)
res_6m_out       <- run_pathway_enrichment_diff_exp(res_6m,       "6m",       olink_mapped, output_dir, filename_suffix = filename_suffix)
res_12m_out      <- run_pathway_enrichment_diff_exp(res_12m,      "12m",      olink_mapped, output_dir, filename_suffix = filename_suffix)

all_results  <- Filter(Negate(is.null), list(Baseline = res_Baseline_out, `6m` = res_6m_out, `12m` = res_12m_out))
reg_list     <- lapply(all_results, `[[`, "reg")
summary_list <- lapply(all_results, `[[`, "summary")

}


# Only run when executed directly (not when sourced from another script)
if (sys.nframe() == 0L) {
  run_proteomics_top_vs_bottom_diff_exp_enrichment(subset_HCT_type = "ALL")
}
