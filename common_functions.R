library(dplyr)
library(tidyr)
library(ggplot2)
library(fgsea)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(stringr)
library(AnnotationDbi)
library(viridis)
library(ggrepel)
library(patchwork)
library(gridExtra)
library(openxlsx)
library(ggVennDiagram)
library(showtext)
library(sysfonts)

# Load centralised configuration (user-editable values live in config.R)
source(here::here("config.R"))

set.seed(random_seed)

# Load font from Google Fonts and activate showtext for all plots
font_add_google(font_style, font_style)
showtext_auto()
showtext_opts(dpi = 600)

# Override theme_bw so all calls inherit font_style as base_family by default
theme_bw <- function(base_size = 11, base_family = font_style, base_line_size = base_size / 22,
                     base_rect_size = base_size / 22) {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family,
                    base_line_size = base_line_size, base_rect_size = base_rect_size)
}

vo2_label <- function(text) {
  # plotmath label helper for VO2peak
  # - Turn a plain label string into a plotmath expression in which "VO2peak" is
  #   typeset as VO with a "2peak" subscript (or "VO2" as VO with a "2" subscript).
  if (!grepl("VO2", text, fixed = TRUE)) return(text)

  one_line <- function(s) {
    m <- regexpr("VO2(peak)?", s)
    if (m == -1L) return(s)
    hit  <- regmatches(s, m)
    pre  <- substr(s, 1, m - 1)
    post <- substr(s, m + attr(m, "match.length"), nchar(s))
    out  <- if (identical(hit, "VO2peak")) quote(VO[2 * peak]) else quote(VO[2])
    if (nzchar(pre))  out <- bquote(.(pre) * .(out))
    if (nzchar(post)) out <- bquote(.(out) * .(post))
    out
  }

  parts <- lapply(strsplit(text, "\n", fixed = TRUE)[[1]], one_line)
  Reduce(function(a, b) bquote(atop(.(a), .(b))), parts, right = TRUE)
}


###
# Centralised path helpers — use these instead of defining input_folder / output_folder in every script
###
get_input_folder  <- function() here::here(input_dir_name)
get_output_folder <- function() here::here(output_dir_name)

append_report_line <- function(report_file, text = "") {
  write(text, file = report_file, append = TRUE)
}

assign_gobp_group <- function(path_name) {
  nm <- toupper(path_name)
  dplyr::case_when(
    grepl(paste0("IMMUNE|INFLAMM|LEUKOCYTE|LYMPHOCYTE|NEUTROPHIL|",
                 "MACROPHAGE|T CELL|B CELL|INTERLEUKIN|INTERFERON|",
                 "CYTOKINE|IMMUNOGLOBULIN|MONONUCLEAR CELL DIFFERENTIATION|",
                 "CELL ACTIVATION|ACUTE PHASE RESPONSE|DEFENSE"), nm)               ~ "Immune / Inflammation",
    grepl("METABOLI|MITOCHONDR|GLYCOLYS|OXIDAT|LIPID|FATTY ACID|CHOLESTEROL", nm)   ~ "Metabolism",
    grepl("CELL CYCLE|MITOTIC|DIVISION|DNA REPLICATION|PROLIFERAT", nm)             ~ "Cell Cycle / Proliferation",
    grepl("APOPTO|CELL DEATH|NECRO", nm)                                            ~ "Cell Death / Stress",
    grepl("ECM|EXTRACELLULAR MATRIX|ADHESION|INTEGRIN|FOCAL ADHESION", nm)          ~ "ECM / Adhesion",
    grepl(paste0("ANGIOGEN|VASCULATURE|TUBE DEVELOPMENT|MORPHOGENESIS",
                 "|REMODELING|CARTILAGE DEVELOPMENT|DEVELOPMENTAL PROCESS"), nm)    ~ "Morphogenesis / Remodeling",
    grepl("SIGNALING|SIGNALLING", nm)                                               ~ "Signaling",
    grepl("MUSCLE|SARCOMERE|MYOFIBRIL", nm)                                         ~ "Muscle",
    TRUE                                                                            ~ "Other"
  )
}

build_feature_summary <- function(reg_list,
                                         file_path) {
  if (length(reg_list) == 0) return(NULL)
  
  df <- bind_rows(
    lapply(names(reg_list), function(lbl) {
      reg_list[[lbl]] %>% mutate(Contrast = lbl)
    })
  )
  
  write.csv(
    df,
    file.path(file_path),
    row.names = FALSE
  )
}

ensure_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Directory created: ", dir_path)
  }
}

generate_correlation_plot <- function(df, sample_meta, y_var, filename,
                                      id_name = "feature_label",
                                      feature_meta = NULL,
                                      label_col = "compound_id_library",
                                      visit_filter = NULL,
                                      n_top_labels = 20,
                                      dedup_labels = TRUE,
                                      method = "pearson",
                                      text_size = 14, point_size = 2,
                                      width = 10, height = 7, title=NULL,
                                      y_label = NULL) {
  ###
  # df:            Wide-format data frame (features x samples), with an id column and numeric sample columns.
  # sample_meta:   Sample metadata data frame containing 'sample_id', 'visit', and y_var column.
  # y_var:         Name of the column in sample_meta to correlate feature abundance against (e.g. "VO2peak_Baseline").
  # filename:      Full output file path for the saved PNG.
  # id_name:       Column name in df used as the feature identifier.
  # feature_meta:  Optional data frame with feature metadata; must contain id_name and label_col columns.
  #                When provided, label_col values are used as point labels for the top correlates.
  # label_col:     Column in feature_meta to use for labeling (default "compound_id_library").
  #                Falls back to id_name if feature_meta is NULL or label_col is not found.
  # visit_filter:  Optional vector of visit labels to restrict samples (e.g. c("Baseline")). NULL = all samples.
  # n_top_labels:  Number of top features (by |r|) to annotate on the plot.
  # dedup_labels:  If TRUE (default), when multiple features share the same label_col value (e.g. the same gene
  #                SYMBOL measured in multiple Olink panels), only the feature with the highest |r| among them
  #                is labelled. All individual data points are still plotted.
  # method:        Correlation method passed to cor(): "pearson" or "spearman".
  # text_size:     Base font size for the plot.
  # point_size:    Size of points.
  # width/height:  Plot dimensions in inches.
  # title:        Optional plot title (default NULL = no title).
  ###

  # Identify sample columns (everything except the id column)
  sample_cols <- setdiff(colnames(df), id_name)

  # Pivot to long format and join outcome_var values
  meta_cols <- intersect(c("sample_id", "visit", y_var), colnames(sample_meta))
  long_df <- df %>%
    tidyr::pivot_longer(
      cols      = all_of(sample_cols),
      names_to  = "sample_id",
      values_to = "abundance"
    ) %>%
    dplyr::left_join(sample_meta %>% dplyr::select(all_of(meta_cols)), by = "sample_id")

  # Optionally restrict to specific visit(s)
  if (!is.null(visit_filter) && "visit" %in% colnames(long_df)) {
    long_df <- long_df %>% dplyr::filter(visit %in% visit_filter)
  }

  long_df <- long_df %>%
    dplyr::filter(!is.na(.data[[y_var]]), !is.na(abundance), is.finite(abundance))

  if (nrow(long_df) == 0) {
    warning("generate_correlation_plot: no valid rows after filtering for y_var = '", y_var, "'. Skipping.")
    return(invisible(NULL))
  }

  # Compute per-metabolite correlation with y_var
  cor_df <- long_df %>%
    dplyr::group_by(.data[[id_name]]) %>%
    dplyr::summarise(
      r = cor(abundance, .data[[y_var]], method = method, use = "complete.obs"),
      n = sum(!is.na(abundance) & !is.na(.data[[y_var]])),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(r)) %>%
    dplyr::arrange(r) %>%
    dplyr::mutate(rank = dplyr::row_number())

  # Join feature metadata to get informative label column
  use_label_col <- id_name  # fallback
  if (!is.null(feature_meta) && label_col %in% colnames(feature_meta) && id_name %in% colnames(feature_meta)) {
    cor_df <- cor_df %>%
      dplyr::left_join(
        feature_meta %>% dplyr::select(all_of(c(id_name, label_col))),
        by = id_name
      ) %>%
      # Use compound_id_library when available, fall back to feature id
      dplyr::mutate(
        .display_label = dplyr::if_else(
          !is.na(.data[[label_col]]) & trimws(.data[[label_col]]) != "",
          .data[[label_col]],
          .data[[id_name]]
        )
      )
    use_label_col <- ".display_label"
  }

  # Select top features by |r| for labeling, with optional deduplication
  if (dedup_labels && use_label_col != id_name) {
    # For each unique display label, retain only the row with the highest |r|.
    # This prevents the same annotation (e.g. gene SYMBOL across Olink panels) from
    # appearing as multiple labels; all individual points are still plotted.
    top_ids <- cor_df %>%
      dplyr::group_by(.data[[use_label_col]]) %>%
      dplyr::slice_max(abs(r), n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::slice_max(abs(r), n = n_top_labels, with_ties = FALSE) %>%
      dplyr::pull(.data[[id_name]])
    top_idx <- which(cor_df[[id_name]] %in% top_ids)
  } else {
    top_idx <- order(abs(cor_df$r), decreasing = TRUE)[seq_len(min(n_top_labels, nrow(cor_df)))]
  }
  r_range    <- max(abs(cor_df$r), na.rm = TRUE)
  rank_range <- max(cor_df$rank, na.rm = TRUE)
  x_pad      <- rank_range * 0.15
  cor_df <- cor_df %>%
    dplyr::mutate(
      label         = dplyr::if_else(dplyr::row_number() %in% top_idx, .data[[use_label_col]], ""),
      label_nudge_y = dplyr::if_else(label != "", sign(r) * r_range * 0.22, 0),
      label_nudge_x = dplyr::if_else(
        label != "",
        dplyr::if_else(rank > median(rank), rank_range * 0.06, -rank_range * 0.06),
        0
      )
    )

  # Visit label for the plot title
  visit_label <- if (!is.null(visit_filter)) paste0(" (", paste(visit_filter, collapse = "/"), ")") else ""
  plot_title <- if (!is.null(title)) title else paste0("Feature Correlation with ", vo2_label(y_var), visit_label)

  p <- ggplot(cor_df, aes(x = rank, y = r)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
    geom_point(aes(color = r), size = point_size, alpha = 0.85) +
    scale_color_gradient2(
      low      = "#2166AC",
      mid      = "#F7F7F7",
      high     = "#B2182B",
      midpoint = 0,
      name     = paste0(method, " r")
    ) +
    ggrepel::geom_text_repel(
      aes(label = label),
      nudge_y            = cor_df$label_nudge_y,
      nudge_x            = cor_df$label_nudge_x,
      size               = text_size / ggplot2::.pt,
      alpha              = 0.8,
      max.overlaps       = Inf,
      min.segment.length = 0,
      force              = 18,
      force_pull         = 0.2,
      box.padding        = 0.7,
      point.padding      = 0.4,
      max.iter           = 20000,
      segment.size       = 0.5,
      segment.color      = "black",
      segment.alpha      = 0.7,
      bg.color           = "white",
      bg.r               = 0.1,
      xlim               = c(-x_pad, rank_range + x_pad),
      na.rm              = TRUE
    ) +
    labs(
      title    = plot_title,
      subtitle = paste0("n = ", nrow(cor_df), " features | top ", min(n_top_labels, nrow(cor_df)), " labeled by |r|"),
      x        = "Rank (low to high r)",
      y        = if (!is.null(y_label)) y_label else paste0(tools::toTitleCase(method), " r")
    ) +
    theme_bw(base_size = text_size) +
    coord_cartesian(xlim = c(-x_pad, rank_range + x_pad), clip = "off") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title    = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "black"),
      # set font size of axis titles and text and tick labels
      axis.title.x = element_text(size = text_size, color = "black"),
      axis.title.y = element_text(size = text_size, color = "black"),
      axis.text.x = element_text(size = text_size, color = "black"),
      axis.text.y = element_text(size = text_size, color = "black"),
      legend.position = "right"
    )

  ensure_dir(dirname(filename))
  ggsave(filename, plot = p, width = width, height = height, dpi = 600)
  message("Correlation plot saved: ", filename)

  # Save correlation coefficient results to csv in the same directory as the plot
  cor_output_file <- file.path(dirname(filename), paste0(tools::file_path_sans_ext(basename(filename)), "_correlations.csv"))
  write.csv(cor_df, cor_output_file, row.names = FALSE)
  message("Correlation coefficients saved: ", cor_output_file)

  return(invisible(cor_df))
}


generate_pca_plot <- function(data_df, sample_meta, color_by, filename,
                              id_colname = "OlinkID", text_size = 14,
                              point_size = 3, width = 8, height = 6,
                              dramatic_color_scale = FALSE,
                              colors = NULL, legend_title = NULL) {
  # Transpose so rows = samples, columns = proteins
  pca_data <- data_df %>%
    dplyr::select(-all_of(id_colname)) %>%
    as.matrix() %>%
    t()

  pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)

  pca_df <- data.frame(
    sample_id = rownames(pca_result$x),
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2]
  ) %>%
    dplyr::left_join(sample_meta, by = "sample_id")

  if (!color_by %in% colnames(pca_df)) {
    stop("color_by not found in PCA metadata: ", color_by)
  }
  pc1_var <- round((pca_result$sdev[1]^2 / sum(pca_result$sdev^2)) * 100, 1)
  pc2_var <- round((pca_result$sdev[2]^2 / sum(pca_result$sdev^2)) * 100, 1)

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[color_by]])) +
    geom_point(size = point_size, alpha = 0.8) +
    labs(x = paste0("Principal Component 1 (", pc1_var, "%)"), y = paste0("Principal Component 2 (", pc2_var, "%)"), color = if (!is.null(legend_title)) legend_title else color_by) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.title.x = element_text(size = text_size, margin = margin(t = 15)),
      axis.title.y = element_text(size = text_size, margin = margin(r = 15)),
      axis.text.x = element_text(size = text_size),
      axis.text.y = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      legend.title = element_text(size = text_size),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  # Add dramatic color scale for continuous variables
  if (is.numeric(pca_df[[color_by]]) && dramatic_color_scale) {
    p <- p + scale_color_gradient2(
      low = "#2166AC",      # Dark blue
      mid = "#F7F7F7",      # Light gray/white
      high = "#B2182B",     # Dark red
      midpoint = median(pca_df[[color_by]], na.rm = TRUE),
      na.value = "gray50"
    )
  }

  # Add custom colors if provided (for categorical variables)
  if (!is.null(colors) && is.character(colors)) {
    p <- p + scale_color_manual(values = colors)
  }

  ggsave(filename, plot = p, width = width, height = height, dpi = 600)
  message("PCA plot saved: ", filename)
}

generate_sig_venn_diagram <- function(res_list, title, output_dir, filename_suffix = "", beta_dir = "any", pval_col = "FDR", pval_cutoff = 0.05, protein_mapping = NULL, font_size = 14, id_col = "OlinkID") {
  # Filter each contrast for significant features by FDR cutoff and optional beta direction,
  # then (for proteins only!) deduplicate to unique gene SYMBOLs (keep lowest p per SYMBOL when multiple OlinkIDs map to same gene)
  use_symbols <- !is.null(protein_mapping) && "SYMBOL" %in% colnames(protein_mapping)

  sig_df_list <- lapply(res_list, function(res) {
    if (is.null(res) || nrow(res) == 0) return(NULL)

    # Step 1: filter by FDR only (determines Venn region membership, same for all beta_dir variants)
    filtered <- res %>% filter(!!rlang::sym(pval_col) < pval_cutoff & !is.na(!!rlang::sym(pval_col)))

    if (use_symbols) {
      filtered <- filtered %>%
        left_join(protein_mapping %>% dplyr::select(OlinkID, SYMBOL) %>% distinct(), by = "OlinkID") %>%
        filter(!is.na(SYMBOL) & SYMBOL != "") %>%
        group_by(SYMBOL) %>%
        arrange(p, .by_group = TRUE) %>%
        dplyr::filter(dplyr::row_number() == 1L) %>%
        ungroup()
    }

    # Step 2: apply beta direction filter AFTER deduplication, so that Venn region
    # membership (which timepoints a gene is significant at) is identical across all
    # three diagrams. Without this, a gene with a sign flip across timepoints (e.g.
    # negative at Baseline, positive at 6m) would appear in the Baseline+6m intersection
    # in "any" but in the Baseline-only region in "negative" — making the diagrams
    # geometrically incomparable and producing unexplained count discrepancies.
    if (beta_dir == "positive") filtered <- filtered %>% filter(beta > 0)
    else if (beta_dir == "negative") filtered <- filtered %>% filter(beta < 0)

    filtered
  })
  sig_lists <- lapply(sig_df_list, function(df) {
    if (is.null(df) || nrow(df) == 0) return(character(0))
    if (use_symbols && "SYMBOL" %in% colnames(df)) unique(df$SYMBOL)
    else unique(df[[id_col]])
  })

  # Skip if no significant features found in any set
  if (all(sapply(sig_lists, length) == 0)) {
    message("No significant features (", pval_col, " < ", pval_cutoff, ", beta_dir = ", beta_dir, ") found for: ", title)
    return(invisible(NULL))
  }

  beta_label <- switch(beta_dir,
    "any"      = "Any Beta",
    "positive" = "Positive Beta",
    "negative" = "Negative Beta"
  )
  plot_title <- if (is.null(title)) {
    paste0(pval_col, " < ", pval_cutoff, " | ", beta_label)
  } else {
    paste0(title, "\n", pval_col, " < ", pval_cutoff, " | ", beta_label)
  }

  # Build Venn diagram: all regions white; all-overlap region is light green with bold count
  # Note: ggVennDiagram >= 1.4 dropped sf dependency; accessors return plain tibbles with X/Y columns.
  n_sets <- length(sig_lists)
  venn_obj       <- ggVennDiagram::Venn(sig_lists)
  venn_data      <- ggVennDiagram::process_data(venn_obj)
  region_edge_df <- ggVennDiagram::venn_regionedge(venn_data)
  setedge_df     <- ggVennDiagram::venn_setedge(venn_data)
  setlabel_df    <- ggVennDiagram::venn_setlabel(venn_data)
  regionlabel_df <- ggVennDiagram::venn_regionlabel(venn_data)

  # All-overlap region id: "1/2" for 2 sets, "1/2/3" for 3 sets, etc.
  all_overlap_id <- paste(seq_len(n_sets), collapse = "/")

  region_edge_df <- region_edge_df %>%
    dplyr::mutate(fill_color = ifelse(id == all_overlap_id, "#a4caa6", "white"))

  regionlabel_df <- regionlabel_df %>%
    dplyr::mutate(label_face = ifelse(id == all_overlap_id, "bold", "plain"))

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(ggplot2::aes(x = X, y = Y, fill = fill_color, group = id),
                          data = region_edge_df, show.legend = FALSE) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_path(ggplot2::aes(x = X, y = Y, group = id),
                       data = setedge_df, color = "black", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(x = X, y = Y, label = name),
                       data = setlabel_df, size = font_size / ggplot2::.pt) +
    ggplot2::geom_text(ggplot2::aes(x = X, y = Y, label = count, fontface = label_face),
                       data = regionlabel_df, size = font_size / ggplot2::.pt) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(hjust = 0.5, size = font_size + 2, face = "bold"),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin      = ggplot2::margin(t = 20, r = 90, b = 20, l = 90, unit = "pt"),
      legend.position  = "none"
    )

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(plot_title)
  }

  # Save plot
  beta_file_label <- switch(beta_dir,
    "any"      = "AnyBeta",
    "positive" = "PosBeta",
    "negative" = "NegBeta"
  )
  filename_tag <- gsub("[^A-Za-z0-9_]", "_", if (is.null(title)) paste(names(res_list), collapse = "_") else title)
  out_file <- file.path(output_dir, paste0("VennDiagram_", filename_tag, "_", beta_file_label, filename_suffix, ".png"))
  ggplot2::ggsave(out_file, p, width = 7, height = 7, dpi = 600, bg = "white")
  message("Saved Venn diagram: ", basename(out_file))

  # Generate Excel summary of significant proteins
  all_sig_ids <- unique(unlist(sig_lists))  # SYMBOLs if use_symbols, else OlinkIDs
  if (length(all_sig_ids) > 0) {
    # Build annotation table
    if (use_symbols) {
      annot <- protein_mapping %>%
        filter(SYMBOL %in% all_sig_ids) %>%
        dplyr::select(any_of(c("SYMBOL", "OlinkID", "UniProt", "ENTREZID"))) %>%
        group_by(SYMBOL) %>%
        dplyr::filter(dplyr::row_number() == 1L) %>%
        ungroup() %>%
        distinct()

      # Fetch gene descriptions via ENTREZID
      if ("ENTREZID" %in% colnames(annot) && any(!is.na(annot$ENTREZID))) {
        gene_desc <- get_gene_descriptions(as.character(annot$ENTREZID))
        annot <- annot %>%
          mutate(
            ENTREZID    = as.character(ENTREZID),
            Description = gene_desc[ENTREZID],
            Description = ifelse(is.na(Description), "", Description)
          ) %>%
          dplyr::select(any_of(c("SYMBOL", "Description", "OlinkID", "UniProt", "ENTREZID")))
      }
    } else {
      annot <- tibble(x = all_sig_ids)
      names(annot) <- id_col
    }

    # Append per-contrast sig flag + beta/p/FDR columns (SYMBOL-deduplicated values)
    join_col <- if (use_symbols) "SYMBOL" else id_col
    for (cname in names(sig_df_list)) {
      df <- sig_df_list[[cname]]
      sig_ids_this <- sig_lists[[cname]]

      if (!is.null(df) && nrow(df) > 0) {
        id_col_vals <- if (use_symbols && "SYMBOL" %in% colnames(df)) df$SYMBOL else df[[id_col]]
        contrast_cols <- df %>%
          mutate(.join_id = id_col_vals) %>%
          filter(.join_id %in% all_sig_ids) %>%
          dplyr::select(.join_id, beta, p, FDR) %>%
          dplyr::rename(
            !!join_col                    := .join_id,
            !!paste0(cname, "_beta") := beta,
            !!paste0(cname, "_p")    := p,
            !!paste0(cname, "_FDR")  := FDR
          ) %>%
          mutate(!!paste0(cname, "_sig") := .data[[join_col]] %in% sig_ids_this)
        annot <- annot %>% left_join(contrast_cols, by = join_col)
      } else {
        annot <- annot %>%
          mutate(
            !!paste0(cname, "_beta") := NA_real_,
            !!paste0(cname, "_p")    := NA_real_,
            !!paste0(cname, "_FDR")  := NA_real_,
            !!paste0(cname, "_sig")  := FALSE
          )
      }
    }

    # Sort: features significant in more contrasts first
    sig_cols <- paste0(names(res_list), "_sig")
    valid_sig_cols <- intersect(sig_cols, colnames(annot))
    if (length(valid_sig_cols) > 0) {
      annot <- annot %>%
        mutate(n_sig_contrasts = rowSums(dplyr::select(., all_of(valid_sig_cols)), na.rm = TRUE)) %>%
        arrange(desc(n_sig_contrasts)) %>%
        dplyr::select(-n_sig_contrasts)
    }

    # Second sheet: only features significant in ALL contrasts
    if (length(valid_sig_cols) > 0) {
      annot_sig_all <- annot %>%
        filter(rowSums(dplyr::select(., all_of(valid_sig_cols)), na.rm = TRUE) == length(valid_sig_cols))
    } else {
      annot_sig_all <- annot
    }

    xlsx_file <- file.path(output_dir, paste0("VennSummary_", filename_tag, "_", beta_file_label, filename_suffix, ".xlsx"))
    write.xlsx(list("all_sig" = annot, "sig_all_contrasts" = annot_sig_all), xlsx_file, rowNames = FALSE)
    message("Saved Venn summary Excel: ", basename(xlsx_file))
  }

  invisible(p)
}


generate_volcano_plot <- function(data, pval_col, beta_col, contrast,
                         p_thresh=0.1, beta_thresh=0.0, output_dir="output", feature_type="Protein", font_size=20, line_width=3, tick_length=10, point_size=4, label_col=NULL, n_label_limit = 50, filename_suffix="", y_top_padding=0.1, y_max_cap=NULL) {
    ###
    # Generate a volcano plot for feature association analysis.
    # Parameters:
    # - data: DataFrame containing the association data
    # - pval_col: Column name for p-values
    # - beta_col: Column name for beta coefficients
    # - contrast: String for contrast name
    # - p_thresh: Threshold for p-value significance
    # - beta_thresh: Threshold for beta coefficient significance
    # - output_dir: Directory to save the plot
    # - feature_type: Type of feature being plotted (e.g., "Gene", "Protein", "Compound")
    # - font_size: Font size for plot text
    # - line_width: Line width for plot
    # - tick_length: Tick length for plot
    # - point_size: Size of points in the plot
    # - label_col: Column name for labels (e.g., "SYMBOL"). If NULL, no labels are added.
    # - n_label_limit: Maximum number of labels to display
    # - y_max_cap: Optional numeric value to cap the y-axis at for visualization
    ###
    
    # Calculate -log10 p-values, with optional y-axis capping
    data <- data %>%
      mutate(
        neg_log10_p = -log10(pmax(.data[[pval_col]], .Machine$double.xmin)),
        beta_val = .data[[beta_col]]
      )
    
    if (!is.null(y_max_cap)) {
      data <- data %>%
        mutate(neg_log10_p = pmin(neg_log10_p, y_max_cap))
    }
    
    data <- data %>%
      filter(is.finite(neg_log10_p), is.finite(beta_val))
    
    if (nrow(data) == 0) {
      message("No valid data points for volcano plot")
      return(NULL)
    }

    pval_display_label <- switch(pval_col,
      "FDR"  = "p-adj",
      "padj" = "p-adj",
      pval_col
    )
    y_axis_label <- paste0("-log10(", pval_display_label, ")")

    # Create data point categories
    data <- data %>%
      mutate(
        Significance = case_when(
          .data[[pval_col]] < p_thresh & abs(beta_val) > beta_thresh & beta_val > 0 ~ "Up",
          .data[[pval_col]] < p_thresh & abs(beta_val) > beta_thresh & beta_val < 0 ~ "Down",
          TRUE ~ "Not Sig."
        ),
        Significance = factor(Significance, 
                            levels = c("Up", "Down", "Not Sig."))
      )
    
    # Count features in each category
    counts <- data %>%
      dplyr::count(Significance) %>%
      arrange(Significance)
    
    n_sig_up <- counts %>% filter(Significance == "Up") %>% pull(n) %>% sum()
    n_sig_down <- counts %>% filter(Significance == "Down") %>% pull(n) %>% sum()
    n_not_sig <- counts %>% filter(Significance == "Not Sig.") %>% pull(n) %>% sum()
    
    if (length(n_sig_up) == 0) n_sig_up <- 0
    if (length(n_sig_down) == 0) n_sig_down <- 0
    if (length(n_not_sig) == 0) n_not_sig <- 0
    
    # Define colors
    colors <- c(
      "Up" = "#B2182B",
      "Down" = "#2166AC", 
      "Not Sig." = "lightgray"
    )
    
    # Plot points by category
    title <- paste0(feature_type, " Volcano Plot - ", contrast)
    p <- ggplot(data, aes(x = beta_val, y = neg_log10_p, color = Significance)) +
      geom_point(alpha = 0.6, size = point_size) +
      scale_color_manual(values = colors, drop = FALSE) +
      scale_x_continuous(expand = expansion(mult = c(0.08, 0.18))) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
      
      # Add threshold lines
      geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", 
                 color = "black", linewidth = line_width / 3) +
      geom_vline(xintercept = c(-beta_thresh, beta_thresh), linetype = "dashed",
                 color = "black", linewidth = line_width / 3) +
      coord_cartesian(clip = "off") +
      
      # Customize Plot: labels
      labs(
        title = title,
        x = "Beta",
        y = paste0("-log10(", pval_col, ")"),
        color = NULL
      ) +
      
      # Customize Plot: theme
      theme_bw(base_size = font_size) +
      theme(
        plot.title = element_text(hjust = 0.5, size = font_size),
        axis.title = element_text(size = font_size),
        axis.text = element_text(size = font_size),
        legend.title = element_text(size = font_size),
        legend.text = element_text(size = font_size),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = line_width / 2),
        axis.ticks = element_line(linewidth = line_width / 3),
        axis.ticks.length = unit(tick_length, "pt"),
        plot.margin = margin(t = 12, r = 36, b = 12, l = 12)
      )
    
    # Add labels for top n significant points by pval_col if label_col is provided
    if (!is.null(label_col) && label_col %in% colnames(data)) {
      data_sig <- data %>%
        filter(Significance != "Not Sig.") %>%
        arrange(.data[[pval_col]]) %>%  # Sort by p-value
        slice_head(n = n_label_limit) %>%  # Select top n_label_limit
        mutate(label = .data[[label_col]])
      
      if (nrow(data_sig) > 0) {
        p <- p + ggrepel::geom_text_repel(
          data = data_sig,
          aes(label = label),
          color = "black",
          size = font_size * 0.3,
          box.padding = 0.7,
          point.padding = 0.5,
          max.overlaps = nrow(data_sig),
          force = 2,
          force_pull = 0.2,
          min.segment.length = 0,
          max.iter = 20000,
          max.time = 2,
          seed = 42,
          segment.color = "gray50",
          segment.size = 0.3,
          bg.color = "white",
          bg.r = 0.15
        )
      }
    }
    
    # Save the plot
    output_file <- file.path(output_dir, paste0(feature_type, "_Volcano_Plot_", contrast, filename_suffix,
                                                ".png"))
    
    ggsave(
      filename = output_file,
      plot = p,
      width = 12,
      height = 10,
      dpi = 600
    )
    
    return(p)
}

get_gene_descriptions <- function(entrez_ids) {
  # Try to get gene descriptions from org.Hs.eg.db
  descriptions <- tryCatch({
    # Get gene names (full names/descriptions)
    gene_names <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = as.character(entrez_ids),
      columns = c("GENENAME"),
      keytype = "ENTREZID"
    )
    
    # Create a named vector for easy matching
    desc_vec <- gene_names$GENENAME
    names(desc_vec) <- gene_names$ENTREZID
    desc_vec
  }, error = function(e) {
    message("Warning: Could not fetch gene descriptions: ", e$message)
    return(NULL)
  })
  
  if (!is.null(descriptions)) {
    return(descriptions)
  } else {
    # Return empty descriptions if fetch failed
    empty_desc <- rep(NA_character_, length(entrez_ids))
    names(empty_desc) <- as.character(entrez_ids)
    return(empty_desc)
  }
}

get_msigdb_set <- function(collection, subcollection = NULL, cache_folder) {
  # Create filename
  if (is.null(subcollection)) {
    cache_file <- file.path(cache_folder, paste0("msigdb_", collection, ".rds"))
    set_name <- paste(collection)
  } else {
    safe_subname <- gsub(":", "_", subcollection)
    cache_file <- file.path(cache_folder, paste0("msigdb_", collection, "_", safe_subname, ".rds"))
    set_name <- paste(collection, subcollection, sep = "_")
  }
  
  # Check if file exists
  if (file.exists(cache_file)) {
    message("Loading ", set_name, " from cache: ", cache_file)
    msig_data <- readRDS(cache_file)
  } else {
    message("Downloading ", set_name, " from MSigDB...")
    if (is.null(subcollection)) {
      msig_data <- safe_msig(species = "Homo sapiens", collection = collection)
    } else {
      msig_data <- safe_msig(species = "Homo sapiens", collection = collection, 
                            subcollection = subcollection)
    }
    
    if (!is.null(msig_data)) {
      message("Saving ", set_name, " to cache: ", cache_file)
      saveRDS(msig_data, cache_file)
    } else {
      warning("Failed to download ", set_name)
    }
  }
  
  return(msig_data)
}

impute_numeric_mean <- function(x) {
  if (all(is.na(x))) return(x)
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}

impute_factor_mode <- function(x) {
  if (all(is.na(x))) return(x)
  tab <- table(x, useNA = "no")
  if (length(tab) == 0) return(x)
  x[is.na(x)] <- names(which.max(tab))
  x
}

make_t2g <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    tibble(term = character(), gene = character())
  } else {
    # Handle different possible column names in msigdbr versions
    gene_col <- if ("ncbi_gene" %in% names(df)) {
      "ncbi_gene"
    } else if ("entrez_gene" %in% names(df)) {
      "entrez_gene"
    } else {
      stop("Cannot find entrez gene column in msigdbr output. Available columns: ", 
           paste(names(df), collapse = ", "))
    }
    df %>% transmute(term = gs_name, gene = as.character(.data[[gene_col]]))
  }
}

msig_to_list <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(list())
  # Handle different possible column names in msigdbr versions
  gene_col <- if ("ncbi_gene" %in% names(df)) {
    "ncbi_gene"
  } else if ("entrez_gene" %in% names(df)) {
    "entrez_gene"
  } else {
    stop("Cannot find entrez gene column in msigdbr output. Available columns: ", 
         paste(names(df), collapse = ", "))
  }
  tmp <- df %>% transmute(term = gs_name, gene = as.character(.data[[gene_col]]))
  split(tmp$gene, tmp$term)
}

normalize_pathway_name <- function(DB, Pathway) {
  x <- toupper(Pathway)
  x <- gsub("^HALLMARK[_ ]", "", x)
  x <- gsub("^REACTOME[_ ]", "", x)
  x <- gsub("^GO[_ ]*BP[_ ]", "", x)  # remove GOBP prefix
  x <- gsub("^KEGG[_ ]", "", x)
  x <- gsub(" SIGNALING PATHWAY$", "", x)
  x <- gsub(" PATHWAY$", "", x)
  x <- gsub(" PATHWAYS$", "", x)
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

normalize_pathway_name_publish <- function(DB, Pathway) {
  # Upper-case first letter and manually defined acronyms/phrases.
  # Otherwise, lower-case all letters.

  # --- Edit this vector to add/remove always-uppercase terms ---
  uppercase_terms <- c("IL", "IFN", "TNF", "TGF", "VEGF", "MAPK",
                       "PI3K", "AKT", "MTOR", "JAK", "STAT",
                       "DNA", "RNA", "MAPK", "TOR")
  lowercase_terms <- c("and", "of", "in", "on", "with", "by", "for", "to", "from", "at", "as", "via")
  phrases_to_substitute <- c(
    "INTERLEUKIN 6" = "Interleukin-6",
    "INTERLEUKIN 1" = "Interleukin-1",
    "INTERLEUKIN 2" = "Interleukin-2",
    "INTERLEUKIN 10" = "Interleukin-10",
    "INTERLEUKIN 12" = "Interleukin-12",
    "INTERLEUKIN 17" = "Interleukin-17",
    "CELL CELL" = "Cell-Cell",
    "INSULIN LIKE" = "Insulin-like",
    " MEDIATED" = "-Mediated"
  )

  x <- toupper(Pathway)
  x <- gsub("^HALLMARK[_ ]", "", x)
  x <- gsub("^REACTOME[_ ]", "", x)
  x <- gsub("^GO[_ ]*BP[_ ]", "", x)  # remove GOBP prefix
  x <- gsub("^KEGG[_ ]", "", x)
  x <- gsub(" SIGNALING PATHWAY$", "", x)
  x <- gsub(" PATHWAY$", "", x)
  x <- gsub(" PATHWAYS$", "", x)
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  x <- tolower(x)
  # Restore uppercase for each listed term's words (space-separated)
  pattern_upper <- paste0("\\b(", paste(uppercase_terms, collapse = "|"), ")\\b")
  x <- gsub(pattern_upper, "\\U\\1", x, perl = TRUE, ignore.case = TRUE)
  # Capitalize the first letter of each word in the result
  x <- gsub("\\b([a-z])", "\\U\\1", x, perl = TRUE)
  # Handle phrases to substitute
  for (phrase in names(phrases_to_substitute)) {
    x <- gsub(phrase, phrases_to_substitute[[phrase]], x, ignore.case = TRUE)
  }
  # Ensure lowercase for common short words (case-insensitive match after title-casing)
  pattern_lower <- paste0("\\b(", paste(lowercase_terms, collapse = "|"), ")\\b")
  x <- gsub(pattern_lower, "\\L\\1", x, perl = TRUE, ignore.case = TRUE)
  # Always capitalize the first character of the full string
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

plot_beta_distribution <- function(reg_result, output_dir, contrast_label, filename_suffix="") {
  # Plot distribution of beta coefficients
  p <- ggplot(reg_result, aes(x = beta)) +
    geom_histogram(bins = 30, fill = "#2166AC", color = "black") +
    theme_minimal() +
    labs(x = "\u03B2", y = "Protein Count")
  
  ggsave(file.path(output_dir, paste0("Beta_Distribution_", contrast_label, filename_suffix, ".png")), p, width = 6, height = 4)
}

plot_distribution_histogram <- function(data, title, xlab, filename,
                                       text_size = 14,
                                       fill_color = "#ffffff",
                                       border_color = "#525252",
                                       width = 8,
                                       height = 6, bin_num = 100, vertical_cutoff_labels = NULL,
                                       region_colors = NULL, y_label = "Count") {
  # Create histogram plot with ggplot2
  # Remove NA values before plotting
  data <- data[!is.na(data)]

  if (length(data) == 0) {
    stop("No valid data to plot after removing NAs")
  }

  use_region_colors <- !is.null(region_colors) &&
    !is.null(vertical_cutoff_labels) &&
    nrow(vertical_cutoff_labels) >= 2

  if (use_region_colors) {
    sorted_cuts <- sort(vertical_cutoff_labels$cutoff)
    lower_cut   <- sorted_cuts[1]
    upper_cut   <- sorted_cuts[length(sorted_cuts)]

    hist_obj   <- hist(data, breaks = bin_num, plot = FALSE)
    bin_widths <- diff(hist_obj$breaks)
    bin_df <- data.frame(
      mid    = hist_obj$mids,
      count  = hist_obj$counts,
      bwidth = bin_widths,
      region = dplyr::case_when(
        hist_obj$mids < lower_cut ~ "below",
        hist_obj$mids > upper_cut ~ "above",
        TRUE                      ~ "between"
      )
    )
    bin_df$region <- factor(bin_df$region, levels = c("below", "between", "above"))

    rc <- if (!is.null(names(region_colors))) region_colors else
      c("below" = region_colors[1], "between" = region_colors[2], "above" = region_colors[3])

    p <- ggplot(bin_df, aes(x = mid, y = count, fill = region, width = bwidth)) +
      geom_col(color = border_color, alpha = 0.8) +
      scale_fill_manual(values = rc, guide = "none") +
      labs(x = xlab, y = y_label)
  } else {
    p <- ggplot(data.frame(value = data), aes(x = value)) +
      geom_histogram(bins = bin_num, fill = fill_color, color = border_color, alpha = 0.8) +
      labs(x = xlab, y = y_label)
  }

  p <- p +
    theme_bw(base_size = text_size) +
    theme(
      plot.title       = element_blank(),
      axis.title.x     = element_text(size = text_size, margin = margin(t = 15)),
      axis.title.y     = element_text(size = text_size, margin = margin(r = 15)),
      axis.text.x      = element_text(size = text_size, color = "black"),
      axis.text.y      = element_text(size = text_size, color = "black"),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.ticks        = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(5, "pt"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(20, 20, 20, 20)
    )

  # Add dashed vertical lines for cutoffs if provided
  if (!is.null(vertical_cutoff_labels) && nrow(vertical_cutoff_labels) > 0) {
    p <- p + geom_vline(
      data = vertical_cutoff_labels,
      aes(xintercept = cutoff),
      linetype = "dashed", color = "#000000", linewidth = 1.5
    )
  }

  # Save plot
  ggsave(filename, plot = p, width = width, height = height, dpi = 600)
  message("Plot saved: ", filename)
}

plot_gsea_combined <- function(gsea_sig, outdir, label, font_size=20, max_num_pathways=20, filename_suffix="", select_pathways=NULL) {
  if (is.null(gsea_sig) || nrow(gsea_sig) == 0) return(NULL)
  
  # Handle both 'padj' (from fgsea) and 'FDR' column names
  fdr_col <- if ("padj" %in% names(gsea_sig)) "padj" else if ("FDR" %in% names(gsea_sig)) "FDR" else NULL
  pval_col <- if ("pval" %in% names(gsea_sig)) "pval" else if ("PValue" %in% names(gsea_sig)) "PValue" else NULL
  pathway_col <- if ("pathway" %in% names(gsea_sig)) "pathway" else if ("Pathway" %in% names(gsea_sig)) "Pathway" else NULL
  nes_col <- if ("NES" %in% names(gsea_sig)) "NES" else if ("Stat" %in% names(gsea_sig)) "Stat" else NULL
  
  if (is.null(fdr_col) || is.null(pval_col) || is.null(pathway_col) || is.null(nes_col)) {
    message("  [GSEA GOBP] Missing required columns for ", label)
    return(NULL)
  }
  
  # Filter by FDR < 0.1 first
  df_filtered <- gsea_sig %>%
    filter(.data[[fdr_col]] < 0.1)
  
  # Check if any pathways pass FDR threshold
  if (nrow(df_filtered) == 0) {
    message("  [GSEA GOBP] No pathways with FDR < 0.1 for ", label)
    return(NULL)
  }
  
  df <- df_filtered %>%
    mutate(
      PathwayClean = normalize_pathway_name_publish("GOBP", .data[[pathway_col]]),
      GOBP_Group   = assign_gobp_group(PathwayClean),
      NES          = as.numeric(.data[[nes_col]]),
      PValue       = .data[[pval_col]],
      FDR          = .data[[fdr_col]]
    ) %>%
    filter(!is.na(GOBP_Group), !is.na(NES), is.finite(NES)) %>%  # Filter out missing/invalid values
    group_by(PathwayClean) %>%
    summarise(
      PValue     = min(PValue, na.rm = TRUE),
      FDR        = min(FDR, na.rm = TRUE),
      NES        = NES[which.min(PValue)],
      GOBP_Group = dplyr::first(GOBP_Group),
      .groups = "drop"
    ) %>%
    filter(is.finite(PValue), is.finite(NES)) %>%  # Remove any Inf/NA values from aggregation
    arrange(PValue)

  if (nrow(df) == 0) {
    message("  [GSEA GOBP] No pathways remain after GOBP group/NES filtering for ", label)
    return(NULL)
  }

  if (!is.null(select_pathways)) {
    df <- df %>% filter(PathwayClean %in% select_pathways)
    if (nrow(df) == 0) {
      message("  [GSEA GOBP] No selected pathways found in results for ", label)
      return(NULL)
    }
  } else {
    df <- df %>% slice_head(n = max_num_pathways) # Limit to top max_num_pathways # of pathways
  }
  
  df$PathwayClean <- factor(df$PathwayClean,
                            levels = df$PathwayClean[order(df$NES)])
  
  # dynamic plot height — allow as short as 3in for few pathways
  n_paths <- nrow(df)
  plot_height <- max(3, min(24, n_paths * 0.65))
  
  # Plot GSEA bar plot
  p <- ggplot(df,
              aes(x = PathwayClean, y = NES, fill = GOBP_Group)) +
    geom_col(width = 0.55) +
    coord_flip(expand = FALSE, clip = "off") +
    scale_fill_manual(values = gobp_palette, drop = FALSE) +
    scale_y_continuous(breaks = function(x) { b <- pretty(x, n = 5); b[b == floor(b)] }) +
    
    theme_bw(base_size = font_size) +
    theme(
      plot.margin  = margin(15, 20, 15, 40),
      axis.text.y  = element_text(size = font_size, color="black"),
      axis.text.x  = element_text(size = font_size, color="black"),
      legend.text  = element_text(size = font_size, color="black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    
    scale_x_discrete(
      labels = function(x) stringr::str_wrap(x, width = 60)
    ) +
    
    labs(
      x     = NULL,
      y     = "NES",
      fill  = "GO:BP Category"
    )
  
  # Calculate width: label area (font-scaled) + fixed panel area for bars/axes.
  # str_wrap caps displayed label width at 60 chars, so use that as the ceiling.
  effective_label_length <- min(max(nchar(as.character(df$PathwayClean)), na.rm = TRUE), 60)
  label_area_width <- effective_label_length * (font_size / 72) * 0.6  # char width ≈ 0.6× font height
  panel_area_width <- 6  # fixed inches for bars, NES axis, and legend
  plot_width <- max(14, min(30, label_area_width + panel_area_width))

  ggsave(
    filename = file.path(outdir, paste0("GSEA_GOBP_Grouped_", label, "_FDRadj", filename_suffix, ".png")),
    plot     = p,
    width    = plot_width,
    height   = plot_height,
    dpi = 600
  )
}

plot_annotation_time_series_heatmap <- function(summary_list, timepoints, dir_filename, padj_cutoff = 0.05, font_size = 16, enrichment_type = "GSEA", max_annotations = 20, sig_cutoff = NULL, title = NULL) {
  # sig_cutoff controls significance dots on the heatmap; defaults to padj_cutoff when NULL
  sig_display_cutoff <- if (is.null(sig_cutoff)) padj_cutoff else sig_cutoff
  # Select classes that have p-adj<padj_cutoff at any of the specified timepoints
  # Collect all significant pathways across timepoints
  sig_pathways <- character(0)
  sig_pathway_records <- list()
  
  for (tp in timepoints) {
    if (!tp %in% names(summary_list) || is.null(summary_list[[tp]][[enrichment_type]])) {
      message(paste0("Warning: Missing ", enrichment_type, " data for timepoint: ", tp))
      next
    }
    
    enrichment_data <- summary_list[[tp]][[enrichment_type]]
    sig_in_tp_df <- if (is.null(padj_cutoff)) {
      enrichment_data %>% dplyr::select(pathway, padj)
    } else {
      enrichment_data %>% filter(padj < padj_cutoff) %>% dplyr::select(pathway, padj)
    }

    sig_in_tp <- sig_in_tp_df %>% pull(pathway)

    if (nrow(sig_in_tp_df) > 0) {
      sig_pathway_records[[length(sig_pathway_records) + 1]] <- sig_in_tp_df
    }
    
    sig_pathways <- union(sig_pathways, sig_in_tp)
  }
  
  if (length(sig_pathways) == 0) {
    message(paste0("No significant ", enrichment_type, " classes found across timepoints"))
    return(NULL)
  }

  # If too many significant pathways, keep those with the lowest min padj across selected timepoints
  if (length(sig_pathway_records) > 0 && length(sig_pathways) > max_annotations) {
    pathway_priority <- bind_rows(sig_pathway_records) %>%
      group_by(pathway) %>%
      summarise(min_padj = min(padj, na.rm = TRUE), .groups = "drop") %>%
      arrange(min_padj)

    sig_pathways <- pathway_priority %>%
      slice_head(n = max_annotations) %>%
      pull(pathway)

    message(
      paste0(
        "More than ", max_annotations, " significant ", enrichment_type,
        " classes found; plotting top ", max_annotations,
        " by minimum p-adj across selected timepoints"
      )
    )
  }
  
  # Create matrix of NES values
  nes_matrix <- matrix(NA, nrow = length(sig_pathways), ncol = length(timepoints),
                       dimnames = list(sig_pathways, timepoints))
  padj_matrix <- matrix(NA_real_, nrow = length(sig_pathways), ncol = length(timepoints),
                        dimnames = list(sig_pathways, timepoints))
  
  for (tp in timepoints) {
    if (!tp %in% names(summary_list) || is.null(summary_list[[tp]][[enrichment_type]])) next
    
    enrichment_data <- summary_list[[tp]][[enrichment_type]]
    for (pathway in sig_pathways) {
      pathway_data <- enrichment_data %>% filter(pathway == !!pathway)
      if (nrow(pathway_data) > 0) {
        nes_matrix[pathway, tp] <- pathway_data$NES[1]
        padj_matrix[pathway, tp] <- pathway_data$padj[1]
      }
    }
  }
  
  # Order pathways by hierarchical clustering across timepoints
  if (nrow(nes_matrix) > 1) {
    # Replace NAs with column means for distance calculation
    nes_for_clustering <- nes_matrix
    for (j in seq_len(ncol(nes_for_clustering))) {
      col_vals <- nes_for_clustering[, j]
      col_mean <- mean(col_vals, na.rm = TRUE)
      if (!is.finite(col_mean)) {
        col_mean <- 0
      }
      col_vals[is.na(col_vals)] <- col_mean
      nes_for_clustering[, j] <- col_vals
    }

    pathway_order <- hclust(dist(nes_for_clustering), method = "complete")$order
  } else {
    pathway_order <- 1
  }

  nes_matrix <- nes_matrix[pathway_order, , drop = FALSE]
  padj_matrix <- padj_matrix[pathway_order, , drop = FALSE]
  sig_pathways <- sig_pathways[pathway_order]
  
  # Clean up pathway names for display
  clean_pathway_names <- ifelse(nchar(sig_pathways) > 50, 
                                paste0(substr(sig_pathways, 1, 47), "..."), 
                                sig_pathways)
  clean_pathway_names <- make.unique(clean_pathway_names, sep = " ")
  rownames(nes_matrix) <- clean_pathway_names
  rownames(padj_matrix) <- clean_pathway_names
  
  # Convert to long format for ggplot
  heatmap_df <- as.data.frame(nes_matrix) %>%
    tibble::rownames_to_column("Pathway") %>%
    tidyr::pivot_longer(cols = -Pathway, names_to = "Timepoint", values_to = "NES") %>%
    mutate(
      Timepoint = factor(Timepoint, levels = timepoints),
      Pathway = factor(Pathway, levels = rev(clean_pathway_names))
    )

  padj_df <- as.data.frame(padj_matrix) %>%
    tibble::rownames_to_column("Pathway") %>%
    tidyr::pivot_longer(cols = -Pathway, names_to = "Timepoint", values_to = "padj") %>%
    mutate(
      Timepoint = factor(Timepoint, levels = timepoints),
      Pathway = factor(Pathway, levels = rev(clean_pathway_names))
    )

  heatmap_df <- heatmap_df %>%
    left_join(padj_df, by = c("Pathway", "Timepoint")) %>%
    mutate(Significant = !is.null(sig_display_cutoff) & !is.na(padj) & padj < sig_display_cutoff)
  
  # Create heatmap
  sig_dot_label <- paste0("p-adj < ", sig_display_cutoff)

  # Set color meaning based on GSEA (NES) vs ORA (direction * -log10(p-adj))
  if (grepl("GSEA", enrichment_type)) {
    color_meaning <- "NES"
  } else if (grepl("ORA", enrichment_type)) {
    color_meaning <- "Direction * -log10(p-adj)"
  } else {
    color_meaning <- "Statistic"
  }

  # Plot heatmap
  p <- ggplot(heatmap_df, aes(x = Timepoint, y = Pathway, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "#2166AC", # Dark blue
      mid = "white",
      high = "#B2182B", # Dark red
      midpoint = 0,
      na.value = "gray80",
      name = color_meaning
    ) +
    labs(
      title = title,
      x = "Timepoint",
      y = "Class"
    ) +
    guides(fill = guide_colorbar(order = 1)) +
    theme_minimal(base_size = font_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = font_size, color = "black"),
      axis.text.y = element_text(size = font_size, color = "black"),
      axis.title.x = element_text(size = font_size, color = "black"),
      axis.title.y = element_text(size = font_size, color = "black"),
      plot.title = element_text(hjust = 0.5, size = font_size, color = "black"),
      legend.title = element_text(size = font_size, color = "black"),
      legend.text = element_text(size = font_size, color = "black"),
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black"),
      legend.position = "right"
    )

  if (nrow(subset(heatmap_df, Significant)) > 0) {
    p <- p +
      geom_point(
        data = subset(heatmap_df, Significant),
        aes(x = Timepoint, y = Pathway, shape = sig_dot_label),
        inherit.aes = FALSE,
        size = 2,
        color = "black"
      ) +
      scale_shape_manual(
        values = setNames(16, sig_dot_label),
        name = "Significance"
      ) +
      guides(shape = guide_legend(order = 2, override.aes = list(size = 3, color = "black")))
  }
  
  # Dynamic height based on number of pathways
  plot_height <- max(8, min(20, length(sig_pathways) * 0.4))
  
  # Save plot
  ggsave(
    filename = dir_filename,
    plot = p,
    width = 10,
    height = plot_height,
    dpi = 600
  )
  
  message("Time series heatmap saved with ", length(sig_pathways),
          " significant pathways for ", enrichment_type,
          " with timepoints: ", paste(timepoints, collapse = ", "))
}

plot_protein_time_series_heatmap <- function(results_list,
                                      proteins,
                                      timepoints,
                                      outputdir,
                                      filename = "Protein_Time_Series_Heatmap.png", font_size=20, filename_suffix="", sig_cutoff = 0.05, metric_col = "beta", pval_col = NULL, sort_type = "hierarchical clustering") {
  # results_list: list of data frames per timepoint contrast
  # Use columns SYMBOL and beta (or value given by metric_col)
  
  # Check for timepoints as contrasts in results_list
  missing_timepoints <- setdiff(timepoints, names(results_list))
  if (length(missing_timepoints) > 0) {
    message("Warning: Missing timepoints in results_list: ", paste(missing_timepoints, collapse = ", "))
  }

  # Check for specified proteins in each results_list dataframe
  all_missing_proteins <- c()
  for (tp in timepoints) {
    if (!tp %in% names(results_list)) {
      message("Warning: Timepoint contrast '", tp, "' not found in results_list.")
      next
    }
    df <- results_list[[tp]]
    missing_proteins <- setdiff(proteins, df$SYMBOL)
    all_missing_proteins <- union(all_missing_proteins, missing_proteins)
    if (length(missing_proteins) > 0) {
      message("Warning: Missing proteins for timepoint '", tp, "': ", paste(missing_proteins, collapse = ", "))
    }
  }

  # Remove proteins missing in all timepoints
  proteins <- setdiff(proteins, all_missing_proteins)
  
  if (length(proteins) == 0) {
    message("Error: No valid proteins found across timepoints.")
    return(NULL)
  }
  
  get_sig_col <- function(df, preferred_col = NULL) {
    if (!is.null(preferred_col) && preferred_col %in% names(df)) return(preferred_col)
    candidates <- c("FDR", "padj", "p")
    found <- candidates[candidates %in% names(df)]
    if (length(found) == 0) return(NA_character_)
    found[1]
  }

  # Extract metric values and significance values for each protein at each timepoint
  metric_matrix <- matrix(NA, nrow = length(proteins), ncol = length(timepoints),
                        dimnames = list(proteins, timepoints))
  sig_matrix <- matrix(NA_real_, nrow = length(proteins), ncol = length(timepoints),
                       dimnames = list(proteins, timepoints))
  
  for (tp in timepoints) {
    if (!tp %in% names(results_list)) next
    
    df <- results_list[[tp]]
    sig_col <- get_sig_col(df, preferred_col = pval_col)
    if (is.na(sig_col)) {
      message("Warning: No significance column found for timepoint '", tp, "'. Expected one of: FDR, padj, p")
    }

    for (protein in proteins) {
      protein_data <- df %>% filter(SYMBOL == protein)
      if (nrow(protein_data) > 0) {
        metric_matrix[protein, tp] <- protein_data[[metric_col]][1]
        if (!is.na(sig_col)) {
          sig_matrix[protein, tp] <- as.numeric(protein_data[[sig_col]][1])
        }
      }
    }
  }
  
  # Remove proteins with all NA values
  valid_rows <- rowSums(!is.na(metric_matrix)) > 0
  if (sum(valid_rows) == 0) {
    message("Error: No proteins with valid ", metric_col, " values.")
    return(NULL)
  }
  metric_matrix <- metric_matrix[valid_rows, , drop = FALSE]
  sig_matrix <- sig_matrix[valid_rows, , drop = FALSE]

  if (sort_type == "hierarchical clustering") {
    # Order proteins by hierarchical clustering across selected timepoints
    if (nrow(metric_matrix) > 1) {
      metric_for_clustering <- metric_matrix
      for (j in seq_len(ncol(metric_for_clustering))) {
        col_vals <- metric_for_clustering[, j]
        col_mean <- mean(col_vals, na.rm = TRUE)
        if (!is.finite(col_mean)) {
          col_mean <- 0
        }
        col_vals[is.na(col_vals)] <- col_mean
        metric_for_clustering[, j] <- col_vals
      }

      protein_order <- hclust(dist(metric_for_clustering), method = "complete")$order
    } else {
      protein_order <- 1
    }
  # else if sort_type is a timepoint, order by lowest to highest NES for that timepoint
  } else if (sort_type %in% timepoints) {
    if (nrow(metric_matrix) > 1) {
      protein_order <- order(metric_matrix[, sort_type], decreasing = TRUE, na.last = TRUE)
    } else {
      protein_order <- 1
    }
  } else if (sort_type == "input") {
    # Preserve the order of the proteins argument as supplied
    protein_order <- seq_len(nrow(metric_matrix))
  } else {
    stop("Warning: Invalid sort_type '", sort_type, "' for time series heatmap.")
  }

  metric_matrix <- metric_matrix[protein_order, , drop = FALSE]
  sig_matrix <- sig_matrix[protein_order, , drop = FALSE]
  proteins_ordered <- rownames(metric_matrix)
  
  # Convert to long format for ggplot
  heatmap_df <- as.data.frame(metric_matrix) %>%
    tibble::rownames_to_column("Protein") %>%
    tidyr::pivot_longer(cols = -Protein, names_to = "Timepoint", values_to = metric_col) %>%
    mutate(Timepoint = factor(Timepoint, levels = timepoints),
           Protein = factor(Protein, levels = rev(proteins_ordered)))

  sig_df <- as.data.frame(sig_matrix) %>%
    tibble::rownames_to_column("Protein") %>%
    tidyr::pivot_longer(cols = -Protein, names_to = "Timepoint", values_to = "SigVal") %>%
    mutate(
      Timepoint = factor(Timepoint, levels = timepoints),
      Protein = factor(Protein, levels = rev(proteins_ordered))
    )

  heatmap_df <- heatmap_df %>%
    left_join(sig_df, by = c("Protein", "Timepoint")) %>%
    mutate(Significant = if (is.null(sig_cutoff)) FALSE else (!is.na(SigVal) & SigVal < sig_cutoff))
  
  # Create heatmap
  sig_col_display <- if (is.null(pval_col) || pval_col == "FDR") "p-adj" else pval_col
  sig_dot_label <- paste0(sig_col_display, " < ", sig_cutoff)

  p <- ggplot(heatmap_df, aes(x = Timepoint, y = Protein, fill = !!sym(metric_col))) +
    geom_tile(color = "white", linewidth = 0.5) +
     scale_fill_gradient2(
      low = "#2166AC", # Dark blue
      mid = "white",
      high = "#B2182B", # Dark red
      midpoint = 0,
      na.value = "gray80",
      name = if (metric_col == "beta") "\u03b2" else metric_col
    ) +
    theme_minimal(base_size = font_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = font_size, color = "black"),
      axis.text.y = element_text(size = font_size, color = "black"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  if (nrow(subset(heatmap_df, Significant)) > 0) {
    p <- p +
      geom_point(
        data = subset(heatmap_df, Significant),
        aes(x = Timepoint, y = Protein, shape = sig_dot_label),
        inherit.aes = FALSE,
        size = 2,
        color = "black"
      ) +
      scale_shape_manual(
        values = setNames(16, sig_dot_label),
        name = "Significance"
      )
  }
  
  # Dynamic height based on number of proteins
  plot_height <- max(6, min(20, length(proteins_ordered) * 0.4))
  
  ggsave(
    filename = file.path(outputdir, gsub(".png$", paste0(filename_suffix, ".png"), filename)),
    plot = p,
    width = 8,
    height = plot_height,
    dpi = 600
  )
  
  message("Time series heatmap saved to: ", file.path(outputdir, gsub(".png$", paste0(filename_suffix, ".png"), filename)))
}

pvalue_histogram <- function(pvalues, # vector of p-values
                             outputdir,
                             label,
                             filename_suffix="",
                             b=0.05, # bin width
                             alpha = 0.05, # significance threshold
                             font_size = 16
                             ) {
  # Filter out NA values and invalid p-values
  pvalues <- pvalues[!is.na(pvalues) & is.finite(pvalues)]
  pvalues <- pvalues[pvalues > 0 & pvalues <= 1]
  
  if (length(pvalues) == 0) {
    message("No valid p-values for histogram: ", label)
    return(NULL)
  }
  
  stopifnot(all(is.numeric(pvalues)) & all(pvalues > 0) & all(pvalues <= 1));
  stopifnot(length(b) == 1 & is.numeric(b) & b >= 0 & b <= 0.2);

  p.df <- data.frame(p = pvalues);
  m <- sum(!is.na(pvalues));
  signal.cutoff <- qbinom(
      p = 1 - alpha,
      size = m,
      prob = b
      );
  qc.cutoff <- qbinom(
      p = 1 - alpha / (1 / b),
      size = m,
      prob = b
      );
  
  p <- ggplot(p.df, aes(x = p)) +
    geom_histogram(binwidth = b, boundary = 0, fill = "skyblue", color = "black") +
    # Option 1:
    # geom_hline(yintercept = c(signal.cutoff, qc.cutoff), 
    #            color = c("#5a6bca", "#1d2f96"), linewidth = 1) +
    # annotate("text", x = 0.98, y = signal.cutoff, label = "Signal cutoff", 
    #          vjust = -0.5, hjust = 1, color = "#5a6bca", size = font_size * 0.25) +
    # annotate("text", x = 0.98, y = qc.cutoff, label = "QC cutoff (Bonferroni)", 
    #          vjust = -0.5, hjust = 1, color = "#1d2f96", size = font_size * 0.25) +
    # Option 2:
    geom_hline(yintercept = qc.cutoff, color = "#1d2f96", linewidth = 1) +
    annotate("text", x = 0.98, y = qc.cutoff, label = "QC cutoff (Bonferroni)", 
             vjust = -0.5, hjust = 1, color = "#1d2f96", size = font_size * 0.25) +
    labs(
      x = "P-values",
      y = "Frequency"
    ) +
    theme_minimal(base_size = font_size)
  ggsave(
    filename = file.path(outputdir, paste0("Pvalue_Histogram_", label, ".png")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 600
  )
}

run_pathway_enrichment <- function(res_df, label, mapping,
                                   base_dir, filename_suffix="", regr_type="Linear", ora_p_col = "p", ora_p_cutoff = 0.05,
                                   use_t_statistic = FALSE) {
  
  if (is.null(res_df)) return(NULL)
  outdir <- file.path(base_dir, label) # create subfolder per contrast
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  message("\n=== PATHWAY ANALYSIS: ", label, " ===")
  
  # Use OlinkID if present; fallback to Protein for backward compatibility
  id_col <- if ("OlinkID" %in% names(res_df)) "OlinkID" else "Protein"
  
  reg <- res_df %>%
    left_join(mapping, by = setNames("OlinkID", id_col)) %>%
    filter(!is.na(ENTREZID)) %>%
    mutate(
      ENTREZID = as.character(ENTREZID),
      beta     = as.numeric(beta),
      p        = as.numeric(p),
      FDR      = as.numeric(FDR)
    ) %>%
    group_by(ENTREZID) %>%
    summarise(
      beta    = mean(beta, na.rm = TRUE),
      p       = suppressWarnings(min(p, na.rm = TRUE)),
      FDR     = suppressWarnings(min(FDR, na.rm = TRUE)),
      SYMBOL  = dplyr::first(SYMBOL),
      UniProt = dplyr::first(UniProt),
      OlinkID = dplyr::first(OlinkID),
      .groups = "drop"
    ) %>%
    filter(is.finite(beta), is.finite(p))
  
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
  
  reg_list[[label]] <<- reg
  
  # Sort reg by ascending p-value and reorder columns for better readability
  reg <- reg %>% 
    arrange(p) %>%
    dplyr::select(SYMBOL, Description, ENTREZID, UniProt, OlinkID, beta, p, FDR)
  
  # Write Excel file with two sheets: all results and significant results (FDR < 0.1)
  reg_sig <- reg %>% 
    filter(FDR < 0.1) %>%
    mutate(
      beta = signif(beta, 2),
      p = signif(p, 2),
      FDR = signif(FDR, 2)
    ) %>%
    dplyr::select(SYMBOL, Description, beta, p, FDR)
  
  write.xlsx(
    list(
      all = reg,
      sig = reg_sig
    ),
    file = file.path(outdir, paste0(regr_type, "_Protein_Results_Annotated_", label, filename_suffix, ".xlsx"))
  )
  
  # Create ranked gene list for GSEA: sign(beta) * -log10(p)
  # Positive ranks = upregulated genes with low p-values
  # Negative ranks = downregulated genes with low p-values
  ranks <- if (use_t_statistic) {
    setNames(reg$statistic, reg$ENTREZID)
  } else {
    setNames(sign(reg$beta) * -log10(pmax(reg$p, .Machine$double.xmin)), reg$ENTREZID)
  }
  
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
  
  genes_up   <- reg %>% filter(!!sym(ora_p_col) < ora_p_cutoff, beta > 0) %>% pull(ENTREZID)
  genes_down <- reg %>% filter(!!sym(ora_p_col) < ora_p_cutoff, beta < 0) %>% pull(ENTREZID)
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
            eps         = 1e-10,   # Set minimum p-value threshold to prevent underestimation
            nproc       = 1 # control random state
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
  
  if (length(all_summaries) > 0) {
    final <- bind_rows(all_summaries)
    summary_list[[label]] <<- final
    
    # Summary_Pathways_<label>.xlsx saves all significant pathways (GSEA and ORA)
    # Sort by ascending PValue
    final <- final %>% arrange(PValue)
    
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
    
    gsea_gobp_sig <- final %>%
      filter(Source == "GSEA", DB == "GOBP", PValue < 0.05, !is.na(GOBP_Group))
  }
  
  message("✓ Completed contrast: ", label)
}

safe_msig <- function(...) {
  out <- try(msigdbr(...), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out
}

extract_gsea_enrich <- function(enrich_all, db) {
  # NES is stored in the Stat column; pval from PValue
  enrich_all %>%
    filter(Source == "GSEA", DB == db) %>%
    transmute(pathway = Pathway, pval = PValue, padj = FDR, NES = Stat, Genes = Genes)
}

extract_ora_enrich <- function(enrich_all, db) {
  # Use signed -log10(FDR) as pseudo-NES (positive = UP, negative = DOWN)
  enrich_all %>%
    filter(Source == "ORA", DB == db) %>%
    mutate(NES = ifelse(Direction == "UP", 1, -1) * pmin(-log10(pmax(FDR, 1e-300)), 10)) %>%
    group_by(Pathway) %>%
    slice_min(order_by = FDR, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(pathway = Pathway, padj = FDR, NES)
}


run_protein_time_series_heatmap_generation <- function(results_list, tps_for_focus,
                                                       top_n_list, outputdir, filename_suffix = "",
                                                       timepoints_label = NULL, metric_col = "beta", pval_col = "p", pval_cutoff = 0.05) {
  tp_tag <- if (!is.null(timepoints_label)) paste0("_", timepoints_label) else ""
  for (tp_focus in tps_for_focus) {
    results_tp_focus <- results_list[[tp_focus]]
    for (top_n in top_n_list) {
      heatmap_proteins_topn <- results_tp_focus %>%
        arrange(.data[[pval_col]]) %>%
        slice_head(n = top_n) %>%
        pull(SYMBOL)
      heatmap_proteins_sorted <- results_tp_focus %>%
        filter(SYMBOL %in% heatmap_proteins_topn) %>%
        arrange(desc(.data[[metric_col]])) %>%
        pull(SYMBOL)
      plot_protein_time_series_heatmap(
        results_list = results_list,
        proteins = heatmap_proteins_sorted,
        timepoints = tps_for_focus,
        outputdir = outputdir,
        filename = paste0("Protein_Time_Series_Heatmap_", tp_focus, "_", top_n, tp_tag, ".png"),
        filename_suffix = filename_suffix,
        metric_col = metric_col,
        sig_cutoff = pval_cutoff,
        pval_col = pval_col
      )
    }
  }
}


run_tp_focus_annotation_time_series_heatmap_generation <- function(annotation_sets, annotation_subdirs,
                                                                    tps_for_focus, timepoints,
                                                                    n_top_annotations, output_dir,
                                                                    filename_suffix = "") {
  for (enrich_name in names(annotation_sets)) {
    enrich_list <- annotation_sets[[enrich_name]]
    subdir <- annotation_subdirs[[enrich_name]]

    for (tp_focus in tps_for_focus) {
      focus_data <- enrich_list[[tp_focus]][[enrich_name]]
      if (is.null(focus_data) || nrow(focus_data) == 0) next

      top_pathways <- focus_data %>%
        arrange(padj) %>%
        slice_head(n = n_top_annotations) %>%
        pull(pathway)

      if (length(top_pathways) == 0) next

      filtered_list <- enrich_list
      for (tp in timepoints) {
        if (!is.null(filtered_list[[tp]][[enrich_name]])) {
          filtered_list[[tp]][[enrich_name]] <- filtered_list[[tp]][[enrich_name]] %>%
            filter(pathway %in% top_pathways)
        }
      }

      plot_annotation_time_series_heatmap(
        filtered_list,
        timepoints,
        file.path(output_dir, subdir,
                  paste0("Time_Series_Heatmap_", enrich_name, "_Top", n_top_annotations, "_", tp_focus, filename_suffix, ".png")),
        enrichment_type = enrich_name,
        padj_cutoff = 1.0,
        sig_cutoff = 0.05,
        max_annotations = n_top_annotations,
        font_size = 12
      )
    }
  }
}

fit_one_feature <- function(df, outcome_var, covars,
                            measurement_col    = "feature_value",
                            impute_measurement = FALSE,
                            return_statistic   = FALSE) {
  na_ret <- if (return_statistic) tibble(beta = NA_real_, p = NA_real_, statistic = NA_real_)
             else                 tibble(beta = NA_real_, p = NA_real_)

  # Check minimum sample size
  if (nrow(df) < 3) return(na_ret)

  # Create outcome and feature value variables
  df <- df %>%
    mutate(
      OUTCOME     = .data[[outcome_var]],
      feature_use = .data[[measurement_col]]
    )

  if (impute_measurement) {
    # Impute missing measurement values with per-feature min (proteomics behaviour)
    if (any(is.na(df$feature_use))) {
      df$feature_use[is.na(df$feature_use)] <- min(df$feature_use, na.rm = TRUE)
    }
  } else {
    # Measurement must be pre-imputed upstream (metabolomics behaviour)
    if (any(is.na(df$feature_use))) {
      stop("Missing measurement values in fit_one_feature(). Check data processing.")
    }
  }

  # Check for zero variance
  if (stats::sd(df$feature_use) == 0) return(na_ret)

  # Define covariates to keep - remove zero-variance ones for this feature's subset
  keep <- covars[sapply(df[covars], function(x) {
    if (is.factor(x) || is.character(x)) return(length(unique(na.omit(x))) > 1)
    stats::sd(x, na.rm = TRUE) > 0
  })]

  # Create formula with kept covariates
  form <- as.formula(
    if (length(keep) == 0) "OUTCOME ~ feature_use"
    else paste("OUTCOME ~ feature_use +", paste(keep, collapse = " + "))
  )

  # Fit linear regression model
  fit <- try(stats::lm(form, data = df), silent = TRUE)
  if (inherits(fit, "try-error")) return(na_ret)

  # Extract coefficient, p-value and optionally t-statistic for feature_use
  cf <- summary(fit)$coefficients
  if (!"feature_use" %in% rownames(cf)) return(na_ret)

  if (return_statistic) {
    tibble(
      beta      = cf["feature_use", "Estimate"],
      p         = cf["feature_use", "Pr(>|t|)"],
      statistic = cf["feature_use", "t value"]
    )
  } else {
    tibble(
      beta = cf["feature_use", "Estimate"],
      p    = cf["feature_use", "Pr(>|t|)"]
    )
  }
}

fit_one_feature_delta <- function(df, outcome_var, covars, delta_col, return_statistic = FALSE) {
  na_ret <- if (return_statistic) tibble(beta = NA_real_, p = NA_real_, statistic = NA_real_)
             else                 tibble(beta = NA_real_, p = NA_real_)

  # Check minimum sample size
  if (nrow(df) < 3) return(na_ret)

  # Create outcome and delta variables
  df <- df %>%
    mutate(
      OUTCOME   = .data[[outcome_var]],
      delta_use = .data[[delta_col]]
    )

  # Remove rows where delta is NA (no Baseline or follow-up measurement)
  df <- df %>% filter(!is.na(delta_use))
  if (nrow(df) < 3) return(na_ret)

  # Check for zero variance in delta
  if (stats::sd(df$delta_use) == 0) return(na_ret)

  # Define covariates to keep - remove zero-variance ones
  keep <- covars[sapply(df[covars], function(x) {
    if (is.factor(x) || is.character(x)) return(length(unique(na.omit(x))) > 1)
    stats::sd(x, na.rm = TRUE) > 0
  })]

  # Create formula with kept covariates
  form <- as.formula(
    if (length(keep) == 0) "OUTCOME ~ delta_use"
    else paste("OUTCOME ~ delta_use +", paste(keep, collapse = " + "))
  )

  # Fit linear regression model
  fit <- try(stats::lm(form, data = df), silent = TRUE)
  if (inherits(fit, "try-error")) return(na_ret)

  # Extract coefficient, p-value and optionally t-statistic for delta_use
  cf <- summary(fit)$coefficients
  if (!"delta_use" %in% rownames(cf)) return(na_ret)

  if (return_statistic) {
    tibble(
      beta      = cf["delta_use", "Estimate"],
      p         = cf["delta_use", "Pr(>|t|)"],
      statistic = cf["delta_use", "t value"]
    )
  } else {
    tibble(
      beta = cf["delta_use", "Estimate"],
      p    = cf["delta_use", "Pr(>|t|)"]
    )
  }
}

fit_one_protein <- function(df, outcome_var, covars, return_statistic = TRUE) {
  # Thin wrapper around fit_one_feature() for proteomics NPX data.
  # Imputes missing NPX with per-protein min.
  fit_one_feature(df, outcome_var, covars,
                  measurement_col    = "NPX_mean",
                  impute_measurement = TRUE,
                  return_statistic   = return_statistic)
}

run_regression_contrast <- function(df, patient_meta, label, outcome_var, covars,
                                    id_col, measurement_col,
                                    visit_filter     = NULL,
                                    impute_covars    = FALSE,
                                    return_statistic = FALSE) {
  # Generic per-feature linear regression contrast.
  # id_col          : grouping identifier column (ie: "OlinkID" or "feature_label")
  # measurement_col : predictor column in df (ie: "NPX_mean" or "normalized_abundance")
  # impute_covars   : TRUE → impute missing covariates (proteomics);
  #                   FALSE → stop if missing covariates (metabolomics)
  # return_statistic: TRUE → also return t-statistic column (metabolomics)

  message("Running regression for: ", label)

  # Filter by visit if specified
  if (!is.null(visit_filter)) {
    d <- df %>% filter(visit == visit_filter)
  } else {
    d <- df
  }

  if (nrow(d) == 0) return(NULL)

  # Join patient metadata and prepare data
  d <- d %>%
    left_join(patient_meta, by = "PTID") %>%
    mutate(across(any_of(c("gender", "Diabetes", "Hypertension")), factor)) %>%
    filter(!is.na(.data[[outcome_var]]))

  # Convert any remaining character covariates to factors
  for (cv in covars) {
    if (cv %in% colnames(d) && is.character(d[[cv]])) d[[cv]] <- factor(d[[cv]])
  }

  if (nrow(d) == 0) return(NULL)

  # Pre-filter covariates: remove those with >50% missing
  valid_covars <- covars[sapply(d[covars], function(x) mean(is.na(x)) <= 0.5)]

  if (impute_covars) {
    for (cv in valid_covars) {
      if (is.factor(d[[cv]])) {
        d[[cv]] <- impute_factor_mode(d[[cv]])
      } else {
        d[[cv]] <- impute_numeric_mean(d[[cv]])
      }
    }
  } else {
    if (length(valid_covars) > 0 && any(sapply(d[valid_covars], function(x) any(is.na(x))))) {
      stop("Missing covariate values in run_regression_contrast(). Check data processing.")
    }
  }

  # Standardise measurement column name for fit_one_feature()
  names(d)[names(d) == measurement_col] <- "feature_value"

  # Group by id_col and fit models
  result <- d %>%
    group_by(!!sym(id_col)) %>%
    group_modify(~fit_one_feature(.x, outcome_var, valid_covars,
                                  measurement_col    = "feature_value",
                                  impute_measurement = impute_covars,
                                  return_statistic   = return_statistic)) %>%
    ungroup()

  # Backward-compatible Protein alias for proteomics downstream scripts
  if (id_col == "OlinkID") result <- result %>% mutate(Protein = OlinkID)

  result <- result %>% mutate(contrast = label)

  # Add FDR correction (Benjamini-Hochberg)
  if (nrow(result) > 0) {
    result <- result %>% mutate(FDR = p.adjust(p, method = "BH"))
  }

  result %>% arrange(p)
}

run_regression_contrast_delta <- function(deltas_df, patient_meta, label, outcome_var, covars,
                                          id_col, delta_col,
                                          return_statistic = FALSE) {
  # Generic delta (double-delta) linear regression contrast.
  # id_col   : grouping identifier column ("OlinkID" or "feature_label")
  # delta_col: column in deltas_df holding the pre-computed delta values
  # Covariates are always imputed (both pipelines impute in delta models).

  message("Running delta regression for: ", label)

  d <- deltas_df %>%
    left_join(patient_meta, by = "PTID") %>%
    mutate(across(any_of(c("gender", "Diabetes", "Hypertension")), factor)) %>%
    filter(!is.na(.data[[outcome_var]]))

  for (cv in covars) {
    if (cv %in% colnames(d) && is.character(d[[cv]])) d[[cv]] <- factor(d[[cv]])
  }

  if (nrow(d) == 0) return(NULL)

  valid_covars <- covars[sapply(d[covars], function(x) mean(is.na(x)) <= 0.5)]

  for (cv in valid_covars) {
    if (is.factor(d[[cv]])) {
      d[[cv]] <- impute_factor_mode(d[[cv]])
    } else {
      d[[cv]] <- impute_numeric_mean(d[[cv]])
    }
  }

  result <- d %>%
    group_by(!!sym(id_col)) %>%
    group_modify(~fit_one_feature_delta(.x, outcome_var, valid_covars,
                                        delta_col        = delta_col,
                                        return_statistic = return_statistic)) %>%
    ungroup()

  # Backward-compatible Protein alias for proteomics downstream scripts
  if (id_col == "OlinkID") result <- result %>% mutate(Protein = OlinkID)

  result <- result %>% mutate(contrast = label)

  if (nrow(result) > 0) {
    result <- result %>% mutate(FDR = p.adjust(p, method = "BH"))
  }

  result %>% arrange(p)
}

remove_ptids_without_transplant_type <- function(patient_metadata_df, data_df, ptid_visit_colnames, report_file, id_colname="OlinkID") {
  patients_before <- nrow(patient_metadata_df)
  patient_metadata_df <- patient_metadata_df %>%
    dplyr::filter(!is.na(TRANSPLANT_TYPE) & TRANSPLANT_TYPE != "" & trimws(TRANSPLANT_TYPE) != "")
  patients_after <- nrow(patient_metadata_df)
  if (patients_before != patients_after) {
    message(paste0("Removed ", patients_before - patients_after, " patients without TRANSPLANT_TYPE (",
                  round((patients_before - patients_after) / patients_before * 100, 2), "% of patients)."))
  }
  append_report_line(report_file, "")
  append_report_line(report_file, "--- Patient Metadata Filtering (missing TRANSPLANT_TYPE) ---")
  append_report_line(report_file, paste0("Patients before: ", patients_before, ", after: ", patients_after, ", removed: ", patients_before - patients_after))
  # Check only ALLO and AUTO remain
  remaining_HCT_types <- unique(patient_metadata_df$TRANSPLANT_TYPE)
  if (!all(remaining_HCT_types %in% c("ALLO", "AUTO"))) {
    stop(paste0("Found unexpected transplant types in patient metadata: ", paste(remaining_HCT_types, collapse = ", ")))
  }

  # Filter data_df sample columns to only include patients present in metadata
  valid_ptids <- unique(patient_metadata_df$PTID)
  keep_cols <- c(id_colname, ptid_visit_colnames[
    sub("^X", "", sub("_.*$", "", ptid_visit_colnames)) %in% as.character(valid_ptids)
  ])
  samples_removed <- length(ptid_visit_colnames) - (length(keep_cols) - 1)
  if (samples_removed > 0) {
    message(paste0("Removed ", samples_removed, " samples from patients without TRANSPLANT_TYPE."))
    data_df <- data_df[, keep_cols, drop = FALSE]
    ptid_visit_colnames <- setdiff(colnames(data_df), id_colname)
  }
  append_report_line(report_file, "For metadata and data:")
  append_report_line(report_file, paste0("Samples removed (no TRANSPLANT_TYPE): ", samples_removed))
  append_report_line(report_file, paste0("Samples remaining after patient filtering: ", length(ptid_visit_colnames)))
  return(list(patient_metadata_df = patient_metadata_df, data_df = data_df, ptid_visit_colnames = ptid_visit_colnames))
}

generate_sample_wise_feature_distribution_plots <- function(data_df, output_dir, filename_stem, filename_suffix = "", id_colname = "OlinkID", sample_colnames = NULL, timepoint_labels = c("Baseline", "6m", "12m"), font_size = 10, bin_size = 50, y_axis_label = "Normalized Abundance", single_plot = TRUE) {
  ###
  # data_df:         Data frame with rows as features and columns as samples (PTID_visit), plus an ID column.
  # output_dir:      Directory to save plots.
  # filename_stem:   Stem for plot filenames (e.g. "protein_distribution_samples_all_timepoints").
  # filename_suffix: Suffix appended to output plot filenames (e.g. "_ALLO", "_prefiltering", or "").
  # id_colname:      Name of the column in data_df that contains feature IDs (default "OlinkID").
  # sample_colnames: Optional vector of column names that correspond to sample measurements.
  #                  If NULL (default), derived as setdiff(colnames(data_df), id_colname).
  # timepoint_labels: Vector of timepoint labels to include in plots (default c("Baseline", "6m", "12m")).
  # font_size:       Base font size for plots.
  # bin_size:        Number of samples per plot when single_plot = FALSE (default 50).
  # y_axis_label:    Label for the y-axis.
  # single_plot:     If TRUE (default), create one combined faceted plot across all timepoints
  #                  instead of separate binned per-timepoint plots.
  ###
  if (is.null(sample_colnames)) {
    sample_colnames <- setdiff(colnames(data_df), id_colname)
  }

  if (length(sample_colnames) == 0) {
    message("Warning: No sample columns available for sample-wise feature distribution plots.")
    return(invisible(NULL))
  }

  sample_distribution_df <- data_df %>%
    dplyr::select(all_of(c(id_colname, sample_colnames))) %>%
    tidyr::pivot_longer(
      cols = all_of(sample_colnames),
      names_to = "sample_id",
      values_to = "feature_abundance"
    ) %>%
    tidyr::separate("sample_id", into = c("PTID", "visit"), sep = "_", remove = FALSE, extra = "merge", fill = "right")

  if (single_plot) {
    plot_df <- sample_distribution_df %>%
      dplyr::filter(.data$visit %in% timepoint_labels) %>%
      dplyr::mutate(
        visit = factor(.data$visit, levels = timepoint_labels),
        PTID_num = as.integer(sub("^X", "", .data$PTID))
      )

    if (nrow(plot_df) == 0) {
      message("Warning: No rows available to plot sample-wise feature distributions.")
      return(invisible(NULL))
    }

    ordered_samples <- plot_df %>%
      dplyr::distinct(.data$visit, .data$sample_id, .data$PTID_num) %>%
      dplyr::arrange(.data$visit, .data$PTID_num, .data$sample_id) %>%
      dplyr::pull(.data$sample_id)

    plot_df <- plot_df %>%
      dplyr::mutate(sample_id = factor(.data$sample_id, levels = ordered_samples))

    n_timepoints <- length(unique(plot_df$visit))
    plot_height <- max(7, 3.5 * n_timepoints)

    p <- ggplot(plot_df, aes(x = .data$sample_id, y = .data$feature_abundance)) +
      geom_boxplot(outlier.size = 0.2, linewidth = 0.3, na.rm = TRUE) +
      facet_wrap(~visit, scales = "free_x", ncol = 1) +
      labs(
        x = "Sample",
        y = y_axis_label
      ) +
      theme_bw(base_size = font_size) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = font_size * 0.5),
        plot.title = element_blank(),
        panel.grid = element_blank()
      )

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(
      filename = file.path(output_dir, paste0(filename_stem, filename_suffix, ".png")),
      plot = p,
      width = 18,
      height = plot_height,
      dpi = 600
    )

    return(invisible(NULL))
  }

  for (current_timepoint in timepoint_labels) {
    timepoint_df <- sample_distribution_df %>%
      dplyr::filter(.data$visit == current_timepoint)

    if (nrow(timepoint_df) == 0) {
      next
    }

    timepoint_folder <- file.path(output_dir, current_timepoint)
    dir.create(timepoint_folder, recursive = TRUE, showWarnings = FALSE)

    ordered_samples <- timepoint_df %>%
      dplyr::distinct(.data$sample_id, .data$PTID) %>%
      dplyr::mutate(PTID_num = as.integer(sub("^X", "", .data$PTID))) %>%
      dplyr::arrange(.data$PTID_num, .data$sample_id) %>%
      dplyr::pull(.data$sample_id)

    sample_bins <- split(ordered_samples, ceiling(seq_along(ordered_samples) / bin_size))

    for (bin_idx in seq_along(sample_bins)) {
      bin_samples <- sample_bins[[bin_idx]]

      if (length(bin_samples) == 0) {
        next
      }

      bin_df <- timepoint_df %>%
        dplyr::filter(.data$sample_id %in% bin_samples) %>%
        dplyr::mutate(sample_id = factor(.data$sample_id, levels = bin_samples))

      sample_start_idx <- (bin_idx - 1) * bin_size + 1
      sample_end_idx <- sample_start_idx + length(bin_samples) - 1

      p <- ggplot(bin_df, aes(x = .data$sample_id, y = .data$feature_abundance)) +
        geom_boxplot(outlier.size = 0.2, linewidth = 0.3, na.rm = TRUE) +
        labs(
          x = "Sample",
          y = y_axis_label
        ) +
        theme_bw(base_size = font_size) +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = font_size * 0.5),
          plot.title = element_blank(),
          panel.grid = element_blank()
        )

      ggsave(
        filename = file.path(timepoint_folder, paste0(filename_stem, "_", sprintf("%03d", sample_start_idx), "_", sprintf("%03d", sample_end_idx), filename_suffix, ".png")),
        plot = p,
        width = 16,
        height = 7,
        dpi = 600
      )
    }
  }
}

generate_sample_median_variability_plots <- function(data_df, output_dir, filename_stem, filename_suffix = "", id_colname = "OlinkID", timepoint_labels = c("Baseline", "6m", "12m"), font_size = 12, y_axis_label = "Sample Median NPX", add_stats = FALSE) {
  # data_df: rows = proteins, cols = samples (PTID_visit)
  # add_stats: if TRUE, run Friedman test (overall) + pairwise paired Wilcoxon (BH-corrected)
  #            and annotate the distribution plot with significance brackets
  sample_data_cols <- setdiff(colnames(data_df), id_colname)
  if (length(sample_data_cols) == 0) {
    message("Warning: No sample columns available for sample-median variability plots.")
    return(invisible(NULL))
  }
  data_matrix <- as.matrix(data_df[, sample_data_cols])
  sample_medians <- matrixStats::colMedians(data_matrix, na.rm = TRUE)
  visit_labels <- sub("^X?[^_]+_", "", sample_data_cols)
  sample_median_df <- tibble::tibble(
    sample_id = sample_data_cols,
    visit = visit_labels,
    sample_median = sample_medians
  ) %>%
    dplyr::filter(.data$visit %in% timepoint_labels, is.finite(.data$sample_median)) %>%
    dplyr::mutate(
      visit = factor(.data$visit, levels = timepoint_labels),
      PTID  = sub("^X?([^_]+)_.*", "\\1", .data$sample_id)
    )
  if (nrow(sample_median_df) == 0) {
    message("Warning: No finite sample medians available for timepoint variability plots.")
    return(invisible(NULL))
  }
  summary_df <- sample_median_df %>%
    dplyr::group_by(.data$visit) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median = stats::median(.data$sample_median, na.rm = TRUE),
      q1 = stats::quantile(.data$sample_median, probs = 0.25, na.rm = TRUE),
      q3 = stats::quantile(.data$sample_median, probs = 0.75, na.rm = TRUE),
      sd = stats::sd(.data$sample_median, na.rm = TRUE),
      .groups = "drop"
    )
  p_distribution <- ggplot(sample_median_df, aes(x = .data$visit, y = .data$sample_median)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.4, fill = "#CFE3F2") +
    geom_jitter(width = 0.15, height = 0, alpha = 0.8, size = 1.6, color = "#2E5C8A") +
    labs(x = "Timepoint", y = y_axis_label) +
    theme_bw(base_size = font_size) +
    theme(plot.title = element_blank(), panel.grid = element_blank())
  p_summary <- ggplot(summary_df, aes(x = .data$visit, y = .data$median)) +
    geom_pointrange(aes(ymin = .data$q1, ymax = .data$q3), color = "#7A1F59", linewidth = 0.6) +
    geom_text(aes(label = paste0("n=", .data$n, ", SD=", round(.data$sd, 3))), vjust = -0.8, size = (font_size * 0.8) / ggplot2::.pt) +
    labs(x = "Timepoint", y = y_axis_label) +
    theme_bw(base_size = font_size) +
    theme(plot.title = element_blank(), panel.grid = element_blank())

  if (add_stats && dplyr::n_distinct(sample_median_df$visit) >= 2) {
    visit_levels <- levels(sample_median_df$visit)

    # --- Friedman test (paired, repeated measures across timepoints) ---
    # Requires patients with measurements at all timepoints
    friedman_p <- NA_real_
    complete_ptids <- sample_median_df %>%
      dplyr::group_by(.data$PTID) %>%
      dplyr::filter(dplyr::n_distinct(.data$visit) == length(visit_levels)) %>%
      dplyr::pull(.data$PTID) %>%
      unique()
    wide_df <- sample_median_df %>%
      dplyr::filter(.data$PTID %in% complete_ptids) %>%
      dplyr::select("PTID", "visit", "sample_median") %>%
      tidyr::pivot_wider(names_from = "visit", values_from = "sample_median") %>%
      tidyr::drop_na()
    if (nrow(wide_df) >= 3 && length(visit_levels) >= 3) {
      friedman_p <- friedman.test(as.matrix(wide_df[, visit_levels]))$p.value
    }

    # --- Pairwise paired Wilcoxon signed-rank tests (BH-corrected) ---
    visit_pairs <- utils::combn(visit_levels, 2, simplify = FALSE)
    pairwise_results <- dplyr::bind_rows(lapply(visit_pairs, function(pair) {
      common_ptids <- intersect(
        sample_median_df %>% dplyr::filter(.data$visit == pair[1]) %>% dplyr::pull(.data$PTID),
        sample_median_df %>% dplyr::filter(.data$visit == pair[2]) %>% dplyr::pull(.data$PTID)
      )
      v1 <- sample_median_df %>% dplyr::filter(.data$visit == pair[1], .data$PTID %in% common_ptids) %>% dplyr::arrange(.data$PTID) %>% dplyr::pull(.data$sample_median)
      v2 <- sample_median_df %>% dplyr::filter(.data$visit == pair[2], .data$PTID %in% common_ptids) %>% dplyr::arrange(.data$PTID) %>% dplyr::pull(.data$sample_median)
      p_val <- if (length(v1) >= 3) wilcox.test(v1, v2, paired = TRUE, exact = FALSE)$p.value else NA_real_
      tibble::tibble(group1 = pair[1], group2 = pair[2], p_raw = p_val, n_pairs = length(v1))
    }))
    pairwise_results <- pairwise_results %>%
      dplyr::mutate(
        p_adj = stats::p.adjust(.data$p_raw, method = "BH"),
        sig_label = dplyr::case_when(
          is.na(.data$p_adj)    ~ "NA",
          .data$p_adj < 0.001   ~ "***",
          .data$p_adj < 0.01    ~ "**",
          .data$p_adj < 0.05    ~ "*",
          TRUE                  ~ "ns"
        )
      )

    # --- Add Friedman p as caption ---
    friedman_caption <- if (!is.na(friedman_p)) {
      paste0("Friedman test: p = ", signif(friedman_p, 3),
             " (n = ", nrow(wide_df), " complete-case patients)")
    } else {
      "Friedman test: insufficient data (needs \u22653 complete-case patients across all timepoints)"
    }
    p_distribution <- p_distribution +
      labs(caption = friedman_caption) +
      theme(plot.caption = element_text(size = font_size - 5, hjust = 0.5, color = "gray40"))

    # --- Add significance brackets to distribution plot ---
    y_vals  <- sample_median_df$sample_median
    y_max   <- max(y_vals, na.rm = TRUE)
    y_range <- diff(range(y_vals, na.rm = TRUE))
    y_step  <- y_range * 0.12
    bracket_height <- y_step * 0.25

    for (i in seq_len(nrow(pairwise_results))) {
      x1 <- match(pairwise_results$group1[i], visit_levels)
      x2 <- match(pairwise_results$group2[i], visit_levels)
      y_bracket <- y_max + y_step * i
      label_text <- paste0(pairwise_results$sig_label[i],
                           "\np.adj=", signif(pairwise_results$p_adj[i], 2),
                           " (n=", pairwise_results$n_pairs[i], ")")
      p_distribution <- p_distribution +
        annotate("segment", x = x1, xend = x2, y = y_bracket, yend = y_bracket, linewidth = 0.4) +
        annotate("segment", x = x1, xend = x1, y = y_bracket, yend = y_bracket - bracket_height, linewidth = 0.4) +
        annotate("segment", x = x2, xend = x2, y = y_bracket, yend = y_bracket - bracket_height, linewidth = 0.4) +
        annotate("text", x = (x1 + x2) / 2, y = y_bracket + bracket_height * 0.6,
                 label = label_text, size = (font_size * 0.35) / ggplot2::.pt, vjust = 0, lineheight = 0.9)
    }
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(output_dir, paste0(filename_stem, "_distribution", filename_suffix, ".png")), p_distribution, width = 9, height = 7, dpi = 600)
  ggsave(file.path(output_dir, paste0(filename_stem, "_summary", filename_suffix, ".png")), p_summary, width = 9, height = 6, dpi = 600)
  invisible(NULL)
}

data_processing_generate_outcome_var_deltas <- function(vo2_data_df, outcome_var_colname="cpet_vo2_adjusted_num"){
  # Ensure visits vector is available (from config.R)
  # Map numeric visits (0, 1, 2) to strings ("Baseline", "6m", "12m") if they are numeric
  vo2_data_df <- vo2_data_df %>%
    mutate(visit = as.character(visit)) %>%
    mutate(visit = case_when(
      visit == "0" ~ "Baseline",
      visit == "1" ~ "6m",
      visit == "2" ~ "12m",
      TRUE ~ visit
    ))

  # Transform vo2_data_df from long to wide format
  vo2_wide <- vo2_data_df %>%
    dplyr::select(PTID, visit, !!sym(outcome_var_colname)) %>%
    filter(!is.na(PTID) & !is.na(visit) & !is.na(!!sym(outcome_var_colname))) %>%
    group_by(PTID, visit) %>%
    summarise(VO2peak = mean(!!sym(outcome_var_colname), na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = visit, values_from = VO2peak) %>%
    # Calculate deltas for later visits relative to first visit (Baseline)
    mutate(
      vo2_delta_6m_Baseline  = if_else(!is.na(Baseline) & !is.na(`6m`), `6m` - Baseline, NA_real_),
      vo2_delta_12m_Baseline = if_else(!is.na(Baseline) & !is.na(`12m`), `12m` - Baseline, NA_real_)
    ) %>%
    dplyr::select(
      PTID,
      VO2peak_Baseline = Baseline,
      VO2peak_6m = `6m`,
      VO2peak_12m = `12m`,
      delta_VO2peak_6m_Baseline = vo2_delta_6m_Baseline,
      delta_VO2peak_12m_Baseline = vo2_delta_12m_Baseline
    )
    return(vo2_wide)
}

impute_covariates_patient_metadata <- function(patient_metadata_df, report_file, covars_num, covars_cat) {
  # Use median for numerical and mode for categorical
  # Add to report which patients needed imputation
  
  # Impute numerical covariates with median value if NA
  covars_num <- covars_num[covars_num %in% colnames(patient_metadata_df)]
  append_report_line(report_file, paste0("Imputing numerical covariates with median value if NA: ", paste(covars_num, collapse = ", ")))
  for (v in covars_num) {
    missing_idx <- which(is.na(patient_metadata_df[[v]]))
    if (length(missing_idx) > 0) {
      ptids <- patient_metadata_df$PTID[missing_idx]
      append_report_line(report_file, paste0("  - ", v, ": Imputed for PTID# ", paste(ptids, collapse = ", ")))
      med <- median(patient_metadata_df[[v]], na.rm = TRUE)
      patient_metadata_df[[v]][missing_idx] <- med
    } else {
      append_report_line(report_file, paste0("  - ", v, ": No missing values to impute."))
    }
  }

  # Impute categorical covariates with mode if NA
  covars_cat <- covars_cat[covars_cat %in% colnames(patient_metadata_df)]
  append_report_line(report_file, paste0("Imputing categorical covariates with mode if NA: ", paste(covars_cat, collapse = ", ")))
  for (v in covars_cat) {
    missing_idx <- which(is.na(patient_metadata_df[[v]]))
    if (length(missing_idx) > 0) {
      ptids <- patient_metadata_df$PTID[missing_idx]
      append_report_line(report_file, paste0("  - ", v, ": Imputed for PTID# ", paste(ptids, collapse = ", ")))
      mode_val <- names(sort(table(patient_metadata_df[[v]], useNA = "no"), decreasing = TRUE))[1]
      patient_metadata_df[[v]][missing_idx] <- mode_val
    } else {
      append_report_line(report_file, paste0("  - ", v, ": No missing values to impute."))
    }
  }

  return(patient_metadata_df)
}

join_patient_metadata_with_bioimpedance <- function(patient_metadata_df, bioimp_data_df)
{
  bioimp_data_df <- bioimp_data_df %>%
    dplyr::select(PTID, BFM, LBM, SMM, BMI, PBF) %>%
    mutate(PTID = as.integer(PTID))
  colnames(bioimp_data_df)[colnames(bioimp_data_df) == "BMI"] <- "bioimp_BMI"
  bioimp_data_df$PSMM <- round(bioimp_data_df$SMM / (bioimp_data_df$BFM + bioimp_data_df$LBM) * 100, 1)
  patient_metadata_df <- patient_metadata_df %>%
    dplyr::left_join(bioimp_data_df, by = "PTID")
  # Note, some patients lack bioimpedance data; downstream, drop NAs to handle this.

  return(patient_metadata_df)
}

build_data_for_pca <- function(data_df, patient_metadata_df_final, id_colname="OlinkID") {
  # Build sample metadata from PTID_visit column names and join TRANSPLANT_TYPE
  sample_cols <- setdiff(names(data_df), id_colname)
  sample_meta <- tibble::tibble(sample_id = sample_cols) %>%
    tidyr::separate(sample_id, into = c("PTID", "visit"), sep = "_", remove = FALSE, extra = "merge", fill = "right") %>%
    mutate(PTID = sub("^X", "", PTID)) %>%  # Remove leading X if present
    mutate(PTID = as.integer(PTID)) %>%
    dplyr::left_join(
      patient_metadata_df_final %>% dplyr::select(PTID, TRANSPLANT_TYPE, gender, VO2peak_Baseline, VO2peak_6m, VO2peak_12m, delta_VO2peak_6m_Baseline, delta_VO2peak_12m_Baseline, BMI, PBF, PSMM),
      by = "PTID"
    ) %>%
    dplyr::mutate(
      gender = dplyr::recode(gender, "M" = "Male", "F" = "Female"),
      VO2peak_visit_specific = dplyr::case_when(
        visit == "Baseline" ~ VO2peak_Baseline,
        visit == "6m"       ~ VO2peak_6m,
        visit == "12m"      ~ VO2peak_12m,
        TRUE ~ NA_real_
      ),
      delta_VO2peak_visit_specific = dplyr::case_when(
        visit == "Baseline" ~ 0,
        visit == "6m"       ~ delta_VO2peak_6m_Baseline,
        visit == "12m"      ~ delta_VO2peak_12m_Baseline,
        TRUE ~ NA_real_
      )
    )

    return(sample_meta)
}


run_all_pca_plot_generation <- function(data_df, sample_meta, output_dir, filename_suffix, subset_HCT_type, font_size, outcome_var = "VO2peak", id_colname = "OlinkID") {

  # Color points by visit
  generate_pca_plot(
    data_df = data_df,
    sample_meta = sample_meta,
    color_by = "visit",
    filename = file.path(output_dir, paste0("proteomics_pca_postfiltering_by_visit", filename_suffix, ".png")),
    text_size = font_size*0.6,
    point_size = 3,
    colors = timepoint_palette,
    legend_title = "Timepoint",
    id_colname = id_colname
  )

  # Color points by TRANSPLANT_TYPE
  if (subset_HCT_type == "ALL") {
  generate_pca_plot(
    data_df = data_df,
    sample_meta = sample_meta,
    color_by = "TRANSPLANT_TYPE",
    filename = file.path(output_dir, paste0("proteomics_pca_postfiltering_by_transplant_type", filename_suffix, ".png")),
    text_size = font_size*0.6,
    point_size = 3,
    colors = c("ALLO" = "#5fa2b6", "AUTO" = "#0a3468"),
    legend_title = "HCT Type",
    id_colname = id_colname
  )
  }

  # Color points by PTID
  generate_pca_plot(
    data_df = data_df,
    sample_meta = sample_meta,
    color_by = "PTID",
    filename = file.path(output_dir, paste0("proteomics_pca_postfiltering_by_ptid", filename_suffix, ".png")),
    text_size = font_size*0.6,
    point_size = 3,
    legend_title = "Patient ID",
    id_colname = id_colname
  )

  # Color points by outcome_var_Baseline
  generate_pca_plot(
    data_df = data_df,
    sample_meta = sample_meta,
    color_by = paste0(outcome_var, "_Baseline"),
    filename = file.path(output_dir, paste0("proteomics_pca_postfiltering_by_", tolower(outcome_var), "_baseline", filename_suffix, ".png")),
    text_size = font_size*0.6,
    point_size = 3,
    dramatic_color_scale = TRUE,
    id_colname = id_colname
  )

  # Color points by visit-specific outcome_var
  generate_pca_plot(
    data_df = data_df,
    sample_meta = sample_meta,
    color_by = paste0(outcome_var, "_visit_specific"),
    filename = file.path(output_dir, paste0("proteomics_pca_postfiltering_by_", tolower(outcome_var), "_visit", filename_suffix, ".png")),
    text_size = font_size*0.6,
    point_size = 3,
    dramatic_color_scale = TRUE,
    id_colname = id_colname
  )

  # Color points by visit-specific outcome_var delta
  generate_pca_plot(
    data_df = data_df,
    sample_meta = sample_meta,
    color_by = paste0("delta_", outcome_var, "_visit_specific"),
    filename = file.path(output_dir, paste0("proteomics_pca_postfiltering_by_delta_", tolower(outcome_var), "_visit", filename_suffix, ".png")),
    text_size = font_size*0.6,
    point_size = 3,
    dramatic_color_scale = TRUE,
    id_colname = id_colname
  )

  # Color points by gender
  generate_pca_plot(
    data_df = data_df,
    sample_meta = sample_meta,
    color_by = "gender",
    filename = file.path(output_dir, paste0("proteomics_pca_postfiltering_by_gender", filename_suffix, ".png")),
    text_size = font_size*0.6,
    point_size = 3,
    legend_title = "Sex",
    id_colname = id_colname
  )

  # Color by patient BMI
  generate_pca_plot(
    data_df = data_df,
    sample_meta = sample_meta,
    color_by = "BMI",
    filename = file.path(output_dir, paste0("proteomics_pca_postfiltering_by_bmi", filename_suffix, ".png")),
    text_size = font_size*0.6,
    point_size = 3,
    dramatic_color_scale = TRUE,
    legend_title = "BMI",
    id_colname = id_colname
  )

  if ("PBF" %in% colnames(sample_meta)) {
    # Color by patient PBF
    generate_pca_plot(
      data_df = data_df,
      sample_meta = sample_meta,
      color_by = "PBF",
      filename = file.path(output_dir, paste0("proteomics_pca_postfiltering_by_pbf", filename_suffix, ".png")),
      text_size = font_size*0.6,
      point_size = 3,
      dramatic_color_scale = TRUE,
      legend_title = "PBF",
      id_colname = id_colname
    )
  }
  
  if ("PSMM" %in% colnames(sample_meta)) {
    # Color by patient PSMM
    generate_pca_plot(
      data_df = data_df,
      sample_meta = sample_meta,
      color_by = "PSMM",
      filename = file.path(output_dir, paste0("proteomics_pca_postfiltering_by_psmm", filename_suffix, ".png")),
      text_size = font_size*0.6,
      point_size = 3,
      dramatic_color_scale = TRUE,
      legend_title = "PSMM",
      id_colname = id_colname
    )
  }
}
