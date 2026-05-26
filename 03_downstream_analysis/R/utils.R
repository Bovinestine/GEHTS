# =============================================================================
# utils.R — Shared plotting utilities for GEHTS downstream analysis
#
# Functions that appear in multiple figure scripts are defined here once.
# Source after config.R:
#   source(here::here("03_downstream_analysis", "R", "utils.R"))
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(ggsci)
  library(ggrepel)
  library(gridExtra)
  library(cowplot)
  library(ggsignif)
})

# -----------------------------------------------------------------------------
# 1. Heatmap of SCTransform residuals (Figures 2b, 3b)
# -----------------------------------------------------------------------------

#' Gene Expression Heatmap Categorized by Phenotype
#'
#' Extracts SCTransform Pearson residuals from a Seurat object, filters by
#' expression quantile, and generates a pheatmap with anabolic/catabolic/
#' housekeeping row annotations. Used for both primary-cell and drug-screen
#' heatmaps by changing the grouping metadata columns.
#'
#' @param seurat_object A processed Seurat object with SCT assay.
#' @param upper_threshold Upper quantile for gene expression filtering (default 0).
#' @param bottom_threshold Lower quantile for gene expression filtering (default 1).
#' @param assay Assay name (default "SCT").
#' @param primary_metadata Column for primary cell grouping (default "drug_condition").
#' @param secondary_metadata Column for secondary grouping (default "file").
#' @param genes_anabolic Anabolic gene panel (default: GENES_ANABOLIC from config).
#' @param genes_catabolic Catabolic gene panel (default: GENES_CATABOLIC from config).
#' @param genes_housekeeping Housekeeping gene panel (default: GENES_HOUSEKEEPING from config).
#' @return A pheatmap object.
#' @export
create_heatmap_by_expression_ps <- function(seurat_object,
                                            upper_threshold = 0,
                                            bottom_threshold = 1,
                                            assay = "SCT",
                                            primary_metadata = "drug_condition",
                                            secondary_metadata = "file",
                                            genes_anabolic = GENES_ANABOLIC,
                                            genes_catabolic = GENES_CATABOLIC,
                                            genes_housekeeping = GENES_HOUSEKEEPING) {

  # Use counts for quantile-based gene filtering (residuals have mean ~0)
  count_data <- GetAssayData(seurat_object, assay = assay, slot = "counts")
  count_data <- count_data[, colSums(is.finite(count_data)) > 0]
  avg_exp <- rowMeans(count_data, na.rm = TRUE)

  q_bottom <- quantile(avg_exp, 1 - bottom_threshold)
  q_top    <- quantile(avg_exp, 1 - upper_threshold)
  sel_genes <- names(avg_exp[avg_exp >= q_bottom & avg_exp <= q_top])

  sel_anabolic <- intersect(sel_genes, genes_anabolic)
  sel_catabolic <- intersect(sel_genes, genes_catabolic)
  sel_housek   <- intersect(sel_genes, genes_housekeeping)

  sorted_genes <- c(
    sel_anabolic[order(avg_exp[sel_anabolic], decreasing = TRUE)],
    sel_catabolic[order(avg_exp[sel_catabolic], decreasing = TRUE)],
    sel_housek[order(avg_exp[sel_housek], decreasing = TRUE)]
  )

  cell_scores <- data.frame(
    cell_name = rownames(seurat_object@meta.data),
    primary   = seurat_object@meta.data[[primary_metadata]],
    secondary = seurat_object@meta.data[[secondary_metadata]],
    stringsAsFactors = FALSE
  ) %>% arrange(primary, secondary)

  # Pearson residuals are already centered — use directly
  scale_data <- GetAssayData(seurat_object, assay = assay, slot = "scale.data")
  avail_genes <- intersect(sorted_genes, rownames(scale_data))

  if (length(avail_genes) < length(sorted_genes)) {
    missing <- setdiff(sorted_genes, avail_genes)
    warning("Genes missing from scale.data (skipped): ",
            paste(missing, collapse = ", "))
  }

  mat <- as.matrix(scale_data[avail_genes, cell_scores$cell_name])
  mat[is.na(mat)] <- 0
  # Cap extreme residuals so outliers don't wash out the color scale
  mat <- pmax(pmin(mat, LOG2FC_CAP), -LOG2FC_CAP)

  max_val <- max(abs(mat), na.rm = TRUE)
  breaks  <- seq(-max_val, max_val, length.out = 101)
  colors  <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

  ann_col <- data.frame(
    Primary   = factor(cell_scores$primary,   levels = unique(cell_scores$primary)),
    Secondary = factor(cell_scores$secondary, levels = unique(cell_scores$secondary))
  )
  rownames(ann_col) <- cell_scores$cell_name

  final_ana  <- intersect(avail_genes, sel_anabolic)
  final_cat  <- intersect(avail_genes, sel_catabolic)
  final_hk   <- intersect(avail_genes, sel_housek)

  ann_row <- data.frame(
    Gene_Type = factor(
      rep(c("Anabolic", "Catabolic", "Housekeeping"),
          c(length(final_ana), length(final_cat), length(final_hk))),
      levels = c("Anabolic", "Catabolic", "Housekeeping")
    )
  )
  rownames(ann_row) <- avail_genes

  n_primary   <- length(unique(cell_scores$primary))
  n_secondary <- length(unique(cell_scores$secondary))

  ann_colors <- list(
    Gene_Type = c(Anabolic    = COL_ANABOLIC,
                  Catabolic   = COL_CATABOLIC,
                  Housekeeping = COL_HOUSEKEEP),
    Primary   = setNames(colorRampPalette(pal_jco()(10))(n_primary),
                         unique(cell_scores$primary)),
    Secondary = setNames(colorRampPalette(PALETTE_NPG)(n_secondary),
                         unique(cell_scores$secondary))
  )

  gap_rows <- c(length(final_ana), length(final_ana) + length(final_cat))

  pheatmap::pheatmap(
    mat               = mat,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    show_rownames     = TRUE,
    show_colnames     = FALSE,
    annotation_col    = ann_col,
    annotation_row    = ann_row,
    annotation_colors = ann_colors,
    color             = colors,
    breaks            = breaks,
    main              = "SCTransform Residuals",
    gaps_row          = gap_rows
  )
}

# -----------------------------------------------------------------------------
# 2. Dimension reduction plot — t-SNE or UMAP (Figures 2e, 3c)
# -----------------------------------------------------------------------------

#' Plot t-SNE or UMAP from a Seurat Object
#'
#' Extracts embeddings and colors points by a metadata column. Applies the NPG
#' palette with "control" locked to blue and "inflammatory" locked to red.
#'
#' @param seurat_object Processed Seurat object with the requested reduction.
#' @param reduction Reduction name: "tsne" or "umap" (default "tsne").
#' @param group_by Metadata column to color by (default "drug_condition").
#' @param font_size Base font size (default 9).
#' @return A ggplot2 object.
#' @export
plot_dim_reduction <- function(seurat_object,
                               reduction = "tsne",
                               group_by  = "drug_condition",
                               font_size = 9) {

  if (!reduction %in% names(seurat_object@reductions)) {
    stop("Reduction '", reduction, "' not found. Run the appropriate Seurat function first.")
  }
  if (!group_by %in% colnames(seurat_object@meta.data)) {
    stop("Column '", group_by, "' not found in metadata.")
  }

  coords <- as.data.frame(Embeddings(seurat_object, reduction = reduction))
  axis_labels <- toupper(sub("([a-z]+)(.*)", "\\1", reduction))
  colnames(coords) <- paste0(axis_labels, "_", 1:2)
  coords[[group_by]] <- as.character(seurat_object@meta.data[rownames(coords), group_by])

  locked <- c(COND_CONTROL = COL_CONTROL, COND_DISEASE = COL_DISEASE)
  names(locked) <- c(COND_CONTROL, COND_DISEASE)

  other_conds <- setdiff(unique(coords[[group_by]]), names(locked))
  extra_cols  <- if (length(other_conds) > 0) {
    avail <- setdiff(PALETTE_NPG, locked)
    if (length(other_conds) > length(avail)) {
      avail <- colorRampPalette(avail)(length(other_conds))
    }
    setNames(avail[seq_along(other_conds)], other_conds)
  } else character(0)

  color_map <- c(locked, extra_cols)

  ggplot(coords, aes(x = .data[[colnames(coords)[1]]],
                     y = .data[[colnames(coords)[2]]],
                     color = .data[[group_by]])) +
    geom_point(alpha = 0.8, size = 1) +
    scale_color_manual(values = color_map) +
    theme_classic(base_size = font_size) +
    labs(x = paste(axis_labels, "1"), y = paste(axis_labels, "2"),
         color = "Condition") +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          legend.position = "right")
}

# -----------------------------------------------------------------------------
# 3. GO module score calculator (Figures 4e, 5e, 5f)
# -----------------------------------------------------------------------------

#' Compute Mean GO Module Scores Per Cell
#'
#' Adds each GO module as a per-cell metadata column in a Seurat object.
#'
#' @param seurat_obj A Seurat object.
#' @param module_list Named list of gene vectors (default: GO_MODULES from config).
#' @param assay Assay to use (default "SCT").
#' @param slot Data slot (default "data").
#' @return The Seurat object with module scores added to metadata.
#' @export
compute_go_scores <- function(seurat_obj,
                              module_list = GO_MODULES,
                              assay = "SCT",
                              slot  = "data") {

  mat <- GetAssayData(seurat_obj, assay = assay, slot = slot)

  calc_mean <- function(genes) {
    valid <- intersect(genes, rownames(mat))
    if (length(valid) == 0) return(rep(0, ncol(mat)))
    if (length(valid) == 1) return(as.numeric(mat[valid, ]))
    colMeans(as.matrix(mat[valid, , drop = FALSE]))
  }

  for (mod in names(module_list)) {
    seurat_obj[[mod]] <- calc_mean(module_list[[mod]])
  }

  seurat_obj
}

# -----------------------------------------------------------------------------
# 4. Bliss synergy helpers (Figures 4c, 5d)
# -----------------------------------------------------------------------------

#' Bliss Independence Predicted Efficacy
#' @param e1 Efficacy of drug A (0–1).
#' @param e2 Efficacy of drug B (0–1).
#' @return Bliss-predicted combination efficacy.
bliss_predict <- function(e1, e2) e1 + e2 - (e1 * e2)

#' Split "&"-Delimited Combination Names Into Drug Pair Columns
#' @param prob_vec Named numeric vector of combination efficacy scores.
#' @param dose Character label for the dose (used as a column).
#' @return A data.frame with Drug1Name, Drug2Name, MeanProbability, Dose columns.
split_combo_names <- function(prob_vec, dose) {
  data.frame(
    Combination     = names(prob_vec),
    MeanProbability = as.numeric(prob_vec),
    Dose            = dose,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Drug1Name = sapply(strsplit(Combination, "&"), `[`, 1),
      Drug2Name = sapply(strsplit(Combination, "&"), `[`, 2)
    )
}

# -----------------------------------------------------------------------------
# 5. Faceted gene-expression boxplots with common legend (Figs 2 & 3 supplements)
# -----------------------------------------------------------------------------

#' Publication-Ready Faceted Boxplots with Common Legend
#'
#' Generates one boxplot per gene with Wilcoxon significance brackets, arranged
#' in a grid with a shared legend at the top and a shared Y-axis label.
#'
#' @param seurat_object A Seurat object.
#' @param genes Character vector of genes to plot.
#' @param datatype Data slot (default "counts").
#' @param control_label Metadata string for healthy baseline (default COND_CONTROL).
#' @param treatment_label Metadata string for disease baseline (default COND_DISEASE).
#' @param extra_label1 Third condition label (default "label3").
#' @param extra_label2 Fourth condition label (default "label4").
#' @param color_palette Named color vector. If NULL uses the default four-color palette.
#' @param nrow Number of rows in the gene grid (default 5).
#' @param base_font_size Base font size (default 6).
#' @param line_thickness Axis/boxplot line width (default 0.5).
#' @param legend_gap Fractional height reserved for the legend (default 0.05).
#' @param y_title Shared Y-axis title.
#' @return A combined cowplot object.
#' @export
draw_boxplots_for_genes_with_common_legend <- function(
    seurat_object,
    genes,
    datatype        = "counts",
    control_label   = COND_CONTROL,
    treatment_label = COND_DISEASE,
    extra_label1    = "label3",
    extra_label2    = "label4",
    color_palette   = NULL,
    nrow            = 5,
    base_font_size  = 6,
    line_thickness  = 0.5,
    legend_gap      = 0.05,
    y_title         = "Normalized gene expression (a.u.)") {

  if (is.null(color_palette)) {
    color_palette <- setNames(
      c("#4DBBD559", "#E64B3599", "#A6CEE3", "#1F78B4"),
      c("Basal", "IL-1β", extra_label1, extra_label2)
    )
  }

  meta <- seurat_object@meta.data
  valid_conds  <- c(control_label, treatment_label, extra_label1, extra_label2)
  cells_to_use <- rownames(meta)[meta$drug_condition %in% valid_conds]
  level_order  <- c("Basal", "IL-1β", extra_label1, extra_label2)

  meta_sub <- meta[cells_to_use, ] %>%
    mutate(ConditionLabel = case_when(
      drug_condition == control_label   ~ "Basal",
      drug_condition == treatment_label ~ "IL-1β",
      drug_condition == extra_label1    ~ extra_label1,
      drug_condition == extra_label2    ~ extra_label2,
      TRUE ~ as.character(drug_condition)
    )) %>%
    mutate(ConditionLabel = factor(ConditionLabel, levels = level_order))

  plot_list <- lapply(genes, function(gene) {
    expr <- as.vector(GetAssayData(seurat_object,
                                   slot = datatype)[gene, cells_to_use, drop = TRUE])
    df <- data.frame(Expression = expr, Condition = meta_sub$ConditionLabel)

    ymax <- max(df$Expression, na.rm = TRUE)
    if (!is.finite(ymax) || ymax == 0) ymax <- 1
    # 8% of ymax per significance bracket — derived from chip well spacing geometry
    step <- ymax * 0.08

    p <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
      geom_boxplot(outlier.shape = NA, linewidth = line_thickness, alpha = 0.85) +
      scale_fill_manual(values = color_palette) +
      labs(title = gene) +
      theme_classic() +
      theme(text          = element_text(size = base_font_size),
            axis.line     = element_line(linewidth = line_thickness, color = "black"),
            axis.text.x   = element_blank(),
            axis.ticks.x  = element_blank(),
            axis.text.y   = element_text(size = base_font_size, color = "black"),
            axis.title    = element_blank(),
            plot.title    = element_text(face = "bold.italic", hjust = 0.5,
                                         size = base_font_size + 2),
            legend.position = "none",
            plot.margin   = ggplot2::margin(10, 5, 2, 5))

    group_counts <- table(df$Condition[!is.na(df$Expression)])
    valid_groups <- names(group_counts)[group_counts >= 2]

    if (length(valid_groups) >= 2) {
      comps <- combn(valid_groups, 2, simplify = FALSE)
      p <- p + geom_signif(
        comparisons      = comps,
        map_signif_level = TRUE,
        y_position       = ymax + step * seq_along(comps),
        textsize         = base_font_size * 0.35,
        tip_length       = 0.01,
        size             = 0.3,
        vjust            = 0.5
      ) + coord_cartesian(ylim = c(0, ymax + step * (length(comps) + 1)))
    } else {
      p <- p + coord_cartesian(ylim = c(0, ymax * 1.1))
    }
    p
  })
  names(plot_list) <- genes

  dummy_df  <- data.frame(Condition  = factor(level_order, levels = level_order),
                           Expression = rep(0.5, 4))
  dummy_plt <- ggplot(dummy_df, aes(Condition, Expression, fill = Condition)) +
    geom_boxplot() +
    scale_fill_manual(values = color_palette) +
    theme_classic() +
    theme(legend.position = "top",
          legend.text     = element_text(size = base_font_size + 2, face = "bold"),
          legend.title    = element_blank(),
          legend.key.size = unit(0.4, "cm"),
          legend.margin   = ggplot2::margin(0, 0, 0, 0))

  get_legend <- function(plt) {
    tmp <- ggplot_gtable(ggplot_build(plt))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    tmp$grobs[[leg]]
  }
  common_legend <- get_legend(dummy_plt)

  ncol_val   <- ceiling(length(genes) / nrow)
  grid_plots <- arrangeGrob(grobs = plot_list, ncol = ncol_val, nrow = nrow)
  arranged   <- grid.arrange(common_legend, grid_plots, ncol = 1,
                              heights = c(legend_gap, 1))

  ggdraw() +
    draw_plot(arranged, x = 0.05, y = 0, width = 0.95, height = 1) +
    draw_label(y_title, x = 0.02, y = 0.5, angle = 90, vjust = 0.5,
               size = base_font_size + 2, fontface = "bold")
}

# -----------------------------------------------------------------------------
# 6. Radar chart scale calculator (Figure 4e)
# -----------------------------------------------------------------------------

#' Calculate Universal Radar Chart Scale Across All Conditions
#'
#' Scans every condition in the Seurat object to find the absolute min/max
#' expression per GO module, so all radar charts share an identical axis scale.
#'
#' @param seurat_obj A Seurat object.
#' @param module_list Named list of gene modules (default: GO_MODULES from config).
#' @param group_by Metadata column containing conditions (default "drug_condition").
#' @param assay Assay name (default "SCT").
#' @param slot Slot name (default "counts").
#' @return A list with elements `min` and `max` (named numeric vectors).
#' @export
calculate_radar_scale <- function(seurat_obj,
                                  module_list = GO_MODULES,
                                  group_by    = "drug_condition",
                                  assay       = "SCT",
                                  slot        = "counts") {

  if (!inherits(seurat_obj, "Seurat")) stop("'seurat_obj' must be a Seurat object.")
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop("Column '", group_by, "' not found in metadata.")
  }

  DefaultAssay(seurat_obj) <- assay
  mat  <- GetAssayData(seurat_obj, slot = slot)
  meta <- seurat_obj@meta.data

  all_metrics <- lapply(unique(meta[[group_by]]), function(cond) {
    cells <- rownames(meta[meta[[group_by]] == cond, ])
    if (length(cells) == 0) return(NULL)
    sapply(module_list, function(genes) {
      avail <- intersect(genes, rownames(mat))
      if (length(avail) > 0) mean(rowMeans(mat[avail, cells, drop = FALSE])) else NA
    })
  })
  all_metrics <- do.call(rbind, Filter(Negate(is.null), all_metrics))

  list(
    min = apply(all_metrics, 2, min, na.rm = TRUE),
    max = apply(all_metrics, 2, max, na.rm = TRUE)
  )
}

# -----------------------------------------------------------------------------
# 6. Per-drug highlighted UMAP grid (Figures 4d, Extended Data)
# -----------------------------------------------------------------------------

plot_highlighted_umaps <- function(seurat_obj,
                                   plot_indices = NULL,
                                   output_pdf,
                                   ncol       = 2,
                                   pdf_width  = 6,
                                   pdf_height = 9) {

  umap_coords <- Embeddings(seurat_obj, reduction = "umap")
  umap_df     <- data.frame(umap_coords,
                             drug_name = seurat_obj@meta.data$drug_name,
                             dose      = seurat_obj@meta.data$dose,
                             target    = seurat_obj@meta.data$target)
  colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")

  drugs     <- unique(umap_df$drug_name)
  all_plots <- lapply(drugs, function(drug) {
    ggplot(umap_df, aes(UMAP_1, UMAP_2,
                        color = (drug_name == drug),
                        alpha = (drug_name == drug))) +
      geom_point(size = 1) +
      scale_alpha_manual(values = c("FALSE" = 0.1, "TRUE" = 0.6)) +
      scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "turquoise3")) +
      labs(title = drug) +
      theme_classic() +
      theme(legend.position = "none",
            plot.title      = element_text(face = "bold", size = 10))
  })
  names(all_plots) <- drugs

  to_print <- if (is.null(plot_indices)) all_plots else all_plots[plot_indices]

  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  pdf(output_pdf, width = pdf_width, height = pdf_height)
  grid.arrange(grobs = to_print, ncol = ncol)
  dev.off()
  invisible(output_pdf)
}
