# =============================================================================
# figure2_supplement.R — Supplementary figures for Figure 2
#
# Subfigures:
#   Extended Data Fig. 5a: Class-summed gene-count boxplot (anabolic/catabolic)
#   Extended Data Fig. 5b: Edge-effect boxplot by chip region
#   Extended Data Fig. 5c: Spatial heatmap of chip uniformity
#   Extended Data Fig. 6: Grid of gene-expression boxplots (all 30 genes)
#   Extended Data Fig. 7: CLR-normalized heatmap — Primary cells (mac.sct,
#     control wells) vs. ATDC5 (atdc5.pub, all 4 37C samples)
#   Supplementary (unnumbered): Lollipop chart — rna.pub (bulk RNA-seq) vs.
#     ATDC5, the rna.pub-based companion to Figure 2f (which uses mac.sct as
#     the primary side); see build_rna_pub_from_featurecounts.R /
#     build_atdc5_from_gse269585.R for data provenance.
#
# Prerequisites: config.R, utils.R, data_loading.R
# Inputs: mac.sct, rna.pub, atdc5.pub (loaded by data_loading.R)
# Outputs: output/supplementary/Sfig{5a,5b,5c,6,7}_*.pdf, Sfig_rna_vs_atdc5_lollipop.pdf
#   (Fig. 6 also emits three comparison variants for the flat/zero-inflated-gene
#   axis issue — Sfig6_boxplots_allgenes_zoomed.pdf (whisker-zoomed),
#   Sfig6_boxplots_allgenes_jitter.pdf (full range + jittered points), and
#   Sfig6_violins_allgenes.pdf (violin + inner quartile boxplot); see
#   draw_boxplots_for_genes_with_common_legend_{zoomed,jitter}() and
#   draw_violins_for_genes_with_common_legend() in utils.R)
#
# NOTE: Extended Data Fig. 7 previously referenced primary.sct/atdc5.sct,
# which are not produced by any script in this repo. Fixed to use mac.sct
# (control-condition wells only, matching Figure 2f) and atdc5.pub (all
# samples -- atdc5.pub already contains only the 4 37C/baseline replicates,
# see build_atdc5_from_gse269585.R). Deliberately CLR-normalized straight
# from raw counts, NOT run through SCTransform first, for the same reason
# Figure 2f and the rna.pub/ATDC5 lollipop above skip it: this is the
# heatmap counterpart of that same primary-vs-ATDC5 relative-expression
# comparison, so both sides should be normalized identically for a fair
# comparison, and atdc5.pub is only 4 bulk-RNA-seq-like samples -- too few
# for SCTransform's per-gene mean-variance regularization (designed for
# thousands of single cells) to be meaningful anyway.
# =============================================================================

source(here::here("03_downstream_analysis", "R", "config.R"))
source(here::here("03_downstream_analysis", "R", "utils.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(viridis)
})

set.seed(SEED)

# =============================================================================
# Extended Data Figure 5a — Class-summed gene-count boxplot
# (Anabolic / Catabolic panels, Basal vs. IL-1β; same structure as Extended
# Data Fig. 6, but each panel sums counts across a gene set rather than
# plotting a single gene.)
# =============================================================================

sfig_5a_genesum <- draw_boxplots_for_gene_sums_with_common_legend(
  mac.sct,
  gene_sets      = list(Anabolic = GENES_ANABOLIC, Catabolic = GENES_CATABOLIC),
  nrow           = 1,
  base_font_size = 8,
  legend_gap     = 0.08
)
ggsave(file.path(OUTPUT_SUP, "Sfig5a_genesum_boxplot.pdf"),
       plot = sfig_5a_genesum, width = 4, height = 4, dpi = 300)

# =============================================================================
# Extended Data Figure 5b — Edge-effect boxplot
# Extended Data Figure 5c — Spatial heatmap of chip uniformity
# =============================================================================

#' Map Well Number to 10×10 Grid Coordinates (Snake Pattern)
#'
#' The GEHTS chip uses a 10×10 well grid loaded in a snake (boustrophedon)
#' pattern: odd rows go left-to-right, even rows go right-to-left.
#'
#' @param well_no Well number (integer 1–100) or a string like "no42".
#' @return A list with elements Row, Col, Region ("Center"/"Edge"/"Corner").
get_spatial_coords <- function(well_no) {
  w_clean <- gsub("no|\\s", "", as.character(well_no), ignore.case = TRUE)
  w <- suppressWarnings(as.numeric(w_clean))
  if (is.na(w)) return(list(Row = NA, Col = NA, Region = NA))

  row_idx <- ceiling(w / 10)
  col_idx <- if (row_idx %% 2 != 0) {
    (w - 1) %% 10 + 1          # odd rows: left to right
  } else {
    10 - ((w - 1) %% 10)       # even rows: right to left
  }

  region <- if ((row_idx %in% c(1, 10)) && (col_idx %in% c(1, 10))) "Corner"
             else if (row_idx %in% c(1, 10) || col_idx %in% c(1, 10)) "Edge"
             else "Center"

  list(Row = row_idx, Col = col_idx, Region = region)
}

meta_df <- mac.sct@meta.data

if (!"Total_Counts" %in% colnames(meta_df)) {
  if ("nCount_RNA" %in% colnames(meta_df)) {
    meta_df$Total_Counts <- meta_df$nCount_RNA
  } else {
    stop("Column 'nCount_RNA' not found in Seurat metadata.")
  }
}

spatial_list <- lapply(meta_df$well_no, get_spatial_coords)
spatial_df   <- dplyr::bind_rows(spatial_list)
plot_data    <- cbind(meta_df, spatial_df)
plot_data    <- plot_data[!is.na(plot_data$Row), ]
plot_data$ChipID <- paste(plot_data$file, plot_data$drug_condition, sep = "_")
plot_data$Region <- factor(plot_data$Region, levels = c("Center", "Edge", "Corner"))

sfig_5b_edge <- ggplot(plot_data, aes(x = Region, y = Total_Counts, fill = Region)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.4) +
  facet_wrap(~ChipID, scales = "free_y") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c(Center = COL_CENTER, Edge = COL_EDGE,
                                Corner = COL_CORNER)) +
  labs(title    = "Edge Effects on Gene Detection Sensitivity",
       subtitle  = "Total Counts by Chip Region",
       y = "Total Gene Counts",
       x = "") +
  theme(legend.position    = "none",
        strip.background   = element_rect(fill = "white", color = "black"),
        strip.text         = element_text(face = "bold"))

ggsave(file.path(OUTPUT_SUP, "Sfig5b_edge_effect_boxplot.pdf"),
       plot = sfig_5b_edge, width = 8, height = 6, dpi = 300)

sfig_5c_spatial <- ggplot(plot_data, aes(x = Col, y = Row, fill = Total_Counts)) +
  geom_tile(color = "white", lwd = 0.2, width = 1, height = 1) +
  facet_wrap(~ChipID, ncol = 4) +
  scale_y_reverse(breaks = 1:10, limits = c(10.5, 0.5)) +
  scale_x_continuous(breaks = 1:10, limits = c(0.5, 10.5)) +
  scale_fill_viridis_c(option = "viridis", name = "Counts") +
  labs(title = "Spatial Uniformity (Raw Counts)", x = "Column", y = "Row") +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid  = element_blank(),
        strip.text  = element_text(face = "bold", size = 10),
        axis.text   = element_text(size = 8))

ggsave(file.path(OUTPUT_SUP, "Sfig5c_spatial_uniformity_heatmap.pdf"),
       plot = sfig_5c_spatial, width = 12, height = 10, dpi = 300)

# =============================================================================
# Extended Data Figure 6 — Boxplots of all 30 genes across 4 conditions
# =============================================================================

sfig_6_boxplots <- draw_boxplots_for_genes_with_common_legend(
  mac.sct,
  genes          = GENES_ALL,
  nrow           = 4,
  base_font_size = 6,
  legend_gap     = 0.03
)
ggsave(file.path(OUTPUT_SUP, "Sfig6_boxplots_allgenes.pdf"),
       plot = sfig_6_boxplots, width = 5, height = 5, dpi = 300)

# Whisker-zoomed variant for comparison -- some genes (e.g. Ucma, Il6,
# Adamts5, Cxcl5) have extreme outlier cells that stretch the raw-max ylim
# enough to flatten a real, significant IQR box near zero even though the
# outlier points themselves are already hidden (outlier.shape = NA). This
# version zooms to the upper Tukey whisker instead; kept as a separate output
# so both can be compared before deciding which to keep for Ext. Data Fig. 6.
sfig_6_boxplots_zoomed <- draw_boxplots_for_genes_with_common_legend_zoomed(
  mac.sct,
  genes          = GENES_ALL,
  nrow           = 4,
  base_font_size = 6,
  legend_gap     = 0.03
)
ggsave(file.path(OUTPUT_SUP, "Sfig6_boxplots_allgenes_zoomed.pdf"),
       plot = sfig_6_boxplots_zoomed, width = 5, height = 5, dpi = 300)

# Jittered-points variant for comparison -- for zero-inflated genes (Q1 =
# median = Q3 = 0: Ucma, Il6, Adamts5, Cxcl5) the box collapses to a point at
# zero regardless of axis scaling, so whisker-zooming just crops the very
# cells that explain a significant result. This keeps the full-range axis
# (like the original) and overlays individual cells so the shift in
# detection rate / elevated counts is visible by eye instead of hidden.
sfig_6_boxplots_jitter <- draw_boxplots_for_genes_with_common_legend_jitter(
  mac.sct,
  genes          = GENES_ALL,
  nrow           = 4,
  base_font_size = 6,
  legend_gap     = 0.03
)
ggsave(file.path(OUTPUT_SUP, "Sfig6_boxplots_allgenes_jitter.pdf"),
       plot = sfig_6_boxplots_jitter, width = 5, height = 5, dpi = 300)

# Violin variant for comparison -- density shape (violin width) still encodes
# a distributional shift for zero-inflated genes even when the box/quartiles
# collapse to a point at zero; thin inner boxplot keeps the quartiles visible.
sfig_6_violins <- draw_violins_for_genes_with_common_legend(
  mac.sct,
  genes          = GENES_ALL,
  nrow           = 4,
  base_font_size = 6,
  legend_gap     = 0.03
)
ggsave(file.path(OUTPUT_SUP, "Sfig6_violins_allgenes.pdf"),
       plot = sfig_6_violins, width = 5, height = 5, dpi = 300)

# =============================================================================
# Supplementary — Lollipop chart: rna.pub (bulk RNA-seq) vs. ATDC5
# =============================================================================
# Companion to Figure 2f, which uses mac.sct (GE-HTS chip) as the primary
# side; this version uses rna.pub (bulk RNA-seq, control replicates) instead,
# so the comparison is reported on both platforms. plot_gene_divergence() is
# defined in utils.R (shared with figure2.R's Figure 2f call). Placed before
# Extended Data Fig. 7 below, which is already broken (see NOTE in the header
# comment) and would otherwise stop this script before this section runs.

rna.pub.ctrl.clr <- NormalizeData(subset(rna.pub, drug_condition == COND_CONTROL),
                                  assay = "RNA", normalization.method = "CLR",
                                  verbose = FALSE)
atdc5.clr.sup     <- NormalizeData(atdc5.pub, assay = "RNA",
                                   normalization.method = "CLR", verbose = FALSE)

sfig_rna_atdc5_lollipop <- plot_gene_divergence(
  GetAssayData(rna.pub.ctrl.clr, assay = "RNA", layer = "data"),
  GetAssayData(atdc5.clr.sup,    assay = "RNA", layer = "data")
)
ggsave(file.path(OUTPUT_SUP, "Sfig_rna_vs_atdc5_lollipop.pdf"),
       plot = sfig_rna_atdc5_lollipop, width = 6, height = 5, dpi = 300)

# =============================================================================
# Extended Data Figure 7 — CLR heatmap: Primary cells vs. ATDC5
# =============================================================================

plot_clr_heatmap <- function(seu_primary, seu_atdc5, genes_of_interest,
                             assay      = "RNA",
                             layer      = "counts",
                             pdf_file   = NULL,
                             pdf_width  = 8,
                             pdf_height = 8) {
  valid_genes <- intersect(genes_of_interest, rownames(seu_primary[[assay]]))
  valid_genes <- intersect(valid_genes,        rownames(seu_atdc5[[assay]]))

  counts_primary <- GetAssayData(seu_primary, assay = assay, layer = layer)[valid_genes, , drop = FALSE]
  counts_atdc5   <- GetAssayData(seu_atdc5,   assay = assay, layer = layer)[valid_genes, , drop = FALSE]

  clr_norm <- function(x) { geo_mean <- exp(mean(log(x + 1))); log((x + 1) / geo_mean) }
  norm_combined <- cbind(apply(counts_primary, 2, clr_norm),
                         apply(counts_atdc5,   2, clr_norm))

  annotation_df <- data.frame(
    CellType = c(rep("Primary", ncol(counts_primary)), rep("ATDC5", ncol(counts_atdc5)))
  )
  rownames(annotation_df) <- colnames(norm_combined)

  pheatmap(
    mat               = norm_combined,
    color             = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    breaks            = seq(-2, 2, length.out = 100),
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    annotation_col    = annotation_df,
    annotation_colors = list(CellType = c(Primary = "#00A087", ATDC5 = "#F39B7F")),
    show_colnames     = FALSE,
    main              = paste0("CLR Relative Expression (", length(valid_genes), "-gene panel)"),
    fontsize_row      = 10,
    filename          = pdf_file,
    width             = pdf_width,
    height            = pdf_height
  )
}

# Primary side: mac.sct's control-condition wells only, matching Figure 2f /
# the Sfig lollipop above -- recomputed locally here (rather than reusing
# figure2.R's mac.sct.ctrl) so this script's Ext. Data Fig. 7 section stays
# self-contained if run on its own.
mac.sct.ctrl.sfig7 <- subset(mac.sct, drug_condition == COND_CONTROL)

plot_clr_heatmap(
  mac.sct.ctrl.sfig7,
  atdc5.pub,
  genes_of_interest = GENES_ALL,
  pdf_file          = file.path(OUTPUT_SUP, "Sfig7_clr_heatmap_primary_atdc5.pdf"),
  pdf_width         = 8,
  pdf_height        = 8
)

save_session_info("figure2_supplement")

