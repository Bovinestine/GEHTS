# =============================================================================
# figure2_supplement.R — Supplementary figures for Figure 2
#
# Subfigures:
#   Extended Data Fig. 5a: Edge-effect boxplot by chip region
#   Extended Data Fig. 5b: Spatial heatmap of chip uniformity
#   Extended Data Fig. 6: Grid of gene-expression boxplots (all 30 genes)
#   Extended Data Fig. 7: CLR-normalized heatmap — Primary cells vs. ATDC5
#
# Prerequisites: config.R, utils.R, data_loading.R
# Inputs: mac.sct, primary.sct, atdc5.sct (loaded by data_loading.R)
# Outputs: output/supplementary/Sfig{5a,5b,6,7}_*.pdf
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
# Extended Data Figure 5a — Edge-effect boxplot
# Extended Data Figure 5b — Spatial heatmap of chip uniformity
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

sfig_5a_edge <- ggplot(plot_data, aes(x = Region, y = Total_Counts, fill = Region)) +
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

ggsave(file.path(OUTPUT_SUP, "Sfig5a_edge_effect_boxplot.pdf"),
       plot = sfig_5a_edge, width = 8, height = 6, dpi = 300)

sfig_5b_spatial <- ggplot(plot_data, aes(x = Col, y = Row, fill = Total_Counts)) +
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

ggsave(file.path(OUTPUT_SUP, "Sfig5b_spatial_uniformity_heatmap.pdf"),
       plot = sfig_5b_spatial, width = 12, height = 10, dpi = 300)

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

plot_clr_heatmap(
  primary.sct,
  atdc5.sct,
  genes_of_interest = GENES_ALL,
  pdf_file          = file.path(OUTPUT_SUP, "Sfig7_clr_heatmap_primary_atdc5.pdf"),
  pdf_width         = 8,
  pdf_height        = 8
)

save_session_info("figure2_supplement")

