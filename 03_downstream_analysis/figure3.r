# GE-HTS project pipeline for in situ sequencing data analysis
# author: Nathan Wooseok Lee
# Date of start: 241224
# Date of update: 250711, 250827 (wide gene expression), 251223 (gene group), 260524 (revise)
# env: conda activate seurat4 / issr

### Figure 3b (GEheatmap)
# ==============================================================================
# Script Name: figure3b_heatmap.R
# Description: Generates a gene expression heatmap (Figure 3b) categorizing 
#              anabolic, inflammatory (catabolic), and housekeeping genes.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(pheatmap)
  library(ggsci) # Added for NPG and JCO color palettes
})

# ------------------------------------------------------------------------------
# 1. Define Gene Categories
# ------------------------------------------------------------------------------
anabolic <- c('Acan','Sox9','Col2a1','Matn1','Matn3','Ucma','Ccnd3','Gadd45g','Pth1r','Gm26633','Col27a1')
inflammatory <- c('Mmp3','Mmp13','Il6', 'Il17b','Adamts5','Igfbp3','Ccl2','Cxcl5','Cxcl1','Fosl2','Tlr2','Tnfrsf1b')
housekeeping <- c('Hprt','Actb','Gapdh','B2m','Ubc','Ppia','Rpl23')
# ------------------------------------------------------------------------------
# 2. Main Heatmap Function
# ------------------------------------------------------------------------------

#' Create Gene Expression Heatmap Categorized by Phenotype
#'
#' Extracts assay data from a Seurat object, filters genes by expression quantiles,
#' groups cells by specified metadata, and generates an annotated pheatmap using 
#' publication-ready ggsci color palettes.
#' @param seurat_object A processed Seurat object.
#' @param upper_threshold Upper quantile limit for filtering genes (default: 0).
#' @param bottom_threshold Lower quantile limit for filtering genes (default: 1).
#' @param assay Name of the assay to pull data from (default: 'SCT').
#' @param primary_metadata Column name in meta.data for primary grouping (e.g., 'dose').
#' @param secondary_metadata Column name in meta.data for secondary grouping (e.g., 'drug_name').
#'
#' @return The heatmap plot.
#' @export
create_heatmap_by_expression_ps <- function(seurat_object, 
                                            upper_threshold = 0, 
                                            bottom_threshold = 1, 
                                            assay = 'SCT', 
                                            primary_metadata = 'dose', 
                                            secondary_metadata = 'drug_name') {
  
  # --- Step 1: Extract COUNTS for Sorting & Filtering ---
  # We MUST use counts here because Pearson residuals (scale.data) have a mean 
  # of ~0, which would break the absolute expression quantile filtering.
  count_data <- GetAssayData(seurat_object, assay = assay, slot = 'counts')
  
  count_data <- count_data[, colSums(is.finite(count_data)) > 0]
  average_expression <- rowMeans(count_data, na.rm = TRUE)
  
  q_bottom <- quantile(average_expression, 1 - bottom_threshold)
  q_top    <- quantile(average_expression, 1 - upper_threshold)
  selected_genes <- names(average_expression[average_expression >= q_bottom & 
                                               average_expression <= q_top])
  
  # --- Step 2: Categorize and Sort Genes ---
  selected_anabolic <- intersect(selected_genes, anabolic)
  selected_inflam   <- intersect(selected_genes, inflammatory)
  selected_housek   <- intersect(selected_genes, housekeeping)
  
  sorted_selected_genes <- c(
    selected_anabolic[order(average_expression[selected_anabolic], decreasing = TRUE)],
    selected_inflam[order(average_expression[selected_inflam], decreasing = TRUE)],
    selected_housek[order(average_expression[selected_housek], decreasing = TRUE)]
  )
  
  # --- Step 3: Prepare Cell Metadata and Grouping ---
  cell_scores <- data.frame(
    cell_name = rownames(seurat_object@meta.data),
    primary   = seurat_object@meta.data[[primary_metadata]],
    secondary = seurat_object@meta.data[[secondary_metadata]],
    stringsAsFactors = FALSE
  )
  
  cell_scores <- cell_scores %>% arrange(primary, secondary)
  
  # --- NEW LOGIC: Extract SCTransform Residuals for Plotting ---
  scale_data <- GetAssayData(seurat_object, assay = assay, slot = 'scale.data')
  
  # Safety check: SCTransform by default only stores residuals for highly variable genes.
  # We must intersect to drop any genes that lack computed residuals.
  available_genes <- intersect(sorted_selected_genes, rownames(scale_data))
  
  if (length(available_genes) < length(sorted_selected_genes)) {
    missing_genes <- setdiff(sorted_selected_genes, available_genes)
    warning(paste("The following genes are missing from scale.data and will be skipped:", 
                  paste(missing_genes, collapse = ", ")))
  }
  
  data_for_heatmap <- as.matrix(scale_data[available_genes, cell_scores$cell_name])
  
  # We NO LONGER run t(scale(t())) or log1p() because Pearson residuals 
  # are already centered, regularized, and normalized!
  
  # Handle potential NaNs 
  data_for_heatmap[is.na(data_for_heatmap)] <- 0
  
  # Cap extreme residuals to prevent outliers from washing out the color scale
  data_for_heatmap[data_for_heatmap > 3] <- 3
  data_for_heatmap[data_for_heatmap < -3] <- -3
  
  # --- Step 4: Configure Visuals & Annotations ---
  
  # Create perfectly symmetric breaks around 0 for the diverging palette
  max_val <- max(abs(data_for_heatmap), na.rm = TRUE)
  breaks <- seq(-max_val, max_val, length.out = 101)
  color_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
  
  # Column (Cell) Annotations
  annotation_col <- data.frame(
    Primary   = factor(cell_scores$primary, levels = unique(cell_scores$primary)),
    Secondary = factor(cell_scores$secondary, levels = unique(cell_scores$secondary))
  )
  rownames(annotation_col) <- cell_scores$cell_name
  
  # Re-calculate lengths for annotations in case some genes were skipped
  final_anabolic <- intersect(available_genes, selected_anabolic)
  final_inflam   <- intersect(available_genes, selected_inflam)
  final_housek   <- intersect(available_genes, selected_housek)
  
  # Row (Gene) Annotations
  annotation_row <- data.frame(
    Gene_Type = factor(rep(c("Anabolic", "Catabolic", "Housekeeping"), 
                           c(length(final_anabolic), length(final_inflam), length(final_housek))),
                       levels = c("Anabolic", "Catabolic", "Housekeeping"))
  )
  rownames(annotation_row) <- available_genes
  
  # ggsci Color Palettes for Annotations
  n_primary <- length(unique(cell_scores$primary))
  primary_colors <- colorRampPalette(pal_jco()(10))(n_primary)
  names(primary_colors) <- unique(cell_scores$primary)
  
  n_secondary <- length(unique(cell_scores$secondary))
  secondary_colors <- colorRampPalette(pal_npg("nrc")(10))(n_secondary)
  names(secondary_colors) <- unique(cell_scores$secondary)
  
  ann_colors <- list(
    Gene_Type = c(Anabolic = "forestgreen", Catabolic = "orange", Housekeeping = "grey"),
    Primary   = primary_colors,
    Secondary = secondary_colors
  )
  
  gap_row_positions <- c(length(final_anabolic), 
                         length(final_anabolic) + length(final_inflam))
  
  # --- Step 5: Generate Plot ---
  
  p <- pheatmap::pheatmap(
    mat               = data_for_heatmap,
    cluster_rows      = FALSE, 
    cluster_cols      = FALSE,
    show_rownames     = TRUE,
    show_colnames     = FALSE,
    annotation_col    = annotation_col,
    annotation_row    = annotation_row,
    annotation_colors = ann_colors,
    color             = color_palette,
    breaks            = breaks,
    main              = "SCTransform Residuals",
    gaps_row          = gap_row_positions
  )
  
  return(p)
}

# Exceution
fig_3b <- create_heatmap_by_expression_ps(sin.sct) 
ggsave("./Figure3b_heatmap_zscore.pdf", plot = fig_3b, width =6, height = 5, units = "in", dpi = 300)

# ==============================================================================
# Script Name: figure3c_main_umap.R
# Description: Generates and saves the main global UMAP for a specific drug dose.
#              Uses an interpolated NPG palette for consistency.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggsci)
})

#' Generate and Save Main UMAP by Dose with NPG Colors
#'
#' @param seurat_obj A Seurat object containing 'dose' and 'drug_name' in meta.data.
#' @param target_dose Numeric value of the dose to subset (default is 10).
#' @param pcs_to_use Numeric vector specifying which principal components to use.
#' @param seed Numeric value of the seed to set.seed() (defualt is 123).
#'
#' @return The processed Seurat object (invisibly) for downstream use.
#' @export
plot_main_umap <- function(seurat_obj, 
                           target_dose = 10, 
                           pcs_to_use = 1:10,
                           seed = 123
                           ) {
  
  message("Subsetting data for dose: ", target_dose)
  seurat_sub <- subset(seurat_obj, subset = dose == target_dose)
  
  message("Running PCA and UMAP...")
  set.seed(seed)
  seurat_sub <- RunPCA(seurat_sub, assay = "SCT", npcs = max(pcs_to_use), verbose = FALSE)
  seurat_sub <- RunUMAP(seurat_sub, reduction = "pca", assay = "SCT", dims = pcs_to_use, verbose = FALSE)
  
  # --- Generate Consistent NPG Color Dictionary ---
  drug_names <- sort(unique(as.character(seurat_sub$drug_name)))
  npg_colors <- colorRampPalette(pal_npg("nrc")(10))(length(drug_names))
  names(npg_colors) <- drug_names
  
  # --- Create Plot ---
  p <- DimPlot(seurat_sub, reduction = "umap", group.by = 'drug_name') +
    scale_color_manual(values = npg_colors) + # Apply the NPG dictionary
    labs(title = paste("Global UMAP (Dose =", target_dose, "\u03BCM)")) +
    theme_classic()
    
  return(p)
}

fig_3c_umap <- plot_main_umap(sin.sct)
ggsave('D:/Research_Repository/in_situ_team/Analysis/figures/figure3/fig3c_UMAP_newColor.pdf', plot = fig_3c_umap, width = 6.5, height = 5, units = "in", dpi = 300)

# ==============================================================================
# Script Name: figure3d_highlight_umaps.R
# Description: Generates individual UMAPs highlighting specific drug conditions (figure 3d).
# ==============================================================================

library(Seurat)
library(ggplot2)
library(gridExtra)

# --- Internal Helper Functions ---

#' Extract UMAP Data for ggplot
#' @keywords internal
.extract_umap_data <- function(seurat_obj) {
  # Use Seurat's Embeddings function for safer extraction across versions
  umap_coords <- Embeddings(seurat_obj, reduction = "umap")
  meta_data <- seurat_obj@meta.data[, c("drug_name", "dose", "target")]
  
  # Combine and ensure standard column names
  df <- data.frame(umap_coords, meta_data)
  colnames(df)[1:2] <- c("UMAP_1", "UMAP_2") 
  return(df)
}

#' Create List of Highlighted UMAPs
#' @keywords internal
.make_drug_umap_list <- function(umap_data) {
  drugs <- unique(umap_data$drug_name)
  
  plot_list <- lapply(drugs, function(condition) {
    ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = (drug_name == condition))) +
      geom_point(aes(alpha = (drug_name == condition)), size = 1) +
      scale_alpha_manual(values = c("FALSE" = 0.1, "TRUE" = 0.6)) +
      scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "turquoise3")) +
      labs(title = paste("Highlighted:", condition)) +
      theme_classic() +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", size = 10))
  })
  
  # Name the list elements for easy targeted subsetting later
  names(plot_list) <- drugs
  return(plot_list)
}

# --- Main Exported Function ---

#' Save Grid of Highlighted UMAPs
#'
#' Extracts UMAP coordinates from a processed Seurat object, generates highlighted 
#' plots for each drug, and saves a specified subset of them in a grid layout.
#'
#' @param seurat_obj A Seurat object that has already been run through RunUMAP().
#' @param plot_indices Numeric vector or character vector of drugs to plot 
#'        (e.g., 1:6, or c("Pamapimod", "XAV")).
#' @param output_pdf Character string specifying the output PDF path.
#' @param ncol Integer specifying the number of columns in the grid layout.
#' @param pdf_width Numeric width of the output PDF.
#' @param pdf_height Numeric height of the output PDF.
#'
#' @export
plot_highlighted_umaps <- function(seurat_obj, 
                                   plot_indices = NULL, 
                                   output_pdf, 
                                   ncol = 2,
                                   pdf_width = 6, 
                                   pdf_height = 9) {
  
  # 1. Extract data and generate all plots
  umap_df <- .extract_umap_data(seurat_obj)
  all_plots <- .make_drug_umap_list(umap_df)
  
  # 2. Subset the requested plots
  if (is.null(plot_indices)) {
    plots_to_print <- all_plots # Print all if none specified
  } else {
    plots_to_print <- all_plots[plot_indices]
  }
  
  # 3. Ensure directory exists and save
  output_dir <- dirname(output_pdf)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  pdf(output_pdf, width = pdf_width, height = pdf_height)
  grid.arrange(grobs = plots_to_print, ncol = ncol)
  dev.off()
  
  message("Figure 3D subset saved to: ", output_pdf)
}


source("scripts/figure3d_highlight_umaps.R")

# 2. Generate Figure 3D permutations using the object returned above
# Supplementary 1-6
plot_highlighted_umaps(
  seurat_obj = sin.sct10,
  plot_indices = 1:6,
  output_pdf = "./figures/figure3/umap_single10_1to6.pdf"
)

# Supplementary 7-12
plot_highlighted_umaps(
  seurat_obj = sin.sct10,
  plot_indices = 7:12,
  output_pdf = "./figures/figure3/umap_single10_7to12.pdf"
)

# Specific Highlight (Pamapimod and XAV)
plot_highlighted_umaps(
  seurat_obj = sin.sct10,
  plot_indices = c(1, 13), # Alternatively, c("Pamapimod", "XAV") if names match
  ncol = 1,
  pdf_width = 3,
  pdf_height = 6,
  output_pdf = "./figures/figure3/umap_single10noHarmonyPamapimodXav.pdf"
)

