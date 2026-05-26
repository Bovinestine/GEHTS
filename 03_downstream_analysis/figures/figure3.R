# =============================================================================
# figure3.R — Single-drug screen analysis (Figure 3)
#
# Subfigures:
#   3b: Gene expression heatmap sorted by dose (sin.sct)
#   3c: Global UMAP at dose 10 μM colored by drug name
#
# Prerequisites: config.R, utils.R, data_loading.R
# Inputs: sin.sct (loaded by data_loading.R)
# Outputs: output/main/Figure3{b,c}_*.pdf
# Note: sin.sct10 (UMAP-embedded dose-10 subset) is kept in the environment
#       for use by figure3_supplement.R (Extended Data) and figure4.R (Fig. 4d)
# =============================================================================

source(here::here("03_downstream_analysis", "R", "config.R"))
source(here::here("03_downstream_analysis", "R", "utils.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(ggsci)
})

# =============================================================================
# Figure 3b — Heatmap of drug-perturbed single-cell expression
# =============================================================================
fig_3b <- create_heatmap_by_expression_ps(
  sin.sct,
  primary_metadata   = "dose",
  secondary_metadata = "drug_name"
)
pdf(file.path(OUTPUT_MAIN, "Figure3b_heatmap_SCTransform.pdf"), width = 6, height = 5)
print(fig_3b)
dev.off()

# =============================================================================
# Figure 3c — Global UMAP at dose 10 μM
# =============================================================================

#' Generate Global UMAP for a Specific Dose
#'
#' Subsets the Seurat object to the target dose, runs PCA + UMAP with a fixed
#' seed, and returns a DimPlot with NPG colors per drug.
#'
#' @param seurat_obj A Seurat object with "dose" and "drug_name" metadata.
#' @param target_dose Numeric dose to subset (default 10).
#' @param pcs_to_use Principal components to use for UMAP (default 1:10).
#' @param seed RNG seed (default: SEED from config).
#' @return A ggplot2 object.
#' @export
plot_main_umap <- function(seurat_obj,
                           target_dose = 10,
                           pcs_to_use  = 1:10,
                           seed        = SEED) {

  sub_obj <- subset(seurat_obj, subset = dose == target_dose)
  set.seed(seed)
  sub_obj <- RunPCA(sub_obj, assay = "SCT",
                    npcs = max(pcs_to_use), verbose = FALSE)
  sub_obj <- RunUMAP(sub_obj, reduction = "pca", assay = "SCT",
                     dims = pcs_to_use, verbose = FALSE)

  drug_names <- sort(unique(as.character(sub_obj$drug_name)))
  npg_colors <- colorRampPalette(PALETTE_NPG)(length(drug_names))
  names(npg_colors) <- drug_names

  DimPlot(sub_obj, reduction = "umap", group.by = "drug_name") +
    scale_color_manual(values = npg_colors) +
    labs(title = paste0("Global UMAP (Dose = ", target_dose, " μM)")) +
    theme_classic()
}

fig_3c <- plot_main_umap(sin.sct, target_dose = 10)

sin.sct10 <- subset(sin.sct, subset = dose == 10)
set.seed(SEED)
sin.sct10 <- RunPCA(sin.sct10, assay = "SCT", npcs = 10, verbose = FALSE)
sin.sct10 <- RunUMAP(sin.sct10, reduction = "pca", assay = "SCT",
                     dims = 1:10, verbose = FALSE)

ggsave(file.path(OUTPUT_MAIN, "Figure3c_UMAP_dose10.pdf"),
       plot = fig_3c, width = 6.5, height = 5, dpi = 300)

save_session_info("figure3")
