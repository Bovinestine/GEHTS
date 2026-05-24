# Project name: GEHTS-chip
# Author: Nathan Wooseok Lee
# conda env: seurat4
# date: 260524


####################
library(ggplot2)
library(dplyr)
library(tidyr)

#' Generate Supplementary 6-Panel Violin Plots (Multi-Target)
#'
#' @param seurat_obj A Seurat object.
#' @param target_conditions Character vector. Exact metadata strings for the drugs to plot.
#' @param group_by Character. Metadata column containing conditions.
#' @param ref_control Character. Metadata string for the healthy baseline.
#' @param ref_disease Character. Metadata string for the disease baseline.
#' @param target_labels Character vector. Clean display labels for the x-axis.
#' @param module_list Named list of gene modules. Must match the radar chart.
#' @param assay Character. Default "SCT".
#' @param slot Character. Default "data".
#'
#' @return A ggplot object containing the faceted violin plots.
#' @export
plot_supplementary_violins_multi <- function(seurat_obj, 
                                             target_conditions, 
                                             group_by = "drug_condition",
                                             ref_control = "control", 
                                             ref_disease = "inflammatory",
                                             target_labels = NULL,
                                             assay = "SCT", 
                                             slot = "data",
                                             module_list = list(
                                               Chondrocyte_Development  = c("Matn1", "Acan", "Sox9", "Col27a1"),
                                               Cartilage_Condensation   = c("Col2a1", "Acan", "Sox9"),
                                               NO_Synthase_Upregulation = c("Tlr2", "Ccl2"),
                                               ECM_Disassembly          = c("Adamts5", "Mmp13"),
                                               Inflammatory_response    = c("Tlr2", "Ccl2", "Il17b", "Tnfrsf1b", "Cxcl1", "Il6", "Cxcl5", "Fosl2"),
                                               Tissue_homeostasis       = c("Sox9", "Col2a1", "Pth1r", "Fosl2")
                                             )) {
  
  if (is.null(target_labels)) target_labels <- target_conditions
  
  # 1. Extract Data
  DefaultAssay(seurat_obj) <- assay
  mat <- GetAssayData(seurat_obj, slot = slot)
  meta <- seurat_obj@meta.data
  
  # Subset to only the requested conditions
  conds_to_plot <- c(ref_control, ref_disease, target_conditions)
  cells_to_keep <- rownames(meta[meta[[group_by]] %in% conds_to_plot, ])
  
  mat_sub <- mat[, cells_to_keep]
  meta_sub <- meta[cells_to_keep, ]
  
  # 2. Calculate Per-Cell Module Scores
  score_list <- lapply(names(module_list), function(mod_name) {
    genes <- intersect(module_list[[mod_name]], rownames(mat_sub))
    
    if (length(genes) > 0) {
      scores <- colMeans(as.matrix(mat_sub[genes, , drop = FALSE]), na.rm = TRUE)
    } else {
      scores <- rep(NA, length(cells_to_keep))
    }
    
    data.frame(
      Cell = cells_to_keep,
      Condition = meta_sub[[group_by]],
      Module = gsub("_", " ", mod_name), 
      Score = scores
    )
  })
  
  # 3. Combine and Format Data
  df_long <- do.call(rbind, score_list)
  
  # Create a named mapping vector to elegantly replace original conditions with clean labels
  plot_labels <- c("Basal", "IL-1β", target_labels)
  cond_mapping <- setNames(plot_labels, conds_to_plot)
  
  df_long <- df_long %>%
    mutate(
      Condition_Clean = cond_mapping[as.character(Condition)],
      # Lock factor levels to maintain the left-to-right order (Basal -> IL-1B -> Drug 1 -> Drug 2)
      Condition_Clean = factor(Condition_Clean, levels = plot_labels)
    ) %>%
    filter(!is.na(Score)) # Drop missing values
  
  # 4. Define Dynamic Colors (Supporting up to 6 total conditions)
  base_colors <- c("#3C5488FF", "#E64B35FF", "#00A087FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF")
  condition_colors <- setNames(base_colors[1:length(plot_labels)], plot_labels)
  
  # 5. Generate Faceted Violin Plot
  p_violins <- ggplot(df_long, aes(x = Condition_Clean, y = Score, fill = Condition_Clean)) +
    
    geom_violin(trim = TRUE, alpha = 0.8, scale = "width", color = "black", linewidth = 0.5) +
    
    geom_boxplot(width = 0.15, fill = "white", color = "black", 
                 outlier.shape = NA, linewidth = 0.5, alpha = 0.8) +
    
    facet_wrap(~ Module, scales = "free_y", ncol = 3) +
    
    scale_fill_manual(values = condition_colors) +
    
    theme_bw() +
    labs(
      title = "Single-Cell Module Distributions (Benchmarking)",
      x = "", 
      y = "Mean Gene Expression (SCT)"
    ) +
    theme(
      legend.position = "none", 
      # 45-degree angled text prevents the x-axis labels from overlapping now that there are 4+ violins
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(face = "bold", size = 12, margin = ggplot2::margin(r = 10)),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      
      strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1),
      strip.text = element_text(face = "bold", size = 10),
      panel.border = element_rect(color = "black", linewidth = 1),
      panel.grid.major.x = element_blank()
    )
  
  return(p_violins)
}

sup.fig_6moduleViolin <- plot_supplementary_violins_multi(
  seurat_obj = merged.sct, 
  target_conditions = c("rapamycin_10", "ALK5 inhibitor IV_10"), 
  target_labels = c("Rapamycin (10 μM)", "ALK5 Inhibitor (10 μM)")
)
ggsave("./supfig4d_6moduleViolin_rapamycin_10_260503.pdf", plot= sup.fig_6moduleViolin, width = 7, height = 5.5, units = "in", dpi = 300)

