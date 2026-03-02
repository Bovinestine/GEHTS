# Project name: GEHTS-chip
# Author: Nathan Wooseok Lee
# conda env: seurat4
# date: 260302


library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(igraph)
library(ggraph)
library(cowplot)
library(svglite)


mac.sct$drug_name <- mac.sct$drug_condition
merged.sct <- merge(mac.sct, sin.sct)

# Figure 4b 
fig_4b_dumbbel <- plot_automated_coexpression(merged.sct)
svglite("./figures/figure4/fig4_dumbbel_Coexpression_analysis_rapamycin_10.svg", width = 9, height = 5.5)
print(fig_dumbbel)
dev.off()
# rapamycin differential co-expression analysis
# for main figure, we draw a straightforward Dumbbell plot.

#' Automated Directional Co-expression Dumbbell Plot
#'
#' Automatically computes the top 10 most divergent Pearson correlations between 
#' two single-cell populations and plots their trajectory.
#'
#' @param seurat_obj A Seurat object containing single-cell expression data.
#' @param group_by Character. Metadata column containing condition labels. Default is "drug_condition".
#' @param cond_ref Character. Exact metadata string for the reference condition.
#' @param cond_treat Character. Exact metadata string for the treatment condition.
#' @param label_ref Character. Legend display name. Defaults to cond_ref.
#' @param label_treat Character. Legend display name. Defaults to cond_treat.
#' @param top_n Integer. Number of top divergent pairs to extract automatically. Default is 10.
#' @param assay Character. The Seurat assay to use. Default is "SCT".
#' @param slot Character. The data slot to use for correlation math. Default is "data".
#' @param color_ref Character. Hex code for reference dot. Default is NPG Red.
#' @param color_treat Character. Hex code for treatment dot. Default is NPG Blue.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @import Seurat
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @importFrom stats cor

plot_automated_coexpression <- function(seurat_obj, 
                                        group_by = "drug_condition",
                                        cond_ref = "inflammatory", 
                                        cond_treat = "rapamycin_10",
                                        label_ref = cond_ref,       # Auto-inherits parameter input
                                        label_treat = cond_treat,   # Auto-inherits parameter input
                                        top_n = 10,
                                        assay = "SCT", 
                                        slot = "data",
                                        color_ref = "#E64B35FF",
                                        color_treat = "#4DBBD5FF") {
  
  # ---------------------------------------------------------------------------
  # 1. INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!inherits(seurat_obj, "Seurat")) stop("Error: 'seurat_obj' must be a Seurat object.")
  if (!group_by %in% colnames(seurat_obj@meta.data)) stop(paste("Error: Column", group_by, "not found."))
  
  meta_data <- seurat_obj@meta.data
  if (!cond_ref %in% meta_data[[group_by]] || !cond_treat %in% meta_data[[group_by]]) {
    stop("Error: Specified conditions not found within the group_by metadata column.")
  }
  
  # ---------------------------------------------------------------------------
  # 2. DATA EXTRACTION
  # ---------------------------------------------------------------------------
  DefaultAssay(seurat_obj) <- assay
  matrix_data <- GetAssayData(seurat_obj, slot = slot)
  
  cells_ref <- rownames(meta_data[meta_data[[group_by]] == cond_ref, ])
  cells_treat <- rownames(meta_data[meta_data[[group_by]] == cond_treat, ])
  
  # Convert to standard dense matrices for cor()
  mat_ref <- as.matrix(matrix_data[, cells_ref])
  mat_treat <- as.matrix(matrix_data[, cells_treat])
  
  # Filter out genes with zero variance in either condition (prevents NA in Pearson)
  var_ref <- apply(mat_ref, 1, var)
  var_treat <- apply(mat_treat, 1, var)
  valid_genes <- names(which(var_ref > 0 & var_treat > 0))
  
  mat_ref <- mat_ref[valid_genes, ]
  mat_treat <- mat_treat[valid_genes, ]
  
  # ---------------------------------------------------------------------------
  # 3. VECTORIZED CORRELATION MATH (Highly Scalable)
  # ---------------------------------------------------------------------------
  # cor() expects cells as rows, genes as columns, so we transpose (t)
  cor_ref_mat <- cor(t(mat_ref), method = "pearson")
  cor_treat_mat <- cor(t(mat_treat), method = "pearson")
  
  # Extract the upper triangle to avoid duplicates (Gene A vs B) and self (Gene A vs A)
  upper_idx <- which(upper.tri(cor_ref_mat), arr.ind = TRUE)
  
  df_all <- data.frame(
    Gene1 = rownames(cor_ref_mat)[upper_idx[, 1]],
    Gene2 = colnames(cor_ref_mat)[upper_idx[, 2]],
    R_Ref = cor_ref_mat[upper_idx],
    R_Treat = cor_treat_mat[upper_idx],
    stringsAsFactors = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # 4. AUTOMATED TOP 'N' DISCOVERY
  # ---------------------------------------------------------------------------
  df_plot <- df_all %>%
    drop_na() %>%
    mutate(
      GenePair = paste0(Gene1, " - ", Gene2),
      Delta_r = abs(R_Treat - R_Ref),
      Category = paste("Top", top_n, "Divergent Interactions") # Unified category
    ) %>%
    arrange(desc(Delta_r)) %>%
    slice_head(n = top_n) %>%
    arrange(Delta_r) # Re-sort ascending so the largest shifts plot at the top visually
  
  # Lock factor levels
  df_plot$GenePair <- factor(df_plot$GenePair, levels = df_plot$GenePair)
  df_plot$Category <- factor(df_plot$Category, levels = unique(df_plot$Category))
  
  # Dynamic Axis Scaling
  max_val <- max(abs(c(df_plot$R_Ref, df_plot$R_Treat)), na.rm = TRUE)
  axis_limit <- ceiling(max_val * 10) / 10 + 0.1 
  
  # ---------------------------------------------------------------------------
  # 5. VISUALIZATION
  # ---------------------------------------------------------------------------
  color_mapping <- setNames(c(color_ref, color_treat), c(label_ref, label_treat))
  
  p_dumbbell <- ggplot(df_plot) +
    
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    
    geom_segment(aes(x = R_Ref, xend = R_Treat, y = GenePair, yend = GenePair),
                 color = "grey40", linewidth = 1.2, 
                 arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
    
    geom_point(aes(x = R_Ref, y = GenePair, fill = !!label_ref), 
               size = 4, shape = 21, color = "black", stroke = 1) +
    
    geom_point(aes(x = R_Treat, y = GenePair, fill = !!label_treat), 
               size = 4, shape = 21, color = "black", stroke = 1) +
    
    scale_fill_manual(values = color_mapping) +
    # facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_continuous(limits = c(-axis_limit, axis_limit), 
                       breaks = round(seq(-axis_limit, axis_limit, length.out = 5), 2)) +
    
    theme_bw() +
    labs(
      title = "Data-Driven Differential Co-expression",
      subtitle = paste0("Unbiased extraction of highest magnitude transcriptomic shifts (", label_ref, " vs. ", label_treat, ")"),
      x = "Pearson Correlation Coefficient (r)",
      y = ""
    ) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 11, face = "bold"),
      
      axis.text.y = element_text(size = 12, face = "bold.italic", color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
      
      strip.text.y.left = element_text(angle = 0, size = 11, face = "bold"),
      strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1),
      strip.placement = "outside",
      
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 1)
    )
  
  return(p_dumbbell)
}



###
# radar chart for single drugs
# Load necessary libraries
library(Seurat)
library(fmsb)
library(dplyr)
library(scales)

# ==============================================================================
# 1. SCALE CALCULATION FUNCTION
# ==============================================================================

#' Calculate Universal Radar Scale
#'
#' Scans all conditions within a Seurat object to calculate the absolute 
#' minimum and maximum expression values for predefined biological modules. 
#' Ensures all subsequent radar charts share an identical, comparable axis scale.
#'
#' @param seurat_obj A Seurat object.
#' @param group_by Character. Metadata column containing conditions. Default is "drug_condition".
#' @param assay Character. Assays to pull data from. Default is "SCT".
#' @param slot Character. Slot to pull data from. Default is "counts".
#' @param module_list Named list of character vectors containing gene symbols for each biological module.
#'
#' @return A list containing 'min' and 'max' numeric vectors for the modules.
#' @export
calculate_radar_scale <- function(seurat_obj, 
                                  group_by = "drug_condition",
                                  assay = "SCT", 
                                  slot = "counts",
                                  module_list = list(
                                    Chondrocyte_Development  = c("Matn1", "Acan", "Sox9", "Col27a1"),
                                    Cartilage_Condensation   = c("Col2a1", "Acan", "Sox9"),
                                    NO_Synthase_Upregulation = c("Tlr2", "Ccl2"),
                                    ECM_Disassembly          = c("Adamts5", "Mmp13"),
                                    Inflammatory_response    = c("Tlr2", "Ccl2", "Il17b", "Tnfrsf1b", "Cxcl1", "Il6", "Cxcl5", "Fosl2"),
                                    Tissue_homeostasis       = c("Sox9", "Col2a1", "Pth1r", "Fosl2")
                                  )) {
  
  if (!inherits(seurat_obj, "Seurat")) stop("Error: 'seurat_obj' must be a Seurat object.")
  if (!group_by %in% colnames(seurat_obj@meta.data)) stop(paste("Error: Column", group_by, "not found in metadata."))
  
  cat("Calculating universal scale across all conditions... This may take a moment.\n")
  DefaultAssay(seurat_obj) <- assay
  mat <- GetAssayData(seurat_obj, slot = slot)
  meta <- seurat_obj@meta.data
  all_conditions <- unique(meta[[group_by]])
  
  all_metrics <- list()
  for (cond in all_conditions) {
    cells <- rownames(meta[meta[[group_by]] == cond, ])
    if (length(cells) == 0) next
    
    metrics <- sapply(module_list, function(genes) {
      avail_genes <- intersect(genes, rownames(mat))
      if (length(avail_genes) > 0) mean(rowMeans(mat[avail_genes, cells, drop = FALSE])) else NA
    })
    all_metrics[[cond]] <- metrics
  }
  
  metrics_df <- do.call(rbind, all_metrics)
  
  # Calculate Universal Min and Max
  min_vals <- apply(metrics_df, 2, min, na.rm = TRUE)
  max_vals <- apply(metrics_df, 2, max, na.rm = TRUE)
  
  cat("Scale calculation complete!\n")
  return(list(min = min_vals, max = max_vals))
}

# ==============================================================================
# 2. PLOTTING FUNCTION
# ==============================================================================

#' Plot Phenotypic Shift Radar Chart
#'
#' Generates a base-R radar chart comparing a target drug condition against 
#' healthy and diseased baselines, mapped to a pre-calculated universal scale.
#'
#' @param seurat_obj A Seurat object.
#' @param target_condition Character. Exact metadata string for the drug to plot.
#' @param radar_scale List. The output from `calculate_radar_scale()`.
#' @param group_by Character. Metadata column containing conditions. Default is "drug_condition".
#' @param ref_control Character. Metadata string for the healthy baseline. Default is "control".
#' @param ref_disease Character. Metadata string for the disease baseline. Default is "inflammatory".
#' @param target_label Character. Clean display label for the drug in the legend.
#' @param module_list Named list of gene modules. Must match the one used for scaling.
#' @param assay Character. Default "SCT".
#' @param slot Character. Default "counts".
#'
#' @return NULL (Draws to the active base R graphics device).
#' @export
plot_radar_fast <- function(seurat_obj, 
                            target_condition, 
                            radar_scale, 
                            group_by = "drug_condition",
                            ref_control = "control", 
                            ref_disease = "inflammatory",
                            target_label = NULL,
                            assay = "SCT", 
                            slot = "counts",
                            module_list = list(
                              Chondrocyte_Development  = c("Matn1", "Acan", "Sox9", "Col27a1"),
                              Cartilage_Condensation   = c("Col2a1", "Acan", "Sox9"),
                              NO_Synthase_Upregulation = c("Tlr2", "Ccl2"),
                              ECM_Disassembly          = c("Adamts5", "Mmp13"),
                              Inflammatory_response    = c("Tlr2", "Ccl2", "Il17b", "Tnfrsf1b", "Cxcl1", "Il6", "Cxcl5", "Fosl2"),
                              Tissue_homeostasis       = c("Sox9", "Col2a1", "Pth1r", "Fosl2")
                            )) {
  
  # 1. Validation
  if (!target_condition %in% seurat_obj@meta.data[[group_by]]) {
    stop(paste("Error: target_condition", target_condition, "not found in metadata."))
  }
  if (is.null(target_label)) target_label <- target_condition
  
  DefaultAssay(seurat_obj) <- assay
  mat <- GetAssayData(seurat_obj, slot = slot)
  meta <- seurat_obj@meta.data
  
  # Internal helper to calculate metrics
  get_metrics <- function(cond) {
    cells <- rownames(meta[meta[[group_by]] == cond, ])
    if (length(cells) == 0) stop(paste("No cells found for condition:", cond))
    
    sapply(module_list, function(genes) {
      avail_genes <- intersect(genes, rownames(mat))
      if (length(avail_genes) > 0) mean(rowMeans(mat[avail_genes, cells, drop = FALSE])) else NA
    })
  }
  
  # 2. Calculate metrics ONLY for the 3 needed lines
  conds_to_plot <- c(ref_control, ref_disease, target_condition)
  metrics_mat <- t(sapply(conds_to_plot, get_metrics))
  rownames(metrics_mat) <- conds_to_plot
  
  # 3. Normalize using the Pre-Calculated Universal Scale
  min_vals <- radar_scale$min
  max_vals <- radar_scale$max
  
  diff_vals <- max_vals - min_vals
  diff_vals[diff_vals == 0] <- 1 # Prevent division by zero
  
  norm_df <- sweep(metrics_mat, 2, min_vals, FUN = "-")
  norm_df <- sweep(norm_df, 2, diff_vals, FUN = "/")
  
  # 4. Prepare dataframe for fmsb::radarchart
  plot_df <- rbind(
    rep(1, ncol(norm_df)), # Row 1 must be Max limits
    rep(0, ncol(norm_df)), # Row 2 must be Min limits
    norm_df                # Rows 3-5 are the actual data
  )
  plot_df <- as.data.frame(plot_df)
  colnames(plot_df) <- gsub("_", "\n", colnames(plot_df))
  
  # 5. Generate the Plot
  radarchart(plot_df, axistype = 1,
             pcol = c("#3C5488FF", "#E64B35FF", "#00A087FF"), # NPG Blue, Red, Teal
             pfcol = c(alpha("#3C5488FF", 0.1), alpha("#E64B35FF", 0.1), alpha("#00A087FF", 0.5)),
             plwd = c(2, 2, 3), plty = c(2, 2, 1),
             cglcol = "grey70", cglty = 1, axislabcol = "grey50", caxislabels = seq(0, 1, 0.25),
             cglwd = 0.8, vlcex = 0.9,
             title = paste("Phenotypic Shift:", target_label))
  
  # 6. Add Legend
  legend("topright", legend = c("Control", "Inflammatory", target_label),
         col = c("#3C5488FF", "#E64B35FF", "#00A087FF"),
         lty = c(2, 2, 1), lwd = c(2, 2, 3), bty = "n", pch = 20, pt.cex = 2,
         inset = c(-0.1, 0))
}

# ==============================================================================
# EXPORT WORKFLOW
# ==============================================================================

# 1. Run the heavy calculation ONCE 
universal_scale <- calculate_radar_scale(merged.sct)

# 2. Setup external device layout for side-by-side exporting
library(svglite)
svglite("./figures/figure4/Radar_Comparison_Rank1_vs_Rank13.svg", width = 12, height = 5.5)

# Save default par settings so we can restore them later
old_par <- par(no.readonly = TRUE)

# Set grid to 1 row, 2 columns, with right margin expanded for the legend
par(mfrow = c(1, 2), mar = c(2, 2, 3, 5), xpd = TRUE)

# 3. Plot Rank 1 (Rapamycin)
plot_radar_fast(merged.sct, 
                target_condition = "rapamycin_10", 
                target_label = "Rapamycin (10 μM)",
                radar_scale = universal_scale)

# 4. Plot Rank 13 (ALK5 inhibitor)
plot_radar_fast(merged.sct, 
                target_condition = "ALK5 inhibitor IV_10", 
                target_label = "ALK5 Inhibitor (10 μM)",
                radar_scale = universal_scale)

# 5. Safely close device and restore graphical parameters
dev.off()
par(old_par)


# Figure 4d _ P38 pathway inhibitors


plot_decoupling_boxplots <- function(seurat_obj, 
                                     group_by = "drug_condition",
                                     cond_ctrl = "control",
                                     cond_inf = "inflammatory",
                                     cond_A = "pamapimod_10",
                                     cond_B = "SB203580_10",
                                     label_ctrl = "Control",
                                     label_inf = "IL-1β",
                                     label_A = "Pamapimod",
                                     label_B = "SB203580",
                                     assay = "SCT", 
                                     slot = "data") {
  
  # 1. Extract Data
  DefaultAssay(seurat_obj) <- assay
  mat <- GetAssayData(seurat_obj, slot = slot)[c("Sox9", "Col2a1"), , drop = FALSE]
  
  df <- data.frame(
    Cell = colnames(seurat_obj),
    Sox9 = as.numeric(mat["Sox9", ]),
    Col2a1 = as.numeric(mat["Col2a1", ]),
    Condition_Raw = seurat_obj@meta.data[[group_by]],
    stringsAsFactors = FALSE
  )
  
  # 2. Filter and Map Labels
  target_conditions <- c(cond_ctrl, cond_inf, cond_A, cond_B)
  df <- df %>% filter(Condition_Raw %in% target_conditions)
  
  label_mapping <- setNames(c(label_ctrl, label_inf, label_A, label_B), target_conditions)
  df$Condition <- factor(label_mapping[df$Condition_Raw], 
                         levels = c(label_ctrl, label_inf, label_A, label_B))
  
  # 3. Reshape for ggplot
  df_melt <- melt(df, id.vars = c("Cell", "Condition_Raw", "Condition"), 
                  variable.name = "Gene", value.name = "Expression")
  
  # 4. Colors
  npg_colors <- setNames(
    c("#4DBBD5FF", "#E64B35FF", "#8491B4FF", "#00A087FF"), 
    c(label_ctrl, label_inf, label_A, label_B)
  )
  
  # 5. Helper Function for Individual Box Plots
  create_box <- function(gene_name, plot_title) {
    ggplot(df_melt %>% filter(Gene == gene_name), aes(x = Condition, y = Expression, fill = Condition)) +
      # Adding jittered points behind the boxplot handles the 50-200 cell count beautifully
      geom_jitter(color = "grey60", size = 1, alpha = 0.5, width = 0.2) +
      geom_boxplot(outlier.shape = NA, alpha = 0.9, color = "black", linewidth = 0.6, width = 0.6) +
      scale_fill_manual(values = npg_colors) +
      theme_classic() +
      labs(title = plot_title, y = "Expression (Log-Normalized)", x = "") +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold.italic", size = 14, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title.y = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      )
  }
  
  # 6. Generate Panels
  p_sox9 <- create_box("Sox9", "Sox9 (Master Regulator)")
  p_col2 <- create_box("Col2a1", "Col2a1 (Functional Matrix)")
  
  # 7. Combine
  final_plot <- p_sox9 + p_col2 + 
    plot_annotation(
      title = "Resolution of Transcriptional Decoupling",
      theme = theme(plot.title = element_text(size = 15, face = "bold"))
    )
  
  return(final_plot)
}

# --- Exceution Exampltion ---
fig_4d_box <- plot_decoupling_boxplots(merged.sct)
print(fig_box)


