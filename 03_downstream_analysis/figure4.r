# Project name: GEHTS-chip
# Author: Nathan Wooseok Lee
# conda env: seurat4
# date: 260524

### Figure 4a (Efficacy Score visualization)
# ==============================================================================
# Description: Ver.1 Function to plot and save a bar plot of drug efficacy probabilities.
# ==============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(randomForest)

#' Generate and Save Drug Efficacy Bar Plot
#'
#' @param prob_sng10_raw A list or named vector of ALL raw probabilities at dose 10.
#' @param prob_sng0.1_raw A list or named vector of ALL raw probabilities at dose 0.1.
#' @return The grouped bar plot with error bars.
#' @export
plot_efficacy_barplot <- function(prob_sng10_raw, prob_sng0.1_raw) {
  
  # 1. Convert raw lists to data frames
  df_10 <- stack(prob_sng10_raw)
  colnames(df_10) <- c("Probability", "Drug")
  df_10$Dose <- "10"
  
  df_01 <- stack(prob_sng0.1_raw)
  colnames(df_01) <- c("Probability", "Drug")
  df_01$Dose <- "0.1"
  
  # 2. Combine and format
  raw_data <- rbind(df_10, df_01)
  raw_data$Dose <- factor(raw_data$Dose, levels = c("0.1", "10"))
  
  # 3. Filter out baseline conditions
  exclude_terms <- c('control', 'inflammatory')
  raw_data <- raw_data[!raw_data$Drug %in% exclude_terms, ]
  
  # 4. Calculate Mean and Standard Error (SEM) for the error bars
  summary_data <- raw_data %>%
    group_by(Drug, Dose) %>%
    summarise(
      Mean_Probability = mean(Probability, na.rm = TRUE),
      SE_Probability = sd(Probability, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # 5. Create the Bar Plot
  bar_plot <- ggplot(summary_data, aes(x = Drug, y = Mean_Probability, fill = Dose)) +
    # Draw the bars
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
             color = "black", linewidth = 0.3, alpha = 0.9) +
    
    # Draw the error bars (Mean +/- SEM)
    geom_errorbar(aes(ymin = Mean_Probability - SE_Probability, 
                      ymax = Mean_Probability + SE_Probability),
                  position = position_dodge(width = 0.8), 
                  width = 0.25, alpha = 0.7, linewidth = 0.5) +
    
    # Apply the discrete blue colors we used previously
    scale_fill_manual(values = c("0.1" = "#A6CEE3", "10" = "#1F78B4")) + 
    
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.grid.major.x = element_blank() # Keep x-axis clean
    ) +
    labs(fill = "Dose (\u03BCM)", 
         x = "Drug", 
         y = "Mean Efficacy Probability")
  
  return(bar_plot)
}

# Exceution
# refer following code files:
# Scoringmethods.r
# PredictionScore.r
rf_model 

# 1. Run the new extraction function to get the raw lists
prob_list_rf_raw <- predict_all_doses_raw(rf_model, sin.sct, cmb.sct)

# 2. Extract the individual doses
prob_sng10_rf_raw <- prob_list_rf_raw[[1]]
prob_sng0.1_rf_raw <- prob_list_rf_raw[[2]]

fig_4b_barplot <- plot_efficacy_barplot(prob_sng10_rf_raw, prob_sng0.1_rf_raw)

# 4. View the plot
ggsave("./Figure4_efficacyBarplot260517.pdf", plot = fig_4b_barplot, width = 7, height = 5.2, units = "in", dpi = 300)

### Figure 4b merged scatter plot
# ==============================================================================
# Script Name: phenotypic_screening.R
# Description: Functions to generate a 2D phenotypic screening plot (Figure 4b).
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggsci)
})

#' Generate and Save Phenotypic Screening Plot
#' 
#' @param seurat_obj A merged Seurat object with 'drug_name', 'dose', and 'drug_condition' in meta.data
#' @param genes_anabolic Vector of anabolic genes
#' @param genes_inflammatory Vector of inflammatory/catabolic genes
#' @param topn Number of anntation display of drug_condition
#' @return The ggplot object
plot_phenotypic_screening <- function(
    seurat_obj, 
    genes_anabolic, 
    genes_inflammatory,
    topn = 5) {
    
  # 1. Calculate coordinates and ratios
  DefaultAssay(seurat_obj) <- "SCT"
  exp_matrix <- GetAssayData(seurat_obj, slot = "data")
  
  # Calculate mean expression for the given gene sets per cell
  seurat_obj$Anabolic_Score <- colMeans(as.matrix(exp_matrix[genes_anabolic, , drop = FALSE]))
  seurat_obj$Catabolic_Score <- colMeans(as.matrix(exp_matrix[genes_inflammatory, , drop = FALSE]))
  
  # Summarize metrics grouping by Drug and Dose
  plot_data <- seurat_obj@meta.data %>%
    group_by(drug_name, dose, drug_condition) %>%
    summarise(
      mean_ana  = mean(Anabolic_Score, na.rm = TRUE),
      se_ana    = sd(Anabolic_Score, na.rm = TRUE) / sqrt(n()),
      mean_cata = mean(Catabolic_Score, na.rm = TRUE),
      se_cata   = sd(Catabolic_Score, na.rm = TRUE) / sqrt(n()),
      .groups   = 'drop'
    ) %>%
    mutate(
      AC_Ratio = mean_ana / mean_cata
    )
  
  # 2. Format Data & Identify Top/Bottom Hits for Labels
  top_hits    <- plot_data %>% top_n(topn, AC_Ratio) %>% pull(drug_condition)
  bottom_hits <- plot_data %>% top_n(topn, -AC_Ratio) %>% pull(drug_condition)
  controls    <- grep("control|inflammatory", plot_data$drug_condition, 
                      value = TRUE, ignore.case = TRUE)
  
  plot_data <- plot_data %>%
    mutate(
      Label_Text = ifelse(drug_condition %in% c(top_hits, bottom_hits, controls), 
                          drug_name, NA),
      
      # --- NEW: Categorize points into the 4 distinct groups for coloring ---
      Group = case_when(
        grepl("control", drug_condition, ignore.case = TRUE) ~ "Basal",
        grepl("inflammatory", drug_condition, ignore.case = TRUE) ~ "IL-1β",
        dose == 10 | dose == "10" ~ "10 µM",
        dose == 0.1 | dose == "0.1" ~ "0.1 µM",
        TRUE ~ "Other"
      ),
      # Set factor levels to control the order in the legend
      Group = factor(Group, levels = c("Basal", "IL-1β", "0.1 µM", "10 µM"))
    )  
  
  # 3. Define Discrete Colors
  condition_colors <- c(
    "Basal" = "#3C5488FF", # NPG Blue
    "IL-1β"           = "#E64B35FF", # NPG Red
    "0.1 µM"          = "#A6CEE3", # Sky blue 
    "10 µM"           = "#1f78b4"  # Steel blue
  )
  
  # 4. Generate Plot
  p <- ggplot(plot_data, aes(x = mean_cata, y = mean_ana)) +
    # Error bars
    geom_errorbar(aes(ymin = mean_ana - se_ana, ymax = mean_ana + se_ana), 
                  color = "grey80", width = 0) +
    geom_errorbarh(aes(xmin = mean_cata - se_cata, xmax = mean_cata + se_cata), 
                   color = "grey80", height = 0) +
    
    # Quadrant lines
    geom_vline(xintercept = mean(range(plot_data$mean_cata, na.rm = TRUE)), 
               linetype = "dotted", color = "grey60") +
    geom_hline(yintercept = mean(range(plot_data$mean_ana, na.rm = TRUE)), 
               linetype = "dotted", color = "grey60") +
    
    # Lines connecting the low and high doses of the same drug
    geom_line(aes(group = drug_name), color = "grey60", 
              linewidth = 0.5, alpha = 0.6) +
    
    # --- CHANGED: Points now use the discrete 'Group' column for fill color ---
    geom_point(aes(fill = Group), shape = 21, size = 6, alpha = 0.7, color = "black") +
    scale_fill_manual(name = "Condition", values = condition_colors) +
    
    # Text labels
    geom_text_repel(aes(label = Label_Text), size = 3.5, fontface = "bold",
                    box.padding = 0.6, max.overlaps = 50, min.segment.length = 0) + 
    
    theme_bw() +
    labs(x = "Catabolic gene sum (a.u.)", y = "Anabolic gene sum (a.u.)") +
    theme(
      panel.grid.minor = element_blank(), 
      axis.title = element_text(face = "bold"),
      legend.position = "right", 
      legend.background = element_rect(fill = "white", color = "grey90")
    )
  
  return(p)
}

mac_sin.sct <- merge(sin.sct, mac.sct)

fig4b_plot <- plot_phenotypic_screening(mac_sin.sct, anabolic, inflammatory, 13)
print(fig4b_plot)
ggsave("./Figure4_anaCata_2dScatter_260502.pdf", plot = fig4b_plot, width = 7, height = 5.2, units = "in", dpi = 300)




#####################################
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

 
### Figure 4d (differential co-expression analysis)
# ==============================================================================
# Description: Function to plot and save a dumbbell plot
# ==============================================================================
#' Automated Directional Co-expression Dumbbell Plot (with Healthy Target)
#'
#' Automatically computes the top 10 most divergent Pearson correlations between 
#' two single-cell populations and plots their trajectory, referencing a healthy baseline.
#'
#' @param seurat_obj A Seurat object containing single-cell expression data.
#' @param group_by Character. Metadata column containing condition labels. Default is "drug_condition".
#' @param cond_ref Character. Exact metadata string for the reference (disease) condition.
#' @param cond_treat Character. Exact metadata string for the treatment condition.
#' @param cond_healthy Character. Exact metadata string for the healthy/control condition.
#' @param label_ref Character. Legend display name. Defaults to cond_ref.
#' @param label_treat Character. Legend display name. Defaults to cond_treat.
#' @param label_healthy Character. Legend display name. Defaults to cond_healthy.
#' @param top_n Integer. Number of top divergent pairs to extract automatically. Default is 10.
#' @param assay Character. The Seurat assay to use. Default is "SCT".
#' @param slot Character. The data slot to use for correlation math. Default is "data".
#' @param color_ref Character. Hex code for reference dot. Default is NPG Red.
#' @param color_treat Character. Hex code for treatment dot. Default is NPG Blue.
#' @param color_healthy Character. Hex code for healthy target. Default is NPG Green.
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
                                        cond_healthy = "control",           # NEW
                                        label_ref = cond_ref,
                                        label_treat = cond_treat,
                                        label_healthy = cond_healthy,       # NEW
                                        top_n = 10,
                                        assay = "SCT", 
                                        slot = "data",
                                        color_ref = "#E64B35FF",
                                        color_treat = "#4DBBD5FF",
                                        color_healthy = "#00A087FF") {      # NEW (NPG Green)
  
  # ---------------------------------------------------------------------------
  # 1. INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!inherits(seurat_obj, "Seurat")) stop("Error: 'seurat_obj' must be a Seurat object.")
  if (!group_by %in% colnames(seurat_obj@meta.data)) stop(paste("Error: Column", group_by, "not found."))
  
  meta_data <- seurat_obj@meta.data
  if (!all(c(cond_ref, cond_treat, cond_healthy) %in% meta_data[[group_by]])) {
    stop("Error: One or more specified conditions not found within the group_by metadata column.")
  }
  
  # ---------------------------------------------------------------------------
  # 2. DATA EXTRACTION
  # ---------------------------------------------------------------------------
  DefaultAssay(seurat_obj) <- assay
  matrix_data <- GetAssayData(seurat_obj, slot = slot)
  
  cells_ref     <- rownames(meta_data[meta_data[[group_by]] == cond_ref, ])
  cells_treat   <- rownames(meta_data[meta_data[[group_by]] == cond_treat, ])
  cells_healthy <- rownames(meta_data[meta_data[[group_by]] == cond_healthy, ])
  
  # Convert to standard dense matrices for cor()
  mat_ref     <- as.matrix(matrix_data[, cells_ref])
  mat_treat   <- as.matrix(matrix_data[, cells_treat])
  mat_healthy <- as.matrix(matrix_data[, cells_healthy])
  
  # Filter out genes with zero variance in ANY condition (prevents NA in Pearson)
  var_ref     <- apply(mat_ref, 1, var)
  var_treat   <- apply(mat_treat, 1, var)
  var_healthy <- apply(mat_healthy, 1, var)
  valid_genes <- names(which(var_ref > 0 & var_treat > 0 & var_healthy > 0))
  
  mat_ref     <- mat_ref[valid_genes, ]
  mat_treat   <- mat_treat[valid_genes, ]
  mat_healthy <- mat_healthy[valid_genes, ]
  
  # ---------------------------------------------------------------------------
  # 3. VECTORIZED CORRELATION MATH
  # ---------------------------------------------------------------------------
  cor_ref_mat     <- cor(t(mat_ref), method = "pearson")
  cor_treat_mat   <- cor(t(mat_treat), method = "pearson")
  cor_healthy_mat <- cor(t(mat_healthy), method = "pearson")
  
  # Extract the upper triangle
  upper_idx <- which(upper.tri(cor_ref_mat), arr.ind = TRUE)
  
  df_all <- data.frame(
    Gene1     = rownames(cor_ref_mat)[upper_idx[, 1]],
    Gene2     = colnames(cor_ref_mat)[upper_idx[, 2]],
    R_Ref     = cor_ref_mat[upper_idx],
    R_Treat   = cor_treat_mat[upper_idx],
    R_Healthy = cor_healthy_mat[upper_idx],
    stringsAsFactors = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # 4. AUTOMATED TOP 'N' DISCOVERY
  # ---------------------------------------------------------------------------
  df_plot <- df_all %>%
    drop_na() %>%
    mutate(
      GenePair = paste0(Gene1, " - ", Gene2),
      Delta_r  = abs(R_Treat - R_Ref), # Still ranking by how much the drug moved the needle
      Category = paste("Top", top_n, "Divergent Interactions") 
    ) %>%
    arrange(desc(Delta_r)) %>%
    slice_head(n = top_n) %>%
    arrange(Delta_r) 
  
  # Lock factor levels
  df_plot$GenePair <- factor(df_plot$GenePair, levels = df_plot$GenePair)
  df_plot$Category <- factor(df_plot$Category, levels = unique(df_plot$Category))
  
  # Dynamic Axis Scaling (Now includes healthy values)
  max_val <- max(abs(c(df_plot$R_Ref, df_plot$R_Treat, df_plot$R_Healthy)), na.rm = TRUE)
  axis_limit <- ceiling(max_val * 10) / 10 + 0.1 
  
  # ---------------------------------------------------------------------------
  # 5. VISUALIZATION
  # ---------------------------------------------------------------------------
  color_mapping <- setNames(
    c(color_ref, color_treat, color_healthy), 
    c(label_ref, label_treat, label_healthy)
  )
  
  p_dumbbell <- ggplot(df_plot) +
    
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    
    # 1. Draw the arrow from Disease -> Treatment
    geom_segment(aes(x = R_Ref, xend = R_Treat, y = GenePair, yend = GenePair),
                 color = "grey40", linewidth = 1.2, 
                 arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
                 
    # 2. Draw the Healthy Target (Using shape 23 = diamond to distinguish it)
    geom_point(aes(x = R_Healthy, y = GenePair, fill = !!label_healthy), 
               size = 4, shape = 23, color = "black", stroke = 1, alpha = 0.85) +
    
    # 3. Draw the Reference (Disease) dot
    geom_point(aes(x = R_Ref, y = GenePair, fill = !!label_ref), 
               size = 4, shape = 21, color = "black", stroke = 1) +
    
    # 4. Draw the Treatment point on top
    geom_point(aes(x = R_Treat, y = GenePair, fill = !!label_treat), 
               size = 4, shape = 21, color = "black", stroke = 1) +
    
    scale_fill_manual(values = color_mapping) +
    scale_x_continuous(limits = c(-axis_limit, axis_limit), 
                       breaks = round(seq(-axis_limit, axis_limit, length.out = 5), 2)) +
    
    theme_bw() +
    labs(
      title = "Data-Driven Differential Co-expression",
      subtitle = paste0("Trajectory from ", label_ref, " towards ", label_healthy, " via ", label_treat),
      x = "Pearson Correlation Coefficient (r)",
      y = ""
    ) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 11, face = "bold"),
      
      axis.text.y = element_text(size = 12, face = "bold.italic", color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.title.x = element_text(size = 12, face = "bold", margin = ggplot2::margin(t = 10)),
      
      strip.text.y.left = element_text(angle = 0, size = 11, face = "bold"),
      strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1),
      strip.placement = "outside",
      
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 1)
    )
  
  return(p_dumbbell)
}

fig_4d_dumbbel <- plot_automated_coexpression(merged.sct)
ggsave("./fig4_dumbbel_Coexpression_analysis_rapamycin_10_260502.pdf", plot= fig_4d_dumbbel, width = 7, height = 5.5, units = "in", dpi = 300)
dev.off()


### Figure 4e (radar chart)
# ==============================================================================
# Description: Function to plot and save a radar chart of GO modality scores
# ==============================================================================
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

