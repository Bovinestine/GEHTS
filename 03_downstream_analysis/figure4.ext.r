# Project name: GEHTS-chip
# Author: Nathan Wooseok Lee
# conda env: seurat4
# date: 260302

# Load necessary libraries
library(Seurat)
library(dplyr)
library(reshape2)
library(igraph)
library(ggraph)
library(ggplot2)
library(cowplot)

#' Generate Differential Co-expression Network Graph
#'
#' Constructs and plots side-by-side circular network graphs to compare 
#' gene co-expression topologies between two single-cell populations.
#'
#' @param seurat_obj A Seurat object containing single-cell expression data.
#' @param cond1_name Character. Exact metadata string for the reference condition.
#' @param cond2_name Character. Exact metadata string for the treatment condition.
#' @param title1 Character. Title for the first network panel. Default is "State: [cond1_name]".
#' @param title2 Character. Title for the second network panel. Default is "State: [cond2_name]".
#' @param cor_threshold Numeric. Absolute Pearson correlation threshold to draw an edge. Default is 0.35.
#' @param group_by Character. The metadata column containing the condition labels. Default is "drug_condition".
#' @param genes_anabolic Character vector. Genes belonging to the anabolic/matrix module.
#' @param genes_catabolic Character vector. Genes belonging to the catabolic/inflammatory module.
#' @param genes_housekeeping Character vector. Genes belonging to the housekeeping module.
#' @param assay Character. The Seurat assay to use. Default is "SCT".
#' @param slot Character. The data slot to use for correlations. Default is "data".
#'
#' @return A cowplot grid object containing two ggraph network plots and a shared legend.
#' @export
#'
generate_differential_network <- function(seurat_obj, 
                                          cond1_name, 
                                          cond2_name, 
                                          title1 = NULL, 
                                          title2 = NULL, 
                                          cor_threshold = 0.35,
                                          group_by = "drug_condition",
                                          assay = "SCT",
                                          slot = "data",
                                          genes_anabolic = c('Acan','Sox9','Col2a1','Matn1','Matn3','Ucma','Ccnd3','Gadd45g','Pth1r','Gm26633','Col27a1'),
                                          genes_catabolic = c('Mmp3','Mmp13','Il6','Il17b','Adamts5','Igfbp3','Ccl2','Cxcl5','Cxcl1','Fosl2','Tlr2','Tnfrsf1b'),
                                          genes_housekeeping = c('Hprt','Actb','Gapdh','B2m','Ubc','Ppia','Rpl23')) {
  
  # ---------------------------------------------------------------------------
  # 1. INPUT VALIDATION & SETUP
  # ---------------------------------------------------------------------------
  if (!inherits(seurat_obj, "Seurat")) stop("Error: 'seurat_obj' must be a Seurat object.")
  
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Error: Column", group_by, "not found in Seurat metadata."))
  }
  
  meta_conditions <- unique(seurat_obj@meta.data[[group_by]])
  if (!cond1_name %in% meta_conditions || !cond2_name %in% meta_conditions) {
    stop("Error: cond1_name or cond2_name not found within the specified group_by metadata column.")
  }
  
  if (is.null(title1)) title1 <- paste("State:", cond1_name)
  if (is.null(title2)) title2 <- paste("State:", cond2_name)
  
  all_genes <- c(genes_anabolic, genes_catabolic, genes_housekeeping)
  
  # ---------------------------------------------------------------------------
  # 2. INTERNAL HELPER: CALCULATE NETWORK
  # ---------------------------------------------------------------------------
  get_network_data <- function(condition_name) {
    # Isolate cells for the specific condition
    cells_use <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[group_by]] == condition_name, ])
    
    DefaultAssay(seurat_obj) <- assay
    avail_genes <- intersect(all_genes, rownames(GetAssayData(seurat_obj, slot = slot)))
    
    if(length(avail_genes) < 2) stop("Not enough target genes found in the Seurat object to build a network.")
    
    expr_mat <- GetAssayData(seurat_obj, slot = slot)[avail_genes, cells_use]
    expr_mat <- as.matrix(expr_mat)
    
    # Calculate Edges (Pearson Correlation)
    cor_mat <- suppressWarnings(cor(t(expr_mat)))
    cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
    
    edges <- melt(cor_mat, na.rm = TRUE)
    colnames(edges) <- c("from", "to", "weight")
    edges$from <- as.character(edges$from)
    edges$to <- as.character(edges$to)
    
    edges <- edges %>%
      filter(abs(weight) >= cor_threshold) %>%
      mutate(
        sign = ifelse(weight > 0, "Positive", "Negative"),
        abs_weight = abs(weight)
      )
    
    # Calculate Nodes (Mean Expression and Module Assignment)
    nodes <- data.frame(
      name = as.character(rownames(expr_mat)), 
      MeanExpr = as.numeric(rowMeans(expr_mat, na.rm = TRUE))
    ) %>%
      mutate(Module = case_when(
        name %in% genes_anabolic ~ "Anabolic",
        name %in% genes_catabolic ~ "Catabolic",
        name %in% genes_housekeeping ~ "Housekeeping",
        TRUE ~ "Other"
      ))
    
    return(graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE))
  }
  
  # ---------------------------------------------------------------------------
  # 3. INTERNAL HELPER: PLOT NETWORK
  # ---------------------------------------------------------------------------
  plot_network <- function(graph_obj, plot_title) {
    ggraph(graph_obj, layout = 'linear', circular = TRUE) + 
      geom_edge_arc(aes(edge_width = abs_weight, color = sign), alpha = 0.6) +
      scale_edge_width(range = c(0.5, 2.5), guide = "none") +
      scale_edge_color_manual(values = c("Positive" = "#E64B35FF", "Negative" = "#3C5488FF"), 
                              name = "Correlation") +
      geom_node_point(aes(size = MeanExpr, fill = Module), 
                      shape = 21, color = "black", stroke = 1) +
      scale_size(range = c(4, 10), name = "Mean Expression") +
      scale_fill_manual(values = c("Anabolic" = "#E64B35FF", 
                                   "Catabolic" = "#3C5488FF", 
                                   "Housekeeping" = "grey70")) +
      
      # Pushes x and y outward by 15% and justifies text away from the center
      geom_node_text(aes(x = x * 1.15, y = y * 1.15, 
                         label = name, 
                         hjust = ifelse(x > 0, 0, 1)), 
                     size = 3.5, fontface = "bold.italic") +
      
      # Expands the plot boundary so labels do not get cut off
      coord_cartesian(clip = "off") +
      scale_x_continuous(expand = expansion(mult = 0.2)) +
      scale_y_continuous(expand = expansion(mult = 0.2)) +      
      theme_graph() + 
      labs(title = plot_title) +
      theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")
      )
  }
  
  # ---------------------------------------------------------------------------
  # 4. EXECUTE PIPELINE & COMBINE
  # ---------------------------------------------------------------------------
  g1 <- get_network_data(cond1_name)
  g2 <- get_network_data(cond2_name)
  
  p1 <- plot_network(g1, title1)
  p2 <- plot_network(g2, title2)
  
  # Extract legend and combine with cowplot
  shared_legend <- get_legend(p1 + theme(legend.box.margin = margin(10, 0, 0, 0)))
  
  plot_row <- plot_grid(
    p1 + theme(legend.position = "none"),
    p2 + theme(legend.position = "none"),
    labels = c("a", "b"), # Nature journals prefer lowercase letters for panels
    label_size = 16,
    ncol = 2
  )
  
  final_fig <- plot_grid(
    plot_row, 
    shared_legend, 
    ncol = 1, 
    rel_heights = c(1, 0.15)
  )
  
  return(final_fig)
}

# ---- Excuetion example ----

# 1. Generate the Extended Data Figure 4
fig_4ext <- generate_differential_network(
  seurat_obj = merged.sct,
  cond1_name = "inflammatory",
  cond2_name = "rapamycin_10",
  title1 = "IL-1β (Disease)",
  title2 = "Rapamycin (10 μM)",
  cor_threshold = 0.35
)

# 2. Export to SVG for publication layout (e.g., Illustrator)
ggsave(
  filename = "./figures/figure4/fig4_ext_Coexpression_analysis_rapamycin_10.svg", 
  plot = fig_4ext, 
  width = 11,     # Widened slightly to ensure circular labels aren't clipped
  height = 5.5, 
  units = "in", 
  dpi = 300
)


############ p38 inhibitor comparison
#' 2D Single-Cell Phase Portrait (Density Contours)
#'
#' Generates a publication-ready 1x4 faceted 2D density contour plot to visualize 
#' transcriptional decoupling or coupling between two genes at single-cell resolution.
#' Background cells are plotted as transparent points, overlaid with topological 
#' density contours to indicate the population's "center of gravity."
#'
#' @param seurat_obj A Seurat object containing single-cell expression data.
#' @param gene_x Character. Gene symbol for the X-axis (e.g., the Master Regulator).
#' @param gene_y Character. Gene symbol for the Y-axis (e.g., the Effector/Matrix gene).
#' @param group_by Character. Metadata column name containing the condition labels. Default is "drug_condition".
#' @param cond_ctrl Character. Exact metadata string for the healthy/control condition.
#' @param cond_inf Character. Exact metadata string for the diseased/inflammatory condition.
#' @param cond_A Character. Exact metadata string for Treatment A (e.g., Pamapimod).
#' @param cond_B Character. Exact metadata string for Treatment B (e.g., SB203580).
#' @param label_ctrl Character. Display name for Control.
#' @param label_inf Character. Display name for Disease.
#' @param label_A Character. Display name for Treatment A.
#' @param label_B Character. Display name for Treatment B.
#' @param assay Character. The Seurat assay to use. Default is "SCT".
#' @param slot Character. The data slot to use. Default is "data" (log-normalized).
#'
#' @return A ggplot2 object.
#' @export
#'
#' @import Seurat
#' @import dplyr
#' @import ggplot2
# Generate the Main Figure Panel
fig_phase_portrait <- plot_sc_phase_portrait(
  seurat_obj = merged.sct,
  gene_x = "Sox9",
  gene_y = "Col2a1",
  cond_A = "pamapimod_10",
  cond_B = "SB203580_10",
  label_A = "Pamapimod (Decoupled)",
  label_B = "SB203580 (Rescued)"
)

fig_facs <- plot_sc_facs_overlay(
  seurat_obj = merged.sct,
  gene_x = "Sox9",
  gene_y = "Col2a1",
  gate_x = 1.5, # Draws a vertical line at Sox9 = 1.5
  gate_y = 1.0  # Draws a horizontal line at Col2a1 = 1.0
)

# Export for publication
cairo_pdf("./figures/main/Fig3_Phase_Portrait.pdf", width = 12, height = 4)
print(fig_phase_portrait)
dev.off()

plot_sc_facs_overlay <- function(seurat_obj, 
                                 gene_x = "Sox9", 
                                 gene_y = "Col2a1",
                                 group_by = "drug_condition",
                                 cond_ctrl = "control",
                                 cond_inf = "inflammatory",
                                 cond_A = "pamapimod_10",
                                 cond_B = "SB203580_10",
                                 label_ctrl = "Control (Healthy)",
                                 label_inf = "IL-1β (Disease)",
                                 label_A = "Pamapimod (Decoupled)",
                                 label_B = "SB203580 (Rescued)",
                                 gate_x = NULL,  # e.g., 1.5
                                 gate_y = NULL,  # e.g., 1.0
                                 assay = "SCT", 
                                 slot = "data") {
  
  # 1. INPUT VALIDATION
  if (!inherits(seurat_obj, "Seurat")) stop("Error: 'seurat_obj' must be a Seurat object.")
  DefaultAssay(seurat_obj) <- assay
  
  # 2. DATA EXTRACTION
  mat_subset <- GetAssayData(seurat_obj, slot = slot)[c(gene_x, gene_y), , drop = FALSE]
  
  df_plot <- data.frame(
    Expr_X = as.numeric(mat_subset[gene_x, ]),
    Expr_Y = as.numeric(mat_subset[gene_y, ]),
    Condition_Raw = seurat_obj@meta.data[[group_by]],
    stringsAsFactors = FALSE
  )
  
  target_conditions <- c(cond_ctrl, cond_inf, cond_A, cond_B)
  df_plot <- df_plot %>% filter(Condition_Raw %in% target_conditions)
  
  label_mapping <- setNames(
    c(label_ctrl, label_inf, label_A, label_B),
    target_conditions
  )
  df_plot$Condition <- label_mapping[df_plot$Condition_Raw]
  
  # Z-ordering: We want the drugs plotted ON TOP of the control/disease baselines
  df_plot$Condition <- factor(df_plot$Condition, 
                              levels = c(label_ctrl, label_inf, label_A, label_B))
  
  # 3. NPG COLORS
  npg_colors <- setNames(
    c("#4DBBD5FF", "#E64B35FF", "#8491B4FF", "#00A087FF"), # Blue, Red, Grey/Blue, Teal
    c(label_ctrl, label_inf, label_A, label_B)
  )
  
  # 4. VISUALIZATION
  p_facs <- ggplot(df_plot, aes(x = Expr_X, y = Expr_Y, color = Condition)) +
    
    # Optional FACS Quadrant Gating Lines
    {if(!is.null(gate_x)) geom_vline(xintercept = gate_x, linetype = "dashed", color = "grey50", linewidth = 0.8)} +
    {if(!is.null(gate_y)) geom_hline(yintercept = gate_y, linetype = "dashed", color = "grey50", linewidth = 0.8)} +
    
    # B. The Single Cells (Points)
    # alpha = 0.5 makes them semi-transparent so overlaps create natural density
    # size = 1.5 is small enough to avoid hairballs with 200 cells
    geom_point(alpha = 0.5, size = 1.5, shape = 16) +

    # Overlaid Density Contours (Fewer bins = cleaner overlap)
    geom_density_2d(linewidth = 1.2, bins = 5) +
    
    scale_color_manual(values = npg_colors) +
    
    theme_classic() +
    labs(
      title = "FACS-Style Population Trajectories",
      subtitle = "Overlaid contours map the coordinate shift of the master regulator vs. effector matrix",
      x = paste(gene_x, "Expression (Log-Normalized)"),
      y = paste(gene_y, "Expression (Log-Normalized)")
    ) +
    theme(
      # Place legend cleanly inside the plot (like flow cytometry plots) or at the top
      legend.position = "right", 
      legend.title = element_blank(),
      legend.text = element_text(size = 12, face = "bold"),
      
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
      axis.text = element_text(size = 12, color = "black"),
      
      # Make it square to emphasize the Cartesian quadrant space
      aspect.ratio = 1,
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
    )
  
  return(p_facs)
}


plot_sc_phase_portrait <- function(seurat_obj, 
                                   gene_x = "Sox9", 
                                   gene_y = "Col2a1",
                                   group_by = "drug_condition",
                                   cond_ctrl = "control",
                                   cond_inf = "inflammatory",
                                   cond_A = "pamapimod_10",
                                   cond_B = "SB203580_10",
                                   label_ctrl = "Control (Healthy)",
                                   label_inf = "IL-1β (Disease)",
                                   label_A = "Pamapimod (10 μM)",
                                   label_B = "SB203580 (10 μM)",
                                   assay = "SCT", 
                                   slot = "data") {
  
  # ---------------------------------------------------------------------------
  # 1. INPUT VALIDATION (Defensive Programming)
  # ---------------------------------------------------------------------------
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Error: 'seurat_obj' must be a Seurat object.")
  }
  
  DefaultAssay(seurat_obj) <- assay
  available_genes <- rownames(GetAssayData(seurat_obj, slot = slot))
  
  if (!gene_x %in% available_genes || !gene_y %in% available_genes) {
    stop(paste("Error: Gene(s)", gene_x, "and/or", gene_y, "not found in the specified assay/slot."))
  }
  
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Error: Column", group_by, "not found in Seurat metadata."))
  }
  
  target_conditions <- c(cond_ctrl, cond_inf, cond_A, cond_B)
  meta_conditions <- unique(seurat_obj@meta.data[[group_by]])
  missing_conds <- setdiff(target_conditions, meta_conditions)
  
  if (length(missing_conds) > 0) {
    stop(paste("Error: The following conditions were not found in metadata:", 
               paste(missing_conds, collapse = ", ")))
  }

  # ---------------------------------------------------------------------------
  # 2. DATA EXTRACTION & WRANGLING
  # ---------------------------------------------------------------------------
  # Extract specific genes to save memory (avoids pulling the whole matrix)
  mat_subset <- GetAssayData(seurat_obj, slot = slot)[c(gene_x, gene_y), , drop = FALSE]
  
  # Build a lightweight dataframe for ggplot
  df_plot <- data.frame(
    Cell = colnames(seurat_obj),
    Expr_X = as.numeric(mat_subset[gene_x, ]),
    Expr_Y = as.numeric(mat_subset[gene_y, ]),
    Condition_Raw = seurat_obj@meta.data[[group_by]],
    stringsAsFactors = FALSE
  )
  
  # Filter to only the 4 requested conditions
  df_plot <- df_plot %>% filter(Condition_Raw %in% target_conditions)
  
  # Map raw metadata strings to publication-ready labels
  label_mapping <- setNames(
    c(label_ctrl, label_inf, label_A, label_B),
    target_conditions
  )
  df_plot$Condition <- label_mapping[df_plot$Condition_Raw]
  
  # Lock factor levels to enforce the 1x4 logical flow: Control -> Disease -> Drug A -> Drug B
  df_plot$Condition <- factor(df_plot$Condition, 
                              levels = c(label_ctrl, label_inf, label_A, label_B))
  
  # ---------------------------------------------------------------------------
  # 3. VISUALIZATION (Nature Publishing Group Aesthetics)
  # ---------------------------------------------------------------------------
  # Assign distinct NPG colors to each condition for the contours
  npg_colors <- setNames(
    c("#4DBBD5FF", "#E64B35FF", "#8491B4FF", "#00A087FF"), 
    c(label_ctrl, label_inf, label_A, label_B)
  )
  
  p_phase <- ggplot(df_plot, aes(x = Expr_X, y = Expr_Y)) +
    
    # A. Background Single Cells (Highly transparent, neutral color to show volume without clutter)
    geom_point(color = "grey80", alpha = 0.3, size = 0.5, shape = 16) +
    
    # B. Topological Density Contours (Colored by condition, bold for visual impact)
    geom_density_2d(aes(color = Condition), linewidth = 1.0, bins = 6) +
    
    # C. Scales and Faceting
    scale_color_manual(values = npg_colors) +
    facet_wrap(~ Condition, nrow = 1) + 
    
    # D. Formatting
    theme_classic() +
    labs(
      title = "Single-Cell Phase Portraits: Functional Matrix Coordination",
      subtitle = "Density contours reveal decoupling of master regulator and effector matrix gene",
      x = paste(gene_x, "Expression (Log-Normalized)"),
      y = paste(gene_y, "Expression (Log-Normalized)")
    ) +
    theme(
      legend.position = "none", # Legend is redundant because facets are labeled
      
      # Axis titles: Make the gene names bold and italicized using expression() styling
      axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
      axis.text = element_text(size = 10, color = "black"),
      
      # Clean facet headers
      strip.text = element_text(face = "bold", size = 11, color = "black"),
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 1),
      
      # Strict box around the panels
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  return(p_phase)
}

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

fig_p38dumbbell <- plot_decoupling_dumbbell(merged.sct)

plot_decoupling_dumbbell <- function(seurat_obj, 
                                     cond_inf = "inflammatory",
                                     cond_ctrl = "control",
                                     cond_A = "pamapimod_10",
                                     cond_B = "SB203580_10",
                                     label_ctrl = "Control",
                                     label_A = "Pamapimod\n(Decoupled)",
                                     label_B = "SB203580\n(Coordinated Rescue)",
                                     assay = "SCT", 
                                     slot = "counts") { # FIX 1: Restored strictly to "counts"
  
  # 1. Extract Data
  DefaultAssay(seurat_obj) <- assay
  
  # Fetch counts directly (No expm1 needed!)
  mat <- as.matrix(GetAssayData(seurat_obj, slot = slot)[c("Sox9", "Col2a1"), , drop = FALSE])
  meta <- seurat_obj@meta.data
  
  # Identify cells
  cells_inf <- rownames(meta[meta$drug_condition == cond_inf, ])
  cells_ctrl <- rownames(meta[meta$drug_condition == cond_ctrl, ])
  cells_A   <- rownames(meta[meta$drug_condition == cond_A, ])
  cells_B   <- rownames(meta[meta$drug_condition == cond_B, ])
  
  if(length(cells_inf) == 0) stop("Inflammatory cells not found!")
  
  # 2. Calculate Means (Matches the stable math from our earlier boxplots)
  mean_inf <- rowMeans(mat[, cells_inf, drop = FALSE])
  mean_ctrl <- rowMeans(mat[, cells_ctrl, drop = FALSE])
  mean_A   <- rowMeans(mat[, cells_A, drop = FALSE])
  mean_B   <- rowMeans(mat[, cells_B, drop = FALSE])
  
  # 3. Calculate Log2 Fold Change vs Disease Baseline (Standard pseudo-count of 1)
  pseudo <- 1 
  lfc_ctrl <- log2((mean_ctrl + pseudo) / (mean_inf + pseudo))
  lfc_A <- log2((mean_A + pseudo) / (mean_inf + pseudo))
  lfc_B <- log2((mean_B + pseudo) / (mean_inf + pseudo))
  
  # 4. Build Plotting Data Frame
  df_plot <- data.frame(
    Drug = c(label_A, label_B, label_ctrl),
    Sox9_LFC = c(lfc_A["Sox9"], lfc_B["Sox9"], lfc_ctrl["Sox9"]),
    Col2a1_LFC = c(lfc_A["Col2a1"], lfc_B["Col2a1"], lfc_ctrl["Col2a1"])
  )
  
  df_plot$Drug <- factor(df_plot$Drug, levels = c(label_A, label_B, label_ctrl))
  
  df_long <- pivot_longer(df_plot, cols = c(Sox9_LFC, Col2a1_LFC), 
                          names_to = "Gene", values_to = "LFC")
  
  df_long$Gene <- gsub("_LFC", "", df_long$Gene)
  df_long$Gene <- factor(df_long$Gene, levels = c("Sox9", "Col2a1"))
  
  # 5. Visual Aesthetics
  gene_colors <- c("Sox9" = "#00A087FF",     
                   "Col2a1" = "#4DBBD5FF")   
  
  # 6. Generate Plot
  p_dumbbell <- ggplot() +
    
    # The Disease Baseline
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 1) +
    
    # The "Decoupling Gap"
    geom_segment(data = df_plot, 
                 aes(x = Sox9_LFC, xend = Col2a1_LFC, y = Drug, yend = Drug),
                 color = "grey60", linewidth = 2) +
    
    # The Genes
    geom_point(data = df_long, 
               aes(x = LFC, y = Drug, fill = Gene), 
               size = 6, shape = 21, color = "black", stroke = 1.2) +
    
    scale_fill_manual(values = gene_colors) +
    theme_classic() +
    
    # FIX 2: Removed the crashing annotate() function entirely. 
    # Moved the baseline definition to the subtitle for a cleaner look.
    labs(
      title = "Divergent Resolution of Anabolic Decoupling",
      subtitle = "Log2(FC) vs. inflammatory baseline (Dashed line = 0)",
      x = "Log2 Fold Change (vs. Disease)",
      y = ""
    ) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 12, face = "bold.italic"),
      
      axis.text.y = element_text(size = 12, face = "bold", color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 15)),
      
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
    )
  
  return(p_dumbbell)
}

# fig_box <- plot_decoupling_boxplots(merged.sct)
# print(fig_box)


# previous figures

library(ggplot2)
library(dplyr)
library(tidyr)

plot_mechanistic_divergence <- function(seurat_obj, drug_A = "SB203580_10", drug_B = "pamapimod_10", ref = "inflammatory") {
  
  # 1. Extract SCT Data
  mat <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")
  meta <- seurat_obj@meta.data
  
  # Defined Genes of Interest
  genes_use <- c('Acan','Sox9','Col2a1','Mmp13','Il6','Adamts5','Ccl2','Cxcl5','Fosl2','Tlr2')
  genes_use <- intersect(genes_use, rownames(mat))
  
  # 2. Calculate Mean Expression per Group
  mean_ref <- rowMeans(mat[genes_use, rownames(meta[meta$drug_condition == ref, ]), drop=FALSE])
  mean_A <- rowMeans(mat[genes_use, rownames(meta[meta$drug_condition == drug_A, ]), drop=FALSE])
  mean_B <- rowMeans(mat[genes_use, rownames(meta[meta$drug_condition == drug_B, ]), drop=FALSE])
  
  # 3. Calculate Log2 Fold Change relative to Inflammatory baseline
  # (Adding small pseudo-count to avoid log(0))
  lfc_A <- log2((mean_A + 0.01) / (mean_ref + 0.01))
  lfc_B <- log2((mean_B + 0.01) / (mean_ref + 0.01))
  
  # 4. Prepare DataFrame
  plot_data <- data.frame(
    Gene = rep(genes_use, 2),
    Log2FC = c(lfc_A, lfc_B),
    Drug = rep(c("SB203580 (p38)", "Pamapimod (p38)"), each = length(genes_use))
  )
  
  # Order genes for clean plotting (e.g., sort by SB203580 response)
  gene_order <- plot_data %>% filter(Drug == "SB203580 (p38)") %>% arrange(Log2FC) %>% pull(Gene)
  plot_data$Gene <- factor(plot_data$Gene, levels = gene_order)
  
  # 5. Generate Lollipop Plot
  ggplot(plot_data, aes(x = Gene, y = Log2FC, color = Drug)) +
    # Draw the stems
    geom_segment(aes(x = Gene, xend = Gene, y = 0, yend = Log2FC), 
                 position = position_dodge(width = 0.5), linewidth = 1) +
    # Draw the lollipops
    geom_point(position = position_dodge(width = 0.5), size = 4) +
    
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    coord_flip() +
    scale_color_manual(values = c("SB203580 (p38)" = "#F39B7FFF", "Pamapimod (p38)" = "#8491B4FF")) +
    
    theme_classic() +
    theme(
      axis.text.y = element_text(face = "italic", size = 12, color = "black"),
      axis.title = element_text(face = "bold", size = 12),
      legend.position = "top",
      legend.title = element_blank()
    ) +
    labs(
      title = "Intra-class Mechanistic Resolution",
      subtitle = "Log2 Fold Change vs. Inflammatory Baseline",
      y = "Log2(Fold Change)",
      x = "Gene"
    )
}

# Run the plot
# Replace with the exact strings from your metadata!
# fig_3_lollipop <- plot_mechanistic_divergence(merged.sct, drug_A = "SB203580_10", drug_B = "pamapimod_10")
# print(fig_3_lollipop)

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork) # Added for proportional combining

# --- Execution ---
fig_3I <- plot_relative_sc_boxplot(merged.sct)
svglite("./figures/figure3/p38_relativeCompare.svg", width = 6, height = 5)
print(fig_3I)
dev.off()


plot_relative_sc_boxplot <- function(seurat_obj, 
                                     cond_ctrl = "control",
                                     cond_inf = "inflammatory",
                                     cond_A = "pamapimod_10", 
                                     cond_B = "SB203580_10",
                                     assay = "SCT") { 
  
  # 1. Custom Gene Selection 
  target_genes <- c("Sox9", "Col2a1", "Acan", "Matn1", 
                    "Mmp13", "Adamts5", "Il6", "Cxcl1")
  
  # 2. Extract CORRECTED COUNTS 
  DefaultAssay(seurat_obj) <- assay
  mat_counts <- as.matrix(GetAssayData(seurat_obj, slot = "counts"))
  avail_genes <- intersect(target_genes, rownames(mat_counts))
  mat_sub <- mat_counts[avail_genes, , drop = FALSE]
  
  meta <- seurat_obj@meta.data
  
  # 3. Identify Cells
  cells_ctrl <- rownames(meta[meta$drug_condition == cond_ctrl, ])
  cells_inf  <- rownames(meta[meta$drug_condition == cond_inf, ])
  cells_A    <- rownames(meta[meta$drug_condition == cond_A, ])
  cells_B    <- rownames(meta[meta$drug_condition == cond_B, ])
  
  if(length(cells_inf) == 0) stop("Inflammatory reference cells not found!")
  
  # 4. Calculate Inflammatory Baseline
  mean_inf <- rowMeans(mat_sub[, cells_inf, drop = FALSE])
  
  # 5. Calculate Single-Cell Log2 Fold Change 
  pseudo <- 1 
  
  calc_sc_lfc <- function(cells) {
    if(length(cells) == 0) return(NULL)
    cell_mat <- mat_sub[, cells, drop = FALSE]
    res <- sweep(cell_mat + pseudo, MARGIN = 1, STATS = mean_inf + pseudo, FUN = "/")
    return(log2(res))
  }
  
  lfc_ctrl <- calc_sc_lfc(cells_ctrl)
  lfc_A    <- calc_sc_lfc(cells_A)
  lfc_B    <- calc_sc_lfc(cells_B)
  
  # 6. Build the Dataframe safely
  plot_data <- data.frame()
  if(!is.null(lfc_ctrl)) { df_ctrl <- melt(t(lfc_ctrl)); df_ctrl$Condition <- cond_ctrl; plot_data <- rbind(plot_data, df_ctrl) }
  if(!is.null(lfc_A))    { df_A    <- melt(t(lfc_A));    df_A$Condition <- cond_A;       plot_data <- rbind(plot_data, df_A) }
  if(!is.null(lfc_B))    { df_B    <- melt(t(lfc_B));    df_B$Condition <- cond_B;       plot_data <- rbind(plot_data, df_B) }
  colnames(plot_data) <- c("Cell", "Gene", "LFC", "Condition")
  
  # 7. Categorize Genes
  cat_tf <- "Master TF\n(The Paradox)"
  cat_mat <- "Matrix Synthesis\n(The Rescue)"
  cat_inf <- "Inflammatory\n(Shared Target)"
  
  plot_data <- plot_data %>%
    mutate(Category = case_when(
      Gene == "Sox9" ~ cat_tf,
      Gene %in% c("Col2a1", "Acan", "Matn1") ~ cat_mat,
      Gene %in% c("Mmp13", "Adamts5", "Il6", "Cxcl1") ~ cat_inf,
      TRUE ~ "Other"
    ))
  
  plot_data$Category <- factor(plot_data$Category, levels = c(cat_tf, cat_mat, cat_inf))
  plot_data$Condition <- factor(plot_data$Condition, levels = c(cond_ctrl, cond_A, cond_B))
  gene_order <- c("Sox9", "Col2a1", "Acan", "Matn1", "Mmp13", "Adamts5", "Il6", "Cxcl1")
  plot_data$Gene <- factor(plot_data$Gene, levels = intersect(gene_order, avail_genes))
  
  # 8. Dynamic NPG Colors
  npg_colors <- setNames(
    c("#4DBBD5FF", "#8491B4FF", "#E64B35FF"), 
    c(cond_ctrl, cond_A, cond_B)
  )
  
  # 9. Helper function to generate clean individual plots
  create_panel <- function(df_sub, show_y_label = FALSE) {
    p <- ggplot(df_sub, aes(x = Gene, y = LFC, fill = Condition)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
      geom_boxplot(width = 0.7, outlier.shape = NA, color = "black", linewidth = 0.5, 
                   position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = npg_colors) +
      facet_wrap(~ Category, scales = "free") +
      theme_classic() +
      theme(
        legend.position = "none", # Legend removed from individual plots!
        axis.text.x = element_text(face = "bold.italic", size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(face = "bold", size = 11, color = "black"),
        strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      )
    
    if (show_y_label) {
      p <- p + ylab("Log2(FC vs. Inflammatory Mean)") + 
        theme(axis.title.y = element_text(face = "bold", size = 12))
    } else {
      p <- p + theme(axis.title.y = element_blank())
    }
    
    return(p)
  }
  
  # 10. Generate the three separate panels
  p1 <- create_panel(plot_data %>% filter(Category == cat_tf), show_y_label = TRUE)
  p2 <- create_panel(plot_data %>% filter(Category == cat_mat), show_y_label = FALSE)
  p3 <- create_panel(plot_data %>% filter(Category == cat_inf), show_y_label = FALSE)
  
  # 11. Extract the Legend using Cowplot
  # Create a dummy plot just to steal its perfectly formatted legend
  dummy_plot <- p1 + theme(
    legend.position = "top", 
    legend.title = element_blank(), 
    legend.text = element_text(size = 12, face = "bold")
  )
  shared_legend <- get_legend(dummy_plot)
  
  # 12. Create Cowplot Title Block
  title_block <- ggdraw() + 
    draw_label("Single-Cell Mechanistic Trajectories", fontface = 'bold', x = 0.01, hjust = 0, size = 15)
  
  subtitle_block <- ggdraw() + 
    draw_label("Log2(FC) of corrected counts vs. Inflammatory Baseline (Dashed Line = 0)", x = 0.01, hjust = 0, size = 11, color = "grey30")
  
  header <- plot_grid(title_block, subtitle_block, ncol = 1, rel_heights = c(1.2, 1))
  
  # 13. Stitch everything together horizontally (rel_widths enforces the 1:3:4 ratio)
  plot_row <- plot_grid(p1, p2, p3, align = "h", axis = "tb", nrow = 1, rel_widths = c(1.4, 3, 4))
  
  # 14. Final Assembly: Header -> Legend -> Plots
  final_plot <- plot_grid(header, shared_legend, plot_row, ncol = 1, rel_heights = c(0.12, 0.08, 1))
  
  return(final_plot)
}