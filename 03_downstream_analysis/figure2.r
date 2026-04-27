# GE-HTS project pipeline for in situ sequencing data analysis
# author: Nathan Wooseok Lee
# functions for figure 2.


'%notin%' <- Negate('%in%')

# Gene categroy 
anabolic <- c('Acan','Sox9','Col2a1','Matn1','Matn3','Ucma','Ccnd3','Gadd45g','Pth1r','Gm26633','Col27a1')
inflammatory <- c('Mmp3','Mmp13','Il6','Adamts5','Igfbp3','Ccl2','Cxcl5','Cxcl1','Fosl2','Tlr2','Tnfrsf1b', 'Il17b')
housekeeping <- c('Hprt','Actb','Gapdh','B2m','Ubc','Ppia','Rpl23')
gene_list <- c(anabolic, inflammatory, housekeeping)



### Figure 2b (gene expression Heatmap)
# object: mac.sct
fig_2b <- create_heatmap_by_expression_ps(mac.sct, primary_metadata = 'drug_condition', secondary_metadata = 'file') 

### Figure 2c
# ==============================================================================
# Script Name: figure2c_transcriptome_correlation.R
# Description: Generates simplified scatter plots (Figure 2c) comparing average  
#              gene expression between two Seurat objects. Features pure black 
#              formatting and strict scientific notation for statistics.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

#' Create Simplified Transcriptome Correlation Scatter Plots
#'
#' Compares average SCT expression between two Seurat objects across shared conditions.
#' Outputs a clean, black-and-white scatter plot with Pearson correlation statistics.
#'
#' @param obj1 First Seurat object.
#' @param obj2 Second Seurat object.
#' @param obj1_name Label for the x-axis (defaults to the variable name).
#' @param obj2_name Label for the y-axis (defaults to the variable name).
#' @param assay Name of the assay to pull data from (default: 'SCT').
#' @param slot Name of the slot to pull data from (default: 'data').
#' @param group_by Metadata column used to average expression (default: 'drug_condition').
#' @param pdf_file Path to save the output PDF.
#' @param pdf_width Width of the output PDF in inches.
#' @param pdf_height Height of the output PDF in inches.
#'
#' @return The filepath of the saved PDF (invisibly).
#' @export
create_correlation_plots <- function(obj1, 
                                     obj2, 
                                     obj1_name = deparse(substitute(obj1)), 
                                     obj2_name = deparse(substitute(obj2)),
                                     assay = "SCT",
                                     slot = "data", 
                                     group_by = "drug_condition", 
                                     pdf_file = "./figures/figure2/fig2c_comparison_plots.pdf", 
                                     pdf_width = 2.5, 
                                     pdf_height = 2.5) {
  
  # --- Step 1: Calculate Average Expression ---
  message("Calculating average expression for both objects...")
  avg_exp1 <- AverageExpression(obj1, return.seurat = FALSE, group.by = group_by, verbose = FALSE)[[assay]]
  avg_exp2 <- AverageExpression(obj2, return.seurat = FALSE, group.by = group_by, verbose = FALSE)[[assay]]
  
  if (slot == 'data') {
    avg_exp1 <- log1p(avg_exp1)
    avg_exp2 <- log1p(avg_exp2)
  } 
  
  # --- Step 2: Harmonize Conditions and Genes ---
  drug_colnames <- intersect(colnames(avg_exp1), colnames(avg_exp2))
  
  if (length(drug_colnames) == 0) {
    stop("Error: No common '", group_by, "' conditions found between the two Seurat objects.")
  }
  
  data_sample1 <- as.matrix(avg_exp1[, drug_colnames, drop = FALSE])
  data_sample2 <- as.matrix(avg_exp2[, drug_colnames, drop = FALSE])
  
  common_genes <- intersect(rownames(data_sample1), rownames(data_sample2))
  data_sample1 <- data_sample1[common_genes, , drop = FALSE]
  data_sample2 <- data_sample2[common_genes, , drop = FALSE]
  
  # --- Step 3: Generate Plots for Each Condition ---
  plot_list <- list()
  
  for (i in seq_along(drug_colnames)) {
    condition_name <- drug_colnames[i]
    exp1_vals <- as.numeric(data_sample1[, i])
    exp2_vals <- as.numeric(data_sample2[, i])
    
    # Pearson Correlation
    test <- cor.test(exp1_vals, exp2_vals, method = "pearson")
    r_val <- test$estimate
    p_val <- test$p.value
    
    # Format p-value to scientific notation with 2 decimal places (e.g., 1.23e-04)
    # R machine minimum is ~2.22e-16. We handle this explicitly for academic accuracy.
    p_text <- ifelse(p_val < 2.22e-16, "< 2.22e-16", formatC(p_val, format = "e", digits = 2))
    
    plot_data <- data.frame(
      Expression_Sample1 = exp1_vals,
      Expression_Sample2 = exp2_vals
    )
    
    # Plotting
    p <- ggplot(plot_data, aes(x = Expression_Sample1, y = Expression_Sample2)) +
      
      # Solid black dots
      geom_point(color = "black", alpha = 0.5, size = 1.0, stroke = 0) +
      
      # Black dashed trendline
      geom_smooth(method = 'glm', color = "black", linewidth = 0.5, linetype = "dashed", se = FALSE) +
      
      # Pure black text annotation for statistics
      annotate("text",
               x = min(plot_data$Expression_Sample1) + (max(plot_data$Expression_Sample1) * 0.05),
               y = max(plot_data$Expression_Sample2) * 0.95,
               label = paste0("Pearson r: ", sprintf("%.2f", r_val), "\np-value: ", p_text),
               color = "black", size = 2.5, hjust = 0) +
      
      # Theming
      labs(x = paste("Avg. Expression in", obj1_name),
           y = paste("Avg. Expression in", obj2_name),
           title = paste("Comparison:", condition_name)) +
      theme_classic(base_size = 7) +
      
      # Force all text to be black
      theme(text = element_text(family = "Helvetica", color = "black"),
            plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = 8),
            axis.text = element_text(color = "black"),
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black")) +
      expand_limits(x = 0, y = 0)
    
    plot_list[[i]] <- p
  }
  
  # --- Step 4: Save Output ---
  output_dir <- dirname(pdf_file)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  pdf(file = pdf_file, width = pdf_width, height = pdf_height)
  invisible(lapply(plot_list, print))
  dev.off()
  
  message("Saved ", length(plot_list), " correlation plots to: ", pdf_file)
  return(invisible(pdf_file))
}


# GEHTS-chip data
mac.sct.ctrl <-subset(mac.sct, drug_condition == 'control')

# raw rna-seq data of MAC (in house) 
gene_list <- c(anabolic, inflammatory, housekeeping)

# prepare RNA-seq data to compare
rna.sct <- SCTransform(rna.pub, vst.flavor='v2',verbose=FALSE, return.only.var.genes = FALSE)
cells_of_interest <- c('control','control.1')
rna.sct <- subset(rna.sct, cells = cells_of_interest)
rna.sct@meta.data$drug_condition <- c('control', 'control')
rna.sct<- subset(rna.sct, features = gene_list)

# MAC uninjured RNA-seq data from Sebastian 2021 paper.
macSebastian$drug_condition <- 'control'

# Example Usage
create_correlation_plots(seurat_object1 = mac.sct.ctrl, seurat_object2 = rna.sct, pdf_file = "./figures/figure2/fig2h_macVSrna_correlation_plots.pdf")
create_correlation_plots(seurat_object1 = mac.sct.ctrl, seurat_object2 = macSebastian, pdf_file = "./figures/figure2/fig2h_macVSsebastian_correlation_plots.pdf")
create_correlation_plots(seurat_object1 = rna.sct, seurat_object2 = macSebastian, pdf_file = "./figures/figure2/fig2h_rnaVSsebastian_correlation_plots.pdf")


### Figure 2d (UMAP plot)
# ==============================================================================
# Script Name: figure2d_tsne_plot.R
# Description: Generates a t-SNE plot highlighting specific drug conditions.
#              Defaults 'control' to NPG Blue and 'inflammatory' to NPG Red, 
#              while dynamically accommodating additional conditions if present.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggsci)
})

#' Generate and Save t-SNE Plot
#'
#' Extracts t-SNE embeddings from a Seurat object and plots them. Applies a 
#' specific NPG color mapping for standard baseline conditions, and auto-fills 
#' colors for any additional conditions found in the metadata.
#'
#' @param seurat_object Processed Seurat object with t-SNE reduction.
#' @param group_by Metadata column to color the points by (default: "drug_condition").
#' @param pdf_file Path to save the output PDF.
#' @param pdf_width Numeric width of the output PDF in inches.
#' @param pdf_height Numeric height of the output PDF in inches.
#' @param font_size Base font size for the plot theme.
#'
#' @return The filepath of the saved PDF (invisibly).
#' @export
generate_tsne_pdf <- function(seurat_object, 
                              group_by = "drug_condition",
                              pdf_file = "./figures/figure2/fig2d_tsne_plot.pdf", 
                              pdf_width = 7, 
                              pdf_height = 5.5, 
                              font_size = 9) {
  
  # --- Step 1: Extract t-SNE coordinates and metadata safely ---
  if (!"tsne" %in% names(seurat_object@reductions)) {
    stop("Error: t-SNE reduction not found in the Seurat object. Run RunTSNE() first.")
  }
  
  tsne_coordinates <- Embeddings(seurat_object, reduction = "tsne")
  tsne_df <- as.data.frame(tsne_coordinates)
  colnames(tsne_df) <- c("tSNE_1", "tSNE_2")
  
  # Safely bind metadata using exact cell names from the embedding
  cells <- rownames(tsne_df)
  if (!group_by %in% colnames(seurat_object@meta.data)) {
    stop(paste("Error: Column", group_by, "not found in metadata."))
  }
  tsne_df[[group_by]] <- as.character(seurat_object@meta.data[cells, group_by])
  
  # Optional: Keep the file mapping commented out but ready for future use
  # tsne_df$file <- seurat_object@meta.data[cells, "file"]
  
  # --- Step 2: Dynamic NPG Color Mapping ---
  npg_palette <- pal_npg("nrc")(10)
  
  # Define the locked conditions
  # [1] is Red, [2] is Blue in standard NPG
  locked_colors <- c("control" = npg_palette[2], "inflammatory" = npg_palette[1])
  
  # Identify any conditions in the data that aren't 'control' or 'inflammatory'
  present_conditions <- unique(tsne_df[[group_by]])
  other_conditions <- setdiff(present_conditions, names(locked_colors))
  
  # If other conditions exist, assign them the remaining NPG colors
  final_color_map <- locked_colors
  if (length(other_conditions) > 0) {
    available_colors <- setdiff(npg_palette, locked_colors)
    
    # If there are more new conditions than remaining colors, interpolate to fit
    if (length(other_conditions) > length(available_colors)) {
      available_colors <- colorRampPalette(available_colors)(length(other_conditions))
    }
    
    names(available_colors) <- other_conditions
    final_color_map <- c(final_color_map, available_colors)
  }
  
  # --- Step 3: Create the t-SNE plot ---
  p <- ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2, color = .data[[group_by]])) +
    geom_point(alpha = 0.8, size = 1) +
    scale_color_manual(values = final_color_map) +
    theme_classic(base_size = font_size) +
    labs(title = "t-SNE Plot", 
         x = "t-SNE 1", 
         y = "t-SNE 2",
         color = "Condition") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "right")
  
  # --- Step 4: Save the plot ---
  output_dir <- dirname(pdf_file)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  pdf(pdf_file, width = pdf_width, height = pdf_height)
  print(p)
  dev.off()
  
  message("Figure 2d saved to: ", pdf_file)
  return(invisible(pdf_file))
}

# Example usage
DefaultAssay(mac.sct) <- "SCT"
set.seed(123)
mac.sct <- RunPCA(mac.sct, assay = "SCT", npcs = 10, verbose = FALSE)
# Run t-SNE using the same PCs (dims 1:10)
mac.sct <- RunTSNE(mac.sct, dims = 1:10, verbose = FALSE)

generate_tsne_pdf(seurat_object = mac.sct, 
                  pdf_file = "./figures/figure2/figdf_tsne_plot.pdf", 
                  pdf_width = 7, 
                  pdf_height = 5.5)

### Figure 2e (volcano plot)
# ==============================================================================
# Script Name: figure2e_volcano_plot.R
# Description: Generates a differential expression volcano plot.
#              Highlights downregulated anabolic genes in NPG Green and 
#              upregulated catabolic genes in JCO Yellow.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(ggsci)
})

#' Generate Annotated Volcano Plot for Differential Expression
#'
#' Runs FindMarkers on a Seurat object and plots the results. Specifically 
#' highlights genes from provided anabolic and catabolic lists that meet 
#' statistical and fold-change thresholds.
#'
#' @param object A processed Seurat object.
#' @param ident.1 The primary identity class to test (e.g., 'inflammatory').
#' @param group.by The metadata column containing the identities.
#' @param assay The assay to use for differential expression (default: 'SCT').
#' @param logfc.threshold The minimum logFC threshold for FindMarkers.
#' @param genes_anabolic Character vector of anabolic genes to highlight (green).
#' @param genes_catabolic Character vector of catabolic genes to highlight (yellow).
#' @param pdf_file Path to save the output PDF.
#' @param pdf_width Width of the output PDF in inches.
#' @param pdf_height Height of the output PDF in inches.
#'
#' @return The filepath of the saved PDF (invisibly).
#' @export
generate_volcano_plot <- function(object, 
                                  ident.1, 
                                  group.by, 
                                  assay = 'SCT', 
                                  logfc.threshold = 0.1,
                                  genes_anabolic,
                                  genes_catabolic,
                                  pdf_file = "./figures/figure2/fig2e_volcano_plot.pdf", 
                                  pdf_width = 8, 
                                  pdf_height = 6) {
  
  # --- Step 1: Run Differential Expression ---
  message("Preparing SCT assay for FindMarkers...")
  object_prep <- PrepSCTFindMarkers(object = object)
  
  message("Running FindMarkers...")
  de_results <- FindMarkers(object_prep, 
                            ident.1 = ident.1, 
                            group.by = group.by, 
                            assay = assay, 
                            logfc.threshold = logfc.threshold)

  # Prepare metrics for plotting
  de_results$logP_adj <- -log10(de_results$p_val_adj)
  de_results$gene <- rownames(de_results)

  # --- Step 2: Classify Genes for Highlighting ---
  # Only highlight if they pass thresholds AND belong to the specific gene lists
  de_results <- de_results %>%
    mutate(gene_type = case_when(
      avg_log2FC >= 0.5 & p_val_adj <= 0.05 & gene %in% genes_catabolic ~ "Catabolic_Up",
      avg_log2FC <= -0.5 & p_val_adj <= 0.05 & gene %in% genes_anabolic ~ "Anabolic_Down",
      TRUE ~ "ns"
    ))

  # --- Step 3: Customize Aesthetics (ggsci palettes) ---
  # Extract Green from NPG and Yellow/Gold from JCO
  npg_green <- pal_npg("nrc")(10)[3]  # #00A087FF (Teal/Green)
  jco_yellow <- pal_jco("default")(10)[2] # #EFC000FF (Yellow/Gold)
  
  cols   <- c("Catabolic_Up" = jco_yellow, "Anabolic_Down" = npg_green, "ns" = "grey80")
  sizes  <- c("Catabolic_Up" = 2, "Anabolic_Down" = 2, "ns" = 1)
  alphas <- c("Catabolic_Up" = 1, "Anabolic_Down" = 1, "ns" = 0.3)

  # --- Step 4: Create the Volcano Plot ---
  # Sort so that highlighted genes are plotted on top of the grey background dots
  de_results <- de_results %>% arrange(factor(gene_type, levels = c("ns", "Anabolic_Down", "Catabolic_Up")))

  p <- ggplot(de_results, aes(x = avg_log2FC, y = logP_adj, label = gene, 
                              fill = gene_type, size = gene_type, alpha = gene_type)) +
    geom_point(shape = 21, color = 'black', stroke = 0.2) +
    
    # Threshold lines
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    
    # Labels (only for the significantly categorized genes)
    geom_text_repel(data = filter(de_results, gene_type != 'ns'), 
                    aes(label = gene), 
                    size = 3,
                    color = "black",
                    box.padding = 0.5,
                    max.overlaps = Inf,
                    show.legend = FALSE) +
    
    # Theming
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),    
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(title = paste("Volcano Plot:", ident.1), 
         x = "Log2 Fold Change", 
         y = "-Log10(Adjusted P-value)",
         fill = "Gene Category") +
    scale_fill_manual(values = cols) +
    scale_size_manual(values = sizes, guide = "none") +
    scale_alpha_manual(values = alphas, guide = "none")
    
    # Note: I removed the hardcoded scale_x_continuous limits to prevent 
    # ggplot from clipping significant data points that fall outside c(-2, 3)

  # --- Step 5: Save Output ---
  output_dir <- dirname(pdf_file)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  pdf(pdf_file, width = pdf_width, height = pdf_height)
  print(p)
  dev.off()

  message("Volcano plot saved to: ", pdf_file)
  return(invisible(pdf_file))
}

# Example usage 
generate_volcano_plot(
  object          = mac.sct, 
  ident.1         = 'inflammatory', 
  group.by        = 'drug_condition', 
  assay           = 'SCT', 
  logfc.threshold = 0.1, 
  genes_anabolic  = anabolic,
  genes_catabolic = inflammatory,
  pdf_file        = "./figures/figure2/fig2e_volcano_plot.pdf", 
  pdf_width       = 4, # Adjusted slightly to accommodate the legend comfortably
  pdf_height      = 2.5
)


### fig 2f
# ==============================================================================
# Script Name: boxplot_expression.R
# Description: Generates publication-ready boxplots for gene expression.
#              Automatically colors 'control' with NPG blue and 'inflammatory'
#              with NPG red. Additional conditions pull from the NPG palette.
#              Conditionally excludes extra labels if left as default.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(cowplot)
  library(ggsignif)
  library(ggsci)
})

#' Draw Boxplots for Genes with a Common Legend
#'
#' @param seurat_object A processed Seurat object.
#' @param genes A character vector of gene names to plot.
#' @param datatype The assay slot to pull data from (default: 'counts').
#' @param control_label The metadata label for the control group.
#' @param treatment_label The metadata label for the treatment group.
#' @param extra_label1 Optional metadata label for a 3rd group.
#' @param extra_label2 Optional metadata label for a 4th group.
#' @param pdf_file Path to save the output PDF.
#' @param pdf_width Numeric width of the output PDF.
#' @param pdf_height Numeric height of the output PDF.
#' @param base_font_size Base font size for ggplot theme.
#' @param line_thickness Line thickness for axes and boxes.
#' @param legend_gap Relative height ratio for the legend vs. plots.
#' @param y_title Label for the shared Y-axis.
#'
#' @return The final ggdraw object (invisibly).
#' @export
draw_boxplots_for_genes_with_common_legend <- function(seurat_object,
                                                       genes,
                                                       datatype = 'counts',
                                                       control_label = "control",
                                                       treatment_label = "inflammatory",
                                                       extra_label1 = "label3",
                                                       extra_label2 = "label4",
                                                       pdf_file = "./figures/figure2/boxplots.pdf",
                                                       pdf_width = 6,
                                                       pdf_height = 6,
                                                       base_font_size = 5,
                                                       line_thickness = 0.5,
                                                       legend_gap = 0.05,
                                                       y_title = "Normalized gene expression (a.u)") {
  
  meta <- seurat_object@meta.data
  
  # --- Step 1: Handle Valid Conditions ---
  # Only include extra labels if they differ from the default placeholder values
  valid_conditions <- c(control_label, treatment_label)
  if (extra_label1 != "label3") valid_conditions <- c(valid_conditions, extra_label1)
  if (extra_label2 != "label4") valid_conditions <- c(valid_conditions, extra_label2)
  
  cells_to_use <- rownames(meta)[meta$drug_condition %in% valid_conditions]
  
  if (length(cells_to_use) == 0) {
    stop("Error: None of the specified conditions were found in 'drug_condition'.")
  }
  
  meta_subset <- meta[cells_to_use, ] %>%
    mutate(ConditionLabel = case_when(
      drug_condition == control_label ~ "Control",
      drug_condition == treatment_label ~ "IL-1B",
      TRUE ~ drug_condition
    ))
  
  # Map actual unique conditions found
  actual_conditions <- unique(meta_subset$ConditionLabel)
  # Ensure specific ordering for factors: Control, IL-1B, then others
  factor_levels <- c("Control", "IL-1B", setdiff(actual_conditions, c("Control", "IL-1B")))
  meta_subset$ConditionLabel <- factor(meta_subset$ConditionLabel, levels = factor_levels)

  # --- Step 2: Setup NPG Color Palette ---
  npg_palette <- pal_npg("nrc")(10)
  # Standard: [2] is Blue (Control), [1] is Red (Inflammatory)
  color_map <- c("Control" = npg_palette[2], "IL-1B" = npg_palette[1])
  
  # Map any extra conditions to remaining NPG colors
  extra_levels <- setdiff(factor_levels, c("Control", "IL-1B"))
  if (length(extra_levels) > 0) {
    available_colors <- setdiff(npg_palette, c(npg_palette[1], npg_palette[2]))
    extra_colors <- available_colors[1:length(extra_levels)]
    names(extra_colors) <- extra_levels
    color_map <- c(color_map, extra_colors)
  }

  # --- Step 3: Generate Individual Boxplots ---
  plot_list <- list()
  
  for (gene in genes) {
    # Extract gene expression safely
    expr <- as.vector(GetAssayData(seurat_object, slot = datatype)[gene, cells_to_use, drop = TRUE])
    df <- data.frame(Expression = expr, Condition = meta_subset$ConditionLabel)
    
    comparisons_list <- combn(as.character(factor_levels), 2, simplify = FALSE)
    ymax <- boxplot.stats(df$Expression)$stats[5]
    if (is.na(ymax) || ymax == 0) ymax <- max(df$Expression, na.rm = TRUE)
    
    p_gene <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
      geom_boxplot(outlier.shape = NA, linewidth = line_thickness) +
      geom_signif(comparisons = comparisons_list,
                  map_signif_level = TRUE,
                  y_position = ymax * (1 + 0.15 * seq_along(comparisons_list)),
                  textsize = base_font_size - 1,
                  tip_length = 0.01,
                  step_increase = 0.05,
                  vjust = 0.5) +
      scale_fill_manual(values = color_map) +
      labs(title = gene) +
      theme_classic() +
      theme(text = element_text(size = base_font_size),
            axis.line = element_line(linewidth = line_thickness),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = base_font_size),
            axis.title = element_blank(),
            plot.title = element_text(face = "bold.italic", hjust = 0.5, size = base_font_size),
            legend.position = "none")
    
    plot_list[[gene]] <- p_gene
  }
  
  # --- Step 4: Generate Common Legend ---
  dummy_df <- data.frame(Condition = factor(factor_levels, levels = factor_levels),
                         Expression = rep(0.5, length(factor_levels)))
  
  dummy_plot <- ggplot(dummy_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot() +
    scale_fill_manual(values = color_map) +
    labs(fill = 'Condition') +
    theme_classic() +
    theme(legend.position = "top",
          legend.text = element_text(size = base_font_size),
          legend.title = element_text(size = base_font_size, face = "bold"),
          legend.margin = margin(0,0,0,0),
          legend.key.size = unit(0.5, "cm"))
  
  common_legend <- get_legend(dummy_plot)
  
  # --- Step 5: Arrange and Save ---
  n_genes <- length(genes)
  ncol <- ceiling(n_genes / 2)
  grid_plots <- arrangeGrob(grobs = plot_list, ncol = ncol, nrow = 2)
  
  arranged_plots <- grid.arrange(common_legend, grid_plots,
                                 ncol = 1,
                                 heights = c(legend_gap, 1))
  
  final_plot <- ggdraw() +
    draw_plot(arranged_plots, x = 0.05, y = 0, width = 0.95, height = 1) +
    draw_label(y_title, x = 0.02, y = 0.5,
               angle = 90, vjust = 0.5, size = base_font_size, fontface = "bold")
  
  output_dir <- dirname(pdf_file)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (!is.null(pdf_file)) {
    pdf(file = pdf_file, width = pdf_width, height = pdf_height)
    print(final_plot)
    dev.off()
    message("Boxplots successfully saved to: ", pdf_file)
  }
  
  return(invisible(final_plot))
}

genes_to_plot <- c("Sox9", "Acan", 
                    "Mmp3", "Mmp13")

draw_boxplots_for_genes_with_common_legend(mac.sct, genes = genes_to_plot,
                        pdf_file = "./figures/figure2/fig2f_boxplots.pdf",
                        pdf_width = 2, pdf_height = 3, base_font_size = 5, legend_gap = 0.03)


### fig 2g (Primary cell VS Cell line)
library(dplyr)
library(ggplot2)

#' Plot Transcriptomic Divergence (Primary vs. ATDC5)
#'
#' Generates an academic-grade lollipop chart comparing relative gene enrichment.
#' Uses a neutral, NPG-inspired color palette for publication.
#'
#' @param norm_primary Numeric matrix/dataframe of CLR-normalized expression for Primary cells.
#' @param norm_atdc5 Numeric matrix/dataframe of CLR-normalized expression for ATDC5 cells.
#' @param target_genes Character vector of representative markers to plot.
#' @return A ggplot2 object.
#' @export
plot_gene_divergence <- function(norm_primary, norm_atdc5, 
                                 target_genes = c("Col2a1", "Acan", "Matn3", "Gadd45g", 
                                                  "Ccnd3", "Il6", "Ccl2", "Cxcl1", "Tlr2", "Rpl23")) {
  
  # ==============================================================================
  # 1. DATA PROCESSING
  # ==============================================================================
  # Calculate difference in means (Primary - ATDC5)
  diff_values <- rowMeans(norm_primary, na.rm = TRUE) - rowMeans(norm_atdc5, na.rm = TRUE)
  
  # Extract and filter target genes
  subset_diff <- diff_values[target_genes]
  subset_diff <- na.omit(subset_diff)
  
  # Create plotting dataframe
  plot_data <- data.frame(
    Gene = names(subset_diff),
    Difference = as.numeric(subset_diff)
  )
  
  # Assign biological categories
  plot_data <- plot_data %>%
    mutate(Category = case_when(
      Gene %in% c("Col2a1", "Acan", "Matn3") ~ "Matrix Identity",
      Gene == "Gadd45g" ~ "Quiescence",
      Gene %in% c("Il6", "Ccl2", "Cxcl1", "Tlr2") ~ "Inflammation",
      Gene == "Ccnd3" ~ "Cell Cycle",
      Gene == "Rpl23" ~ "Metabolism",
      TRUE ~ "Other"
    ))
  
  # Sort genes by difference for ordered plotting
  plot_data <- plot_data %>% arrange(Difference)
  plot_data$Gene <- factor(plot_data$Gene, levels = plot_data$Gene)
  
  # ==============================================================================
  # 2. AESTHETICS & PALETTE (Neutral NPG-Inspired)
  # ==============================================================================
  # Color palette carefully chosen for scientific neutrality and colorblind safety
  npg_neutral_colors <- c(
    "Matrix Identity" = "#00A087", # NPG Teal/Green (Good/Primary)
    "Quiescence"      = "#8491B4", # NPG Grey-Blue (Resting/Primary)
    "Inflammation"    = "#F39B7F", # NPG Muted Orange (Bad/ATDC5)
    "Cell Cycle"      = "#C87A9E", # Soft Dusty Pink (Dividing/ATDC5)
    "Metabolism"      = "#919C4A"  # Soft Olive Green (Active/ATDC5)
  )
  
  # Dynamically calculate annotation positions based on data bounds
  y_max <- max(plot_data$Difference)
  y_min <- min(plot_data$Difference)
  n_genes <- nrow(plot_data)
  
  # ==============================================================================
  # 3. GENERATE FIGURE
  # ==============================================================================
  p <- ggplot(plot_data, aes(x = Gene, y = Difference, color = Category)) +
    
    # Structural elements
    geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
    geom_segment(aes(x = Gene, xend = Gene, y = 0, yend = Difference), linewidth = 1.2) +
    geom_point(size = 4.5) +
    
    # Dynamic Annotations (placed near the top and bottom of the Y-axis)
    annotate("text", x = n_genes - 1, y = y_max * 0.85, 
             label = "Enriched in\nPRIMARY", color = "grey30", 
             fontface = "bold", hjust = 1, size = 4) +
    annotate("text", x = 2, y = y_min * 0.85, 
             label = "Enriched in\nATDC5", color = "grey30", 
             fontface = "bold", hjust = 0, size = 4) +
    
    # Formatting
    coord_flip() +
    scale_color_manual(values = npg_neutral_colors) +
    theme_classic() + # Better for journals than theme_minimal()
    labs(
      title = "Transcriptomic Divergence: Primary vs. ATDC5",
      subtitle = "Relative Enrichment of Representative Markers",
      y = "Relative Enrichment (Primary - ATDC5, CLR Normalized)",
      x = NULL
    ) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      axis.text.y = element_text(size = 11, face = "bold.italic", color = "black"), # Genes often italicized
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey20", margin = margin(b = 15))
    )
  
  return(p)
}

# ==============================================================================
# EXAMPLE USAGE:
# my_plot <- plot_gene_divergence(norm_primary, norm_atdc5)
# print(my_plot)
# ggsave("Transcriptomic_Divergence.pdf", plot = my_plot, width = 8, height = 6, dpi = 300)
# ==============================================================================