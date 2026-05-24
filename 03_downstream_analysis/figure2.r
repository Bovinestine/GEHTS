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
ggsave("./Figure2b_heatmap_zscore.pdf", plot = fig_2b, width =6, height = 5, units = "in", dpi = 300)

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


### fig 2d (Primary cell VS Cell line)
# ==============================================================================
# Description: Generates a lollipop plot highlighting difference between two samples.
# ==============================================================================
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


### Figure 2e (UMAP plot)
# ==============================================================================
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

### Figure 2f (volcano plot)
# ==============================================================================
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
      avg_log2FC >= 0.5 & p_val_adj <= 0.05 & gene %in% genes_catabolic ~ "Catabolic",
      avg_log2FC <= -0.5 & p_val_adj <= 0.05 & gene %in% genes_anabolic ~ "Anabolic",
      TRUE ~ "ns"
    ))

  # --- Step 3: Customize Aesthetics (ggsci palettes) ---
  # Extract Green from NPG and Yellow/Gold from JCO
  npg_green <- pal_npg("nrc")(10)[3]  # #00A087FF (Teal/Green)
  jco_yellow <- pal_jco("default")(10)[2] # #EFC000FF (Yellow/Gold)
  
  cols   <- c("Catabolic" = jco_yellow, "Anabolic" = npg_green, "ns" = "grey80")
  sizes  <- c("Catabolic" = 2, "Anabolic" = 2, "ns" = 1)
  alphas <- c("Catabolic" = 1, "Anabolic" = 1, "ns" = 0.3)

  # --- Step 4: Create the Volcano Plot ---
  # Sort so that highlighted genes are plotted on top of the grey background dots
  de_results <- de_results %>% arrange(factor(gene_type, levels = c("ns", "Anabolic", "Catabolic")))

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
    scale_alpha_manual(values = alphas, guide = "none") +
    scale_x_continuous(breaks = seq(-10, 10, by = 0.5))
    
    # Note: I removed the hardcoded scale_x_continuous limits to prevent 
    # ggplot from clipping significant data points that fall outside c(-2, 3)


  return(p)
}

# Example usage 
sfig_volcano <- generate_volcano_plot(
  object          = mac.sct, 
  ident.1         = 'inflammatory', 
  group.by        = 'drug_condition', 
  assay           = 'SCT', 
  logfc.threshold = 0.1, 
  genes_anabolic  = anabolic,
  genes_catabolic = inflammatory
)
ggsave('./sfig_volcanoplot.pdf', plot = sfig_volcano, width = 6, height = 3, units = "in", dpi = 300)

