# GE-HTS project pipeline for in situ sequencing data analysis
# author: Nathan Wooseok Lee
# functions for figure 2.

library(Seurat)
library(ggplot2) # the latest version 3.5.0 released at 24.02.23 is desirable.
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(ggpubr) # for stat_compare_means
library(scales)  # For gradient color palette
library(dplyr)
library(patchwork)
library(cowplot)
# For reproducibility
set.seed(123)
'%notin%' <- Negate('%in%')

# Gene categroy 
anabolic <- c('Acan','Sox9','Col2a1','Matn1','Matn3','Ucma','Ccnd3','Gadd45g','Pth1r','Gm26633','Col27a1')
inflammatory <- c('Mmp3','Mmp13','Il6','Adamts5','Igfbp3','Ccl2','Cxcl5','Cxcl1','Fosl2','Tlr2','Tnfrsf1b', 'Il17b')
housekeeping <- c('Hprt','Actb','Gapdh','B2m','Ubc','Ppia','Rpl23')
gene_list <- c(anabolic, inflammatory, housekeeping)

# figure 2.c 
sum_genes_and_plot(mac.sct, pdf_width = 3, pdf_height = 2.5, base_font_size = 7, pdf_file='./figures/figure2/figure2c_barplot.pdf')#251230

# Function to calculate sum of gene counts per category for each cell and save plot as PDF
sum_genes_and_plot <- function(seurat_object, 
                               anabolic_genes = anabolic,
                               inflammatory_genes = inflammatory, 
                               control_label = "control",
                               treatment_label = "inflammatory",
                               pdf_file = "./figures/figure2/sum_genes_and_barplot.pdf",
                               pdf_width = 4, 
                               pdf_height = 6,
                               base_font_size = 10,
                               datatype = 'counts',
                               common_x_label = "",
                               common_y_label = "Sum of normalized gene counts",
                               title1 = "Anabolic Gene",
                               title2 = "Catabolic Gene") {
  
  # Add computed sums to the Seurat object's metadata
  seurat_object[["anabolic_sum"]] <- Matrix::colSums(GetAssayData(seurat_object, slot = datatype)[anabolic_genes, , drop = FALSE])
  seurat_object[["inflammatory_sum"]] <- Matrix::colSums(GetAssayData(seurat_object, slot = datatype)[inflammatory_genes, , drop = FALSE])
  
  # Subset cells by condition (using the given metadata field)
  control_cells <- subset(seurat_object, subset = drug_condition == control_label)
  treatment_cells <- subset(seurat_object, subset = drug_condition == treatment_label)
  
  # Change the treatment label if any
  if (treatment_label == 'inflammatory'){
    treatment_label <- "IL-1β"
  }

  # Prepare data frames for each gene group:
  anabolic_df <- data.frame(
    condition = rep(c(control_label, treatment_label),
                    times = c(nrow(control_cells@meta.data), nrow(treatment_cells@meta.data))),
    sum = c(control_cells@meta.data$anabolic_sum, treatment_cells@meta.data$anabolic_sum)
  )
  
  inflammatory_df <- data.frame(
    condition = rep(c(control_label, treatment_label),
                    times = c(nrow(control_cells@meta.data), nrow(treatment_cells@meta.data))),
    sum = c(control_cells@meta.data$inflammatory_sum, treatment_cells@meta.data$inflammatory_sum)
  )
  
  # Determine maximum y values for proper label placement
  # anabolic_ymax <- max(anabolic_df$sum, na.rm = TRUE)
  # inflammatory_ymax <- max(inflammatory_df$sum, na.rm = TRUE)
  anabolic_ymax <- boxplot.stats(anabolic_df$sum)$stats[5]
  inflammatory_ymax <- boxplot.stats(inflammatory_df$sum)$stats[5]

  
  # Create the anabolic gene plot (remove axis titles so common labels can be added later)
  p_anabolic <- ggplot(anabolic_df, aes(x = condition, y = sum, fill = condition)) +
    geom_boxplot(outliers = FALSE) +
    #stat_compare_means(aes(group = condition), label = "p.format", 
    #                   label.y = anabolic_ymax * 1.05, size = base_font_size * 0.6) +
    stat_compare_means(comparisons=list(c(control_label,treatment_label)), label = 'p.signif', 
                        label.y = anabolic_ymax * 1.05, size = base_font_size) +
    labs(title = title1) +
    theme_classic() +
    theme(text = element_text(size = base_font_size),
          # legend.position = "top",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(face = "plain", hjust = 0.5, size = base_font_size + 2))
  
  # Create the inflammatory gene plot (also remove individual axis titles)
  p_inflammatory <- ggplot(inflammatory_df, aes(x = condition, y = sum, fill = condition)) +
    geom_boxplot(outliers = FALSE) +
    #stat_compare_means(aes(group = condition), label = "p.format", 
    #                   label.y = inflammatory_ymax * 1.05, size = base_font_size * 0.6) +
    stat_compare_means(comparisons=list(c(control_label, treatment_label)), label = 'p.signif', 
                        label.y = inflammatory_ymax * 1.05, size = base_font_size) +
    labs(title = title2) +
    theme_classic() +
    theme(text = element_text(size = base_font_size),
          # legend.position = "top",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(face = "plain", hjust = 0.5, size = base_font_size + 2))
  
  p_anabolic$layers[[2]]$aes_params$textsize <- base_font_size
  p_inflammatory$layers[[2]]$aes_params$textsize <- base_font_size

  # Combine the two plots side by side with a single (collected) legend at the top.
  #combined_patch <- (p_anabolic + p_inflammatory + plot_layout(guides = "collect")) &
  #  theme(legend.position = "top", legend.key.height = unit(0.3, "cm"),
  #    legend.margin = margin(t = 2, b = 2))
  combined_patch <- (p_anabolic + p_inflammatory + plot_layout(guides = "collect")) &
    theme(legend.position = "none")
  
  
  
  # Now, add common x and y axis labels using cowplot. Note that you may need to adjust
  # the x and y positions (here 0.5 and 0 for x label; 0 and 0.5 for y label) and vertical/horizontal justification.
  final_plot <- ggdraw(combined_patch) +
    draw_label(common_x_label, x = 0.5, y = 0, vjust = -1.2, size = base_font_size + 2) +
    draw_label(common_y_label, x = 0, y = 0.5, angle = 90, vjust = 1.2, size = base_font_size + 2)
  
  # Save the final figure to a PDF file
  pdf(file = pdf_file, width = pdf_width, height = pdf_height)
  print(final_plot)
  dev.off()
  
  invisible(final_plot)
}


### Figure 2d (grid of gene expression box plots)
genes_to_plot <- c("Sox9", "Acan", 
                    "Mmp3", "Mmp13")

draw_boxplots_for_genes_with_common_legend(mac.sct, genes = genes_to_plot,
                        pdf_file = "./figures/figure2/fig2d_boxplots.pdf",
                        pdf_width = 2, pdf_height = 3, base_font_size = 5, legend_gap = 0.03)

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
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(cowplot)
  library(ggsignif)  # for p-value annotation
  
  meta <- seurat_object@meta.data
  valid_conditions <- c(control_label, treatment_label, extra_label1, extra_label2)
  cells_to_use <- rownames(meta)[ meta$drug_condition %in% valid_conditions ]
  
  meta_subset <- meta[cells_to_use, ] %>%
    mutate(ConditionLabel = case_when(
      drug_condition == control_label ~ "Control",
      drug_condition == treatment_label ~ "IL-1B",
      drug_condition == extra_label1 ~ extra_label1,
      drug_condition == extra_label2 ~ extra_label2,
      TRUE ~ drug_condition
    ))
  
  plot_list <- list()
  
  for (gene in genes) {
    expr <- as.vector(GetAssayData(seurat_object, slot = datatype)[gene, cells_to_use, drop = TRUE])
    df <- data.frame(Expression = expr,
                     Condition = meta_subset$ConditionLabel)
    
    # Get all pairwise comparisons (for 4 groups => 6 pairs)
    groups <- unique(df$Condition)
    comparisons_list <- combn(groups, 2, simplify = FALSE)
    
    # Estimate y-position for significance annotation
    ymax <- boxplot.stats(df$Expression)$stats[5]
    
    p_gene <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
      geom_boxplot(outlier.shape = NA) +
      geom_signif(comparisons = comparisons_list,
                  map_signif_level = TRUE,
                  y_position = ymax * (1 + 0.15 * seq_along(comparisons_list)),
                  textsize = base_font_size - 1,
                  tip_length = 0.01,
                  step_increase = 0.05) +
      labs(title = gene) +
      theme_classic() +
      theme(text = element_text(size = base_font_size),
            axis.line = element_line(linewidth = line_thickness),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = base_font_size),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(face = "plain", hjust = 0.5, size = base_font_size),
            legend.position = "none")
    
    plot_list[[gene]] <- p_gene
  }
  
  dummy_df <- data.frame(Condition = factor(c("Control", "IL-1B", extra_label1, extra_label2),
                                            levels = c("Control", "IL-1B", extra_label1, extra_label2)),
                         Expression = rep(0.5, 4))
  
  dummy_plot <- ggplot(dummy_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot() +
    labs(fill = 'Condition') +
    theme_classic() +
    theme(legend.position = "top",
          legend.text = element_text(size = base_font_size),
          legend.title = element_text(size = base_font_size),
          legend.margin = margin(0,0,0,0))
  
  get_legend <- function(myplot) {
    tmp <- ggplot_gtable(ggplot_build(myplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  common_legend <- get_legend(dummy_plot)
  
  n_genes <- length(genes)
  ncol <- ceiling(n_genes / 2)
  grid_plots <- arrangeGrob(grobs = plot_list, ncol = ncol, nrow = 2)
  
  arranged_plots <- grid.arrange(common_legend, grid_plots,
                                 ncol = 1,
                                 heights = c(legend_gap, 1))
  
  final_plot <- ggdraw() +
    draw_plot(arranged_plots, x = 0.05, y = 0, width = 0.95, height = 1) +
    draw_label(y_title, x = 0.02, y = 0.5,
               angle = 90, vjust = 0.5, size = base_font_size)
  
  if (!is.null(pdf_file)) {
    pdf(file = pdf_file, width = pdf_width, height = pdf_height)
    print(final_plot)
    dev.off()
  }
  
  return(final_plot)
}

### Figure 2e (gene expression Heatmap)
# object: mac.sct
create_heatmap_by_expression_ps(mac.sct, upper_threshold = 0, bottom_threshold = 1, 
                                assay='SCT', slot='counts', primary_metadata = 'drug_condition', secondary_metadata = 'file',
                                pdf_file = './figures/figure2/GEheatmap_251230.pdf', pdf_width = 3, pdf_height = 4.5) 

create_heatmap_by_expression_ps <- function(seurat_object, upper_threshold = 0, bottom_threshold = 1, 
                                           assay = 'SCT', slot = 'counts', 
                                           primary_metadata = 'dose', secondary_metadata = 'drug_name', 
                                           log1p = TRUE, pdf_file = "heatmap_expression.pdf", pdf_width = 5, pdf_height = 5) {
  # Get expression data
  expression_data <- GetAssayData(seurat_object, assay = assay, slot = slot)

  # Ensure there are no non-finite values
  expression_data <- expression_data[, colSums(is.finite(expression_data)) > 0]

  average_expression <- rowMeans(expression_data, na.rm = TRUE)

  # Identify highly expressed and barely expressed genes
  selected_genes <- names(average_expression[average_expression >= quantile(average_expression, 1 - bottom_threshold) & 
                                               average_expression <= quantile(average_expression, 1 - upper_threshold)])
  title <- ""

  # Split selected genes into categories
  selected_anabolic_genes <- intersect(selected_genes, anabolic)
  selected_inflammatory_genes <- intersect(selected_genes, inflammatory)
  selected_housekeeping_genes <- intersect(selected_genes, housekeeping)

  # Combine and sort genes within each category
  sorted_selected_genes <- c(
    selected_anabolic_genes[order(average_expression[selected_anabolic_genes], decreasing = TRUE)],
    selected_inflammatory_genes[order(average_expression[selected_inflammatory_genes], decreasing = TRUE)],
    selected_housekeeping_genes[order(average_expression[selected_housekeeping_genes], decreasing = TRUE)]
  )

  # Prepare cell scores with primary and secondary metadata
  cell_scores <- data.frame(
    cell_name = rownames(seurat_object@meta.data),
    primary = seurat_object@meta.data[[primary_metadata]],
    secondary = seurat_object@meta.data[[secondary_metadata]],
    stringsAsFactors = FALSE
  )

  # Group samples by primary and then secondary metadata
  cell_scores <- cell_scores %>%
    arrange(primary, secondary)

  # Function to create heatmap
  create_heatmap <- function(genes, data, cell_scores, title, log1p, pdf_file) {
    data_for_heatmap <- data[genes, cell_scores$cell_name] # Reorder columns based on cell_scores$cell_name

    if (log1p) {
      # Log-transform the data
      data_for_heatmap <- log1p(data_for_heatmap)
    }

    # Calculate color breaks
    data_range <- range(data_for_heatmap)
    breaks <- seq(data_range[1], data_range[2], length.out = 101)
    color_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(length(breaks) - 1)

    # Annotations
    annotation_data <- data.frame(
      Primary = factor(cell_scores$primary, levels = unique(cell_scores$primary)),
      Secondary = factor(cell_scores$secondary, levels = unique(cell_scores$secondary))
    )
    rownames(annotation_data) <- cell_scores$cell_name

    row_annotation <- data.frame(Gene_Type = factor(
      rep(c("Anabolic", "Catabolic", "Housekeeping"), 
          c(length(selected_anabolic_genes), length(selected_inflammatory_genes), length(selected_housekeeping_genes))),
      levels = c("Anabolic", "Catabolic", "Housekeeping")
    ))
    rownames(row_annotation) <- genes

    # Create custom color palette for drugs
    n_primary <- length(unique(cell_scores$primary))
    primary_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(max(n_primary, 2))
    names(primary_colors) <- unique(cell_scores$primary)

    n_secondary <- length(unique(cell_scores$secondary))
    secondary_colors <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Set2"))(max(n_secondary, 7))
    names(secondary_colors) <- unique(cell_scores$secondary)

    ann_colors <- list(
      Gene_Type = c(Anabolic = "forestgreen", Catabolic = "orange", Housekeeping = "grey"),
      Primary = primary_colors,
      Secondary = secondary_colors
    )

    # Save heatmap as a PDF
    pdf(pdf_file, width = pdf_width, height = pdf_height)
    pheatmap::pheatmap(data_for_heatmap,
                       fontsize = 6,
                       cluster_rows = FALSE, 
                       cluster_cols = FALSE,
                       show_rownames = TRUE,
                       show_colnames = FALSE,
                       annotation_col = annotation_data,
                       annotation_row = row_annotation,
                       annotation_colors = ann_colors,
                       color = color_palette,
                       breaks = breaks,
                       gaps_row = c(length(selected_anabolic_genes), length(selected_anabolic_genes) + length(selected_inflammatory_genes)))  # Add a gap between gene types
    dev.off()
  }

  # Create heatmap for the selected genes
  create_heatmap(sorted_selected_genes, expression_data, cell_scores, title, log1p, pdf_file)

  # Return the path to the saved PDF
  return(pdf_file)
}


### Figure 2f (UMAP plot)
# object: mac.sct
library(ggplot2)

# Define the function
generate_tsne_pdf <- function(seurat_object, pdf_file = "tsne_plot.pdf", 
                              pdf_width = 6, pdf_height = 5, font_size = 9) {
  
  # Step 1: Extract t-SNE coordinates
  # Note: Seurat usually stores t-SNE embeddings under the key "tsne"
  tsne_coordinates <- Embeddings(seurat_object, "tsne")
  
  # Convert to a data frame
  tsne_df <- as.data.frame(tsne_coordinates)
  # Standard Seurat t-SNE columns are usually tSNE_1, tSNE_2
  colnames(tsne_df) <- c("tSNE_1", "tSNE_2")
  
  # Add metadata
  tsne_df$drug_condition <- seurat_object@meta.data$drug_condition
  tsne_df$file <- seurat_object@meta.data$file
  
  # Define shape mapping if you need
  #shape_mapping <- scale_shape_manual(values = c(7, 9, 15, 16, 17, 18, 19))
  
  # Step 2: Create the t-SNE plot
  p <- ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2,
   #shape = file, 
   color = drug_condition)) +
    geom_point() +
    #shape_mapping +
    theme_classic(base_size = font_size) +
    labs(title = "t-SNE Plot", x = "t-SNE 1", y = "t-SNE 2")
  
  # Step 3: Save the plot to a PDF file
  pdf(pdf_file, width = pdf_width, height = pdf_height)
  print(p)
  dev.off()
  
  # Return the PDF file path
  return(pdf_file)
}

# Example usage
DefaultAssay(mac.sct) <- "SCT"
mac.sct <- RunPCA(mac.sct, assay = "SCT", npcs = 10, verbose = FALSE)
# Run t-SNE using the same PCs (dims 1:10)
mac.sct <- RunTSNE(mac.sct, dims = 1:10, verbose = FALSE)

generate_tsne_pdf(seurat_object = mac.sct, 
                  pdf_file = "./figures/figure2/fig2f_tsne_plot.pdf", 
                  pdf_width = 7, 
                  pdf_height = 5.5)

### Figure 2g (volcano plot)
library(dplyr)
library(ggplot2)
library(ggrepel)


# Define the function
generate_volcano_plot <- function(object, ident.1, group.by, assay, logfc.threshold, 
                                  pdf_file = "volcano_plot.pdf", 
                                  pdf_width = 8, pdf_height = 6) {
  # Step 1: Prepare data
  mac_deprep <- PrepSCTFindMarkers(object = object)
  de_results <- FindMarkers(mac_deprep, ident.1 = ident.1, group.by = group.by, 
                            assay = assay, logfmacc.threshold = logfc.threshold)

  # Calculate -log10(p-value adjusted)
  de_results$logP_adj <- -log10(de_results$p_val_adj)
  de_results$gene <- rownames(de_results)

  # Add gene classification
  de_results <- de_results %>%
    mutate(gene_type = case_when(
      avg_log2FC >= 0.5 & p_val_adj <= 0.05 ~ "up",
      avg_log2FC <= -0.5 & p_val_adj <= 0.05 ~ "down",
      TRUE ~ "ns"
    ))

  # Step 2: Customize aesthetics
  cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey")
  sizes <- c("up" = 2, "down" = 2, "ns" = 1)
  alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

  # Step 3: Create the volcano plot
  p <- ggplot(de_results, aes(x = avg_log2FC, y = logP_adj, label = gene, 
                              fill = gene_type, size = gene_type, alpha = gene_type)) +
    geom_point(shape = 21, color = 'black') +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "blue") +
    geom_text_repel(data = filter(de_results, gene_type == 'up'), aes(label = gene), size = 3) +
    geom_text_repel(data = filter(de_results, gene_type == 'down'), aes(label = gene), size = 3) +
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),    
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "right"
    ) +
    labs(title = "Volcano Plot", x = "Log Fold Change", y = "-Log10 P-value") +
    scale_fill_manual(values = cols) +
    scale_size_manual(values = sizes) +
    scale_alpha_manual(values = alphas) +
    scale_x_continuous(breaks = seq(-2, 3, 0.5), limits = c(-2, 3))

  # Step 4: Save the plot to a PDF file
  pdf(pdf_file, width = pdf_width, height = pdf_height)
  print(p)
  dev.off()

  # Return the PDF file path
  return(pdf_file)
}

# Example usage 
generate_volcano_plot(object = mac.sct, ident.1 = 'inflammatory', group.by = 'drug_condition', 
                     assay = 'SCT', logfc.threshold = 0.1, pdf_file = "./figures/figure2/fig2g_volcano_plot.pdf", 
                     pdf_width = 4, pdf_height = 2.5)


### Figure 2.h (Correlation of RNA-seq and ISS)
# raw rna-seq data of MAC (in house) 
gene_list <- c(anabolic, inflammatory, housekeeping)

# prepare RNA-seq data to compare
rna.sct <- SCTransform(rna.pub, vst.flavor='v2',verbose=FALSE, return.only.var.genes = FALSE)
cells_of_interest <- c('control','control.1','inflammatory','inflammatory.1')
rna.sct <- subset(rna.sct, cells = cells_of_interest)
rna.sct@meta.data$drug_condition <- c('control', 'control', 'inflammatory', 'inflammatory')
rna.sct<- subset(rna.sct, features = gene_list)

# MAC uninjured RNA-seq data from Sebastian 2021 paper.
macSebastian$drug_condition <- 'control'

# Example Usage
create_correlation_plots(seurat_object1 = mac.sct, seurat_object2 = rna.sct, pdf_file = "./figures/figure2/fig2h_macVSrna_correlation_plots.pdf")
create_correlation_plots(seurat_object1 = mac.sct, seurat_object2 = macSebastian, pdf_file = "./figures/figure2/fig2h_macVSsebastian_correlation_plots.pdf")
create_correlation_plots(seurat_object1 = rna.sct, seurat_object2 = macSebastian, pdf_file = "./figures/figure2/fig2h_rnaVSsebastian_correlation_plots.pdf")

### Pearson correlation coefficient and its p-value.
create_correlation_plots <- function(seurat_object1, seurat_object2, slot = 'data', group.by = 'drug_condition', pdf_file = "comparison_plots.pdf", pdf_width = 2.5, pdf_height = 2.5) {
  
  # Calculate average expression for each Seurat object
  avg_exp1 <- AverageExpression(seurat_object1, return.seurat = FALSE, group.by = group.by, verbose = FALSE)
  avg_exp2 <- AverageExpression(seurat_object2, return.seurat = FALSE, group.by = group.by, verbose = FALSE)
  avg_exp1 <- avg_exp1[["SCT"]]
  avg_exp2 <- avg_exp2[["SCT"]]
  if (slot == 'data') {
    avg_exp1 <- log1p(avg_exp1)
    avg_exp2 <- log1p(avg_exp2)
  } 
  # Get common drug conditions  
  drug_colnames <- intersect(colnames(avg_exp1), colnames(avg_exp2))
  if (length(drug_colnames) == 0) {
    print("There is no common drug condition between the two seurat objects.")
    if (colnames(avg_exp2) == 'all') {
      # If avg_exp2 is not a matrix, convert it to one:
      if (!is.matrix(avg_exp2)) {
        avg_exp2 <- as.matrix(avg_exp2)
      }
      drug_colnames <- names(table(seurat_object2@meta.data[[group.by]]))
      colnames(avg_exp2) <- drug_colnames
    } else if (colnames(avg_exp1) == 'all') {
      # If avg_exp1 is not a matrix, convert it to one:
      if (!is.matrix(avg_exp1)) {
        avg_exp1 <- as.matrix(avg_exp1)
      }
      drug_colnames <- names(table(seurat_object1@meta.data[[group.by]]))
      colnames(avg_exp1) <- drug_colnames
    }
  }
  
 
  if (length(drug_colnames) > 1) {
    data_sample1 <- avg_exp1[, drug_colnames, drop = FALSE]
    data_sample2 <- avg_exp2[, drug_colnames, drop = FALSE]
  } else {
    data_sample1 <- as.matrix(avg_exp1[, drug_colnames, drop = FALSE])
    data_sample2 <- as.matrix(avg_exp2[, drug_colnames, drop = FALSE])
  }
  
  # Subset the data frames to include only the common genes
  common_genes <- intersect(rownames(data_sample1), rownames(data_sample2))
  data_sample1 <- data_sample1[common_genes, , drop = FALSE]
  data_sample2 <- data_sample2[common_genes, , drop = FALSE]
  
  # Extract the first 4 letters of the Seurat object names
  obj1_name <- substr(deparse(substitute(seurat_object1)), 1, 4)
  obj2_name <- substr(deparse(substitute(seurat_object2)), 1, 4)
  
  # Create plots for each column (drug condition) in the data samples
  plot_list <- list()
  for (i in seq_along(drug_colnames)) {
    # Use cor.test to obtain both correlation coefficient and p-value
    test <- cor.test(as.numeric(data_sample1[, i]), as.numeric(data_sample2[, i]), method = "pearson")
    correlation_coefficient <- test$estimate
    p_value <- test$p.value
    
    plot_data <- data.frame(
      Expression_Sample1 = as.numeric(data_sample1[, i]),
      Expression_Sample2 = as.numeric(data_sample2[, i]),
      Gene = common_genes
    )
    
    # Calculate the slope to fit a line through the origin
    slope <- mean(plot_data$Expression_Sample2) / mean(plot_data$Expression_Sample1)
    
    # Add first principal component for coloring
    plot_data$PC <- predict(prcomp(~Expression_Sample1 + Expression_Sample2, data = plot_data))[, 1]
    
    plot <- ggplot(plot_data, aes(x = Expression_Sample1, y = Expression_Sample2, label = Gene, color = PC)) +
      geom_point(size = 1.5, alpha = 0.9, show.legend = FALSE) +
      scale_color_gradient(low = "#0091ff", high = "#f0650e") +
      geom_text(color = 'black', check_overlap = TRUE, hjust = 1, vjust = 1, size = 1.5) +
      # geom_abline(intercept = 0, slope = slope, color = "black", linetype = "dashed", linewidth = 0.5) +
      geom_smooth(method = 'glm') +
      labs(x = paste("Avg. Expression in", obj1_name),
           y = paste("Avg. Expression in", obj2_name),
           title = paste("Gene Expression Comparison -", drug_colnames[i])) +
      annotate("text",
               x = max(plot_data$Expression_Sample1) * 0.1,
               y = max(plot_data$Expression_Sample2) * 0.9,
               label = paste("Pearson r:", round(correlation_coefficient, 2),
                             "\np-value:", p_value),
               color = "blue", size = 1.5) +
      theme_classic(base_size = 5) +
      theme(text = element_text(family = "Helvetica", size = 5),
            axis.text = element_text(face = "plain"),
            plot.title = element_text(hjust = 0.5, size = 7)) +
      theme(legend.position = "none") +
      expand_limits(x = 0, y = 0)
    
    plot_list[[i]] <- plot
  }
  
  # Save all plots to a single PDF
  pdf(file = pdf_file, width = pdf_width, height = pdf_height)
  for (plot in plot_list) {
    print(plot)
  }
  dev.off()
}
