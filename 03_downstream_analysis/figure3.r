# GE-HTS project pipeline for in situ sequencing data analysis
# author: Nathan Wooseok Lee
# Date of start: 241224
# Date of update: 250711, 250827 (wide gene expression), 251223 (gene group)
# env: conda activate seurat4 / issr

# functions for figure 3.
# Ensure required packages are loaded
library(ggplot2)
library(ggpubr)     # for stat_compare_means()
library(cowplot)    # for ggdraw() and draw_label()
library(dplyr)
library(ggsignif)
library(Seurat)
library(gridExtra)
# For reproducibility
set.seed(123)
# Gene categroy
anabolic <- c('Acan','Sox9','Col2a1','Matn1','Matn3','Ucma','Ccnd3','Gadd45g','Pth1r','Gm26633','Col27a1')
inflammatory <- c('Mmp3','Mmp13','Il6', 'Il17b','Adamts5','Igfbp3','Ccl2','Cxcl5','Cxcl1','Fosl2','Tlr2','Tnfrsf1b')
housekeeping <- c('Hprt','Actb','Gapdh','B2m','Ubc','Ppia','Rpl23')

### Figure 3b (GEheatmap)
create_heatmap_by_expression_ps(sin.sct, upper_threshold = 0, bottom_threshold = 1, 
                                assay='SCT', slot='counts', primary_metadata = 'dose', secondary_metadata = 'drug_name',
                                pdf_file = './figures/figure3/fig3b_GEheatmap.pdf', pdf_width = 7, pdf_height = 5) 

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
  title <- "Zoomed in Expressed Genes"

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
    primary_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(max(n_primary, 13))
    names(primary_colors) <- unique(cell_scores$primary)

    n_secondary <- length(unique(cell_scores$secondary))
    secondary_colors <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Set2"))(max(n_secondary, 8))
    names(secondary_colors) <- unique(cell_scores$secondary)

    ann_colors <- list(
      Gene_Type = c(Anabolic = "forestgreen", Catabolic = "orange", Housekeeping = "grey"),
      Primary = primary_colors,
      Secondary = secondary_colors
    )

    # Save heatmap as a PDF
    pdf(pdf_file, width = pdf_width, height = pdf_height)
    pheatmap::pheatmap(data_for_heatmap,
                       cluster_rows = FALSE, 
                       cluster_cols = FALSE,
                       show_rownames = TRUE,
                       show_colnames = FALSE,
                       annotation_col = annotation_data,
                       annotation_row = row_annotation,
                       annotation_colors = ann_colors,
                       color = color_palette,
                       breaks = breaks,
                       main = title,
                       gaps_row = c(length(selected_anabolic_genes), length(selected_anabolic_genes) + length(selected_inflammatory_genes)))  # Add a gap between gene types
    dev.off()
  }

  # Create heatmap for the selected genes
  create_heatmap(sorted_selected_genes, expression_data, cell_scores, title, log1p, pdf_file)

  # Return the path to the saved PDF
  return(pdf_file)
}

### Figure 3c: UMAP dimensional reduction for showing reproducibility and separability among drug conditions
# targets <- c('mTOR', 'ALK', 'PI3K', 'Hedgehog','JNK', 'AKT', 'p38 MAPK', 'NF-kB', 'TNKS') # total 8 + 1 targets
# Mapping of drugs to their targets
drug_to_target <- c(
  "ALK5 inhibitor IV" = "ALK",
  "BMS-345541" = "NF-kB",
  "CAPE" = "NF-kB",
  "JNK inhibitor V" = "JNK",
  "KU0063794" = "mTOR",
  "LY294002" = "PI3K",
  "MK-2206 dihydrochloride" = "AKT",
  "pamapimod" = "p38 MAPK",
  "rapamycin" = "mTOR",
  "SANT-1" = "Hedgehog",
  "SB203580" = "p38 MAPK",
  "SB431542" = "ALK",
  "SB525334" = "ALK",
  "XAV" = "TNKS"
)
sin.sct$target <- drug_to_target[sin.sct$drug_name]

sin.sct10 <- subset(sin.sct, subset = dose ==10) # this is for main figure
sin.sct10 <- RunPCA(sin.sct10, assay = "SCT", npcs = 10, verbose = FALSE)
sin.sct10 <- RunUMAP(sin.sct10, reduction = "pca", assay = "SCT", dims = 1:10, verbose = FALSE)

pdf('./figures/figure3/fig3c_umap_single10.pdf', width = 6.5, height = 4)
DimPlot(sin.sct10, reduction = "umap", group.by='drug_name') # used for figure 3 c
dev.off()


# Figure 3d: highliting single drug conditions at a time along side with figure 3c.
umap_sin10 <- extractUmapData(sin.sct10, dose = 10)
plots_sin10 <- makeDrugUmapList(umap_sin10)

pdf('./figures/figure3/fig3d_umap_PamapimodXav_251230.pdf', width = 3, height = 6)
grid.arrange(grobs= plots_sin10[c(1,13)], ncol =1) # pamapimod and XAV
dev.off()

# for supplementary data
pdf('./figures/figure3/Sfig3d_umap_single10_1to6_251230.pdf', width = 6, height = 9)
grid.arrange(grobs = plots_sin10[1:6], ncol = 2) # check out the all conditions
dev.off()
pdf('./figures/figure3/Sfig3d_umap_single10_7to12_251230.pdf', width = 6, height = 9)
grid.arrange(grobs = plots_sin10[7:12], ncol = 2) # check out the all conditions
dev.off()
pdf('./figures/figure3/Sfig3d_umap_single10_8to13_251230.pdf', width = 6, height = 9)
grid.arrange(grobs = plots_sin10[8:13], ncol = 2) # check out the all conditions
dev.off()


# functions

extractUmapData <- function(seurat_obj, dose){
    seurat_obj <- seurat_obj[,seurat_obj$dose == dose]
    umap_data <- seurat_obj[["umap"]]@cell.embeddings
    umap <- data.frame(umap_data, seurat_obj@meta.data[, c("drug_name", "dose", "target")])    
    return(umap)
}

makeDrugUmapList <- function(umap_data){
    lapply(unique(umap_data$drug_name), function(condition) {
    ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = (drug_name == condition))) +
        geom_point(aes(alpha = (drug_name == condition))) +
        scale_alpha_manual(values = c("FALSE" = 0.1, "TRUE" = 0.6)) +
        scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "turquoise3")) +
        labs(title = paste("Highlighted:", condition)) +
        theme_classic() +
        theme(legend.position = "none")
    })
}


### Figure 3e and g
# fig.3e
sum_genes_and_plot_wide(mac_sin.sct, extra_label1 = 'XAV_0.1', extra_label2 = 'XAV_10', pdf_width = 2, pdf_height = 2, pdf_file='./figures/figure3/fig3e_boxplot.pdf')
# fig.3g
sum_genes_and_plot_wide(mac_sin.sct, extra_label1 = 'pamapimod_10', extra_label2 = 'XAV_10', pdf_width = 2, pdf_height = 2, pdf_file='./figures/figure3/fig3f_boxplot.pdf')

sum_genes_and_plot_wide <- function(seurat_object, 
                                    anabolic_genes = anabolic,
                                    inflammatory_genes = inflammatory,
                                    control_label = "control",
                                    treatment_label = "inflammatory",
                                    extra_label1 = "label3",
                                    extra_label2 = "label4",
                                    pdf_file = "./figures/figure3/fig3_sum_barplot.pdf",
                                    pdf_width = 6, 
                                    pdf_height = 6,
                                    base_font_size = 5,
                                    line_thickness = 0.5,
                                    datatype = 'counts',
                                    common_x_label = "",
                                    common_y_label = "Sum of normalized gene counts",
                                    title1 = "Anabolic Gene",
                                    title2 = "Catabolic Gene",
                                    comparisons_list = NULL) {
  
  # Compute sums and add them to metadata
  seurat_object[["anabolic_sum"]] <- Matrix::colSums(GetAssayData(seurat_object, slot = datatype)[anabolic_genes, , drop = FALSE])
  seurat_object[["inflammatory_sum"]] <- Matrix::colSums(GetAssayData(seurat_object, slot = datatype)[inflammatory_genes, , drop = FALSE])
  
  # Create vector of conditions to include
  conditions <- c(control_label, treatment_label, extra_label1, extra_label2)
  
  # Subset cells having a drug_condition in the provided conditions
  all_cells <- subset(seurat_object, subset = drug_condition %in% conditions)
    
  # Create data frames for anabolic and inflammatory gene sums
  anabolic_df <- data.frame(
    condition = all_cells@meta.data$drug_condition,
    sum = all_cells@meta.data$anabolic_sum
  )
  
  inflammatory_df <- data.frame(
    condition = all_cells@meta.data$drug_condition,
    sum = all_cells@meta.data$inflammatory_sum
  )

  if (treatment_label == "inflammatory") {
    new_label <- "IL-1β"
    anabolic_df$condition[anabolic_df$condition == treatment_label] <- new_label
    inflammatory_df$condition[inflammatory_df$condition == treatment_label] <- new_label
    treatment_label <- new_label
  }

  # Use the upper whisker (boxplot.stats) to position p-value labels
  anabolic_ymax <- boxplot.stats(anabolic_df$sum)$stats[5]
  inflammatory_ymax <- boxplot.stats(inflammatory_df$sum)$stats[5]
  
  # Define the list of all pairwise comparisons
  if (is.null(comparisons_list)){
    comparisons_list <- list(
      c(treatment_label, extra_label2),
      c(extra_label1, extra_label2)
    )
  } else if (comparisons_list == 'all') {
    comparisons_list <- combn(conditions, 2, simplify = FALSE) 
  }
  
  # Create anabolic gene plot
  p_anabolic <- ggplot(anabolic_df, aes(x = condition, y = sum, fill = condition)) +
    geom_boxplot(outliers = FALSE) +
    geom_signif(comparisons = comparisons_list, map_signif_level = TRUE,
                textsize = 3, y_position = anabolic_ymax, step_increase = 0.1, tip_length = 0.01) +
    labs(title = title1) +
    theme_classic() +
    theme(text = element_text(size = base_font_size),
          axis.line = element_line(linewidth = line_thickness),
          legend.position = "none", # REMOVED LEGEND
          # ENABLED X AXIS TEXT (TICKS)
          axis.text.x = element_text(size = base_font_size, angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(size = base_font_size),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(face = "plain", hjust = 0.5, size = base_font_size+2))
  
  # Create inflammatory gene plot
  p_inflammatory <- ggplot(inflammatory_df, aes(x = condition, y = sum, fill = condition)) +
    geom_boxplot(outliers = FALSE) +
    geom_signif(comparisons = comparisons_list, map_signif_level = TRUE, 
                textsize = 3, y_position = inflammatory_ymax, step_increase = 0.1, tip_length = 0.01) +
    labs(title = title2) +
    theme_classic() +
    theme(text = element_text(size = base_font_size),
          axis.line = element_line(linewidth = line_thickness),
          legend.position = "none", # REMOVED LEGEND
          # ENABLED X AXIS TEXT (TICKS)
          axis.text.x = element_text(size = base_font_size, angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(size = base_font_size),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(face = "plain", hjust = 0.5, size = base_font_size+2))
  
  # Combine the two plots side by side
  # No need to remove legend here as it is removed in theme() above
  combined_plots <- arrangeGrob(p_anabolic, p_inflammatory, ncol = 2)

  # No longer need get_legend or the grid.arrange step for the legend
  
  # Add common x and y axis labels using cowplot
  # Use combined_plots directly
  final_plot <- ggdraw(combined_plots) +
    # Optional: If you want a common X label at the very bottom, uncomment below
    # draw_label(common_x_label, x = 0.5, y = 0, vjust = -0.5, size = base_font_size + 2) +
    draw_label(common_y_label, x = 0, y = 0.5, angle = 90, vjust = 1.2, size = base_font_size + 2)
  
  # Save the final figure to a PDF file
  pdf(file = pdf_file, width = pdf_width, height = pdf_height)
  print(final_plot)
  dev.off()
  
  invisible(final_plot)
}


### Figure 3f and h (10 gene average box plots)

genes_to_plot <- c("Col2a1", "Matn3", "Gadd45g", "Sox9", "Acan", 
                    "Mmp3", "Mmp13", "Cxcl1", "Fosl2", "Tnfrsf1b") 

# fig.3f
draw_boxplots_for_genes_with_common_legend(mac_sin.sct, genes = genes_to_plot, extra_label1 = 'XAV_0.1', extra_label2 = 'XAV_10',
                        pdf_file = "./figures/figure3/fig3f_boxplots.pdf",
                        pdf_width = 4, pdf_height = 2.6, base_font_size = 6, legend_gap = 0.03)

# fig.3h
draw_boxplots_for_genes_with_common_legend(mac_sin.sct, genes = genes_to_plot, extra_label1 = 'pamapimod_10', extra_label2 = 'XAV_10',
                        pdf_file = "./figures/figure3/fig3h_boxplots.pdf",
                        pdf_width = 4, pdf_height = 2.6, base_font_size = 6, legend_gap = 0.03)
# for Supplementary data: p38 inhibitor comparision
draw_boxplots_for_genes_with_common_legend(mac_sin.sct, genes = genes_to_plot, extra_label1 = 'pamapimod_10', extra_label2 = 'SB203580_10',
                        pdf_file = "./figures/figure3/Sfig3_boxplots_p38inhibitors.pdf",
                        pdf_width = 4, pdf_height = 2.6, base_font_size = 6, legend_gap = 0.03)

draw_boxplots_for_genes_with_common_legend <- function(seurat_object,
                                                       genes,
                                                       datatype = 'counts',
                                                       control_label = "control",
                                                       treatment_label = "inflammatory",
                                                       extra_label1 = "label3",
                                                       extra_label2 = "label4",
                                                       pdf_file = "./figures/figure3/boxplots.pdf",
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
      drug_condition == treatment_label ~ "IL-1β",
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
  
  dummy_df <- data.frame(Condition = factor(c("Control", "IL-1β", extra_label1, extra_label2),
                                            levels = c("Control", "IL-1β", extra_label1, extra_label2)),
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




### Figure 3i (Efficacy Score heatmap)
# refer following code files:
# Scoringmethods.r
# PredictionScore.r

# Combine into a single data frame
'%notin%' <- Negate('%in%') # define negate of %in%
prob_sng10_rf
prob_sng0.1_rf

# Example usage
efficacy_heatmap(prob_sng10, prob_sng0.1, pdf_file = "./figures/figure3/efficacy_heatmap.pdf", 
                     pdf_width = 5, pdf_height = 2)


# Define a function to create a heatmap PDF
efficacy_heatmap <- function(prob_sng10, prob_sng0.1, 
                                 pdf_file = "heatmap_mean_probability.pdf", 
                                 pdf_width = 5, pdf_height = 2) {
  # Step 1: Filter out 'control' and 'inflammatory'
  prob_sng10_filtered <- prob_sng10[!names(prob_sng10) %in% c('control', 'inflammatory')]
  prob_sng0.1_filtered <- prob_sng0.1[!names(prob_sng0.1) %in% c('control', 'inflammatory')]

  # Step 2: Prepare heatmap data
  heatmap_data <- data.frame(
    Drug = c(names(prob_sng10_filtered), names(prob_sng0.1_filtered)),
    Dose = c(rep("10", length(prob_sng10_filtered)), rep("0.1", length(prob_sng0.1_filtered))),
    MeanProbability = c(prob_sng10_filtered, prob_sng0.1_filtered)
  )

  # Step 3: Create heatmap plot
  heatmap_plot <- ggplot(heatmap_data, aes(x = Drug, y = Dose, fill = MeanProbability)) +
    geom_tile() +
    scale_fill_viridis(option = 'viridis', na.value = "white") +
    coord_fixed() +  # Ensure square tiles
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(fill = "Mean Probability", x = "Drug", y = "Dose")

  # Step 4: Save the plot to a PDF
  pdf(pdf_file, width = pdf_width, height = pdf_height)
  print(heatmap_plot)
  dev.off()

  # Return the path to the PDF file
  return(pdf_file)
}




