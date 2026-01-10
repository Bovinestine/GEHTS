set.seed(123)
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
# Gene categroy # IL17b is relocated to inflammatory from anabolic
anabolic <- c('Acan','Sox9','Col2a1','Matn1','Matn3','Ucma','Ccnd3','Gadd45g','Pth1r','Gm26633','Col27a1')
inflammatory <- c('Mmp3','Mmp13','Il6','Adamts5','Igfbp3','Ccl2','Cxcl5','Cxcl1','Fosl2','Tlr2','Tnfrsf1b', 'Il17b')
housekeeping <- c('Hprt','Actb','Gapdh','B2m','Ubc','Ppia','Rpl23')
gene_list <- c(anabolic, inflammatory, housekeeping)

mac69.sct <- subset(mac.SCT, subset = file %in% c('230608','230914'))
mac69 <- subset(mac, subset = file %in% c('230608','230914'))

# batch QC test
# Visualize nCount_RNA distribution by batch ("file")
VlnPlot(cmb.sct, 
        features = "nCount_SCT", 
        group.by = "file", 
        pt.size = 0.1, # Adjusts the size of the dots
        log = TRUE) +  # Optional: Log scale helps if there are huge outliers
  NoLegend() +
  ggtitle("Distribution of RNA Counts per Batch")

RidgePlot(mac.sct, 
          features = "nCount_SCT", 
          group.by = "file") +
  ggtitle("nCount_RNA Density by Batch")

# Calculate stats for nCount_RNA grouped by "file"
batch_stats <- mac.SCT@meta.data %>%
  group_by(file) %>%
  summarise(
    n_Cells = n(),
    Mean_nCount = mean(nCount_RNA),
    Median_nCount = median(nCount_RNA),
    SD_nCount = sd(nCount_RNA),
    Min_nCount = min(nCount_RNA),
    Max_nCount = max(nCount_RNA)
  )

# View the table
print(batch_stats)
write.csv(batch_stats, file="./figures/supple/batch_statistics.csv", row.names=FALSE)

# 1. Calculate the % of cells detecting each HK gene per batch
hk_detection_stats <- mac.SCT@meta.data %>%
  bind_cols(t(as.matrix(mac.SCT@assays$SCT@counts[housekeeping, ]))) %>%
  group_by(file) %>%
  summarise(
    # Calculate % of cells with count > 0 for key genes
    Ubc_Detection_Pct = mean(Ubc > 0) * 100,
    Actb_Detection_Pct = mean(Actb > 0) * 100,
    Gapdh_Detection_Pct = mean(Gapdh > 0) * 100,
    
    # Calculate Mean Expression levels
    Mean_Ubc = mean(Ubc),
    Mean_Actb = mean(Actb),
    Mean_Gapdh = mean(Gapdh)
  )

print(hk_detection_stats)
write.csv(hk_detection_stats, file="./figures/supple/hk_detection_statistics.csv", row.names=FALSE)


### Supplementary Figure 5 (Figure 2d (grid of gene expression box plots))
# anabolic
draw_boxplots_for_genes_with_common_legend(mac69.sct, genes = anabolic,
                        pdf_file = "./figures/figure.sup/S.fig.5_boxplots_ana.pdf",
                        pdf_width = 5, pdf_height = 3, base_font_size = 6, legend_gap = 0.03)
# catabolic
draw_boxplots_for_genes_with_common_legend(mac69.sct, genes = inflammatory,
                        pdf_file = "./figures/figure.sup/S.fig.5_boxplots_cata.pdf",
                        pdf_width = 5, pdf_height = 3, base_font_size = 6, legend_gap = 0.03)
# housekeeping
draw_boxplots_for_genes_with_common_legend(mac69.sct, genes = housekeeping,
                        pdf_file = "./figures/figure.sup/S.fig.5_boxplots_hk.pdf",
                        pdf_width = 3.5, pdf_height = 3, base_font_size = 6, legend_gap = 0.03)

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

### test edge effect in a chip
# 260105

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork) # For arranging plots

# --- 1. Define Helper Function for Spatial Mapping (Snake Pattern) ---
get_spatial_coords <- function(well_no) {
  # Ensure well_no is a character string to avoid Factor issues
  w_str <- as.character(well_no)
  
  # aggressive cleaning: remove 'no' and any whitespace
  w_clean <- gsub("no", "", w_str, ignore.case = TRUE)
  w_clean <- gsub(" ", "", w_clean)
  w <- as.numeric(w_clean)
  
  # Check if parsing failed (returns NA)
  if (is.na(w)) {
    return(list(Row = NA, Col = NA, Region = NA))
  }
  
  # Calculate Row (1-10)
  row_idx <- ceiling(w / 10)
  
  # Calculate Col (1-10) - Snake Pattern
  if (row_idx %% 2 != 0) {
    # Odd Rows (1, 3...): Left to Right
    col_idx <- (w - 1) %% 10 + 1
  } else {
    # Even Rows (2, 4...): Right to Left
    col_idx <- 10 - ((w - 1) %% 10)
  }
  
  # Define Region
  if ((row_idx %in% c(1, 10)) & (col_idx %in% c(1, 10))) {
    region <- "Corner"
  } else if (row_idx %in% c(1, 10) | col_idx %in% c(1, 10)) {
    region <- "Edge"
  } else {
    region <- "Center"
  }
  
  return(list(Row = row_idx, Col = col_idx, Region = region))
}

# --- 2. Prepare Data from Seurat Object ---
# Extract metadata
meta_df <- mac.sct@meta.data

# Ensure we have a "Total Counts" metric (Sensitivity)
# Usually 'nCount_RNA' or 'nCount_SCT'. Adjust if your slot is named differently.
if (!"Total_Counts" %in% colnames(meta_df)) {
  if ("nCount_RNA" %in% colnames(meta_df)) {
    meta_df$Total_Counts <- meta_df$nCount_RNA 
  } else {
    # Fallback or error if neither exists
    stop("Column 'nCount_RNA' not found in metadata.")
  }
}

# Apply the spatial function to every row
spatial_list <- lapply(meta_df$well_no, get_spatial_coords)
spatial_df <- dplyr::bind_rows(spatial_list)

# Merge back into the plotting dataframe
plot_data <- cbind(meta_df, spatial_df)

# IMPORTANT: Remove rows where parsing failed (NAs) to prevents errors in plotting
plot_data <- plot_data[!is.na(plot_data$Row), ]

# Create a clear Label for the Chip (e.g., "230608_control")
plot_data$ChipID <- paste(plot_data$file, plot_data$drug_condition, sep = "_")

# Set factor levels for clean plotting order
plot_data$Region <- factor(plot_data$Region, levels = c("Center", "Edge", "Corner"))

# --- 3. Generate Boxplot (Main Figure) ---
# This shows the statistical comparison of sensitivity across regions
p_boxplot <- ggplot(plot_data, aes(x = Region, y = Total_Counts, fill = Region)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.4) + # Adds raw points
  facet_wrap(~ChipID, scales = "free_y") +             # Separate panel per chip
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("Center" = "#66c2a5", "Edge" = "#fc8d62", "Corner" = "#8da0cb")) +
  labs(title = "Assessment of Edge Effects on Gene Detection Sensitivity",
       subtitle = "Comparison of Total Counts across Chip Regions",
       y = "Total Gene Counts",
       x = "") +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(face = "bold"))

# Print and Save Boxplot
print(p_boxplot)
ggsave("./figures/figure.sup/Edge_Effect_Boxplot.pdf", plot = p_boxplot, width = 8, height = 6)


# --- 4. Generate Spatial Heatmaps (Validation/Supplementary) ---
# This reconstructs the visual grid of the chip
if (nrow(plot_data) > 0) {
  
  p_heatmaps <- ggplot(plot_data, aes(x = Col, y = Row, fill = Total_Counts)) +
    # Force tiles to fill the grid square
    geom_tile(color = "white", lwd = 0.2, width = 1, height = 1) + 
    
    # Facet by ChipID (Creates one panel per chip automatically)
    facet_wrap(~ChipID, ncol = 4) + 
    
    # Axis formatting
    scale_y_reverse(breaks = 1:10, limits = c(10.5, 0.5)) + 
    scale_x_continuous(breaks = 1:10, limits = c(0.5, 10.5)) +
    scale_fill_viridis_c(option = "viridis", name = "Counts") +
    
    # Labels and Theme
    labs(title = "Spatial Uniformity Check (Raw Counts)", x = "Column", y = "Row") +
    coord_fixed() + 
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 10), # Chip titles
      axis.text = element_text(size = 8),
      legend.position = "right"
    )
  
  print(p_heatmaps)
  ggsave("./figures/figure.sup/Figure_Spatial_Heatmaps_Combined.pdf", plot = p_heatmaps, width = 12, height = 10)
  message("Saved Heatmaps to ./figures/figure.sup/Figure_Spatial_Heatmaps_Combined.pdf")
  
} else {
  stop("Error: plot_data is empty. Check 'well_no' format in your object.")
}


### Figure 2e Old (UMAP cluster heatmap)
library(ggplot2)
library(pheatmap)

png("./figures/figure2/heatmap_body.png", width = 2000, height = 2000, res = 300) # High resolution
pheatmap(
  distance_matrix_m,
  cluster_rows = FALSE,  # Disable dendrograms
  cluster_cols = FALSE,  # Disable dendrograms
  show_rownames = FALSE, # Remove row names
  show_colnames = FALSE
)
dev.off()


# Combine the PNG heatmap with vector annotations and dendrograms
pdf("./figures/figure2/final_heatmap_combined.pdf", width = 4, height = 4)

# Add the rasterized heatmap as the background
grid::grid.raster(png::readPNG("heatmap_body.png"))

# Add vector annotations and dendrograms
pheatmap(
  distance_matrix_m,
  color = NA,  # Disable heatmap colors again
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = annotation_df,
  annotation_col = annotation_df,
  show_rownames = FALSE,
  show_colnames = FALSE
)

# Close the PDF device
dev.off()

# Combine rasterized heatmap with labels (optional)
pdf("./figures/figure2/final_heatmap.pdf", width = 4, height = 4)
grid::grid.raster(png::readPNG("./figures/figure2/heatmap_rasterized.png"))
# Add additional vector elements, e.g., annotations
dev.off()

# Define the function to generate Heatmap PDF
generate_umap_heatmap_pdf <- function(umap_object, pdf_file = "umap_heatmap.pdf", 
                                 pdf_width = 6, pdf_height = 6) {
  # Step 1: Extract UMAP embeddings and metadata
  umap_1data <- umap_object[["umap"]]@cell.embeddings
  drug_conditions <- umap_object$drug_condition

  # Ensure the order of drug_conditions matches the order of umap_1data
  drug_conditions <- drug_conditions[rownames(umap_1data)]

  # Create an annotation data frame for the drug_condition
  annotation_df <- data.frame(DrugCondition = drug_conditions)

  # Calculate distance matrix
  distance_matrix <- dist(umap_1data)

  # Convert to a matrix for heatmap
  distance_matrix_m <- as.matrix(distance_matrix)

  # Step 2: Define the color palette
  color_palette <- colorRampPalette(colors = c("blue4", 'white', "red4"))(100)

  # Step 3: Create and save the heatmap to a PDF file
  pdf(pdf_file, width = pdf_width, height = pdf_height)
  pheatmap(
    distance_matrix_m, 
    clustering_distance_rows = distance_matrix, 
    clustering_distance_cols = distance_matrix,
    show_rownames = FALSE, 
    show_colnames = FALSE
  )
  dev.off()

  # Return the PDF file path
  return(pdf_file)
}

generate_umap_heatmap_pdf(umap_object = mac.sct.umap, pdf_file = "./figures/figure2/umap_heatmap.pdf", 
                          pdf_width = 4, pdf_height = 4)

