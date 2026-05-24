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

# Gene categroy 
anabolic <- c('Acan','Sox9','Col2a1','Matn1','Matn3','Ucma','Ccnd3','Gadd45g','Pth1r','Gm26633','Col27a1')
inflammatory <- c('Mmp3','Mmp13','Il6','Adamts5','Igfbp3','Ccl2','Cxcl5','Cxcl1','Fosl2','Tlr2','Tnfrsf1b', 'Il17b')
housekeeping <- c('Hprt','Actb','Gapdh','B2m','Ubc','Ppia','Rpl23')
gene_list <- c(anabolic, inflammatory, housekeeping)



### extended data Figure 6 (grid of gene expression box plots)
sfig_boxplot <- draw_boxplots_for_genes_with_common_legend(mac.sct, genes = gene_list, nrow=4,
                        base_font_size = 6, legend_gap = 0.03)
ggsave("./S.fig.6_boxplots_allgenes.pdf", plot=sfig_boxplot, width = 5, height = 5, units = "in", dpi = 300)

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(cowplot)
library(ggsignif)

#' Generate Publication-Ready Faceted Boxplots with Common Legend
#'
#' @param seurat_object A Seurat object containing single-cell data.
#' @param genes Character vector of genes to plot.
#' @param datatype Character. The data slot to use. Default is 'counts'.
#' @param control_label Character. Metadata string for the healthy baseline.
#' @param treatment_label Character. Metadata string for the disease baseline.
#' @param extra_label1 Character. Metadata string for treatment 1.
#' @param extra_label2 Character. Metadata string for treatment 2.
#' @param color_palette Named vector of hex colors for the conditions.
#' @param nrow Number of rows of genes in grid
#' @param base_font_size Base font size for ggplot theme.
#' @param line_thickness Thickness of axis lines and boxplot borders.
#' @param legend_gap Spacing ratio for the top legend.
#' @param y_title Label for the shared Y-axis.
#'
#' @return A combined cowplot object.
#' @export
draw_boxplots_for_genes_with_common_legend <- function(seurat_object,
                                                       genes,
                                                       datatype = 'counts',
                                                       control_label = "control",
                                                       treatment_label = "inflammatory",
                                                       extra_label1 = "label3",
                                                       extra_label2 = "label4",
                                                       color_palette = NULL,
                                                       nrow = 5,
                                                       pdf_file = "./figures/figure2/boxplots.pdf",
                                                       pdf_width = 6,
                                                       pdf_height = 6,
                                                       base_font_size = 6,
                                                       line_thickness = 0.5,
                                                       legend_gap = 0.05,
                                                       y_title = "Normalized gene expression (a.u)") {
  
  # 1. Setup Default Colors if none provided (NPG Palette)
  if (is.null(color_palette)) {
    color_palette <- setNames(
      c("#4DBBD559", "#E64B3599", "#A6CEE3", "#1f78b4"), 
      c("Basal", "IL-1β", extra_label1, extra_label2)
    )
  }
  
  # 2. Extract and Format Metadata
  meta <- seurat_object@meta.data
  valid_conditions <- c(control_label, treatment_label, extra_label1, extra_label2)
  cells_to_use <- rownames(meta)[ meta$drug_condition %in% valid_conditions ]
  
  # Map conditions to clean labels and lock the factor levels
  level_order <- c("Basal", "IL-1β", extra_label1, extra_label2)
  
  meta_subset <- meta[cells_to_use, ] %>%
    mutate(ConditionLabel = case_when(
      drug_condition == control_label ~ "Basal",
      drug_condition == treatment_label ~ "IL-1β",
      drug_condition == extra_label1 ~ extra_label1,
      drug_condition == extra_label2 ~ extra_label2,
      TRUE ~ as.character(drug_condition)
    )) %>%
    mutate(ConditionLabel = factor(ConditionLabel, levels = level_order))
  
  plot_list <- list()
  
  # 3. Generate Individual Gene Plots
  for (gene in genes) {
    expr <- as.vector(GetAssayData(seurat_object, slot = datatype)[gene, cells_to_use, drop = TRUE])
    df <- data.frame(Expression = expr,
                     Condition = meta_subset$ConditionLabel)
    
    # Calculate a safe y-max for significance bars based on absolute max value
    ymax <- max(df$Expression, na.rm = TRUE)
    if (is.infinite(ymax) || is.na(ymax) || ymax == 0) ymax <- 1
    step_size <- ymax * 0.08 # 8% step increase between brackets
    
    p_gene <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
      geom_boxplot(outlier.shape = NA, linewidth = line_thickness, alpha = 0.85) +
      scale_fill_manual(values = color_palette) +
      labs(title = gene) +
      theme_classic() +
      theme(text = element_text(size = base_font_size),
            axis.line = element_line(linewidth = line_thickness, color = "black"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = base_font_size, color = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(face = "bold.italic", hjust = 0.5, size = base_font_size + 2),
            legend.position = "none",
            plot.margin = ggplot2::margin(t = 10, r = 5, b = 2, l = 5))

    # Only test groups that have at least 2 valid observations
    group_counts <- table(df$Condition[!is.na(df$Expression)])
    valid_groups <- names(group_counts)[group_counts >= 2]
    
    if (length(valid_groups) >= 2) {
      comparisons_list <- combn(valid_groups, 2, simplify = FALSE)
      
      p_gene <- p_gene + geom_signif(
        comparisons = comparisons_list,
        map_signif_level = TRUE,
        y_position = ymax + (step_size * seq_along(comparisons_list)),
        textsize = base_font_size * 0.35,
        tip_length = 0.01,
        size = 0.3,  # --- BUG FIX: Reverted from linewidth to size ---
        vjust = 0.5
      )

      # Dynamically expand Y-axis limit so top significance bars aren't cut off
      top_bracket_y <- ymax + (step_size * (length(comparisons_list) + 1))
      p_gene <- p_gene + coord_cartesian(ylim = c(0, top_bracket_y))
    } else {
      # If not enough groups for a statistical test, just plot the boxes
      p_gene <- p_gene + coord_cartesian(ylim = c(0, ymax * 1.1))
    }
    
    plot_list[[gene]] <- p_gene
  }

  # 4. Generate Common Legend
  dummy_df <- data.frame(Condition = factor(level_order, levels = level_order),
                         Expression = rep(0.5, 4))
  
  dummy_plot <- ggplot(dummy_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot() +
    scale_fill_manual(values = color_palette) +
    theme_classic() +
    theme(legend.position = "top",
          legend.text = element_text(size = base_font_size + 2, face = "bold"),
          legend.title = element_blank(),
          legend.key.size = unit(0.4, "cm"),
          legend.margin = ggplot2::margin(0,0,0,0))
  
  # Helper to extract the legend grob
  get_legend <- function(myplot) {
    tmp <- ggplot_gtable(ggplot_build(myplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  common_legend <- get_legend(dummy_plot)
  
  # 5. Arrange the Grid
  n_genes <- length(genes)
  ncol <- ceiling(n_genes / nrow)
  grid_plots <- arrangeGrob(grobs = plot_list, ncol = ncol, nrow = nrow)
  
  arranged_plots <- grid.arrange(common_legend, grid_plots,
                                 ncol = 1,
                                 heights = c(legend_gap, 1))
  
  # 6. Final Assembly with Cowplot
  final_plot <- ggdraw() +
    draw_plot(arranged_plots, x = 0.05, y = 0, width = 0.95, height = 1) +
    draw_label(y_title, x = 0.02, y = 0.5,
               angle = 90, vjust = 0.5, size = base_font_size + 2, fontface = "bold")
  
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
