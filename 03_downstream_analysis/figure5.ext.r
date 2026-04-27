### Gene Ontology Visualization functions.
### 260226 
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggsci)

plot_go_trajectory <- function(go_df) {
  
  # 1. Separate the data by type for layered plotting
  controls <- go_df %>% filter(type == "Control")
  singles <- go_df %>% filter(type == "Single Agent")
  combos <- go_df %>% filter(type == "Combination")
  
  # 2. Create the linking dataframe to draw arrows from Singles to Combos
  links <- singles %>%
    left_join(combos %>% select(group_id, Anabolic_Combo = Anabolic_Score, Catabolic_Combo = Catabolic_Score), 
              by = "group_id")
  
  # 3. Build the Plot
  p <- ggplot() +
    
    # Arrows (colored by group)
    geom_segment(data = links, 
                 aes(x = Catabolic_Score, y = Anabolic_Score, 
                     xend = Catabolic_Combo, yend = Anabolic_Combo, color = group_id),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"), 
                 linetype = "dashed", alpha = 0.6, linewidth = 0.6) +
    
    # Plot Single Agents (mapping shape to 'type')
    geom_point(data = singles, aes(x = Catabolic_Score, y = Anabolic_Score, color = group_id, shape = type), 
               size = 3.5, alpha = 0.8) +
    
    # Plot Combinations (mapping shape to 'type', using fill for the inside color)
    geom_point(data = combos, aes(x = Catabolic_Score, y = Anabolic_Score, fill = group_id, shape = type), 
               size = 5, color = "black", stroke = 1.2) +
    
    # Plot Controls (mapping shape, explicitly defining fill based on condition)
    geom_point(data = controls, aes(x = Catabolic_Score, y = Anabolic_Score, shape = type), 
               fill = ifelse(controls$condition == "Healthy Control", "white", "black"), 
               color = "black", size = 5, stroke = 1.2) +
    
    # Add Text Labels
    geom_text_repel(data = go_df, aes(x = Catabolic_Score, y = Anabolic_Score, label = condition),
                    box.padding = 0.8, point.padding = 0.5,
                    segment.color = "grey50", segment.size = 0.5,
                    fontface = ifelse(go_df$type == "Combination", "bold", "plain"),
                    size = ifelse(go_df$type == "Combination", 4.5, 3.5)) +
    
    # --- AESTHETICS & SCALES ---
    
    # 1. Apply the NPG color palette
    scale_color_npg() + 
    scale_fill_npg() +
    
    # 2. Define exactly which shapes to use
    scale_shape_manual(name = "Condition Type", 
                       values = c("Control" = 23,        # Diamond
                                  "Single Agent" = 15,   # Square
                                  "Combination" = 21)) + # Circle with border
    
    # 3. NBT formatting
    theme_classic(base_size = 14) +
    theme(
      axis.line = element_line(linewidth = 0.5),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      legend.background = element_rect(fill = "transparent", color = NA)
    ) +
    
    # 4. Legend Control: Hide the color/fill legends, style the shape legend
    guides(
      color = "none", 
      fill = "none",
      shape = guide_legend(override.aes = list(size = 4, fill = "grey50", color = "black"))
    ) +
    
    labs(
      x = "Catabolic GO Score Sum",
      y = "Anabolic GO Score Sum"
    )
  
  return(p)
}

fig5e_go_traj <- plot_go_trajectory(go_df)
ggsave("./figures/figure5/Figure5e_GO_trajectory.pdf", plot = fig5e_go_traj, width = 6, height = 5, units = "in", dpi = 300)


#### figure Bubble Plot
library(ggplot2)
library(dplyr)
library(tidyr)

#' Plot GO Term Bubble Plot for Target Combinations
#'
#' @param seurat_obj A Seurat object containing GO_Score1 to GO_Score6 in metadata
#' @param top5_combos A character vector of the 5 combination names to plot
#' @return A publication-ready ggplot object
plot_go_bubble <- function(seurat_obj, top5_combos) {
  
  # 1. Extract metadata
  meta <- seurat_obj@meta.data
  
  # 2. Map cells to their respective plotting groups
  meta <- meta %>%
    mutate(plot_group = case_when(
      drug_condition == "inflammatory" ~ "Inflammatory\n(Disease)",
      drug_condition == "control" ~ "Healthy\nControl",
      std_combi %in% top5_combos ~ std_combi,
      TRUE ~ NA_character_ # Ignore all other single agents and weak combos
    )) %>%
    filter(!is.na(plot_group))
  
  # 3. Define the GO columns and assign publication-ready labels
  go_mapping <- c(
    "GO_Score1" = "Chondrocyte Dev (Anabolic)",
    "GO_Score2" = "Cartilage Condensation (Anabolic)",
    "GO_Score3" = "Tissue Homeostasis (Anabolic)",
    "GO_Score4" = "NO Synthase Upreg (Catabolic)",
    "GO_Score5" = "ECM Disassembly (Catabolic)",
    "GO_Score6" = "Inflammatory Response (Catabolic)"
  )
  
  # 4. Reshape data and calculate DotPlot statistics
  dot_df <- meta %>%
    select(plot_group, all_of(names(go_mapping))) %>%
    pivot_longer(cols = starts_with("GO_Score"), names_to = "go_col", values_to = "score") %>%
    mutate(go_term = go_mapping[go_col]) %>%
    group_by(plot_group, go_term) %>%
    summarise(
      mean_score = mean(score, na.rm = TRUE),
      # Calculate % of cells where the module is actively expressed (score > 0)
      pct_expressed = (sum(score > 0, na.rm = TRUE) / n()) * 100,
      .groups = "drop"
    )
  
  # 5. Scale the mean_score (Z-score) per GO term
  # This is critical: it ensures the color gradient is relative across conditions 
  # for a specific pathway, making it easy to see where rescue happens.
  dot_df <- dot_df %>%
    group_by(go_term) %>%
    mutate(scaled_score = scale(mean_score)[, 1]) %>%
    ungroup()
  
  # 6. Factor ordering for the narrative flow
  # X-Axis: Disease on the left, Combos in the middle, Healthy on the right
  dot_df$plot_group <- factor(dot_df$plot_group, 
                              levels = c("Inflammatory\n(Disease)", top5_combos, "Healthy\nControl"))
  
  # Y-Axis: Group Catabolic terms together, and Anabolic terms together
  catabolic_terms <- c("Inflammatory Response (Catabolic)", "ECM Disassembly (Catabolic)", "NO Synthase Upreg (Catabolic)")
  anabolic_terms  <- c("Tissue Homeostasis (Anabolic)", "Cartilage Condensation (Anabolic)", "Chondrocyte Dev (Anabolic)")
  dot_df$go_term <- factor(dot_df$go_term, levels = rev(c(catabolic_terms, anabolic_terms)))
  
  # 7. Build the ggplot
  p <- ggplot(dot_df, aes(x = plot_group, y = go_term)) +
    
    # Use shape 21 (circle with border) for a crisp, publication look
    geom_point(aes(size = pct_expressed, fill = scaled_score), shape = 21, color = "black", stroke = 0.5) +
    
    # Diverging color scale (Blue = Suppressed, Red = Activated)
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, 
                         name = "Relative\nExpression\n(Z-score)") +
    
    # Scale dot sizes appropriately
    scale_size_continuous(range = c(2, 9), limits = c(0, 100), name = "% of Cells\nActive") +
    
    # Publication Theme
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
      axis.text.y = element_text(face = "bold", color = "black"),
      axis.title = element_blank(), # Remove axis titles as the labels are self-explanatory
      panel.grid.major = element_line(color = "grey85", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10)
    )
  
  return(p)
}

# --- Execution Example ---
# Assuming 'top5_combos' is a vector containing the 5 combination strings from the previous step
fig5e_plot <- plot_go_bubble(merged.sct, top5_combos)
ggsave("./figures/figure5/Figure5e_GO_BubblePlot.pdf", plot = fig5e_plot, width = 6, height = 5, units = "in", dpi = 300)



### STACKED Violin plot
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)

#' Plot Stacked Violins for Targeted Modules
#'
#' @param seurat_obj A Seurat object with GO_Score1 to 6 in metadata
#' @param top5_combos Character vector of the 5 combination names to plot
#' @return A publication-ready ggplot object
plot_stacked_violins <- function(seurat_obj, top5_combos) {
  
  # 1. Extract metadata
  meta <- seurat_obj@meta.data
  
  # 2. Map cells to their respective plotting groups
  meta <- meta %>%
    mutate(plot_group = case_when(
      drug_condition == "inflammatory" ~ "Inflammatory\n(Disease)",
      drug_condition == "control" ~ "Healthy\nControl",
      std_combi %in% top5_combos ~ std_combi,
      TRUE ~ NA_character_ # Drop everything else
    )) %>%
    filter(!is.na(plot_group))
  
  # 3. Define the columns and assign publication-ready module names
  go_mapping <- c(
    "GO_Score1" = "Chondrocyte Dev\n(Anabolic)",
    "GO_Score2" = "Cartilage Condensation\n(Anabolic)",
    "GO_Score3" = "Tissue Homeostasis\n(Anabolic)",
    "GO_Score4" = "NO Synthase Upreg\n(Catabolic)",
    "GO_Score5" = "ECM Disassembly\n(Catabolic)",
    "GO_Score6" = "Inflammatory Response\n(Catabolic)"
  )
  
  # 4. Reshape data into long format for ggplot faceting
  vln_df <- meta %>%
    select(plot_group, all_of(names(go_mapping))) %>%
    pivot_longer(cols = starts_with("GO_Score"), names_to = "go_col", values_to = "score") %>%
    mutate(module = go_mapping[go_col])
  
  # 5. Factor ordering for the narrative flow
  # X-Axis: Disease -> Combos -> Healthy
  vln_df$plot_group <- factor(vln_df$plot_group, 
                              levels = c("Inflammatory\n(Disease)", top5_combos, "Healthy\nControl"))
  
  # Y-Axis Facets: Catabolic on top, Anabolic on bottom
  catabolic_terms <- c("Inflammatory Response\n(Catabolic)", "ECM Disassembly\n(Catabolic)", "NO Synthase Upreg\n(Catabolic)")
  anabolic_terms  <- c("Tissue Homeostasis\n(Anabolic)", "Cartilage Condensation\n(Anabolic)", "Chondrocyte Dev\n(Anabolic)")
  vln_df$module <- factor(vln_df$module, levels = c(catabolic_terms, anabolic_terms))
  
  # 6. Build the Stacked Violin Plot
  p <- ggplot(vln_df, aes(x = plot_group, y = score, fill = plot_group)) +
    
    # scale="width" ensures all violins are visually comparable even if cell counts vary slightly
    geom_violin(scale = "width", trim = TRUE, color = "black", linewidth = 0.4, alpha = 0.9) +
    
    # (Optional) Add a tiny boxplot inside the violin to show the median and IQR
    geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.3) +
    
    # Apply the Nature Publishing Group color palette
    scale_fill_npg() +
    
    # Stack the plots vertically based on the module
    # switch = "y" moves the facet labels to the left side where the Y-axis title usually is
    facet_grid(module ~ ., scales = "free_y", switch = "y") +
    
    # NBT Formatting
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none", # Color is already encoded by the X-axis labels
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
      axis.text.y = element_text(size = 9, color = "grey30"),
      axis.title = element_blank(), # Facet labels act as the Y-axis title
      
      # Clean up the facet strips to look like continuous Y-axis labels
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 11),
      
      # Adjust panel spacing so the stacked violins sit tightly together
      panel.spacing = unit(0.1, "lines"),
      axis.line = element_line(linewidth = 0.5)
    )
  
  return(p)
}

# --- Execution Example ---
fig5e_violins <- plot_stacked_violins(merged.sct, top5_combos)
ggsave("./figures/figure5/Figure5e_StackedViolins.pdf", plot = fig5e_violins, width = 7, height = 8, units = "in", dpi = 300)

