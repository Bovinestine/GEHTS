# Project name: GEHTS-chip
# Author: Nathan Wooseok Lee
# date: 260302

anabolic <- c('Acan','Sox9','Col2a1','Matn1','Matn3','Ucma','Ccnd3','Gadd45g','Pth1r','Gm26633','Col27a1')
inflammatory <- c('Mmp3','Mmp13','Il6', 'Il17b','Adamts5','Igfbp3','Ccl2','Cxcl5','Cxcl1','Fosl2','Tlr2','Tnfrsf1b')
housekeeping <- c('Hprt','Actb','Gapdh','B2m','Ubc','Ppia','Rpl23')

### Dataset processing

# ====================================================================
# 1. Prepare the Control Object (mac.sct)
# ====================================================================
# Assign ground-truth efficacy based on the drug_condition column
mac.sct$efficacy <- ifelse(mac.sct$drug_condition == "control", 1.0, 
                           ifelse(mac.sct$drug_condition == "inflammatory", 0.0, NA))

# To ensure smooth merging, add placeholder columns for the drug names
mac.sct$drug_name1 <- mac.sct$drug_condition
mac.sct$drug_name2 <- mac.sct$drug_condition

# ====================================================================
# 2. Prepare the Combinations Object (cmb.sct)
# ====================================================================
# Create the alphabetized standard combination name in cmb.sct metadata
cmb.sct$std_combi <- paste0(
  pmin(as.character(cmb.sct$drug_name1), as.character(cmb.sct$drug_name2)), 
  "&", 
  pmax(as.character(cmb.sct$drug_name1), as.character(cmb.sct$drug_name2))
)

# Alphabetize the names in the ML efficacy vector (prob_cmb0.1_rf)
eff_names <- names(prob_cmb0.1_rf)
eff_splits <- strsplit(eff_names, "&")
std_eff_names <- paste0(
  pmin(sapply(eff_splits, `[`, 1), sapply(eff_splits, `[`, 2)), 
  "&", 
  pmax(sapply(eff_splits, `[`, 1), sapply(eff_splits, `[`, 2))
)

# Create a named vector and map it directly into the cmb.sct metadata
efficacy_map <- as.numeric(prob_cmb0.1_rf)
names(efficacy_map) <- std_eff_names
cmb.sct$efficacy <- efficacy_map[cmb.sct$std_combi]

# ====================================================================
# 3. Merge the Objects
# ====================================================================
# Merge using add.cell.ids to prevent any cell barcode collisions
merged.sct <- merge(x = mac.sct, y = cmb.sct, add.cell.ids = c("mac", "cmb"))

# Ensure the active assay is SCT
DefaultAssay(merged.sct) <- "SCT"

# For a targeted 30-gene panel, all genes are variable features
VariableFeatures(merged.sct) <- rownames(merged.sct)





### Macroscopic view of the combination pairs. 

# Required dependencies for this function
library(dplyr)
library(ggplot2)
library(ggrepel)

#' Figure 5d Plot Bliss Synergy vs. ML Efficacy for High-Throughput Screening Combinations
#'
#' Plot Bliss Synergy vs. ML Efficacy for High-Throughput Screening Combinations
#'
#' Evaluates machine learning-derived efficacy scores against calculated Bliss 
#' independence to identify highly synergistic, disease-modifying drug combinations.
#' The function isolates the top candidates based on a composite rank-sum of 
#' raw efficacy and synergistic excess.
#'
#' @param cmb_rf A named numeric vector of ML efficacy scores for combinations. 
#' @param sng_rf A named numeric vector of ML efficacy scores for single agents.
#' @param top_n Integer. The number of top synergistic combinations to highlight (default = 5).
#'
#' @return A named list containing two elements:
#'   \item{plot}{A publication-ready ggplot2 object visualizing the synergy landscape.}
#'   \item{top_candidates}{A dataframe containing the exact top_n candidates and their metrics.}
#' @export
plot_bliss_synergy <- function(cmb_rf, sng_rf, top_n = 5) {
  
  # -------------------------------------------------------------------
  # 1. Input Validation
  # -------------------------------------------------------------------
  if (is.null(names(cmb_rf)) || is.null(names(sng_rf))) {
    stop("Input Error: Both cmb_rf and sng_rf must be named numeric vectors.")
  }
  
  # -------------------------------------------------------------------
  # 2. Parse Combinations and Calculate Bliss Synergy
  # -------------------------------------------------------------------
  cmb_names <- names(cmb_rf)
  splits <- strsplit(cmb_names, "&")
  
  drug1 <- trimws(sapply(splits, `[`, 1))
  drug2 <- trimws(sapply(splits, `[`, 2))
  
  eff_obs <- as.numeric(cmb_rf)
  eff_A <- as.numeric(sng_rf[drug1])
  eff_B <- as.numeric(sng_rf[drug2])
  
  bliss_pred <- eff_A + eff_B - (eff_A * eff_B)
  synergy_score <- eff_obs - bliss_pred
  
  # -------------------------------------------------------------------
  # 3. Build Dataset and Handle Replicates/Missing Data
  # -------------------------------------------------------------------
  df <- data.frame(
    # Standardize combination name to ensure no A&B vs B&A duplicates
    std_combi = paste0(pmin(drug1, drug2), "&", pmax(drug1, drug2)),
    drug1 = drug1,
    drug2 = drug2,
    E_obs = eff_obs,
    Bliss_pred = bliss_pred,
    Synergy = synergy_score
  ) %>% 
    filter(!is.na(Synergy)) %>%
    # Group and average in case there are multiple technical replicates
    group_by(std_combi, drug1, drug2) %>%
    summarise(
      E_obs = mean(E_obs, na.rm = TRUE),
      Synergy = mean(Synergy, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (nrow(df) == 0) {
    stop("Execution Error: No valid combinations remained after matching single-agent names.")
  }
  
  # -------------------------------------------------------------------
  # 4. Algorithmic Top Candidate Selection (Rank-Sum Method)
  # -------------------------------------------------------------------
  df <- df %>%
    mutate(
      Rank_E = rank(-E_obs, ties.method = "min"),
      Rank_S = rank(-Synergy, ties.method = "min"),
      Rank_Sum = Rank_E + Rank_S
    ) %>%
    arrange(Rank_Sum, desc(E_obs))
  
  # Isolate the top hits
  highlight_df <- head(df, top_n) %>%
    mutate(display_name = paste0(drug1, " + ", drug2))
  
  # -------------------------------------------------------------------
  # 5. Construct the Volcano-Style Scatter Plot
  # -------------------------------------------------------------------
  p <- ggplot(df, aes(x = E_obs, y = Synergy)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_point(aes(fill = Synergy), shape = 21, size = 3, 
               color = "black", stroke = 0.3, alpha = 0.85) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", 
                         midpoint = 0, name = "Bliss Synergy\nScore") +
    geom_point(data = highlight_df, aes(fill = Synergy), 
               shape = 21, size = 4.5, color = "black", stroke = 1.5) +
    geom_text_repel(data = highlight_df, aes(label = display_name),
                    box.padding = 1.0, point.padding = 0.5,
                    nudge_x = -0.01, nudge_y = 0.03, 
                    segment.color = "black", segment.size = 0.5,
                    max.overlaps = Inf, fontface = "bold", size = 4.0) +
    theme_classic(base_size = 14) +
    theme(
      axis.line = element_line(linewidth = 0.5),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10)
    ) +
    labs(
      x = "ML Efficacy Score",
      y = "Bliss Synergy Score"
    )
  
  # -------------------------------------------------------------------
  # 6. Return Both the Plot and the Underlying Data
  # -------------------------------------------------------------------
  return(list(
    plot = p,
    top_candidates = highlight_df
  ))
}

# --- Execution Example ---
# 1. Run the function
synergy_results <- plot_bliss_synergy(prob_cmb0.1_rf, prob_sng0.1_rf, top_n = 5)

# 2. Extract and save the plot (Figure 5d)
fig5d_plot <- synergy_results$plot
ggsave("./figures/figure5/Figure5d_Efficacy_vs_Synergy.pdf", plot = fig5d_plot, width = 6.5, height = 5, units = "in", dpi = 300)

# 3. Extract the Top 5 Combinations vector (Pass this to your Heatmap function!)
top5_combos <- synergy_results$top_candidates$std_combi

# 4. Extract the parent single agents (For your Mechanistic Barplots)
target_singles <- unique(c(synergy_results$top_candidates$drug1, synergy_results$top_candidates$drug2))


### figure 5e ML efficacy axis VS Biological axis (is the efficacy triage perfect? No)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(ggsci)

#' Plot ML Efficacy vs. Biological Rescue Index
#'
#' Dynamically calculates clinical GO modules using exact mean expression to 
#' construct a Biological Rescue Index (Anabolic - Catabolic expression). 
#' Plots this orthogonal biological metric against machine learning-derived 
#' efficacy scores to isolate true regenerative synergies.
#'
#' @param seurat_obj A Seurat object containing an "SCT" assay and drug condition metadata.
#' @param prob_cmb0.1_rf Named numeric vector of ML Efficacy scores for combinations.
#' @param top5_combos Character vector of the 5 top combination names (std_combi format).
#' @return A publication-ready ggplot2 scatter plot object.
#' @export
plot_biological_alignment <- function(seurat_obj, prob_cmb0.1_rf, top5_combos) {
  
  # ====================================================================
  # 1. Dynamically Compute the 6 Clinical GO Modules
  # ====================================================================
  go_terms_genes <- list(
    GO_Score1 = c("Matn1", "Acan", "Sox9", "Col27a1"), # Chondrocyte Dev
    GO_Score2 = c("Col2a1", "Acan", "Sox9"),           # Cartilage Condensation
    GO_Score3 = c("Sox9", "Col2a1", "Pth1r", "Fosl2"), # Tissue Homeostasis
    GO_Score4 = c("Tlr2", "Ccl2"),                     # NO Synthase Upreg
    GO_Score5 = c("Adamts5", "Mmp13"),                 # ECM Disassembly
    GO_Score6 = c("Tlr2", "Ccl2", "Il17b", "Tnfrsf1b", "Cxcl1", "Il6", "Cxcl5", "Fosl2") # Inflammatory
  )
  
  # Extract the log-normalized data matrix
  expr_matrix <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")
  
  # Helper function to compute mean expression for valid genes
  calc_mean <- function(genes, mat) {
    valid_genes <- intersect(genes, rownames(mat))
    if(length(valid_genes) == 0) return(rep(0, ncol(mat)))
    if(length(valid_genes) == 1) return(mat[valid_genes, ])
    return(colMeans(mat[valid_genes, , drop = FALSE]))
  }
  
  # Compute and assign the 6 GO Scores directly to metadata
  for (mod_name in names(go_terms_genes)) {
    seurat_obj[[mod_name]] <- calc_mean(go_terms_genes[[mod_name]], expr_matrix)
  }
  
  # ====================================================================
  # 2. Extract and Calculate the Biological Rescue Index
  # ====================================================================
  meta <- seurat_obj@meta.data %>%
    mutate(
      Anabolic_Score = GO_Score1 + GO_Score2 + GO_Score3,
      Catabolic_Score = GO_Score4 + GO_Score5 + GO_Score6,
      Rescue_Index = Anabolic_Score - Catabolic_Score,
      
      # Standardize the combination identifier for merging
      plot_id = case_when(
        drug_condition %in% c("control", "inflammatory") ~ drug_condition,
        drug_name1 != drug_name2 ~ paste0(
          pmin(as.character(drug_name1), as.character(drug_name2)), "&",
          pmax(as.character(drug_name1), as.character(drug_name2))
        ),
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(plot_id))
  
  # Aggregate the Rescue Index by condition
  bio_df <- meta %>%
    group_by(plot_id) %>%
    summarise(Rescue_Index = mean(Rescue_Index, na.rm = TRUE), .groups = "drop")
  
  # ====================================================================
  # 3. Standardize and Aggregate the ML Efficacy Scores
  # ====================================================================
  eff_names <- names(prob_cmb0.1_rf)
  splits <- strsplit(eff_names, "&")
  std_eff_names <- paste0(
    pmin(trimws(sapply(splits, `[`, 1)), trimws(sapply(splits, `[`, 2))), "&",
    pmax(trimws(sapply(splits, `[`, 1)), trimws(sapply(splits, `[`, 2)))
  )
  
  eff_df <- data.frame(
    plot_id = std_eff_names,
    Efficacy = as.numeric(prob_cmb0.1_rf)
  ) %>% 
    group_by(plot_id) %>% 
    summarise(Efficacy = mean(Efficacy, na.rm = TRUE), .groups = "drop")
  
  # Add the absolute ground-truth efficacy for the controls
  controls_eff <- data.frame(
    plot_id = c("control", "inflammatory"),
    Efficacy = c(1.0, 0.0) 
  )
  all_eff <- bind_rows(eff_df, controls_eff)
  
  # ====================================================================
  # 4. Merge Datasets and Prepare Visual Layers
  # ====================================================================
  plot_df <- inner_join(bio_df, all_eff, by = "plot_id")
  
  plot_df <- plot_df %>%
    mutate(
      category = case_when(
        plot_id == "control" ~ "Control",
        plot_id == "inflammatory" ~ "IL-1β",
        plot_id %in% top5_combos ~ "Top 5 Hit",
        TRUE ~ "Other Combinations"
      ),
      # Create clean display labels for the repels
      label = case_when(
        plot_id == "control" ~ "Control",
        plot_id == "inflammatory" ~ "IL-1β",
        plot_id %in% top5_combos ~ gsub("&", " + ", plot_id),
        TRUE ~ ""
      )
    )
  
  # Segregate data for layered plotting
  bg_df <- plot_df %>% filter(category == "Other Combinations")
  top_df <- plot_df %>% filter(category == "Top 5 Hit")
  ctrl_df <- plot_df %>% filter(category %in% c("Control", "IL-1β")) %>%
    mutate(fill_color = ifelse(plot_id == "control", "white", "black"))
  
  # ====================================================================
  # 5. Build the Scatter Plot
  # ====================================================================
  p <- ggplot() +
    
    # Add a baseline for biological neutrality (Anabolic = Catabolic)
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.8) +
    
    # Plot background combinations (greyed out)
    geom_point(data = bg_df, aes(x = Efficacy, y = Rescue_Index), 
               color = "grey70", size = 2.5, alpha = 0.6) +
    
    # Plot Top 5 hits (highlighted with NPG colors)
    geom_point(data = top_df, aes(x = Efficacy, y = Rescue_Index, fill = plot_id), 
               shape = 21, size = 5, color = "black", stroke = 1.2) +
    
    # Plot Controls (distinct diamond shapes)
    geom_point(data = ctrl_df, aes(x = Efficacy, y = Rescue_Index), 
               shape = 23, fill = ctrl_df$fill_color, color = "black", size = 5, stroke = 1.2) +
    
    # Add smart labels to the highlighted points
    geom_text_repel(data = bind_rows(top_df, ctrl_df), 
                    aes(x = Efficacy, y = Rescue_Index, label = label),
                    box.padding = 1.2, point.padding = 0.5,
                    segment.color = "black", segment.size = 0.5,
                    fontface = "bold", size = 4.5,
                    max.overlaps = Inf) +
    
    # Thematic and aesthetic formatting
    scale_fill_npg() +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none", # Legend not needed due to direct labeling
      axis.line = element_line(linewidth = 0.5),
      axis.title = element_text(face = "bold")
    ) +
    labs(
      x = "ML Efficacy Score",
      y = bquote(* Sigma * "Anabolic - "* Sigma *"Catabolic Expression")
    )
  
  return(p)
}

# --- Execution Example ---
fig5e_efficacyVSBio_plot <- plot_biological_alignment(merged.sct, prob_cmb0.1_rf, top5_combos)
ggsave("./figures/figure5/Figure5e_BiologicalAlignment.pdf", plot = fig5e_efficacyVSBio_plot, width = 5.5, height = 6, units = "in", dpi = 300)



#' Figure 5f biological function heatmap
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)
library(scales) 

#' Plot Compact Log2FC Heatmap for Efficacy and Safety
#'
#' Dynamically computes 6 clinical GO modules and 3 core cellular state modules 
#' using exact mean expression calculation (bypassing AddModuleScore artifacts). 
#' Evaluates top combinations using Log2FC relative to the disease baseline. 
#' Drops the reference column and uses vertical labels to optimize publication space.
#' 
#' @param seurat_obj A Seurat object containing an "SCT" assay.
#' @param top5_combos A character vector of the top 5 combination names.
#' @return A publication-ready ggplot2 heatmap object.
#' @export
plot_heatmap_log2fc_compact <- function(seurat_obj, top5_combos) {
  
  # ====================================================================
  # 1. Define All 9 Module Gene Sets
  # ====================================================================
  # The 6 Clinical GO Modules
  go_terms_genes <- list(
    GO_Score1 = c("Matn1", "Acan", "Sox9", "Col27a1"),
    GO_Score2 = c("Col2a1", "Acan", "Sox9"),
    GO_Score3 = c("Sox9", "Col2a1", "Pth1r", "Fosl2"),
    GO_Score4 = c("Tlr2", "Ccl2"),
    GO_Score5 = c("Adamts5", "Mmp13"),
    GO_Score6 = c("Tlr2", "Ccl2", "Il17b", "Tnfrsf1b", "Cxcl1", "Il6", "Cxcl5", "Fosl2")
  )
  
  # The 3 Core Cellular State Modules
  state_genes <- list(
    State_Viability = c('Hprt', 'Actb', 'Gapdh', 'B2m', 'Ubc', 'Ppia', 'Rpl23'),
    State_Stress    = c('Gadd45g', 'Igfbp3'),
    State_Prolif    = c('Ccnd3')
  )
  
  # ====================================================================
  # 2. Safely Compute Mean Expression for All Modules
  # ====================================================================
  # Extract the log-normalized data matrix
  expr_matrix <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")
  
  # Helper function to compute mean expression for valid genes
  calc_mean <- function(genes, mat) {
    valid_genes <- intersect(genes, rownames(mat))
    if(length(valid_genes) == 0) return(rep(0, ncol(mat)))
    if(length(valid_genes) == 1) return(mat[valid_genes, ])
    return(colMeans(mat[valid_genes, , drop = FALSE]))
  }
  
  # Compute and assign the 6 GO Scores directly to metadata
  for (mod_name in names(go_terms_genes)) {
    seurat_obj[[mod_name]] <- calc_mean(go_terms_genes[[mod_name]], expr_matrix)
  }
  
  # Compute and assign the 3 Safety Scores directly to metadata
  for (mod_name in names(state_genes)) {
    seurat_obj[[mod_name]] <- calc_mean(state_genes[[mod_name]], expr_matrix)
  }
  
  # ====================================================================
  # 3. Extract and Format Metadata
  # ====================================================================
  meta <- seurat_obj@meta.data %>%
    mutate(plot_group = case_when(
      drug_condition == "inflammatory" ~ "IL-1β",
      drug_condition == "control" ~ "Control",
      std_combi %in% top5_combos ~ std_combi,
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(plot_group))
  
  module_mapping <- c(
    "GO_Score1" = "Chondrocyte Dev",
    "GO_Score2" = "Cartilage Condensation",
    "GO_Score3" = "Tissue Homeostasis",
    "GO_Score4" = "NO Synthase Upreg",
    "GO_Score5" = "ECM Disassembly",
    "GO_Score6" = "Inflammatory Response",
    "State_Prolif"    = "Proliferation",
    "State_Viability" = "Global Viability",
    "State_Stress"    = "Cellular Stress"
  )
  
  # ====================================================================
  # 4. Calculate Log2FC and Filter Reference
  # ====================================================================
  heat_df <- meta %>%
    select(plot_group, all_of(names(module_mapping))) %>%
    pivot_longer(cols = -plot_group, names_to = "raw_col", values_to = "score") %>%
    mutate(module = module_mapping[raw_col]) %>%
    group_by(plot_group, module) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")
  
  pc <- 0.01 # Pseudocount to prevent log2(0)
  
  # Calculate baseline, then compute Log2FC
  heat_df <- heat_df %>%
    group_by(module) %>%
    mutate(
      disease_baseline = mean_score[plot_group == "IL-1β"],
      log2fc = log2((mean_score + pc) / (disease_baseline + pc))
    ) %>%
    ungroup() %>%
    # DROP the IL-1β reference column after calculation to save space
    filter(plot_group != "IL-1β")
  
  # ====================================================================
  # 5. Factor Ordering for Visual Narrative
  # ====================================================================
  # Define X-axis order (Disease column is gone, plot top combos + healthy)
  heat_df$plot_group <- factor(heat_df$plot_group, 
                               levels = c(top5_combos, "Control"))
  
  catabolic_terms <- c("Inflammatory Response", "ECM Disassembly", "NO Synthase Upreg")
  anabolic_terms  <- c("Tissue Homeostasis", "Cartilage Condensation", "Chondrocyte Dev")
  safety_terms    <- c("Proliferation", "Global Viability", "Cellular Stress")
  
  heat_df <- heat_df %>%
    mutate(category = case_when(
      module %in% catabolic_terms ~ "Catabolic",
      module %in% anabolic_terms ~ "Anabolic",
      module %in% safety_terms ~ "Cellular State"
    ))
  
  heat_df$category <- factor(heat_df$category, 
                             levels = c("Catabolic", "Anabolic", "Cellular State"))
  heat_df$module <- factor(heat_df$module, 
                           levels = rev(c(catabolic_terms, anabolic_terms, safety_terms)))
  
  # ====================================================================
  # 6. Build the Heatmap
  # ====================================================================
  npg_red  <- pal_npg("nrc")(10)[1] 
  npg_blue <- pal_npg("nrc")(10)[4] 
  
  # Cap extreme Log2FC values to prevent color washout
  limit_val <- 3
  heat_df$log2fc_capped <- pmax(pmin(heat_df$log2fc, limit_val), -limit_val)
  
  p <- ggplot(heat_df, aes(x = plot_group, y = module, fill = log2fc_capped)) +
    
    geom_tile(color = "white", linewidth = 0.8) +
    
    scale_fill_gradient2(low = npg_blue, mid = "#F7F7F7", high = npg_red, midpoint = 0, 
                         limits = c(-limit_val, limit_val), oob = squish,
                         name = "Treatment Effect\n(Log2FC vs Disease)") +
    
    facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
    
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
      axis.text.y = element_text(face = "bold", color = "black"),
      axis.title = element_blank(),
      
      strip.placement = "outside",
      strip.background = element_rect(fill = "grey90", color = "white"),
      # Changed angle to 90 for vertical text reading bottom-to-top
      strip.text.y.left = element_text(angle = 90, face = "bold", size = 11, color = "black"),
      
      panel.grid = element_blank(), 
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10)
    )
  
  return(p)
}

fig5e_heatmap <- plot_heatmap_log2fc_compact(merged.sct, top5_combos)
ggsave("./figures/figure5/Figure5e_GO_Heatmap.pdf", plot = fig5e_heatmap, width =6, height = 5, units = "in", dpi = 300)

### figure 5 g (single additive vs combi synergy of Modules)
### the code is not available in R.
### Refer to the python code.