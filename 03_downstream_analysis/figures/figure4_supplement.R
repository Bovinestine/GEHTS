# =============================================================================
# figure4_supplement.R — Supplementary figures for Figure 4
#
# Subfigures:
#   Extended Data: 5-fold CV performance comparison across 5 ML models
#   Extended Data: 6-panel GO module score violin plots
#
# Prerequisites: config.R, utils.R, data_loading.R, scoring_methods.R
# Inputs: PRED_DIR/cv_results_per_fold.csv (written by scoring_methods.R),
#         mac.sct, sin.sct (loaded by data_loading.R)
# Outputs: output/supplementary/Sfig_cv_comparison.pdf,
#          output/supplementary/Sfig_GO_module_violins.pdf
# =============================================================================

source(here::here("03_downstream_analysis", "R", "config.R"))
source(here::here("03_downstream_analysis", "R", "utils.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
})

set.seed(SEED)

# =============================================================================
# Extended Data — 5-fold CV performance comparison (5 ML models)
# =============================================================================

cv_per_fold <- read.csv(file.path(PRED_DIR, "cv_results_per_fold.csv"))

cv_long <- melt(as.data.frame(cv_per_fold),
                id.vars       = c("Model", "Fold"),
                variable.name = "Metric",
                value.name    = "Value")

sfig_cv <- ggplot(cv_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  facet_wrap(~Metric, scales = "free_y") +
  labs(title = paste0(N_FOLDS_CV, "-Fold CV Performance by Model"),
       y = "Metric Value", x = "Model") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_SUP, "Sfig_cv_comparison.pdf"),
       plot = sfig_cv, width = 10, height = 6, dpi = 300)

# =============================================================================
# Extended Data — 6-panel GO module violin plots
# =============================================================================

merged.sct <- merge(mac.sct, sin.sct)

plot_supplementary_violins <- function(seurat_obj,
                                       target_conditions,
                                       module_list   = GO_MODULES,
                                       group_by      = "drug_condition",
                                       ref_control   = COND_CONTROL,
                                       ref_disease   = COND_DISEASE,
                                       target_labels = NULL,
                                       assay         = "SCT",
                                       slot          = "data") {

  if (is.null(target_labels)) target_labels <- target_conditions
  DefaultAssay(seurat_obj) <- assay
  mat  <- GetAssayData(seurat_obj, slot = slot)
  meta <- seurat_obj@meta.data

  conds_to_plot <- c(ref_control, ref_disease, target_conditions)
  cells         <- rownames(meta[meta[[group_by]] %in% conds_to_plot, ])
  mat_sub       <- mat[, cells]
  meta_sub      <- meta[cells, ]

  score_list <- lapply(names(module_list), function(mod) {
    genes  <- intersect(module_list[[mod]], rownames(mat_sub))
    scores <- if (length(genes) > 0) colMeans(as.matrix(mat_sub[genes, , drop = FALSE])) else rep(NA, length(cells))
    data.frame(Cell = cells, Condition = meta_sub[[group_by]],
               Module = gsub("_", " ", mod), Score = scores)
  })
  df_long <- do.call(rbind, score_list)

  plot_labels <- c("Basal", "IL-1β", target_labels)
  cond_map    <- setNames(plot_labels, conds_to_plot)
  df_long <- df_long %>%
    mutate(Condition_Clean = factor(cond_map[as.character(Condition)],
                                    levels = plot_labels)) %>%
    filter(!is.na(Score))

  base_cols   <- PALETTE_NPG[c(2, 1, 3, 8, 5, 9)]
  cond_colors <- setNames(base_cols[seq_along(plot_labels)], plot_labels)

  ggplot(df_long, aes(x = Condition_Clean, y = Score, fill = Condition_Clean)) +
    geom_violin(trim = TRUE, alpha = 0.8, scale = "width",
                color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.15, fill = "white", color = "black",
                 outlier.shape = NA, linewidth = 0.5, alpha = 0.8) +
    facet_wrap(~Module, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = cond_colors) +
    theme_bw() +
    labs(title = "GO Module Distributions", x = "", y = "Mean Gene Expression (SCT)") +
    theme(legend.position = "none",
          axis.text.x      = element_text(angle = 45, hjust = 1, face = "bold",
                                          size = 11, color = "black"),
          axis.text.y      = element_text(size = 10, color = "black"),
          axis.title.y     = element_text(face = "bold", size = 12,
                                          margin = margin(r = 10)),
          plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
          strip.background = element_rect(fill = "grey90", color = "black"),
          strip.text       = element_text(face = "bold", size = 10),
          panel.border     = element_rect(color = "black"),
          panel.grid.major.x = element_blank())
}

sfig_violin <- plot_supplementary_violins(
  merged.sct,
  target_conditions = c("rapamycin_10", "ALK5 inhibitor IV_10"),
  target_labels     = c("Rapamycin (10 μM)", "ALK5 Inhibitor (10 μM)")
)
ggsave(file.path(OUTPUT_SUP, "Sfig_GO_module_violins.pdf"),
       plot = sfig_violin, width = 7, height = 5.5, dpi = 300)

save_session_info("figure4_supplement")
