# =============================================================================
# figure4.R — Drug efficacy and phenotypic characterization (Figure 4)
#             + Extended Data (6-module violin plots)
#
# Subfigures:
#   4b: Efficacy bar plot (mean ± SEM per drug, two doses)
#   4c: 2D phenotypic screening scatter (anabolic vs. catabolic)
#   4d: Highlighted UMAPs — pamapimod and XAV939 at dose 10 μM
#   4e: Dumbbell plot of differential co-expression (disease→treatment→healthy)
#   4f: Radar charts of 6 GO module scores per drug
#
# Prerequisites: config.R, utils.R, data_loading.R, prediction_score.R
# Inputs: sin.sct, cmb.sct, mac.sct, prob_sng10_rf, prob_sng0.1_rf,
#         prob_cmb0.1_rf (all from data_loading / prediction_score)
# Outputs: output/main/Figure4*.pdf
# =============================================================================

source(here::here("03_downstream_analysis", "R", "config.R"))
source(here::here("03_downstream_analysis", "R", "utils.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(ggsci)
  library(fmsb)
  library(scales)
  library(svglite)
})

# --- Data preparation ---
sin0.1.sct <- subset(sin.sct, subset = dose == 0.1)
cmb0.1.sct <- subset(cmb.sct, subset = dose1 == 0.1)

# Standardize combination drug name: Drug1&Drug2, swapping XAV to second position
cmb0.1.sct$drug_name <- ifelse(
  startsWith(cmb0.1.sct$drug_name1, "XAV"),
  paste(cmb0.1.sct$drug_name2, cmb0.1.sct$drug_name1, sep = "&"),
  paste(cmb0.1.sct$drug_name1, cmb0.1.sct$drug_name2, sep = "&")
)

mac.sct$drug_name <- mac.sct$drug_condition
mac_sin.sct   <- merge(sin.sct, mac.sct)
mac_sin0.1.sct <- merge(sin0.1.sct, mac.sct)
mac_cmb0.1.sct <- merge(cmb0.1.sct, mac.sct)

merged.sct <- merge(mac.sct, sin.sct)

# =============================================================================
# Figure 4b — Efficacy bar plot
# =============================================================================

#' Grouped Bar Plot of Drug Efficacy Probabilities (Two Doses)
#'
#' @param prob_sng10_raw Named list of per-cell RF predictions at 10 μM.
#' @param prob_sng0.1_raw Named list of per-cell RF predictions at 0.1 μM.
#' @return A ggplot2 bar plot.
#' @export
plot_efficacy_barplot <- function(prob_sng10_raw, prob_sng0.1_raw) {

  bind_dose <- function(preds, dose_label) {
    df <- stack(preds)
    colnames(df) <- c("Probability", "Drug")
    df$Dose <- dose_label
    df
  }
  raw_data <- rbind(bind_dose(prob_sng10_raw, "10"),
                    bind_dose(prob_sng0.1_raw, "0.1"))
  raw_data$Dose <- factor(raw_data$Dose, levels = c("0.1", "10"))
  raw_data <- raw_data[!raw_data$Drug %in% c(COND_CONTROL, COND_DISEASE), ]

  summary_data <- raw_data %>%
    group_by(Drug, Dose) %>%
    summarise(Mean_Probability = mean(Probability, na.rm = TRUE),
              SE_Probability   = sd(Probability, na.rm = TRUE) / sqrt(n()),
              .groups = "drop")

  ggplot(summary_data, aes(x = Drug, y = Mean_Probability, fill = Dose)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             color = "black", linewidth = 0.3, alpha = 0.9) +
    geom_errorbar(aes(ymin = Mean_Probability - SE_Probability,
                      ymax = Mean_Probability + SE_Probability),
                  position = position_dodge(width = 0.8),
                  width = 0.25, alpha = 0.7, linewidth = 0.5) +
    scale_fill_manual(values = c("0.1" = COL_DOSE_LOW, "10" = COL_DOSE_HIGH)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          panel.grid.major.x = element_blank()) +
    labs(fill = "Dose (μM)", x = "Drug", y = "Mean Efficacy Probability")
}

prob_list_raw <- predict_all_doses(rf_model, sin.sct, cmb.sct, save_dir = NULL)
fig_4b <- plot_efficacy_barplot(prob_list_raw$pred_single_10,
                                prob_list_raw$pred_single_0.1)
ggsave(file.path(OUTPUT_MAIN, "Figure4b_efficacy_barplot.pdf"),
       plot = fig_4b, width = 7, height = 5.2, dpi = 300)

# =============================================================================
# Figure 4c — 2D phenotypic screening scatter
# =============================================================================

#' 2D Phenotypic Scatter: Anabolic vs. Catabolic Gene Sum
#'
#' @param seurat_obj Merged Seurat object with drug_name, dose, drug_condition.
#' @param genes_anabolic Anabolic gene panel.
#' @param genes_inflammatory Catabolic gene panel.
#' @param topn Number of top/bottom hits to label (default 5).
#' @return A ggplot2 object.
#' @export
plot_phenotypic_screening <- function(seurat_obj,
                                      genes_anabolic     = GENES_ANABOLIC,
                                      genes_inflammatory = GENES_CATABOLIC,
                                      topn = 5) {

  DefaultAssay(seurat_obj) <- "SCT"
  exp_mat <- GetAssayData(seurat_obj, slot = "data")
  seurat_obj$Anabolic_Score  <- colMeans(as.matrix(exp_mat[genes_anabolic, ,    drop = FALSE]))
  seurat_obj$Catabolic_Score <- colMeans(as.matrix(exp_mat[genes_inflammatory,, drop = FALSE]))

  plot_data <- seurat_obj@meta.data %>%
    group_by(drug_name, dose, drug_condition) %>%
    summarise(mean_ana  = mean(Anabolic_Score,  na.rm = TRUE),
              se_ana    = sd(Anabolic_Score,  na.rm = TRUE) / sqrt(n()),
              mean_cata = mean(Catabolic_Score, na.rm = TRUE),
              se_cata   = sd(Catabolic_Score, na.rm = TRUE) / sqrt(n()),
              .groups = "drop") %>%
    mutate(AC_Ratio = mean_ana / mean_cata)

  top_hits    <- plot_data %>% top_n(topn,  AC_Ratio) %>% pull(drug_condition)
  bottom_hits <- plot_data %>% top_n(topn, -AC_Ratio) %>% pull(drug_condition)
  controls    <- grep("control|inflammatory", plot_data$drug_condition,
                      value = TRUE, ignore.case = TRUE)

  plot_data <- plot_data %>%
    mutate(
      Label_Text = ifelse(drug_condition %in% c(top_hits, bottom_hits, controls),
                          drug_name, NA),
      Group = case_when(
        grepl("control",     drug_condition, ignore.case = TRUE) ~ "Basal",
        grepl("inflammatory", drug_condition, ignore.case = TRUE) ~ "IL-1β",
        dose == 10  | dose == "10"  ~ "10 µM",
        dose == 0.1 | dose == "0.1" ~ "0.1 µM",
        TRUE ~ "Other"
      ),
      Group = factor(Group, levels = c("Basal", "IL-1β", "0.1 µM", "10 µM"))
    )

  cond_colors <- c(Basal = COL_CONTROL, `IL-1β` = COL_DISEASE,
                   `0.1 µM` = COL_DOSE_LOW, `10 µM` = COL_DOSE_HIGH)

  ggplot(plot_data, aes(x = mean_cata, y = mean_ana)) +
    geom_errorbar(aes(ymin = mean_ana - se_ana, ymax = mean_ana + se_ana),
                  color = "grey80", width = 0) +
    geom_errorbarh(aes(xmin = mean_cata - se_cata, xmax = mean_cata + se_cata),
                   color = "grey80", height = 0) +
    geom_vline(xintercept = mean(range(plot_data$mean_cata, na.rm = TRUE)),
               linetype = "dotted", color = "grey60") +
    geom_hline(yintercept = mean(range(plot_data$mean_ana, na.rm = TRUE)),
               linetype = "dotted", color = "grey60") +
    geom_line(aes(group = drug_name), color = "grey60", linewidth = 0.5, alpha = 0.6) +
    geom_point(aes(fill = Group), shape = 21, size = 6, alpha = 0.7, color = "black") +
    scale_fill_manual(name = "Condition", values = cond_colors) +
    geom_text_repel(aes(label = Label_Text), size = 3.5, fontface = "bold",
                    box.padding = 0.6, max.overlaps = 50, min.segment.length = 0) +
    theme_bw() +
    labs(x = "Catabolic gene sum (a.u.)", y = "Anabolic gene sum (a.u.)") +
    theme(panel.grid.minor = element_blank(),
          axis.title       = element_text(face = "bold"),
          legend.position  = "right")
}

fig_4c <- plot_phenotypic_screening(mac_sin.sct, topn = 13)
ggsave(file.path(OUTPUT_MAIN, "Figure4c_phenotypic_scatter.pdf"),
       plot = fig_4c, width = 7, height = 5.2, dpi = 300)

# =============================================================================
# Figure 4d — Highlighted UMAPs: pamapimod (index 1) and XAV939 (index 13)
# Requires sin.sct10 created by figure3.R
# =============================================================================

plot_highlighted_umaps(
  sin.sct10,
  plot_indices = c(1, 13),
  output_pdf   = file.path(OUTPUT_MAIN, "Figure4d_highlighted_UMAP_pamapimod_XAV.pdf"),
  ncol         = 2,
  pdf_width    = 6,
  pdf_height   = 3.5
)

# =============================================================================
# Figure 4e (co-expression) — Dumbbell plot
# =============================================================================

#' Directional Co-expression Dumbbell Plot (Disease → Treatment → Healthy)
#'
#' Automatically selects the top-N most divergent gene pairs and draws
#' arrows showing how Pearson correlation changes across three conditions.
#'
#' @param seurat_obj A Seurat object.
#' @param group_by Metadata column with condition labels (default "drug_condition").
#' @param cond_ref Disease reference condition (default COND_DISEASE).
#' @param cond_treat Treatment condition to evaluate.
#' @param cond_healthy Healthy baseline (default COND_CONTROL).
#' @param top_n Number of top gene pairs to display (default 10).
#' @param assay Assay to use (default "SCT").
#' @param slot Data slot (default "data").
#' @return A ggplot2 dumbbell plot.
#' @export
plot_automated_coexpression <- function(seurat_obj,
                                        group_by    = "drug_condition",
                                        cond_ref    = COND_DISEASE,
                                        cond_treat  = "rapamycin_10",
                                        cond_healthy = COND_CONTROL,
                                        top_n       = 10,
                                        assay       = "SCT",
                                        slot        = "data",
                                        color_ref    = COL_DISEASE,
                                        color_treat  = PALETTE_NPG[4],
                                        color_healthy = PALETTE_NPG[3]) {

  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop("Column '", group_by, "' not found in metadata.")
  }
  all_conds <- c(cond_ref, cond_treat, cond_healthy)
  if (!all(all_conds %in% seurat_obj@meta.data[[group_by]])) {
    stop("One or more conditions not found in '", group_by, "' metadata.")
  }

  DefaultAssay(seurat_obj) <- assay
  mat  <- GetAssayData(seurat_obj, slot = slot)
  meta <- seurat_obj@meta.data

  get_mat <- function(cond) {
    cells <- rownames(meta[meta[[group_by]] == cond, ])
    as.matrix(mat[, cells])
  }
  mat_ref     <- get_mat(cond_ref)
  mat_treat   <- get_mat(cond_treat)
  mat_healthy <- get_mat(cond_healthy)

  valid <- Reduce(intersect, lapply(list(mat_ref, mat_treat, mat_healthy),
                                    function(m) names(which(apply(m, 1, var) > 0))))
  mat_ref     <- mat_ref[valid, ]
  mat_treat   <- mat_treat[valid, ]
  mat_healthy <- mat_healthy[valid, ]

  cor_ref     <- cor(t(mat_ref),    method = "pearson")
  cor_treat   <- cor(t(mat_treat),  method = "pearson")
  cor_healthy <- cor(t(mat_healthy), method = "pearson")

  idx <- which(upper.tri(cor_ref), arr.ind = TRUE)
  df_all <- data.frame(
    Gene1     = rownames(cor_ref)[idx[, 1]],
    Gene2     = colnames(cor_ref)[idx[, 2]],
    R_Ref     = cor_ref[idx],
    R_Treat   = cor_treat[idx],
    R_Healthy = cor_healthy[idx],
    stringsAsFactors = FALSE
  )

  df_plot <- df_all %>%
    drop_na() %>%
    mutate(GenePair = paste0(Gene1, " - ", Gene2),
           Delta_r  = abs(R_Treat - R_Ref)) %>%
    arrange(desc(Delta_r)) %>%
    slice_head(n = top_n) %>%
    arrange(Delta_r)
  df_plot$GenePair <- factor(df_plot$GenePair, levels = df_plot$GenePair)

  axis_limit <- ceiling(max(abs(c(df_plot$R_Ref, df_plot$R_Treat,
                                  df_plot$R_Healthy)), na.rm = TRUE) * 10) / 10 + 0.1

  color_map <- setNames(c(color_ref, color_treat, color_healthy),
                        c(cond_ref, cond_treat, cond_healthy))

  ggplot(df_plot) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_segment(aes(x = R_Ref, xend = R_Treat, y = GenePair, yend = GenePair),
                 color = "grey40", linewidth = 1.2,
                 arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
    geom_point(aes(x = R_Healthy, y = GenePair, fill = !!cond_healthy),
               size = 4, shape = 23, color = "black", stroke = 1, alpha = 0.85) +
    geom_point(aes(x = R_Ref,   y = GenePair, fill = !!cond_ref),
               size = 4, shape = 21, color = "black", stroke = 1) +
    geom_point(aes(x = R_Treat, y = GenePair, fill = !!cond_treat),
               size = 4, shape = 21, color = "black", stroke = 1) +
    scale_fill_manual(values = color_map) +
    scale_x_continuous(limits = c(-axis_limit, axis_limit),
                       breaks = round(seq(-axis_limit, axis_limit, length.out = 5), 2)) +
    theme_bw() +
    labs(title    = "Differential Co-expression",
         subtitle  = paste(cond_ref, "→", cond_healthy, "via", cond_treat),
         x = "Pearson r", y = "") +
    theme(legend.position = "top",
          legend.title    = element_blank(),
          axis.text.y     = element_text(size = 12, face = "bold.italic"),
          axis.text.x     = element_text(size = 11),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank())
}

fig_4e_dumbbell <- plot_automated_coexpression(merged.sct)
ggsave(file.path(OUTPUT_MAIN, "Figure4e_coexpression_dumbbell.pdf"),
       plot = fig_4e_dumbbell, width = 7, height = 5.5, dpi = 300)

# =============================================================================
# Figure 4f — Radar charts of 6 GO module scores
# =============================================================================

#' Plot Phenotypic Shift Radar Chart for One Drug Condition
#'
#' Compares a target drug against healthy (control) and diseased (inflammatory)
#' baselines using a universal scale computed across all conditions.
#'
#' @param seurat_obj A Seurat object.
#' @param target_condition Exact metadata string for the drug to plot.
#' @param radar_scale Output of calculate_radar_scale() (from utils.R).
#' @param module_list Named list of gene modules (default: GO_MODULES).
#' @param group_by Metadata column (default "drug_condition").
#' @param ref_control Healthy baseline label (default COND_CONTROL).
#' @param ref_disease Disease baseline label (default COND_DISEASE).
#' @param target_label Display label for the legend.
#' @param assay Assay (default "SCT").
#' @param slot Data slot (default "counts").
#' @return NULL — draws to the active graphics device.
#' @export
plot_radar_fast <- function(seurat_obj,
                            target_condition,
                            radar_scale,
                            module_list  = GO_MODULES,
                            group_by     = "drug_condition",
                            ref_control  = COND_CONTROL,
                            ref_disease  = COND_DISEASE,
                            target_label = NULL,
                            assay        = "SCT",
                            slot         = "counts") {

  if (!target_condition %in% seurat_obj@meta.data[[group_by]]) {
    stop("Condition '", target_condition, "' not found in metadata.")
  }
  if (is.null(target_label)) target_label <- target_condition

  DefaultAssay(seurat_obj) <- assay
  mat  <- GetAssayData(seurat_obj, slot = slot)
  meta <- seurat_obj@meta.data

  get_metrics <- function(cond) {
    cells <- rownames(meta[meta[[group_by]] == cond, ])
    sapply(module_list, function(genes) {
      avail <- intersect(genes, rownames(mat))
      if (length(avail) > 0) mean(rowMeans(mat[avail, cells, drop = FALSE])) else NA
    })
  }

  conds <- c(ref_control, ref_disease, target_condition)
  met   <- t(sapply(conds, get_metrics))

  min_v <- radar_scale$min
  max_v <- radar_scale$max
  dif_v <- max_v - min_v
  dif_v[dif_v == 0] <- 1

  norm_df <- sweep(sweep(met, 2, min_v, "-"), 2, dif_v, "/")

  plot_df <- as.data.frame(rbind(rep(1, ncol(norm_df)),
                                  rep(0, ncol(norm_df)),
                                  norm_df))
  colnames(plot_df) <- gsub("_", "\n", colnames(plot_df))

  radarchart(plot_df, axistype = 1,
             pcol  = c(COL_CONTROL, COL_DISEASE, PALETTE_NPG[3]),
             pfcol = c(alpha(COL_CONTROL, 0.1), alpha(COL_DISEASE, 0.1),
                       alpha(PALETTE_NPG[3], 0.5)),
             plwd  = c(2, 2, 3), plty = c(2, 2, 1),
             cglcol = "grey70", cglty = 1, axislabcol = "grey50",
             caxislabels = seq(0, 1, 0.25), cglwd = 0.8, vlcex = 0.9,
             title = paste("Phenotypic Shift:", target_label))

  legend("topright",
         legend = c("Control", "Inflammatory", target_label),
         col    = c(COL_CONTROL, COL_DISEASE, PALETTE_NPG[3]),
         lty    = c(2, 2, 1), lwd = c(2, 2, 3), bty = "n",
         pch = 20, pt.cex = 2, inset = c(-0.1, 0))
}

# Compute universal scale across all conditions (run once)
universal_scale <- calculate_radar_scale(merged.sct, module_list = GO_MODULES)
save(universal_scale, file = file.path(DATA_DIR, "uni_scale.Rdata"))

## Radar chart (fmsb base-R graphics) ----
svglite(file.path(OUTPUT_MAIN, "Figure4f_radar_top_vs_bottom.svg"),
        width = 12, height = 5.5)
old_par <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 5), xpd = TRUE)
plot_radar_fast(merged.sct, "rapamycin_10",        universal_scale,
                target_label = "Rapamycin (10 μM)")
plot_radar_fast(merged.sct, "ALK5 inhibitor IV_10", universal_scale,
                target_label = "ALK5 Inhibitor (10 μM)")
dev.off()
par(old_par)

save_session_info("figure4")
