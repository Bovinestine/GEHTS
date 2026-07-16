# =============================================================================
# figure5.R — Drug combination synergy and biological validation (Figure 5)
#             + Extended Data (GO trajectory, bubble plot, stacked violins)
#
# Subfigures:
#   5c: Bliss synergy heatmap (Obs − Pred difference, Bliss independence)
#   5d: Bliss synergy vs. ML efficacy scatter (top 5 highlighted)
#   5e: ML efficacy vs. biological rescue index scatter
#   5f: Compact dot plot, 7 GO modules (Δ residual vs. IL-1β; dot size =
#       Wilcoxon/BH significance)
#   Ext: Extended Figure 5 — 9-module synergy barplot (figure4d.ipynb)
#   Ext: Extended Figure 5g — per-module synergy with significance (figure5g.ipynb)
#
# Prerequisites: config.R, utils.R, data_loading.R, prediction_score.R
# Inputs: cmb.sct, mac.sct, pred_combo_0.1_rf, pred_single_0.1_rf
# Outputs: output/main/Figure5*.pdf, output/supplementary/Sfig5ext_*.pdf
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
  library(scales)
  library(reshape2)
  library(viridis)
})

# pred_combo_0.1_rf / pred_single_0.1_rf (from prediction_score.R) are lists of
# per-cell (per-well) predictions grouped by drug/combo name, not flat
# scalars -- collapse each entry to its mean here, once, so every downstream
# consumer in this file (and generate_bliss_pdf()/split_combo_names() in
# utils.R) can keep treating them as ordinary named numeric vectors.
pred_combo_0.1_rf <- vapply(pred_combo_0.1, mean, numeric(1))
pred_single_0.1_rf <- vapply(pred_single_0.1, mean, numeric(1))

# --- Build merged object with ML efficacy labels ---
mac.sct$efficacy <- ifelse(mac.sct$drug_condition == COND_CONTROL, 1.0,
                     ifelse(mac.sct$drug_condition == COND_DISEASE,  0.0, NA))
mac.sct$drug_name1 <- mac.sct$drug_condition
mac.sct$drug_name2 <- mac.sct$drug_condition

cmb.sct$std_combi <- paste0(
  pmin(as.character(cmb.sct$drug_name1), as.character(cmb.sct$drug_name2)), "&",
  pmax(as.character(cmb.sct$drug_name1), as.character(cmb.sct$drug_name2))
)
eff_splits   <- strsplit(names(pred_combo_0.1_rf), "&")
std_eff_names <- paste0(pmin(sapply(eff_splits, `[`, 1), sapply(eff_splits, `[`, 2)), "&",
                         pmax(sapply(eff_splits, `[`, 1), sapply(eff_splits, `[`, 2)))
efficacy_map <- setNames(as.numeric(pred_combo_0.1_rf), std_eff_names)
cmb.sct$efficacy <- unname(efficacy_map[cmb.sct$std_combi])

merged.sct <- merge(mac.sct, cmb.sct, add.cell.ids = c("mac", "cmb"))
DefaultAssay(merged.sct) <- "SCT"
VariableFeatures(merged.sct) <- rownames(merged.sct)

# =============================================================================
# Figure 5c — Bliss synergy heatmap (Obs − Pred difference)
# =============================================================================

.create_triangular_matrix <- function(data, dose) {
  dose_data <- data %>% filter(Dose == dose, Drug1Name != Drug2Name)
  if (nrow(dose_data) == 0) {
    warning("No combination data for dose = ", dose)
    return(matrix(NA, 0, 0))
  }
  drugs <- sort(unique(c(dose_data$Drug1Name, dose_data$Drug2Name)))
  mat   <- matrix(NA, nrow = length(drugs), ncol = length(drugs),
                  dimnames = list(drugs, drugs))
  for (i in seq_len(nrow(dose_data))) {
    d1 <- dose_data$Drug1Name[i]; d2 <- dose_data$Drug2Name[i]
    if (which(drugs == d1) <= which(drugs == d2)) mat[d1, d2] <- dose_data$MeanProbability[i]
    else mat[d2, d1] <- dose_data$MeanProbability[i]
  }
  mat
}

.create_predicted_matrix <- function(sng_vec, cmb_df) {
  drugs <- sort(unique(c(cmb_df$Drug1Name, cmb_df$Drug2Name)))
  mat   <- matrix(NA, nrow = length(drugs), ncol = length(drugs),
                  dimnames = list(drugs, drugs))
  for (d1 in drugs) for (d2 in drugs) {
    if (which(drugs == d1) < which(drugs == d2)) {
      mat[d1, d2] <- bliss_predict(sng_vec[d1], sng_vec[d2])
    }
  }
  mat
}

.find_shared_limits <- function(matrices) {
  vals <- unlist(lapply(matrices, function(m) m[is.finite(m)]))
  c(min(vals, na.rm = TRUE), max(vals, na.rm = TRUE))
}

.plot_heatmap_diff <- function(mat, title, lim) {
  df    <- melt(mat, na.rm = FALSE, varnames = c("Drug1", "Drug2"))
  drugs <- rownames(mat)
  ggplot(df, aes(x = Drug1, y = Drug2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "#0000CC", high = "#CC0000", mid = "white",
                         midpoint = 0, limits = lim, na.value = NA) +
    labs(title = title, fill = "Synergy") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(angle = 45, hjust = 1))
}

#' Bliss Synergy Difference Heatmap (Observed − Predicted)
#'
#' Computes the Bliss independence prediction from single-agent efficacies,
#' then plots the residual (observed − predicted) as a diverging heatmap.
#' Positive values (red) indicate synergy; negative values (blue) indicate
#' antagonism relative to the independence model.
#'
#' @param prob_sng Named numeric vector of single-agent efficacy scores.
#' @param prob_cmb Named numeric vector of combination efficacy scores.
#' @param dose Character label for the dose (e.g., "0.1").
#' @param pdf_file Output PDF path.
#' @param pdf_width PDF width in inches (default 5).
#' @param pdf_height PDF height in inches (default 4.2).
#' @export
generate_bliss_pdf <- function(prob_sng, prob_cmb, dose,
                               pdf_file   = file.path(OUTPUT_MAIN, "Figure5c_bliss.pdf"),
                               pdf_width  = 5,
                               pdf_height = 4.2) {

  `%notin%` <- Negate(`%in%`)
  prob_sng <- prob_sng[names(prob_sng) %notin% c(COND_CONTROL, COND_DISEASE)]
  prob_cmb <- prob_cmb[names(prob_cmb) %notin% c(COND_CONTROL, COND_DISEASE)]

  cmb_df   <- split_combo_names(prob_cmb, dose)
  obs_mat  <- .create_triangular_matrix(cmb_df, dose)
  pred_mat <- .create_predicted_matrix(prob_sng, cmb_df)
  diff_mat <- obs_mat - pred_mat

  bliss_lim <- .find_shared_limits(list(diff_mat))
  p <- .plot_heatmap_diff(diff_mat,
                          paste("Synergy (Obs − Pred) —", dose, "μM"),
                          bliss_lim)
  pdf(pdf_file, width = pdf_width, height = pdf_height)
  print(p)
  dev.off()
  invisible(pdf_file)
}

generate_bliss_pdf(pred_single_0.1_rf, pred_combo_0.1_rf, dose = "0.1",
                   pdf_file = file.path(OUTPUT_MAIN, "Figure5c_bliss_rf_0.1uM.pdf"))

# =============================================================================
# Figure 5d — Bliss Synergy vs. ML Efficacy
# =============================================================================

#' Scatter Plot: Bliss Synergy vs. ML Efficacy for All Combinations
#'
#' Identifies the top N combinations by rank-sum of efficacy and synergy,
#' highlights them, and returns both the plot and the top-candidate data.frame.
#'
#' @param cmb_rf Named numeric vector of ML efficacy scores for combinations.
#' @param sng_rf Named numeric vector of ML efficacy scores for single agents.
#' @param top_n Number of top hits to highlight (default 5).
#' @return A list with elements plot (ggplot2) and top_candidates (data.frame).
#' @export
plot_bliss_synergy <- function(cmb_rf, sng_rf, top_n = 5) {

  if (is.null(names(cmb_rf)) || is.null(names(sng_rf))) {
    stop("Both cmb_rf and sng_rf must be named numeric vectors.")
  }

  splits  <- strsplit(names(cmb_rf), "&")
  drug1   <- trimws(sapply(splits, `[`, 1))
  drug2   <- trimws(sapply(splits, `[`, 2))
  eff_obs <- as.numeric(cmb_rf)
  eff_A   <- as.numeric(sng_rf[drug1])
  eff_B   <- as.numeric(sng_rf[drug2])

  df <- data.frame(
    std_combi  = paste0(pmin(drug1, drug2), "&", pmax(drug1, drug2)),
    drug1      = drug1,
    drug2      = drug2,
    E_obs      = eff_obs,
    Bliss_pred = bliss_predict(eff_A, eff_B),
    Synergy    = eff_obs - bliss_predict(eff_A, eff_B)
  ) %>%
    filter(!is.na(Synergy), drug1 != drug2) %>%
    group_by(std_combi, drug1, drug2) %>%
    summarise(E_obs   = mean(E_obs,   na.rm = TRUE),
              Synergy = mean(Synergy, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(Rank_E   = rank(-E_obs,   ties.method = "min"),
           Rank_S   = rank(-Synergy, ties.method = "min"),
           Rank_Sum = Rank_E + Rank_S) %>%
    arrange(Rank_Sum, desc(E_obs))

  highlight_df <- head(df, top_n) %>%
    mutate(display_name = paste0(drug1, " + ", drug2))

  p <- ggplot(df, aes(x = E_obs, y = Synergy)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey50", linewidth = 0.8) +
    geom_point(shape = 21, size = 3, fill = "grey70",
               color = "black", stroke = 0.3, alpha = 0.85) +
    geom_point(data = highlight_df, aes(fill = std_combi),
               shape = 21, size = 5, color = "black", stroke = 1.2) +
    scale_fill_npg() +
    geom_text_repel(data = highlight_df, aes(label = display_name),
                    box.padding = 1.0, point.padding = 0.5,
                    segment.color = "black", segment.size = 0.5,
                    fontface = "bold", size = 4.0, max.overlaps = Inf) +
    theme_classic(base_size = 14) +
    labs(x = "ML Efficacy Score", y = "Bliss Synergy Score") +
    theme(axis.title = element_text(face = "bold"),
          legend.position = "none")

  list(plot = p, top_candidates = highlight_df)
}

synergy_results <- plot_bliss_synergy(pred_combo_0.1_rf, pred_single_0.1_rf, top_n = 5)
top5_combos     <- synergy_results$top_candidates$std_combi
target_singles  <- unique(c(synergy_results$top_candidates$drug1,
                             synergy_results$top_candidates$drug2))

ggsave(file.path(OUTPUT_MAIN, "Figure5d_bliss_synergy_scatter.pdf"),
       plot = synergy_results$plot, width = 5, height = 5, dpi = 300)

# =============================================================================
# Figure 5e — ML Efficacy vs. Biological Rescue Index
# =============================================================================

#' Scatter: ML Efficacy vs. Biological Rescue Index (Anabolic − Catabolic)
#'
#' @param seurat_obj Merged Seurat object with combination and control cells.
#' @param prob_cmb Named numeric vector of ML combination efficacy scores.
#' @param top5_combos Character vector of top-5 combination names.
#' @return A ggplot2 object.
#' @export
plot_biological_alignment <- function(seurat_obj, prob_cmb, top5_combos) {

  DefaultAssay(seurat_obj) <- "SCT"
  exp_mat <- GetAssayData(seurat_obj, layer = "data")
  seurat_obj$Anabolic_Score  <- colMeans(as.matrix(exp_mat[GENES_ANABOLIC,  , drop = FALSE]))
  seurat_obj$Catabolic_Score <- colMeans(as.matrix(exp_mat[GENES_CATABOLIC, , drop = FALSE]))

  meta <- seurat_obj@meta.data %>%
    mutate(
      Rescue_Index    = Anabolic_Score - Catabolic_Score,
      plot_id = case_when(
        drug_condition %in% c(COND_CONTROL, COND_DISEASE) ~ drug_condition,
        drug_name1 != drug_name2 ~ paste0(
          pmin(as.character(drug_name1), as.character(drug_name2)), "&",
          pmax(as.character(drug_name1), as.character(drug_name2))
        ),
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(plot_id))

  bio_df <- meta %>%
    group_by(plot_id) %>%
    summarise(Rescue_Index = mean(Rescue_Index, na.rm = TRUE), .groups = "drop")

  eff_splits    <- strsplit(names(prob_cmb), "&")
  std_eff_names <- paste0(pmin(trimws(sapply(eff_splits, `[`, 1)),
                               trimws(sapply(eff_splits, `[`, 2))), "&",
                           pmax(trimws(sapply(eff_splits, `[`, 1)),
                                trimws(sapply(eff_splits, `[`, 2))))
  eff_df <- data.frame(plot_id  = std_eff_names,
                        Efficacy = as.numeric(prob_cmb)) %>%
    group_by(plot_id) %>%
    summarise(Efficacy = mean(Efficacy, na.rm = TRUE), .groups = "drop")

  all_eff  <- bind_rows(eff_df,
                         data.frame(plot_id  = c(COND_CONTROL, COND_DISEASE),
                                    Efficacy = c(1.0, 0.0)))
  plot_df  <- inner_join(bio_df, all_eff, by = "plot_id") %>%
    mutate(
      category = case_when(
        plot_id == COND_CONTROL ~ "Control",
        plot_id == COND_DISEASE ~ "IL-1β",
        plot_id %in% top5_combos ~ "Top 5 Hit",
        TRUE ~ "Other Combinations"
      ),
      label = case_when(
        plot_id == COND_CONTROL ~ "Control",
        plot_id == COND_DISEASE ~ "IL-1β",
        plot_id %in% top5_combos ~ gsub("&", " + ", plot_id),
        TRUE ~ ""
      )
    )

  bg_df   <- plot_df %>% filter(category == "Other Combinations")
  top_df  <- plot_df %>% filter(category == "Top 5 Hit")
  ctrl_df <- plot_df %>% filter(category %in% c("Control", "IL-1β")) %>%
    mutate(fill_color = ifelse(plot_id == COND_CONTROL, "white", "black"))

  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.8) +
    geom_point(data = bg_df, aes(x = Efficacy, y = Rescue_Index),
               color = "grey70", size = 2.5, alpha = 0.6) +
    geom_point(data = top_df, aes(x = Efficacy, y = Rescue_Index, fill = plot_id),
               shape = 21, size = 5, color = "black", stroke = 1.2) +
    geom_point(data = ctrl_df, aes(x = Efficacy, y = Rescue_Index),
               shape = 23, fill = ctrl_df$fill_color,
               color = "black", size = 5, stroke = 1.2) +
    geom_text_repel(data = bind_rows(top_df, ctrl_df),
                    aes(x = Efficacy, y = Rescue_Index, label = label),
                    box.padding = 1.2, point.padding = 0.5,
                    fontface = "bold", size = 4.5, max.overlaps = Inf) +
    scale_fill_npg() +
    theme_classic(base_size = 14) +
    labs(x = "ML Efficacy Score",
         y = expression(Sigma * "Anabolic" ~ "-" ~ Sigma * "Catabolic")) +
    theme(legend.position = "none", axis.title = element_text(face = "bold"))
}

fig5e_plot <- plot_biological_alignment(merged.sct, pred_combo_0.1_rf, top5_combos)
ggsave(file.path(OUTPUT_MAIN, "Figure5e_biological_alignment.pdf"),
       plot = fig5e_plot, width = 5, height = 5, dpi = 300)


# =============================================================================
# Figure 5f — Compact dot plot, all 7 GO modules (residual-based)
#
# plot_dotplot_diff() is the shared, module-list-parameterized function
# defined once in utils.R (also used by Fig 4e) -- see utils.R section 3b.
# Color is the plain difference of module residual means from the IL-1beta
# baseline (log2FC is invalid once scores go negative); size is Wilcoxon/BH
# significance of that shift.
# =============================================================================

# unname() is required here: case_when()'s combo branch pulls from the named
# vector merged.sct$std_combi, so the combined result carries a names
# attribute ("" for the literal-string branches). Seurat's `$<-` matches
# metadata by name when one is present, so an un-stripped vector silently
# drops every "" -named (i.e. Control/IL-1b) entry to NA on assignment.
merged.sct$plot_group <- unname(case_when(
  merged.sct$drug_condition == COND_DISEASE ~ "IL-1β",
  merged.sct$drug_condition == COND_CONTROL ~ "Control",
  merged.sct$std_combi %in% top5_combos ~ merged.sct$std_combi,
  TRUE ~ NA_character_
))
fig5f_seurat_obj      <- merged.sct[, !is.na(merged.sct$plot_group)]
fig5f_plot_conditions <- c(top5_combos, "Control")
fig5f_condition_labels <- setNames(c(gsub("&", " + ", top5_combos), "Control"), fig5f_plot_conditions)

fig5f_dotplot_significance <- plot_dotplot_diff(
  fig5f_seurat_obj,
  plot_conditions   = fig5f_plot_conditions,
  condition_labels  = fig5f_condition_labels,
  group_by          = "plot_group",
  ref_disease       = "IL-1β"
)
ggsave(file.path(OUTPUT_MAIN, "Figure5f_dotplot_significance.pdf"),
       plot = fig5f_dotplot_significance, width = 7, height = 5.5, dpi = 300)

save_session_info("figure5")
