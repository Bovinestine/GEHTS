# =============================================================================
# prediction_score.R — Generate per-cell ML efficacy predictions
#
# Unified prediction interface for all five model types. Single and combination
# drug predictions share the same core function; a `combo` argument controls
# how drug names are assembled from metadata.
#
# Prerequisites: config.R, data_loading.R, scoring_methods.R (models must exist)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(randomForest)
  library(xgboost)
  library(glmnet)
  library(kernlab)
  library(ggplot2)
  library(viridis)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# 1. Core prediction function
# -----------------------------------------------------------------------------

#' Predict Efficacy Scores Per Cell
#'
#' Extracts SCT counts from a Seurat object, aligns features to the model,
#' runs prediction, and groups raw scores by drug name into a named list.
#'
#' @param seurat_obj A Seurat object with SCT assay.
#' @param model A trained model (randomForest, cv.glmnet, xgb.Booster, svm, gausspr).
#' @param combo Logical. If TRUE, assembles drug names from drug_name1 & drug_name2
#'   columns (combination experiments). If FALSE, uses drug_name column.
#' @return A named list of numeric vectors (raw per-cell predictions per drug).
#' @export
predict_efficacy <- function(seurat_obj, model, combo = FALSE) {

  drug_data <- t(as.matrix(GetAssayData(seurat_obj, assay = "SCT", slot = "counts")))

  # Feature alignment
  if (inherits(model, "xgb.Booster")) {
    if (!is.null(model$feature_names)) {
      missing_feats <- setdiff(model$feature_names, colnames(drug_data))
      if (length(missing_feats) > 0) {
        warning(length(missing_feats), " model feature(s) absent in data: ",
                paste(head(missing_feats, 5), collapse = ", "),
                if (length(missing_feats) > 5) "..." else "")
      }
      drug_data <- drug_data[, model$feature_names, drop = FALSE]
    }
    preds <- as.numeric(predict(model, newdata = drug_data))

  } else if (inherits(model, "cv.glmnet")) {
    model_feats <- rownames(coef(model, s = "lambda.min"))[-1]
    shared      <- intersect(model_feats, colnames(drug_data))
    drug_data   <- drug_data[, shared, drop = FALSE]
    preds <- as.numeric(predict(model, newx = drug_data,
                                s = "lambda.min", type = "response"))

  } else if (inherits(model, "randomForest")) {
    preds <- as.numeric(predict(model, newdata = drug_data))

  } else if (inherits(model, "gausspr")) {
    preds <- as.numeric(predict(model, newdata = as.data.frame(drug_data)))

  } else if (inherits(model, "svm")) {
    preds <- as.numeric(predict(model, newdata = drug_data))

  } else {
    stop("Unsupported model class: ", paste(class(model), collapse = ", "))
  }

  # Drug-name grouping
  if (combo) {
    drug_names <- paste0(seurat_obj@meta.data$drug_name1, "&",
                         seurat_obj@meta.data$drug_name2)
  } else {
    drug_names <- seurat_obj@meta.data$drug_name
  }

  split(preds, drug_names)
}

# -----------------------------------------------------------------------------
# 2. Wrapper: predict all doses (single and combo)
# -----------------------------------------------------------------------------

#' Predict Efficacy for All Standard Dose Subsets
#'
#' Subsets the single-drug and combo Seurat objects by dose, runs
#' predict_efficacy() for each, and optionally saves the results.
#'
#' @param model A trained model object.
#' @param single_seurat Seurat object for single-drug experiments.
#' @param combo_seurat Seurat object for combination experiments.
#' @param save_dir Directory for saving .Rdata prediction files (default: PRED_DIR).
#' @return A named list with pred_single_10, pred_single_0.1, pred_combo_0.1.
#' @export
predict_all_doses <- function(model,
                              single_seurat,
                              combo_seurat,
                              save_dir = PRED_DIR) {

  single_10  <- subset(single_seurat, subset = dose == 10)
  single_0.1 <- subset(single_seurat, subset = dose == 0.1)
  combo_0.1  <- subset(combo_seurat,  subset = dose1 == 0.1)

  pred_single_10  <- predict_efficacy(single_10,  model, combo = FALSE)
  pred_single_0.1 <- predict_efficacy(single_0.1, model, combo = FALSE)
  pred_combo_0.1  <- predict_efficacy(combo_0.1,  model, combo = TRUE)

  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    save(pred_single_10,  file = file.path(save_dir, "prob_sng10_rf.Rdata"))
    save(pred_single_0.1, file = file.path(save_dir, "prob_sng0.1_rf.Rdata"))
    save(pred_combo_0.1,  file = file.path(save_dir, "prob_cmb0.1_rf.Rdata"))
    message("Predictions saved to: ", save_dir)
  }

  list(pred_single_10  = pred_single_10,
       pred_single_0.1 = pred_single_0.1,
       pred_combo_0.1  = pred_combo_0.1)
}

# -----------------------------------------------------------------------------
# 3. Heatmap of mean predictions by drug and dose
# -----------------------------------------------------------------------------

#' Heatmap of Mean Prediction Probability by Drug and Dose
#'
#' @param pred_list Output from predict_all_doses().
#' @return A ggplot2 heatmap.
#' @export
plot_prediction_heatmap <- function(pred_list) {
  summarize_preds <- function(preds, dose_label) {
    data.frame(
      Drug            = names(preds),
      Dose            = dose_label,
      MeanProbability = sapply(preds, mean, na.rm = TRUE)
    )
  }

  hmap_data <- rbind(
    summarize_preds(pred_list$pred_single_10,  "10 μM"),
    summarize_preds(pred_list$pred_single_0.1, "0.1 μM")
  )
  hmap_data$Dose <- factor(hmap_data$Dose, levels = c("0.1 μM", "10 μM"))

  ggplot(hmap_data, aes(x = Drug, y = Dose, fill = MeanProbability)) +
    geom_tile() +
    scale_fill_viridis(option = "viridis", na.value = "white",
                       name = "Mean\nProbability") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Drug", y = "Dose")
}

# -----------------------------------------------------------------------------
# 4. Feature importance plot (RF only)
# -----------------------------------------------------------------------------

#' Bar Plot of Random Forest Feature Importance (Top 20 Genes)
#'
#' @param rf_model A randomForest object trained with importance = TRUE.
#' @param genes_anabolic Anabolic gene panel (default: GENES_ANABOLIC).
#' @param genes_catabolic Catabolic gene panel (default: GENES_CATABOLIC).
#' @param genes_housekeeping Housekeeping gene panel (default: GENES_HOUSEKEEPING).
#' @param top_n Number of top genes to show (default 20).
#' @return A ggplot2 bar plot.
#' @export
plot_feature_importance <- function(rf_model,
                                    genes_anabolic    = GENES_ANABOLIC,
                                    genes_catabolic   = GENES_CATABOLIC,
                                    genes_housekeeping = GENES_HOUSEKEEPING,
                                    top_n = 20) {

  imp_mat <- importance(rf_model, type = 1)  # %IncMSE
  imp_df <- data.frame(Gene       = rownames(imp_mat),
                       Importance = as.numeric(imp_mat[, 1])) %>%
    mutate(Group = case_when(
      Gene %in% genes_anabolic     ~ "Anabolic",
      Gene %in% genes_catabolic    ~ "Catabolic",
      Gene %in% genes_housekeeping ~ "Housekeeping",
      TRUE ~ "Other"
    )) %>%
    arrange(desc(Importance)) %>%
    slice_head(n = top_n)

  ggplot(imp_df, aes(x = reorder(Gene, Importance), y = Importance, fill = Group)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.5, width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c(Anabolic     = COL_DISEASE,
                                 Catabolic    = COL_CONTROL,
                                 Housekeeping = "grey70")) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 11, face = "italic", color = "black"),
          axis.text.x = element_text(size = 11, color = "black"),
          axis.title  = element_text(size = 13, face = "bold"),
          legend.position = c(0.8, 0.2),
          legend.background = element_rect(color = "black", fill = "white"),
          legend.title = element_blank()) +
    labs(x = "Gene Symbol", y = "Importance (% Increase in MSE)")
}

# -----------------------------------------------------------------------------
# 5. Main execution
# -----------------------------------------------------------------------------

# Predictions using the final RF model
pred_list_rf <- predict_all_doses(rf_model, sin.sct, cmb.sct)

# Assign to globally named objects expected by figure scripts
prob_sng10_rf  <- pred_list_rf$pred_single_10
prob_sng0.1_rf <- pred_list_rf$pred_single_0.1
prob_cmb0.1_rf <- pred_list_rf$pred_combo_0.1

# Feature importance figure (Extended Data / Supplementary)
p_importance <- plot_feature_importance(rf_final)
ggsave(file.path(OUTPUT_SUP, "feature_importance_RF.pdf"),
       plot = p_importance, width = 6, height = 5, dpi = 300)

save_session_info("prediction_score")
