# GE-HTS project pipeline for in situ sequencing data analysis
# author: Nathan Wooseok Lee
# Date of start: 250630
# Date of update: 250701
# env: conda activate seurat4 / issr

# ------- Function definition -------

predict_and_average_generic <- function(seurat_object, model) {
  # Extract and transpose expression
  drug_data <- GetAssayData(seurat_object, assay = "SCT", slot = "counts")
  drug_data <- t(as.matrix(drug_data))

  # Align features by model type
  if (inherits(model, "xgb.Booster")) {
    # XGBoost requires feature_names match
    if (!is.null(model$feature_names)) {
      drug_data <- drug_data[, model$feature_names, drop = FALSE]
    }
    predicted_scores <- predict(model, newdata = drug_data)

  } else if (inherits(model, "cv.glmnet")) {
    # glmnet: use non-zero features
    model_features <- rownames(coef(model, s = "lambda.min"))[-1]  # drop intercept
    model_features <- intersect(model_features, colnames(drug_data))
    drug_data <- drug_data[, model_features, drop = FALSE]
    predicted_scores <- predict(model, newx = drug_data, s = "lambda.min", type = "response")[,1]

  } else if (inherits(model, "randomForest")) {
    predicted_scores <- predict(model, newdata = drug_data)

  } else if (inherits(model, "gausspr")) {
    predicted_scores <- predict(model, newdata = drug_data)

  } else {
    stop("Unsupported model type for prediction.")
  }

  # Group and average predictions
  drug_names <- seurat_object@meta.data$drug_name
  mean_scores <- tapply(predicted_scores, drug_names, mean, na.rm = TRUE)
  return(mean_scores)
}


predict_and_average_generic_cmb <- function(seurat_object,model) {
  # Extract and transpose expression
  drug_data <- GetAssayData(seurat_object, assay = "SCT", slot = "counts")
  drug_data <- t(as.matrix(drug_data))

  # Align features by model type
  if (inherits(model, "xgb.Booster")) {
    # XGBoost requires feature_names match
    if (!is.null(model$feature_names)) {
      drug_data <- drug_data[, model$feature_names, drop = FALSE]
    }
    predicted_scores <- predict(model, newdata = drug_data)

  } else if (inherits(model, "cv.glmnet")) {
    # glmnet: use non-zero features
    model_features <- rownames(coef(model, s = "lambda.min"))[-1]  # drop intercept
    model_features <- intersect(model_features, colnames(drug_data))
    drug_data <- drug_data[, model_features, drop = FALSE]
    predicted_scores <- predict(model, newx = drug_data, s = "lambda.min", type = "response")[,1]

  } else if (inherits(model, "randomForest")) {
    predicted_scores <- predict(model, newdata = drug_data)

  } else if (inherits(model, "gausspr")) {
    predicted_scores <- predict(model, newdata = drug_data)

  } else {
    stop("Unsupported model type for prediction.")
  }

  # Construct custom combination names
  drug_names1 <- seurat_object@meta.data$drug_name1
  drug_names2 <- seurat_object@meta.data$drug_name2
  combo_names <- paste0(drug_names1, "&", drug_names2)
  names(predicted_scores) <- combo_names

  # Aggregate by custom drug combination name
  mean_scores <- tapply(predicted_scores, combo_names, mean, na.rm = TRUE)

  return(mean_scores)
}

predict_all_doses <- function(model, single_seurat, combo_seurat, save_dir = "./Rdata/prediction") {
  # Safely extract model name like "rf_model" from the function argument
  model_expr <- substitute(model)
  model_name <- deparse(model_expr)

  # Use part before "_model" as the model symbol
  if (grepl("_model", model_name)) {
    model_symbol <- sub("_model.*$", "", basename(model_name))
  } else {
    model_symbol <- "model"
  }

  # Ensure save directory exists
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  # Subset single-drug data by dose
  single_cmb10 <- subset(single_seurat, subset = dose == 10)
  # single_cmb1   <- subset(single_seurat, subset = dose == 1)
  single_cmb0.1 <- subset(single_seurat, subset = dose == 0.1)

  # Subset combo-drug data by dose
  # cmb1   <- subset(combo_seurat, subset = dose1 == 1)
  cmb0.1 <- subset(combo_seurat, subset = dose1 == 0.1)

  # Predict for single-drug
  prob_sng10 <- predict_and_average_generic(single_cmb10, model)
  # prob_sng1  <- predict_and_average_generic(single_cmb1, model)
  prob_sng0.1 <- predict_and_average_generic(single_cmb0.1, model)

  # Predict for combinations
  # prob_cmb1   <- predict_and_average_generic_cmb(cmb1, model)
  prob_cmb0.1 <- predict_and_average_generic_cmb(cmb0.1, model)

  # Save each object using model symbol in filename
  save(prob_sng10, file = file.path(save_dir, sprintf("prob_sng10_%s.Rdata", model_symbol)))
  # save(prob_sng1,  file = file.path(save_dir, sprintf("prob_sng1_%s.Rdata", model_symbol)))
  save(prob_sng0.1, file = file.path(save_dir, sprintf("prob_sng0.1_%s.Rdata", model_symbol)))
  # save(prob_cmb1,   file = file.path(save_dir, sprintf("prob_cmb1_%s.Rdata", model_symbol)))
  save(prob_cmb0.1, file = file.path(save_dir, sprintf("prob_cmb0.1_%s.Rdata", model_symbol)))

  return(list(
    prob_sng10 = prob_sng10,
    #prob_sng1  = prob_sng1,
    prob_sng0.1 = prob_sng0.1,
    #prob_cmb1   = prob_cmb1,
    prob_cmb0.1 = prob_cmb0.1
  ))
}


plot_single_drug_heatmap <- function(prob_list) {
  # Extract dose 10 and dose 0.1 single-drug scores
  prob_sng10 <- prob_list[[1]]
  prob_sng0.1 <- prob_list[[2]]
  
  # Construct heatmap input data frame
  heatmap_data <- data.frame(
    Drug = c(names(prob_sng10), names(prob_sng0.1)),
    Dose = c(rep("10", length(prob_sng10)), rep("0.1", length(prob_sng0.1))),
    MeanProbability = c(prob_sng10, prob_sng0.1)
  )
  
  # Plot the heatmap
  p <- ggplot(heatmap_data, aes(x = Drug, y = Dose, fill = MeanProbability)) +
    geom_tile() +
    scale_fill_viridis(option = 'viridis', na.value = "white") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(fill = "Mean Probability", x = "Drug", y = "Dose")
  
  return(p)
}

# ------ main excute --------
# 260219 
# 1. Train the FINAL model on 100% of the data for the figure
# Setting importance = TRUE is crucial to get %IncMSE, which is the most robust metric
set.seed(123)
rf_final <- randomForest(x = X, y = y, ntree = 500, importance = TRUE)
# 2. Extract Feature Importance
# type = 1 extracts %IncMSE (Percent Increase in Mean Squared Error if the gene is permuted)
imp_matrix <- importance(rf_final, type = 1)
imp_df <- data.frame(
  Gene = rownames(imp_matrix),
  Importance = as.numeric(imp_matrix[, 1])
)

# 3. Annotate Genes by Biological Group
# Using your exact gene lists
imp_df <- imp_df %>%
  mutate(
    Group = case_when(
      Gene %in% anabolic ~ "Anabolic",
      Gene %in% inflammatory ~ "Inflammatory",
      Gene %in% housekeeping ~ "Housekeeping",
      TRUE ~ "Other"
    )
  ) %>%
  # Sort by importance and select the top 20 for a clean figure
  arrange(desc(Importance)) %>%
  slice_head(n = 20) 

# 4. Generate the Publication Plot
ggplot(imp_df, aes(x = reorder(Gene, Importance), y = Importance, fill = Group)) +
  
  # Horizontal bars (width=0.7 leaves nice spacing, size=0.5 adds crisp borders)
  geom_bar(stat = "identity", color = "black", size = 0.5, width = 0.7) +
  
  # Flip coordinates so gene names are easily readable on the Y-axis
  coord_flip() +
  
  # Color mapping matching the NPG palette from the previous scatter plot
  scale_fill_manual(
    values = c(
      "Anabolic" = "#E64B35FF",     # NPG Red
      "Inflammatory" = "#3C5488FF", # NPG Dark Blue
      "Housekeeping" = "grey70"     # Neutral Grey
    )
  ) +
  
  # Clean, minimalist theme
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 11, face = "italic", color = "black"), # Italicize gene symbols!
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = c(0.8, 0.2), # Inset legend to save space
    legend.background = element_rect(color = "black", fill = "white", size = 0.5),
    legend.title = element_blank()
  ) +
  labs(
    title = "Figure 3G: ML Feature Importance",
    x = "Gene Symbol",
    y = "Importance (% Increase in MSE)"
  )


# 251230
# Load the best-performing Random Forest model, when the training was done with mac.sct69
rf_model <- readRDS("./Rdata/prediction/Random Forest_model.rds")

# Subset drug-perturbed samples
# 251230


# 251230 function calling
prob_list_rf <- predict_all_doses(rf_model, sin.sct, cmb.sct)
prob_sng10_rf <-prob_list_rf[[1]]
prob_sng0.1_rf <-prob_list_rf[[2]]
prob_cmb0.1_rf <-prob_list_rf[[3]]

# save probability
save(prob_sng10_rf, file = "./Rdata/prediction/prob_sng10_rf_260219.Rdata")
#save(prob_sng1_rf, file = "./Rdata/prediction/prob_sng1_rf.Rdata")
save(prob_sng0.1_rf, file = "./Rdata/prediction/prob_sng0.1_rf_260219.Rdata")

#save(prob_cmb1_rf, file = "./Rdata/prediction/prob_cmb1_rf.Rdata")
save(prob_cmb0.1_rf, file = "./Rdata/prediction/prob_cmb0.1_rf_260219.Rdata")



# 250701 function calling
prob_list_rf <- predict_all_doses(rf_model, single.SCT, cmb.SCT)
prob_list_xgb <- predict_all_doses(xgb_model, single.SCT, cmb.SCT)
prob_list_svr <- predict_all_doses(svr_model, single.SCT, cmb.SCT)
prob_list_en <- predict_all_doses(en_model, single.SCT, cmb.SCT)
prob_list_gpr <- predict_all_doses(gpr_model, single.SCT, cmb.SCT)

# ------ visualization ----------
# for various models
hm_sp_rf <- plot_single_drug_heatmap(prob_list_rf)
print(hm_sp_rf)

hm_sp_xgb <- plot_single_drug_heatmap(prob_list_xgb)
print(hm_sp_xgb)

hm_sp_svr <- plot_single_drug_heatmap(prob_list_svr)
print(hm_sp_svr)

hm_sp_en <- plot_single_drug_heatmap(prob_list_en)
print(hm_sp_en)

hm_sp_gpr <- plot_single_drug_heatmap(prob_list_gpr)
print(hm_sp_gpr)
# Combine into a single data frame
# all doses
heatmap_data <- data.frame(
  Drug = c(names(prob_sng10_rf), names(prob_sng1_rf), names(prob_sng0.1_rf)),
  Dose = c(rep("10", length(prob_sng10_rf)), rep("1", length(prob_sng1_rf)), rep("0.1", length(prob_sng0.1_rf))),
  MeanProbability = c(prob_sng10_rf, prob_sng1_rf, prob_sng0.1_rf)
)
# without 1 uM dose
# 250630 This figure is used for Figure.3i
heatmap_data_d1 <- data.frame(
  Drug = c(names(prob_sng10_rf), names(prob_sng0.1_rf)),
  Dose = c(rep("10", length(prob_sng10_rf)),rep("0.1", length(prob_sng0.1_rf))),
  MeanProbability = c(prob_sng10_rf, prob_sng0.1_rf)
)

# Draw heatmap 
# change the heatmap dataset
ggplot(heatmap_data_d1, aes(x = Drug, y = Dose, fill = MeanProbability)) +
  geom_tile() +
  # scale_fill_gradient(low = "blue", high = "red") +
  scale_fill_viridis(option='viridis', na.value = "white") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "Mean Probability", x = "Drug", y = "Dose")
  


# previously using all mac.SCT GPR was the best model in terms of RMSE.
# Load trained GPR model
gpr_model <- readRDS("./Rdata/prediction/GPR_model.rds")

# Subset the perturbation data (e.g., 10μM dose)
single_cmb10 <- subset(single.SCT, subset = dose == 10)

# Get mean efficacy score per drug
mean_scores_10uM <- predict_and_average_gpr(single_cmb10, gpr_model)

print(mean_scores_10uM)


##### test with old Regression model
# Predict probabilities from logistic regression model
pred_prob <- as.vector(predict(cfit, newx = testSet_x, s = "lambda.min", type = "response"))

# Calculate regression-like metrics
r2_val <- cor(testSet_y, pred_prob)^2
rmse_val <- rmse(testSet_y, pred_prob)
mae_val <- mae(testSet_y, pred_prob)
pearson_val <- cor(testSet_y, pred_prob)

cat("\n Elastic Net Model as Regression:\n")
cat(sprintf("  R²       = %.4f\n", r2_val))
cat(sprintf("  RMSE     = %.4f\n", rmse_val))
cat(sprintf("  MAE      = %.4f\n", mae_val))
cat(sprintf("  Pearson  = %.4f\n", pearson_val))


### ---------- obsolete functions ------------ ###


predict_and_average_rf <- function(seurat_object, rf_model) {
  # Extract expression data (genes x cells), transpose to samples x genes
  drug_data <- GetAssayData(seurat_object, assay = "SCT", slot = "counts")
  drug_data <- t(as.matrix(drug_data))
  # Align columns to match model's expected feature names (XGBoost)
  if (!is.null(rf_model$feature_names)) {
    common_features <- intersect(rf_model$feature_names, colnames(drug_data))
    drug_data <- drug_data[, rf_model$feature_names, drop = FALSE]  # exact order
  }
  # Predict efficacy scores
  predicted_scores <- predict(rf_model, newdata = drug_data)

  # Group by drug_name and compute mean
  drug_names <- seurat_object@meta.data$drug_name
  names(predicted_scores) <- drug_names
  mean_scores <- tapply(predicted_scores, drug_names, mean, na.rm = TRUE)
  

  return(mean_scores)
}

predict_and_average_rfcmb <- function(seurat_object, rf_model) {
  # Extract expression data (genes x cells), transpose to samples x genes
  drug_data <- GetAssayData(seurat_object, assay = "SCT", slot = "counts")
  drug_data <- t(as.matrix(drug_data))
  if (!is.null(rf_model$feature_names)) {
    common_features <- intersect(rf_model$feature_names, colnames(drug_data))
    drug_data <- drug_data[, rf_model$feature_names, drop = FALSE]  # exact order
  }
  # Predict efficacy scores
  predicted_scores <- predict(rf_model, newdata = drug_data)

  # Construct custom combination names
  drug_names1 <- seurat_object@meta.data$drug_name1
  drug_names2 <- seurat_object@meta.data$drug_name2
  combo_names <- paste0(drug_names1, "&", drug_names2)
  names(predicted_scores) <- combo_names

  # Aggregate by custom drug combination name
  mean_scores <- tapply(predicted_scores, combo_names, mean, na.rm = TRUE)

  return(mean_scores)
}


predict_and_average_gpr <- function(seurat_object, gpr_model) {
  # Extract expression matrix (genes x cells)
  drug_data <- GetAssayData(seurat_object, assay = "SCT", slot = 'counts')
  drug_data <- t(as.matrix(drug_data))  # Samples as rows

  # Predict efficacy scores (continuous)
  predicted_scores <- predict(gpr_model, newdata = as.data.frame(drug_data))

  # Group-wise mean score per drug_name
  drug_names <- seurat_object@meta.data$drug_name
  mean_scores <- tapply(predicted_scores, drug_names, mean, na.rm = TRUE)

  return(mean_scores)
}
