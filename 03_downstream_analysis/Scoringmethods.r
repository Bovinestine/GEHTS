# GE-HTS project pipeline for in situ sequencing data analysis
# author: Nathan Wooseok Lee

# Load required libraries
library(Seurat)
library(glmnet)
library(randomForest)
library(xgboost)
library(e1071)
library(kernlab)
library(caret)
library(Metrics)
library(dplyr)

# -------- Prepare Data --------
prepare_data <- function(seurat_obj, gene_list, score_col = "efficacy_score") {
  expr_data <- GetAssayData(seurat_obj, assay = 'SCT', slot = "counts")[gene_list, ]
  expr_data <- t(as.matrix(expr_data))
  # Extract labels
  labels <- seurat_obj$drug_condition
  # Ensure labels are binary
  unique_labels <- unique(labels)
  if (length(unique_labels) != 2) {
    stop("There should be exactly two unique drug conditions")
  }
  # Convert labels to 0 and 1
  labels <- 2 - as.integer(factor(labels))
  if (any(is.na(labels))) stop("Efficacy scores contain NA values.")
  return(list(X = expr_data, y = labels))
}

# -------- Cross-validation function --------
cross_validate_model <- function(X, y, model_type = "rf", folds = 5, return_model = FALSE) {
  set.seed(123)
  #folds_idx <- createFolds(y, k = folds, list = TRUE) # if y is a factor
  folds_idx <- createFolds(1:length(y), k = folds, list = TRUE) #simply split the row indices of y, regardless of whether y is numeric (as in regression) or factor (as in classification)

  
  metrics <- data.frame(RMSE = numeric(folds), MAE = numeric(folds),
                        R2 = numeric(folds), Pearson = numeric(folds))
  
  for (i in 1:folds) {
    test_idx <- folds_idx[[i]]
    X_train <- X[-test_idx, ]
    y_train <- y[-test_idx]
    X_test <- X[test_idx, ]
    y_test <- y[test_idx]
    
    if (model_type == "glmnet") {
      model <- cv.glmnet(X_train, y_train, family = "gaussian", alpha = 0.5)
      preds <- predict(model, newx = X_test, s = "lambda.min")
    } else if (model_type == "rf") {
      model <- randomForest(x = X_train, y = y_train, ntree = 500)
      preds <- predict(model, newdata = X_test)
    } else if (model_type == "xgb") {
      dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
      dtest <- xgb.DMatrix(data = as.matrix(X_test))
      params <- list(objective = "reg:squarederror", eval_metric = "rmse", eta = 0.1, max_depth = 4)
      model <- xgb.train(params = params, data = dtrain, nrounds = 100, verbose = 0)
      preds <- predict(model, newdata = dtest)
    } else if (model_type == "svm") {
      model <- svm(x = X_train, y = y_train, type = "eps-regression", kernel = "radial")
      preds <- predict(model, newdata = X_test)
    } else if (model_type == "gpr") {
      df_train <- as.data.frame(X_train)
      df_train$y <- y_train
      model <- gausspr(y ~ ., data = df_train, kernel = "rbfdot")
      preds <- predict(model, newdata = as.data.frame(X_test))
    } else {
      stop("Unknown model type")
    }
    
    metrics$RMSE[i] <- rmse(y_test, preds)
    metrics$MAE[i] <- mae(y_test, preds)
    metrics$R2[i] <- cor(y_test, preds)^2
    metrics$Pearson[i] <- cor(y_test, preds)
  }
  
  summary <- colMeans(metrics) %>%
    round(4)
  sd_vals <- apply(metrics, 2, sd) %>%
    round(4)
  if (return_model) return(model)
  return(list(mean = summary, sd = sd_vals, per_fold = metrics))
}

# -------- Run Cross-Validation on All Models --------
anabolic <- c('Acan','Sox9','Col2a1','Matn1','Matn3','Ucma','Ccnd3','Gadd45g','Pth1r','Gm26633','Col27a1') # 251230
inflammatory <- c('Mmp3','Mmp13','Il6','Il17b','Adamts5','Igfbp3','Ccl2','Cxcl5','Cxcl1','Fosl2','Tlr2','Tnfrsf1b')
housekeeping <- c('Hprt','Actb','Gapdh','B2m','Ubc','Ppia','Rpl23')
all_genes <- c(anabolic, inflammatory, housekeeping)

data <- prepare_data(mac.sct, all_genes) 

X <- data$X
y <- data$y

model_list <- c("glmnet", "rf", "xgb", "svm", "gpr")
model_names <- c("Elastic Net", "Random Forest", "XGBoost", "SVR", "GPR")
cv_results <- list()

for (i in seq_along(model_list)) {
  cat("\nRunning 5-fold CV for:", model_names[i], "\n")
  res <- cross_validate_model(X, y, model_type = model_list[i], folds = 5)
  cv_results[[model_names[i]]] <- res
}

# -------- Display Summary Results --------
for (name in names(cv_results)) {
  cat("\n Results for", name, ":\n")
  print(data.frame(Mean = cv_results[[name]]$mean, SD = cv_results[[name]]$sd))
}

# -------- Find Best Model by RMSE --------
rmse_values <- sapply(cv_results, function(x) x$mean["RMSE"])
best_model <- names(which.min(rmse_values))
cat("\n Best model by RMSE (lower is better):", best_model, "\n")

# -------- Save the outputs ---------
# Create output dir
output_dir <- "./Rdata/prediction/"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Save models after training
for (i in seq_along(model_list)) {
  cat("\nRunning 5-fold CV for:", model_names[i], "\n")
  model <- cross_validate_model(X, y, model_type = model_list[i], folds = 5, return_model = TRUE)
  # cv_results[[model_names[i]]] <- res
  saveRDS(model, file = file.path(output_dir, paste0(model_names[i], "_model.rds")))
}


# Save CV metircs as CSV
# Per-fold results
cv_summary_df <- do.call(rbind, lapply(names(cv_results), function(name) {
  df <- cv_results[[name]]$per_fold
  df$Model <- name
  df$Fold <- seq_len(nrow(df))
  return(df)
}))
write.csv(cv_summary_df, file.path(output_dir, "cv_results_per_fold.csv"), row.names = FALSE)

# Summary (mean ± SD)
cv_summary_stats <- do.call(rbind, lapply(names(cv_results), function(name) {
  cbind(Model = name,
        t(cv_results[[name]]$mean),
        t(cv_results[[name]]$sd))
}))
colnames(cv_summary_stats)[2:9] <- c("RMSE_mean", "MAE_mean", "R2_mean", "Pearson_mean",
                                     "RMSE_sd", "MAE_sd", "R2_sd", "Pearson_sd")
write.csv(cv_summary_stats, file.path(output_dir, "cv_results_summary.csv"), row.names = FALSE)

# Visualization code
library(tidyverse)
library(reshape2)

# -------- Compile CV Results into DataFrame --------
cv_summary_df <- as.data.frame(cv_summary_df)

# -------- Tidy Format for ggplot2 --------
cv_long <- melt(cv_summary_df, id.vars = c("Model", "Fold"),
                variable.name = "Metric", value.name = "Value")

# -------- Boxplot Comparison --------
p <- ggplot(cv_long, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  facet_wrap(~ Metric, scales = "free_y") +
  labs(title = "5-Fold CV Performance by Model",
       y = "Metric Value", x = "Model") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "cv_comparison_plot.pdf"), plot = p, width = 10, height = 6)
# or you can save as pdf file from the R window popup.
