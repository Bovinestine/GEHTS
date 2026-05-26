# =============================================================================
# scoring_methods.R — ML model training and cross-validation
#
# Trains five regression models (Elastic Net, RF, XGBoost, SVR, GPR) on the
# primary MAC SCT data using 5-fold cross-validation. Saves models and CV
# metrics to PRED_DIR. The final RF model is trained on 100% of data for
# downstream prediction (see prediction_score.R).
#
# Prerequisites: config.R, data_loading.R
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(glmnet)
  library(randomForest)
  library(xgboost)
  library(e1071)
  library(kernlab)
  library(caret)
  library(Metrics)
  library(dplyr)
  library(ggplot2)
})

set.seed(SEED)

# -----------------------------------------------------------------------------
# 1. Data preparation
# -----------------------------------------------------------------------------

#' Extract Expression Matrix and Binary Efficacy Labels from Seurat Object
#'
#' Assigns 1 to "control" cells and 0 to "inflammatory" cells.
#'
#' @param seurat_obj A Seurat object. Must contain exactly two unique drug_condition values.
#' @param gene_list Character vector of genes to use as features.
#' @return A list with elements X (expression matrix) and y (binary label vector).
#' @export
prepare_data <- function(seurat_obj, gene_list) {
  expr <- GetAssayData(seurat_obj, assay = "SCT", slot = "counts")[gene_list, ]
  X    <- t(as.matrix(expr))

  labels <- seurat_obj$drug_condition
  if (length(unique(labels)) != 2) {
    stop("Expected exactly 2 unique drug_condition values; got: ",
         paste(unique(labels), collapse = ", "))
  }

  # Explicit conversion: control=1, inflammatory=0
  y <- ifelse(labels == COND_CONTROL, 1L, 0L)
  if (any(is.na(y))) stop("NA values found in efficacy labels.")

  list(X = X, y = y)
}

# -----------------------------------------------------------------------------
# 2. Cross-validation
# -----------------------------------------------------------------------------

#' 5-Fold Cross-Validation for a Single Model Type
#'
#' @param X Numeric matrix (cells × genes).
#' @param y Numeric vector of binary labels (0/1).
#' @param model_type One of "glmnet", "rf", "xgb", "svm", "gpr".
#' @param folds Number of CV folds (default: N_FOLDS_CV from config).
#' @param return_model If TRUE, returns the model trained on the last fold instead of metrics.
#' @return A list with elements mean, sd (per-metric), and per_fold (data.frame).
#' @export
cross_validate_model <- function(X, y,
                                 model_type   = "rf",
                                 folds        = N_FOLDS_CV,
                                 return_model = FALSE) {

  set.seed(SEED)
  folds_idx <- createFolds(seq_along(y), k = folds, list = TRUE)

  metrics <- data.frame(RMSE    = numeric(folds),
                        MAE     = numeric(folds),
                        R2      = numeric(folds),
                        Pearson = numeric(folds))

  for (i in seq_len(folds)) {
    test_idx <- folds_idx[[i]]
    X_train  <- X[-test_idx, ]
    y_train  <- y[-test_idx]
    X_test   <- X[test_idx, ]
    y_test   <- y[test_idx]

    set.seed(SEED + i)  # per-fold seed for full reproducibility

    if (model_type == "glmnet") {
      # Elastic Net: alpha=0.5 (equal L1/L2 mix)
      model <- cv.glmnet(X_train, y_train, family = "gaussian", alpha = 0.5)
      preds <- as.numeric(predict(model, newx = X_test, s = "lambda.min"))

    } else if (model_type == "rf") {
      # Random Forest: 500 trees, default mtry = floor(sqrt(ncol(X)))
      model <- randomForest(x = X_train, y = y_train, ntree = 500)
      preds <- as.numeric(predict(model, newdata = X_test))

    } else if (model_type == "xgb") {
      # XGBoost: eta=0.1, max_depth=4, 100 rounds
      dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
      dtest  <- xgb.DMatrix(data = as.matrix(X_test))
      params <- list(objective  = "reg:squarederror",
                     eval_metric = "rmse",
                     eta        = 0.1,
                     max_depth  = 4)
      model  <- xgb.train(params = params, data = dtrain,
                          nrounds = 100, verbose = 0)
      preds  <- as.numeric(predict(model, newdata = dtest))

    } else if (model_type == "svm") {
      # Support Vector Regression: RBF kernel
      model <- svm(x = X_train, y = y_train,
                   type = "eps-regression", kernel = "radial")
      preds <- as.numeric(predict(model, newdata = X_test))

    } else if (model_type == "gpr") {
      # Gaussian Process Regression: RBF kernel via kernlab
      df_train   <- as.data.frame(X_train)
      df_train$y <- y_train
      model <- gausspr(y ~ ., data = df_train, kernel = "rbfdot")
      preds <- as.numeric(predict(model, newdata = as.data.frame(X_test)))

    } else {
      stop("Unknown model_type: '", model_type, "'. ",
           "Choose from: glmnet, rf, xgb, svm, gpr.")
    }

    metrics$RMSE[i]    <- rmse(y_test, preds)
    metrics$MAE[i]     <- mae(y_test, preds)
    metrics$R2[i]      <- cor(y_test, preds)^2
    metrics$Pearson[i] <- cor(y_test, preds)
  }

  if (return_model) return(model)

  list(
    mean     = round(colMeans(metrics), 4),
    sd       = round(apply(metrics, 2, sd), 4),
    per_fold = metrics
  )
}

# -----------------------------------------------------------------------------
# 3. Run CV on all five models
# -----------------------------------------------------------------------------

data_for_cv <- prepare_data(mac.sct, GENES_ALL)
X <- data_for_cv$X
y <- data_for_cv$y

model_types <- c("glmnet", "rf",             "xgb",      "svm", "gpr")
model_names <- c("Elastic Net", "Random Forest", "XGBoost", "SVR", "GPR")

cv_results <- setNames(
  lapply(seq_along(model_types), function(i) {
    cat("\nRunning", N_FOLDS_CV, "-fold CV for:", model_names[i], "\n")
    cross_validate_model(X, y, model_type = model_types[i], folds = N_FOLDS_CV)
  }),
  model_names
)

# -----------------------------------------------------------------------------
# 4. Report and select best model
# -----------------------------------------------------------------------------
for (nm in names(cv_results)) {
  cat("\nResults for", nm, ":\n")
  print(data.frame(Mean = cv_results[[nm]]$mean, SD = cv_results[[nm]]$sd))
}

rmse_vals  <- sapply(cv_results, function(x) x$mean["RMSE"])
best_model_name <- names(which.min(rmse_vals))
cat("\nBest model by RMSE:", best_model_name, "\n")

# -----------------------------------------------------------------------------
# 5. Save models and CV metrics
# -----------------------------------------------------------------------------
if (!dir.exists(PRED_DIR)) dir.create(PRED_DIR, recursive = TRUE)

for (i in seq_along(model_types)) {
  cat("\nSaving model:", model_names[i], "\n")
  m <- cross_validate_model(X, y, model_type = model_types[i],
                            folds = N_FOLDS_CV, return_model = TRUE)
  saveRDS(m, file.path(PRED_DIR, paste0(model_names[i], "_model.rds")))
}

# Per-fold CSV
cv_per_fold <- do.call(rbind, lapply(names(cv_results), function(nm) {
  df <- cv_results[[nm]]$per_fold
  df$Model <- nm
  df$Fold  <- seq_len(nrow(df))
  df
}))
write.csv(cv_per_fold,
          file.path(PRED_DIR, "cv_results_per_fold.csv"),
          row.names = FALSE)

# Summary CSV (mean ± SD)
cv_summary <- do.call(rbind, lapply(names(cv_results), function(nm) {
  cbind(Model = nm,
        t(cv_results[[nm]]$mean),
        t(cv_results[[nm]]$sd))
}))
colnames(cv_summary)[2:9] <- c("RMSE_mean", "MAE_mean", "R2_mean", "Pearson_mean",
                                "RMSE_sd",   "MAE_sd",   "R2_sd",   "Pearson_sd")
write.csv(cv_summary,
          file.path(PRED_DIR, "cv_results_summary.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. Train final RF model on 100% of data for downstream prediction
# -----------------------------------------------------------------------------
set.seed(SEED)
rf_final <- randomForest(x = X, y = y, ntree = 500, importance = TRUE)
saveRDS(rf_final, file.path(PRED_DIR, "Random Forest_model.rds"))
message("Final RF model saved.")

save_session_info("scoring_methods")
