# GE-HTS project pipeline for in situ sequencing data analysis
# author: Nathan Wooseok Lee
# Date of start: 250728
# Date of update: 
# env: conda activate seurat4

platforms <- read.csv("./Transcriptomic_Screening_Platform_Comparison.csv")
rownames(platforms) <- platforms$Platform


plot_platform_comparison <- function(data) {
  library(ggplot2)
  library(scales)
  
  # Check required columns exist
  required_cols <- c("Platform", "Cell_Input_per_Condition", "Cost_per_Condition", 
                     "Target_Genes", "Combo_Demonstrated")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }
  # Ensure numeric values
  data$Cell_Input <- as.numeric(data$Cell_Input)
  data$Cost_per_Condition <- as.numeric(data$Cost_per_Condition)
  data$Target_Genes <- as.numeric(data$Target_Genes)

  # Add log-scaled size
  data$Log_Target_Genes <- log10(data$Target_Genes)

  # Plot
  p<-ggplot(data, aes(
    x = Cell_Input_per_Condition,
    y = Cost_per_Condition,
    size = Log_Target_Genes,
    color = Combo_Demonstrated,
    label = Platform
  )) +
    geom_point(alpha = 0.8) +
    geom_text(vjust = -1.1, size = 3) +
    scale_x_log10(labels = comma) +
    scale_y_continuous(labels = dollar) +
    scale_color_manual(values = c("No" = "blue", "Yes" = "darkred")) +
    scale_size_continuous(range = c(2, 8), name = "log(Target Genes)") +
    labs(
      x = "Cells per Condition (log scale)",
      y = "Cost per Condition (USD)",
      size = "Target Genes",
      color = "Combo Demonstrated",
      title = "Comparison of Transcriptomic Screening Platforms"
    ) +
    theme_minimal(base_size = 14)
    print(p)
}

plot_platform_comparison(platforms)