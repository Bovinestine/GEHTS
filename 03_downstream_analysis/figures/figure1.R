# =============================================================================
# figure1.R — Platform comparison bubble plot (Figure 1)
#
# Input:  PLATFORM_CSV (defined in config.R)
# Output: output/main/Figure1_platform_comparison.pdf
# =============================================================================

source(here::here("03_downstream_analysis", "R", "config.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

#' Bubble Plot Comparing Transcriptomic Screening Platforms
#'
#' X-axis: cells per condition (log scale).
#' Y-axis: cost per condition (USD).
#' Bubble size: log10(target genes).
#' Color: whether drug combination screening has been demonstrated.
#'
#' @param data A data.frame with columns Platform, Cell_Input_per_Condition,
#'   Cost_per_Condition, Target_Genes, Combo_Demonstrated.
#' @return A ggplot2 object.
#' @export
plot_platform_comparison <- function(data) {

  required_cols <- c("Platform", "Cell_Input_per_Condition", "Cost_per_Condition",
                     "Target_Genes", "Combo_Demonstrated")
  missing <- setdiff(required_cols, colnames(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  data$Cell_Input_per_Condition <- as.numeric(data$Cell_Input_per_Condition)
  data$Cost_per_Condition       <- as.numeric(data$Cost_per_Condition)
  data$Target_Genes             <- as.numeric(data$Target_Genes)
  data$Log_Target_Genes         <- log10(data$Target_Genes)

  ggplot(data, aes(x     = Cell_Input_per_Condition,
                   y     = Cost_per_Condition,
                   size  = Log_Target_Genes,
                   color = Combo_Demonstrated,
                   label = Platform)) +
    geom_point(alpha = 0.8) +
    geom_text(vjust = -1.1, size = 3) +
    scale_x_log10(labels = comma) +
    scale_y_continuous(labels = dollar) +
    scale_color_manual(values = c("No" = COL_CONTROL, "Yes" = COL_DISEASE)) +
    scale_size_continuous(range = c(2, 8), name = "log₁₀(Target Genes)") +
    labs(x     = "Cells per Condition (log scale)",
         y     = "Cost per Condition (USD)",
         color = "Combo Demonstrated") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right")
}

# --- Execution ---
platforms <- read.csv(PLATFORM_CSV)

fig1_plot <- plot_platform_comparison(platforms)
ggsave(file.path(OUTPUT_MAIN, "Figure1_platform_comparison.pdf"),
       plot = fig1_plot, width = 7, height = 5, dpi = 300)

save_session_info("figure1")
