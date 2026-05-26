# =============================================================================
# figure3_supplement.R — Supplementary figure for Figure 3
#
# Subfigures:
#   Extended Data Fig. 8: Highlighted UMAP grid for remaining 11 drug conditions
#   Extended Data Fig. 9: Boxplots of 10 key genes (5 anabolic + 5 catabolic)
#                         across control, inflammatory, pamapimod_10, sb203580_10
#
# Prerequisites: config.R, utils.R, data_loading.R, figure3.R
# Inputs: sin.sct (loaded by data_loading.R),
#         sin.sct10 (UMAP-embedded subset created by figure3.R)
# Outputs: output/supplementary/Sfig_fig3ext_boxplots.pdf,
#          output/supplementary/Sfig_fig3ext8_highlighted_UMAPs.pdf
# =============================================================================

source(here::here("03_downstream_analysis", "R", "config.R"))
source(here::here("03_downstream_analysis", "R", "utils.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

set.seed(SEED)



# =============================================================================
# Extended Data Fig. 8 — Highlighted UMAPs for remaining drug conditions
# Indices 2:12 of sin.sct10 (all drugs except pamapimod [1] and XAV939 [13],
# which appear in Figure 4d)
# =============================================================================

plot_highlighted_umaps(
  sin.sct10,
  plot_indices = 2:12,
  output_pdf   = file.path(OUTPUT_SUP, "Sfig_fig3ext8_highlighted_UMAPs.pdf"),
  ncol         = 3,
  pdf_width    = 9,
  pdf_height   = 13
)

# =============================================================================
# Extended Data Fig. 9 — Boxplots of 10 key genes (5 anabolic + 5 catabolic)
#                         across control, inflammatory, pamapimod_10, sb203580_10
# =============================================================================

GENES_FIG3_EXT <- c(
  "Col2a1", "Matn3", "Gadd45g", "Sox9", "Acan",
  "Mmp3",   "Mmp13", "Cxcl1",   "Fosl2", "Tnfrsf1b"
)

sfig_fig3_boxplots <- draw_boxplots_for_genes_with_common_legend(
  sin.sct,
  genes         = GENES_FIG3_EXT,
  extra_label1  = "pamapimod_10",
  extra_label2  = "sb203580_10",
  nrow          = 2,
  base_font_size = 6,
  legend_gap    = 0.03
)

ggsave(file.path(OUTPUT_SUP, "Sfig_fig3ext_boxplots.pdf"),
       plot   = sfig_fig3_boxplots,
       width  = 7,
       height = 5,
       dpi    = 300)

save_session_info("figure3_supplement")

