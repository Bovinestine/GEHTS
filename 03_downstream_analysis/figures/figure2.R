# =============================================================================
# figure2.R — Primary MAC characterization (Figure 2)
#
# Subfigures:
#   2b: Gene expression heatmap (SCTransform residuals, mac.sct)
#   2c: Transcriptome correlation scatter plots (MAC vs. RNA-seq vs. Sebastian)
#   2d: t-SNE plot of primary MAC cells
#   2e: Volcano plot (anabolic vs. catabolic, inflammatory vs. control)
#   2f: Lollipop chart (Primary MAC vs. ATDC5 cell line)
#
# Prerequisites: config.R, utils.R, data_loading.R
# Inputs: mac.sct, rna.pub, macSebastian, atdc5.pub (loaded by data_loading.R)
# Outputs: output/main/Figure2{b-f}_*.pdf
# =============================================================================

source(here::here("03_downstream_analysis", "R", "config.R"))
source(here::here("03_downstream_analysis", "R", "utils.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(ggsci)
  library(cowplot)
})

`%notin%` <- Negate(`%in%`)

# =============================================================================
# Figure 2b — Gene expression heatmap
# =============================================================================
fig_2b <- create_heatmap_by_expression_ps(
  mac.sct,
  primary_metadata   = "drug_condition",
  secondary_metadata = "file"
)
pdf(file.path(OUTPUT_MAIN, "Figure2b_heatmap_SCTransform.pdf"), width = 6, height = 5)
print(fig_2b)
dev.off()

# =============================================================================
# Figure 2c — t-SNE plot
# =============================================================================
DefaultAssay(mac.sct) <- "SCT"
set.seed(SEED)
mac.sct <- RunPCA(mac.sct, assay = "SCT", npcs = 10, verbose = FALSE)
mac.sct <- RunTSNE(mac.sct, dims = 1:10, verbose = FALSE)

fig_2c <- plot_dim_reduction(mac.sct, reduction = "tsne",
                              group_by = "drug_condition")
ggsave(file.path(OUTPUT_MAIN, "Figure2c_tsne.pdf"),
       plot = fig_2c, width = 7, height = 5.5, dpi = 300)

# =============================================================================
# Figure 2d — Volcano plot (inflammatory vs. control DE)
# =============================================================================

#' Annotated Volcano Plot for Differential Expression
#'
#' Highlights anabolic and catabolic genes that meet the logFC/adjusted-p thresholds.
#'
#' @param object A Seurat object.
#' @param ident.1 Primary identity class for FindMarkers (e.g., "inflammatory").
#' @param group.by Metadata column containing identities.
#' @param assay Assay for differential expression (default "SCT").
#' @param logfc.threshold Minimum absolute logFC for FindMarkers (default 0.1).
#' @param genes_anabolic Anabolic gene panel.
#' @param genes_catabolic Catabolic gene panel.
#' @return A ggplot2 object.
#' @export
generate_volcano_plot <- function(object,
                                  ident.1,
                                  group.by,
                                  assay           = "SCT",
                                  logfc.threshold = 0.1,
                                  genes_anabolic  = GENES_ANABOLIC,
                                  genes_catabolic = GENES_CATABOLIC) {

  object_prep <- PrepSCTFindMarkers(object = object)
  de <- FindMarkers(object_prep, ident.1 = ident.1, group.by = group.by,
                    assay = assay, logfc.threshold = logfc.threshold)

  de$logP_adj <- -log10(de$p_val_adj)
  de$gene     <- rownames(de)

  de <- de %>%
    mutate(gene_type = case_when(
      avg_log2FC >=  0.5 & p_val_adj <= 0.05 & gene %in% genes_catabolic ~ "Catabolic",
      avg_log2FC <= -0.5 & p_val_adj <= 0.05 & gene %in% genes_anabolic  ~ "Anabolic",
      TRUE ~ "ns"
    )) %>%
    arrange(factor(gene_type, levels = c("ns", "Anabolic", "Catabolic")))

  npg_green  <- pal_npg("nrc")(10)[3]
  jco_yellow <- pal_jco("default")(10)[2]

  cols   <- c(Catabolic = jco_yellow, Anabolic = npg_green, ns = "grey80")
  sizes  <- c(Catabolic = 2, Anabolic = 2, ns = 1)
  alphas <- c(Catabolic = 1, Anabolic = 1, ns = 0.3)

  ggplot(de, aes(x = avg_log2FC, y = logP_adj, label = gene,
                 fill = gene_type, size = gene_type, alpha = gene_type)) +
    geom_point(shape = 21, color = "black", stroke = 0.2) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    geom_text_repel(data = filter(de, gene_type != "ns"),
                    aes(label = gene), size = 3, color = "black",
                    box.padding = 0.5, max.overlaps = Inf,
                    show.legend = FALSE) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position  = "right",
          plot.title       = element_text(face = "bold", hjust = 0.5)) +
    labs(title  = paste("DE:", ident.1, "vs. control"),
         x      = "Log₂ Fold Change",
         y      = "-Log₁₀(Adjusted P-value)",
         fill   = "Gene Category") +
    scale_fill_manual(values = cols) +
    scale_size_manual(values = sizes, guide = "none") +
    scale_alpha_manual(values = alphas, guide = "none")
}

fig_2d <- generate_volcano_plot(mac.sct, ident.1 = COND_DISEASE,
                                 group.by = "drug_condition")
ggsave(file.path(OUTPUT_MAIN, "Figure2d_volcano.pdf"),
       plot = fig_2d, width = 6, height = 3, dpi = 300)


# =============================================================================
# Figure 2e — Transcriptome correlation scatter plots
# =============================================================================

#' Create Transcriptome Correlation Scatter Plots
#'
#' Compares average SCT expression between two Seurat objects across shared
#' conditions and saves one scatter per condition as a PDF.
#'
#' @param obj1 First Seurat object.
#' @param obj2 Second Seurat object.
#' @param obj1_name X-axis label.
#' @param obj2_name Y-axis label.
#' @param assay Assay to use (default "SCT").
#' @param slot Slot to use (default "data").
#' @param group_by Metadata column to average over (default "drug_condition").
#' @param pdf_file Output PDF path.
#' @param pdf_width PDF width in inches.
#' @param pdf_height PDF height in inches.
#' @return The output PDF path (invisibly).
#' @export
create_correlation_plots <- function(obj1, obj2,
                                     obj1_name  = deparse(substitute(obj1)),
                                     obj2_name  = deparse(substitute(obj2)),
                                     assay      = "SCT",
                                     slot       = "data",
                                     group_by   = "drug_condition",
                                     pdf_file   = file.path(OUTPUT_MAIN, "Figure2c_correlation.pdf"),
                                     pdf_width  = 2.5,
                                     pdf_height = 2.5) {

  avg1 <- AverageExpression(obj1, return.seurat = FALSE,
                            group.by = group_by, verbose = FALSE)[[assay]]
  avg2 <- AverageExpression(obj2, return.seurat = FALSE,
                            group.by = group_by, verbose = FALSE)[[assay]]

  if (slot == "data") { avg1 <- log1p(avg1); avg2 <- log1p(avg2) }

  shared_conds <- intersect(colnames(avg1), colnames(avg2))
  if (length(shared_conds) == 0) {
    stop("No shared '", group_by, "' conditions between the two objects.")
  }

  avg1 <- as.matrix(avg1[, shared_conds, drop = FALSE])
  avg2 <- as.matrix(avg2[, shared_conds, drop = FALSE])
  shared_genes <- intersect(rownames(avg1), rownames(avg2))
  avg1 <- avg1[shared_genes, , drop = FALSE]
  avg2 <- avg2[shared_genes, , drop = FALSE]

  plot_list <- lapply(seq_along(shared_conds), function(i) {
    x_val <- as.numeric(avg1[, i])
    y_val <- as.numeric(avg2[, i])
    test  <- cor.test(x_val, y_val, method = "pearson")
    p_txt <- ifelse(test$p.value < 2.22e-16, "< 2.22e-16",
                    formatC(test$p.value, format = "e", digits = 2))

    df <- data.frame(x = x_val, y = y_val)
    ggplot(df, aes(x, y)) +
      geom_point(color = "black", alpha = 1, size = 2, stroke = 0) +
      geom_smooth(method = "glm", color = "blue", linewidth = 1,
                   se = TRUE) +
      annotate("text",
               x = min(df$x) + diff(range(df$x)) * 0.05,
               y = max(df$y) * 0.95,
               label = paste0("Pearson r: ", sprintf("%.2f", test$estimate),
                              "\np: ", p_txt),
               color = "black", size = 2.5, hjust = 0) +
      labs(x = paste("Avg. Expression —", obj1_name),
           y = paste("Avg. Expression —", obj2_name)) +
      theme_classic(base_size = 7) +
      theme(text      = element_text(color = "black"),
            axis.text = element_text(color = "black")) +
      expand_limits(x = 0, y = 0)
  })

  dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)
  pdf(pdf_file, width = pdf_width, height = pdf_height)
  invisible(lapply(plot_list, print))
  dev.off()
  message("Saved ", length(plot_list), " correlation plots to: ", pdf_file)
  invisible(pdf_file)
}

# Prepare RNA-seq object for comparison
# rna.pub: rebuilt from 6 raw featureCounts outputs (Rdata/RNAseq/), n=3
# control + n=3 inflammatory replicates -- see
# 03_downstream_analysis/R/build_rna_pub_from_featurecounts.R. Only the
# control replicates are used here (basal-expression comparison across
# platforms), selected by drug_condition rather than hardcoded cell names
# so this doesn't silently break if replicate count changes again.
# min_cells = 1 (SCTransform's default is 5): rna.pub's "cells" are bulk
# RNA-seq samples, not single cells, so a gene detected in only 1-4 of 6
# samples is real biological signal, not single-cell dropout -- the
# default filter otherwise silently drops genes like Il6 (near-zero in
# control, detected in inflammatory), which is in the 30-gene panel.
rna.sct <- SCTransform(rna.pub, vst.flavor = "v2", verbose = FALSE,
                       return.only.var.genes = FALSE, min_cells = 1)
rna.sct <- subset(rna.sct, drug_condition == COND_CONTROL)
rna.sct <- subset(rna.sct, features = GENES_ALL)

mac.sct.ctrl <- subset(mac.sct, drug_condition == COND_CONTROL)

create_correlation_plots(mac.sct.ctrl, rna.sct,
                         obj1_name = "MAC (GE-HTS)", obj2_name = "MAC (Smart-seq2)",
                         pdf_file  = file.path(OUTPUT_MAIN, "Figure2e_issVSsmart.pdf"))
create_correlation_plots(mac.sct.ctrl, macSebastian,
                         obj1_name = "MAC (GE-HTS)", obj2_name = "MAC (10x Chromium)",
                         pdf_file  = file.path(OUTPUT_MAIN, "Figure2e_issVSchromium.pdf"))
create_correlation_plots(rna.sct, macSebastian,
                         obj1_name = "MAC (Smart-seq2)", obj2_name = "MAC (10x Chromium)",
                         pdf_file  = file.path(OUTPUT_MAIN, "Figure2e_smartVSchromium.pdf"))

# =============================================================================
# Figure 2f — Lollipop chart (Primary MAC vs. ATDC5)
# plot_gene_divergence() is defined in utils.R (shared with the rna.pub-vs-
# ATDC5 supplementary version of this comparison in figure2_supplement.R).
# =============================================================================
# Primary side: mac.sct's control-condition wells (GE-HTS chip, RNA assay
# raw counts) -- mac.sct.ctrl was already computed in the Figure 2c section
# above. ATDC5 side: atdc5.pub's 37C replicates (GSE269585) -- the only
# condition in that dataset confirmed as a valid baseline counterpart (37C =
# physiological/untreated; 40C/44C are a temperature-response axis unrelated
# to this project). Both sides are CLR-normalized per sample/well from raw
# counts, as plot_gene_divergence() expects. (The rna.pub [bulk RNA-seq] vs.
# ATDC5 version of this comparison is a supplementary figure -- see
# figure2_supplement.R.)

mac.sct.ctrl.clr <- NormalizeData(mac.sct.ctrl, assay = "RNA",
                                  normalization.method = "CLR", verbose = FALSE)
atdc5.clr         <- NormalizeData(atdc5.pub, assay = "RNA",
                                   normalization.method = "CLR", verbose = FALSE)

norm_primary <- GetAssayData(mac.sct.ctrl.clr, assay = "RNA", layer = "data")
norm_atdc5   <- GetAssayData(atdc5.clr,        assay = "RNA", layer = "data")

fig_2f <- plot_gene_divergence(norm_primary, norm_atdc5)
ggsave(file.path(OUTPUT_MAIN, "Figure2f_lollipop.pdf"),
       plot = fig_2f, width = 6, height = 5, dpi = 300)

save_session_info("figure2")
