# GE-HTS project pipeline for in situ sequencing data analysis
# author: Nathan Wooseok Lee

library(pheatmap)

# prepare data
sin0.1.sct <- subset(sin.sct, subset = dose == 0.1)
cmb0.1.sct <- subset(cmb.sct, subset = dose1 == 0.1)
cmb0.1.sct$drug_name <- paste(cmb0.1.sct$drug_name1, cmb0.1.sct$drug_name2, sep = '&')
mac.sct$drug_name<-mac.sct$drug_condition

# cmb0.1.sct has weird XAV combi's. Let's switch them.
# Modify the paste function to swap drug_name1 and drug_name2 if drug_name1 starts with "XAV"
cmb0.1.sct$drug_name <- ifelse(
  startsWith(cmb0.1.sct$drug_name1, "XAV"),
  paste(cmb0.1.sct$drug_name2, cmb0.1.sct$drug_name1, sep = '&'),  # Swap drug_name1 and drug_name2
  paste(cmb0.1.sct$drug_name1, cmb0.1.sct$drug_name2, sep = '&')   # Keep original order
)

# merge mac and drug treated seurat objects
mac_sin0.1.sct <- merge(sin0.1.sct, mac.sct)
mac_cmb0.1.sct <- merge(cmb0.1.sct, mac.sct)


### Figure 4b (Heatmap)
create_heatmap_by_expression_n(mac_cmb0.1.sct, prob_cmb0.1_rf, pdf_file = './figures/figure4/heatmapCombi.pdf', log1p=TRUE, pdf_width =10, pdf_height = 10) 

# Draw heatmap of drug condition sorted by efficacy score or other predefined vectors. # 240802 #
create_heatmap_by_expression_n <- function(seurat_object, 
                                           drug_probabilities, # probability of efficacy of drugs: prob_sng10, prob_sng1 , prob_sng0.1, prob_cmb1, prob_cmb0.1
                                           upper_threshold = 0, # The upper range of genes that you want to display (0: highest ~ 1: lowest)
                                           bottom_threshold = 1, 
                                           assay = 'SCT', 
                                           slot = 'counts', 
                                           log1p = TRUE, 
                                           pdf_file = "./figures/heatmap.pdf", 
                                           pdf_width = 4, 
                                           pdf_height = 4) {
  # Get expression data
  expression_data <- GetAssayData(seurat_object, assay = assay, slot = slot)
 
  # Calculate average expression using normalized data
  average_expression <- rowMeans(expression_data, na.rm = TRUE)
  
  # Identify genes within the specified expression range
  selected_genes <- names(average_expression[average_expression >= quantile(average_expression, 1- bottom_threshold) & 
                                             average_expression <= quantile(average_expression, 1- upper_threshold)])
  title <- "Zoomed in Expressed Genes"
  
  # Split selected genes into anabolic and inflammatory categories
  selected_anabolic_genes <- intersect(selected_genes, anabolic)
  selected_inflammatory_genes <- intersect(selected_genes, inflammatory)
  selected_housekeeping_genes <- intersect(selected_genes, housekeeping)

  # Combine and sort genes within each category
  sorted_selected_genes <- c(
    selected_anabolic_genes[order(average_expression[selected_anabolic_genes], decreasing = TRUE)],
    selected_inflammatory_genes[order(average_expression[selected_inflammatory_genes], decreasing = TRUE)],
    selected_housekeeping_genes[order(average_expression[selected_housekeeping_genes], decreasing = TRUE)]
  )

  # Prepare cell scores (unchanged)
  cell_scores <- data.frame(
    cell_name = rownames(seurat_object@meta.data),
    drug_name = ifelse(is.na(seurat_object@meta.data$drug_name) | seurat_object@meta.data$drug_name == "", 
                       seurat_object@meta.data$drug_condition, 
                       seurat_object@meta.data$drug_name),
    stringsAsFactors = FALSE
  )
  
  # Add the scores for each drug_name
  cell_scores$score <- drug_probabilities[cell_scores$drug_name]
  
  # Separate control and inflammatory
  control_scores <- cell_scores[cell_scores$drug_name == 'control', ]
  inflammatory_scores <- cell_scores[cell_scores$drug_name == 'inflammatory', ]
  other_scores <- cell_scores[!(cell_scores$drug_name %in% c('control', 'inflammatory')), ]
  
  # Arrange the other scores by score in descending order
  other_scores <- other_scores %>%
    arrange(desc(score))
  
  # Combine the scores with control first and inflammatory last
  cell_scores <- rbind(control_scores, other_scores, inflammatory_scores)

  # Function to create heatmap for a given set of genes
  create_heatmap <- function(genes, normalized_data, cell_scores, title, log1p) {
    data_for_heatmap <- normalized_data[genes, ]
    # data_for_heatmap <- remove_outliers(data_for_heatmap)

    if (ncol(data_for_heatmap) == 0) {
      stop("No data available after removing outlier samples.")
    }

    if (log1p) {
        # Log-transform the data
        data_for_heatmap <- log1p(data_for_heatmap)
    }
    # Calculate color breaks
    data_range <- range(data_for_heatmap)
    breaks <- seq(data_range[1], data_range[2], length.out = 101)
    color_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(length(breaks) - 1)
    

    # Update cell_scores to match the remaining columns in data_for_heatmap
    cell_scores <- cell_scores[cell_scores$cell_name %in% colnames(data_for_heatmap), ]
    data_for_heatmap_ordered <- data_for_heatmap[, cell_scores$cell_name]
    
    # Ensure there are no non-finite values in the ordered data
    data_for_heatmap_ordered <- data_for_heatmap_ordered[, apply(data_for_heatmap_ordered, 2, function(x) all(is.finite(x)))]
    
    if (ncol(data_for_heatmap_ordered) == 0) {
      stop("No data available after removing non-finite values.")
    }

    annotation_data <- data.frame(Drug = cell_scores$drug_name[match(colnames(data_for_heatmap_ordered), cell_scores$cell_name)])
    rownames(annotation_data) <- colnames(data_for_heatmap_ordered)
    annotation_data$Drug <- factor(annotation_data$Drug, levels = unique(annotation_data$Drug))
    
    row_annotation <- data.frame(Gene_Type = factor(
      rep(c("Anabolic", "Catabolic", "Housekeeping"), c(length(selected_anabolic_genes), length(selected_inflammatory_genes), length(selected_housekeeping_genes))),
      levels = c("Anabolic", "Catabolic", "Housekeeping")
    ))

    
    # Create custom color palette for drugs
    n_drugs <- length(unique(cell_scores$drug_name))
    drug_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(max(n_drugs, 13))
    names(drug_colors) <- unique(cell_scores$drug_name)
    
    
    ann_colors <- list(
      Gene_Type = c(Anabolic = "forestgreen", Catabolic = "orange", Housekeeping = "grey"),
      Drug = drug_colors
    )
    
    rownames(row_annotation) <- genes
    # Save the heatmap as a PDF with specified font size and style
    pdf(file = pdf_file, width = pdf_width, height = pdf_height, family = "Helvetica")
    pheatmap(data_for_heatmap_ordered, 
             cluster_rows = FALSE, 
             cluster_cols = FALSE,
             show_rownames = TRUE,
             show_colnames = FALSE,
             annotation_col = annotation_data,
             annotation_row = row_annotation,
             annotation_colors = ann_colors,
             color = color_palette,
             main = title,
             fontsize = 5,  # Set font size to 5pt
             fontsize_row = 5,  # Row names font size
             fontsize_col = 5,  # Column names font size
             gaps_row = c(length(selected_anabolic_genes), length(selected_anabolic_genes) + length(selected_inflammatory_genes)))
    dev.off()  # Close the PDF device
  }
  
  # Create heatmap for the selected genes
  create_heatmap(sorted_selected_genes, expression_data, cell_scores, title, log1p=log1p)
}



### Figure 4c (Bliss model heatmap)
library(dplyr)
library(reshape2) # melt included
library(viridis)

generate_bliss_pdf(prob_sng0.1_rf, prob_cmb0.1_rf, dose = '0.1', pdf_file = './figures/figure4/bliss_rf_260110.pdf', pdf_width = 15, pdf_height = 4.2) 

# Define a function to process and plot heatmaps as a PDF
generate_bliss_pdf <- function(prob_sng0.1, prob_cmb0.1, dose, pdf_file = "bliss.pdf", 
                                 pdf_width = 15, pdf_height = 8) {
  # Step 1: Filter out 'control' and 'inflammatory'
  '%notin%' <- Negate('%in%')
  if (is.null(names(prob_sng0.1))) {
  stop("Error: prob_sng0.1 has no names.")
  }
  if (is.null(names(prob_cmb0.1))) {
  stop("Error: prob_cmb0.1 has no names.")
  }
  prob_sng0.1_filtered <- prob_sng0.1[names(prob_sng0.1) %notin% c('control', 'inflammatory')]
  prob_cmb0.1_filtered <- prob_cmb0.1[names(prob_cmb0.1) %notin% c('control', 'inflammatory')]

  # Step 2: Convert the combination probabilities into a data frame
  prob_cmb0.1_df <- convert_probvec2df(prob_cmb0.1_filtered, dose)

  # Step 3: Create triangular matrix for observed data
  list_of_matrices <- list()
  list_of_matrices[[dose]] <- create_triangular_matrix(prob_cmb0.1_df, dose)

  # Step 4: Create matrices of predicted combination efficacies
  predicted_matrices <- list()
  predicted_matrices[[dose]] <- create_predicted_matrix(prob_sng0.1_filtered, prob_cmb0.1_df)

  # Step 5: Calculate the differences between observed and predicted matrices
  difference_matrices <- list()
  difference_matrices[[dose]] <- calculate_difference_matrices(list_of_matrices[[1]], predicted_matrices[[1]])
 
  # Step 6: Find shared limits for plotting
  all_scores <- c(list_of_matrices, predicted_matrices)
  score_shared_limits <- find_shared_limits(all_scores)
  bliss_shared_limits <- find_shared_limits(difference_matrices)

  # Step 7: Generate heatmaps
  p01 <- plot_heatmap(list_of_matrices[[dose]], paste("Observed efficacy for dose", dose, "uM"), score_shared_limits)
  p02 <- plot_heatmap(predicted_matrices[[dose]], paste("Predicted efficacy for dose", dose, "uM"), score_shared_limits)
  p03 <- plot_heatmap_dif(difference_matrices[[dose]], paste("Efficacy difference for dose", dose, "uM"), bliss_shared_limits)

  # Step 8: Save all plots into a single PDF
  pdf(pdf_file, width = pdf_width, height = pdf_height)
  print(p01+p02+p03)
  dev.off()

  # Return the PDF file path
  return(pdf_file)
}

# Define a function to separate drug names and group by them
convert_probvec2df <- function(element_vector, dose) {
  # Step 1: Convert the vector into a data frame
  results_df <- data.frame(
    Combination = names(element_vector),
    MeanProbability = as.numeric(element_vector),
    Dose = dose,
    stringsAsFactors = FALSE
  )

  # Step 2: Separate 'Combination' into 'Drug1Name' and 'Drug2Name'
  results_df <- results_df %>%
    mutate(
      Drug1Name = sapply(strsplit(Combination, "&"), `[`, 1),
      Drug2Name = sapply(strsplit(Combination, "&"), `[`, 2)
    )

  # Return the grouped data frame
  return(results_df)
}


# Function to create a triangular matrix for a given dose
create_triangular_matrix <- function(data, dose) {
  # Filter data for the given dose
  dose_data <- data %>% filter(Dose == dose)
  if (nrow(dose_data) == 0) {
    warning(paste("No combination data found for dose =", dose))
    return(matrix(NA, 0, 0))  # or NULL if you'd rather halt downstream
  }

  # Remove rows where Drug1Name is the same as Drug2Name
  dose_data <- dose_data %>% filter(Drug1Name != Drug2Name)
  # Get unique drug names
  unique_drugs <- sort(unique(c(dose_data$Drug1Name, dose_data$Drug2Name)))

  # Create an empty matrix
  prob_matrix <- matrix(NA, nrow = length(unique_drugs), ncol = length(unique_drugs))
  rownames(prob_matrix) <- unique_drugs
  colnames(prob_matrix) <- unique_drugs

  # Populate the matrix
  for (i in 1:nrow(dose_data)) {
    drug1 <- dose_data$Drug1Name[i]
    drug2 <- dose_data$Drug2Name[i]
    prob <- dose_data$MeanProbability[i]
    
    # Ensure matrix is triangular
    if (which(unique_drugs == drug1) <= which(unique_drugs == drug2)) {
      prob_matrix[drug1, drug2] <- prob
    } else {
      prob_matrix[drug2, drug1] <- prob
    }
  }

  return(prob_matrix)
}


# Function to calculate Bliss predicted efficacy
bliss_predict <- function(efficacy1, efficacy2) {
  efficacy1 + efficacy2 - (efficacy1 * efficacy2)
}

# Function to create a matrix of predicted combination efficacies
create_predicted_matrix <- function(single_dose_data, combination_data) {
  unique_drugs <- sort(unique(c(combination_data$Drug1Name, combination_data$Drug2Name)))
  predicted_matrix <- matrix(NA, nrow = length(unique_drugs), ncol = length(unique_drugs))
  rownames(predicted_matrix) <- unique_drugs
  colnames(predicted_matrix) <- unique_drugs
  
  for (drug1 in unique_drugs) {
    for (drug2 in unique_drugs) {
      if (which(unique_drugs == drug1) < which(unique_drugs == drug2)) {
        efficacy1 <- single_dose_data[drug1]
        efficacy2 <- single_dose_data[drug2]
        predicted_matrix[drug1, drug2] <- bliss_predict(efficacy1, efficacy2)
      }
    }
  }
  
  return(predicted_matrix)
}

# Function to calculate the differences between observed and predicted efficacies
calculate_difference_matrices <- function(observed_matrix, predicted_matrix) {
  difference_matrix <- observed_matrix - predicted_matrix
  return(difference_matrix)
}

find_shared_limits <- function(matrices) {
  overall_min <- Inf
  overall_max <- -Inf
  
  for (matrix in matrices) {
    matrix_min <- min(matrix, na.rm = TRUE)
    matrix_max <- max(matrix, na.rm = TRUE)
    
    if (matrix_min < overall_min) {
      overall_min <- matrix_min
    }
    
    if (matrix_max > overall_max) {
      overall_max <- matrix_max
    }
  }
  
  return(c(overall_min, overall_max))
}

# Draw heatmap for a triangle matrix
# need to make the scale the same across matirces
plot_heatmap <- function(matrix, title, value_limit) {
  long_format <- melt(matrix, na.rm = TRUE, varnames = c("Drug1", "Drug2"))
  unique_drugs <- rownames(matrix)
  ggplot(long_format, aes(x = Drug1, y = Drug2, fill = value)) +
    geom_tile() +
    scale_x_discrete(limits = unique_drugs) +
    scale_y_discrete(limits = unique_drugs) +
    scale_fill_viridis(option='viridis', na.value = "white", limits=value_limit) +
    labs(title = title, fill = "Mean value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(angle = 45, hjust = 1))
}

# Draw heatmap for a triangle matrix of which values are diverging at 0.
# need to make the scale the same across matirces
plot_heatmap_dif <- function(matrix, title, value_limit) {
  long_format <- melt(matrix, na.rm = FALSE, varnames = c("Drug1", "Drug2"))
  unique_drugs <- rownames(matrix)
  value_bottom <- value_limit[1]
  value_top <- value_limit[2]
  ggplot(long_format, aes(x = Drug1, y = Drug2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "#0000CC", high = "#CC0000", mid = "white", midpoint = 0, 
                         space = "Lab", na.value = NA, limits = c(value_bottom, value_top)) +
    labs(title = title, fill = "Difference") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.text.y = element_text(angle = 45, hjust = 1))
}


### Figure 4d 
# Convert seurat data to anndata for python
sceasy::convertFormat(mac_cmb0.1.sct, from = "seurat", to = "anndata", outFile = "./Rdata/pub/mac_cmb0.1.sct.h5ad", main_layer = "counts", assay = "SCT")
sceasy::convertFormat(mac_sin0.1.sct, from = "seurat", to = "anndata", outFile = "./Rdata/pub/mac_sin0.1.sct.h5ad", main_layer = "counts", assay = "SCT")

# Run figure4d.ipynb code for plotting.


### Figure 4e (radar chart)
library(scales)
library(fmsb)

# Define Seurat objects and corresponding drug lists
all0.1.sct.list <- list(mac.sct, sin0.1.sct, cmb0.1.sct)

# Get all unique combination strings from the Seurat object metadata
all_combos <- unique(cmb0.1.sct$drug_name)
drug_cmb <- all_combos[sapply(strsplit(all_combos, "&"), function(x) x[1] != x[2])]# Filter out self-combinations (where Drug A == Drug B)
drug_single <- unique(sin0.1.sct$drug_name)
drug_mac  <- c("control", "inflammatory")
drug_lists <- list(drug_mac, drug_single, drug_cmb)


# Define a scale of radar
uni_scale <- calculate_universal_scale_for_radar(all0.1.sct.list, drug_lists, hknorm = FALSE) # make sure you followed all sequence in this code file.
# If you run into error: No cell found, you did not followed the meta.data setting at the begining of the code.

radar_df_mac_cmb0.1<- prep_radar_chart(mac_cmb0.1.sct, c(drug_mac, drug_cmb), uni_scale, hknorm = FALSE)
radar_df_mac_sin0.1<- prep_radar_chart(mac_sin0.1.sct, c(drug_mac, drug_single), uni_scale, hknorm = FALSE)


# The example drug
selected_drug <- "SB431542&SANT-1"

# Concatenate selected drug with control and inflammatory  (selection from top and bottom 10% of efficacy & bliss)
selected_drug_cmb0.1 <- c("SB431542&SANT-1", "rapamycin&KU0063794", "SB525334&XAV", "rapamycin&SANT-1", "MK-2206 dihydrochloride&pamapimod", "KU0063794&BMS-345541", "ALK5 inhibitor IV&MK-2206 dihydrochloride", "KU0063794&ALK5 inhibitor IV")

selected_drug_sin0.1 <- c("ALK5 inhibitor IV", "KU0063794", "JNK inhibitor V", "SANT-1", "KU0063794", "LY294002", "SB525334", "XAV", "SB203580", "SB431542", "pamapimod", "MK-2206 dihydrochloride",
                          "BMS-345541")

selected_drug_mac_cmb0.1 <- c(drug_mac, selected_drug_cmb0.1)
selected_drug_mac_sin0.1 <- c(drug_mac, selected_drug_sin0.1)
# Prep dataframe before drawing radar charts
radar_df_selected <- prep_radar_chart(mac_cmb0.1.sct, selected_drug_mac_cmb0.1, uni_scale, hknorm = FALSE)
radar_df_mac_sin0.1<- prep_radar_chart(mac_sin0.1.sct, selected_drug_mac_sin0.1, uni_scale, hknorm = FALSE)

plot_radar_chart_A_drug(radar_df_selected, selected_drug)

# Save pdf files
plot_radar_chart_A_drug_pdf(radar_df_mac_cmb0.1, selected_drug, 
                            pdf_file = "./figures/figure4/radar_chart_cmb.pdf", width = 4, height = 4) 
plot_radar_chart_A_drug_pdf(radar_df_mac_sin0.1, 'SB431542', 
                            pdf_file = "./figures/figure4/radar_chart_SB431542.pdf", width = 4, height = 4) 
plot_radar_chart_A_drug_pdf(radar_df_mac_sin0.1, 'SANT-1', 
                            pdf_file = "./figures/figure4/radar_chart_SANT-1.pdf", width = 4, height = 4) 


# # for loop for drawing all drugs in dataframe, if you need.
# for (drug in selected_drug_sin0.1) {plot_radar_chart_A_drug(radar_df_mac_sin0.1, drug)}

# Get the top 10% combinations (returns a dataframe)
top_hits <- extract_top_combinations(
    prob_sng = prob_sng0.1_rf, 
    prob_cmb = prob_cmb0.1_rf, 
    dose = '0.1'
)

# Save to CSV if needed for downstream tools
write.csv(top_hits, "./figures/supple/top_synergistic_combinations.csv", row.names = FALSE)

extract_top_combinations <- function(prob_sng, prob_cmb, dose, top_percent = 0.10) {
  
  # --- Step 1: Data Cleaning (Same as your generate_bliss_pdf) ---
  '%notin%' <- Negate('%in%')
  
  # Filter out controls
  prob_sng_clean <- prob_sng[names(prob_sng) %notin% c('control', 'inflammatory')]
  prob_cmb_clean <- prob_cmb[names(prob_cmb) %notin% c('control', 'inflammatory')]
  
  # Convert combination vector to dataframe
  prob_cmb_df <- convert_probvec2df(prob_cmb_clean, dose)
  
  # --- Step 2: Calculate Matrices ---
  # 1. Observed Efficacy Matrix (Upper triangular)
  matrix_observed <- create_triangular_matrix(prob_cmb_df, dose)
  
  # 2. Predicted Efficacy Matrix (Bliss)
  matrix_predicted <- create_predicted_matrix(prob_sng_clean, prob_cmb_df)
  
  # 3. Synergy Matrix (Observed - Predicted)
  matrix_synergy <- calculate_difference_matrices(matrix_observed, matrix_predicted)
  
  # --- Step 3: Flatten Matrices to Data Frame ---
  # We use melt to turn the matrices into lists of drug pairs
  # na.rm = TRUE removes the empty lower triangle
  df_observed <- melt(matrix_observed, na.rm = TRUE, varnames = c("Drug1", "Drug2"), value.name = "Observed_Efficacy")
  df_synergy  <- melt(matrix_synergy,  na.rm = TRUE, varnames = c("Drug1", "Drug2"), value.name = "Bliss_Synergy")
  
  # Merge them by Drug names to ensure alignment
  combined_df <- inner_join(df_observed, df_synergy, by = c("Drug1", "Drug2"))
  
  # --- Step 4: Filter Top 10% Efficacy & Sort by Synergy ---
  
  # Calculate the threshold for the top X% of efficacy
  efficacy_threshold <- quantile(combined_df$Observed_Efficacy, probs = (1 - top_percent))
  
  top_combinations <- combined_df %>%
    # Filter: Keep only rows where Observed Efficacy is in the top 10%
    filter(Observed_Efficacy >= efficacy_threshold) %>%
    # Sort: Highest Bliss Synergy first
    arrange(desc(Bliss_Synergy)) %>%
    mutate(
        Dose = dose,
        Rank = row_number() # Add a rank column for convenience
    )
    
  return(top_combinations)
}


calculate_universal_scale_for_radar <- function(seurat_objs, drug_lists, assay = "SCT", slot = "counts", hknorm = FALSE) {
  # Initialize lists to store the metrics for references

  # Define the genes associated with each GO term
  go_terms_genes <- list(
    Chondrocyte_Development = c("Matn1", "Acan", "Sox9", "Col27a1"),
    Cartilage_Condensation = c("Col2a1", "Acan", "Sox9"),
    NO_Synthase_Upregulation = c("Tlr2", "Ccl2"),
    ECM_Disassembly = c("Adamts5", "Mmp13"),
    Inflammatory_response = c("Tlr2", "Ccl2", "Il17b", "Tnfrsf1b", "Cxcl1", "Il6", "Cxcl5", "Fosl2"),
    Tissue_homeostasis = c("Sox9", "Col2a1", "Pth1r", "Fosl2")
  )

  min_values <- list()
  max_values <- list()
  i <- 0
  for (seurat_obj in seurat_objs) {
    i <- i + 1
    drug_names <- drug_lists[[i]]  # Use the specific drug list for each Seurat object
    
    # Initialize a list to store metrics for all drugs
    all_metrics <- list()

    # Calculate the metric scores for each drug
    for (drug in drug_names) {
      message(paste("Processing drug:", drug, "in Seurat object index", i))  # Debugging message
      
      # Subset the Seurat object by the current drug
      subset_obj_d <- subset(seurat_obj, subset = drug_name == drug)
      
      # Check if the subset has no cells
      if (ncol(subset_obj_d) == 0) {
        warning(paste("No cells found for drug:", drug, "in Seurat object index", i))
        next  # Skip to the next drug if no cells are found
      }

      mat <- GetAssayData(subset_obj_d, assay = assay, slot = slot)
      # Calculate the metric scores for each GO term
      metrics <- sapply(go_terms_genes, function(genes) {
        available_genes <- intersect(genes, rownames(mat))
        if (length(available_genes) > 0) {
          if (hknorm) {
            hkmean <- colMeans(mat[housekeeping, , drop = FALSE])
            mean(rowMeans(mat[available_genes, , drop = FALSE] / hkmean))
          } else {
            mean(rowMeans(mat[available_genes, , drop = FALSE]))
          }
        } else {
          NA
        }
      })
      
      all_metrics[[drug]] <- metrics
    }

    # Convert the list to a data frame
    all_metrics_df <- do.call(rbind, all_metrics)
    # rownames(all_metrics_df) <- drug_names

    # Determine the min and max values across all drugs for each metric
    min_values[[i]] <- apply(all_metrics_df, 2, min, na.rm = TRUE)
    max_values[[i]] <- apply(all_metrics_df, 2, max, na.rm = TRUE)
  }

  all_min_values <- do.call(rbind, min_values)
  all_max_values <- do.call(rbind, max_values)

  # Calculate universal min/max
  uni_min_values <- apply(all_min_values, 2, min, na.rm = TRUE)
  uni_max_values <- apply(all_max_values, 2, max, na.rm = TRUE)

  return(list(min_values = uni_min_values, max_values = uni_max_values))
}

prep_radar_chart <- function(seurat_obj, drug_names = NULL, uni_scale = NULL, assay = "SCT", slot = "counts", hknorm = FALSE) {
  
  # 1. Determine the correct metadata column name once
  if ("drug_name" %in% colnames(seurat_obj@meta.data)) {
    meta_col <- "drug_name"
  } else if ("drug_condition" %in% colnames(seurat_obj@meta.data)) {
    meta_col <- "drug_condition"
  } else {
    stop("Neither 'drug_name' nor 'drug_condition' found in Seurat object's meta.data.")
  }
  
  # 2. If no drug names provided, grab all from the identified column
  if (is.null(drug_names)) {
    drug_names <- unique(seurat_obj@meta.data[[meta_col]])
  }
  
  # Cache available drugs from metadata for fast lookup
  available_drugs_in_meta <- unique(seurat_obj@meta.data[[meta_col]])
  
  # Define the genes associated with each GO term
  go_terms_genes <- list(
    Chondrocyte_Development = c("Matn1", "Acan", "Sox9", "Col27a1"),
    Cartilage_Condensation = c("Col2a1", "Acan", "Sox9"),
    NO_Synthase_Upregulation = c("Tlr2", "Ccl2"),
    ECM_Disassembly = c("Adamts5", "Mmp13"),
    Inflammatory_response = c("Tlr2", "Ccl2", "Il17b", "Tnfrsf1b", "Cxcl1", "Il6", "Cxcl5", "Fosl2"),
    Tissue_homeostasis = c("Sox9", "Col2a1", "Pth1r", "Fosl2")
  )
  
  all_metrics <- list()
  
  for (drug in drug_names) {
    print(paste("Processing requested drug:", drug))
    
    # Handle A&B vs B&A 
    target_drug <- drug # Default to what was requested
    
    # If the drug isn't in metadata exactly as is...
    if (!(drug %in% available_drugs_in_meta) && grepl("&", drug)) {
      # Split by '&' and try reversing
      parts <- strsplit(drug, "&")[[1]]
      if (length(parts) == 2) {
        rev_drug <- paste(parts[2], parts[1], sep = "&")
        
        # If the reverse exists, use it!
        if (rev_drug %in% available_drugs_in_meta) {
          print(paste("  -> Exact match not found. Switching to found reverse pair:", rev_drug))
          target_drug <- rev_drug
        }
      }
    }
    # Subset cells based on the resolved 'target_drug' name
    # We use 'cells' argument which is safer than 'subset' for string variables
    cells_to_keep <- colnames(seurat_obj)[seurat_obj@meta.data[[meta_col]] == target_drug]
    
    if (length(cells_to_keep) == 0) {
      warning(paste("No cells found for drug condition:", drug))
      next
    }
    
    subset_obj_d <- subset(seurat_obj, cells = cells_to_keep)
    
    # Process data
    mat <- GetAssayData(subset_obj_d, assay = assay, slot = slot)
    
    metrics <- sapply(go_terms_genes, function(genes) {
      available_genes <- intersect(genes, rownames(mat))
      if (length(available_genes) > 0) {
        if (hknorm) {
          # Ensure 'housekeeping' is defined in your global environment
          hkmean <- colMeans(mat[housekeeping, , drop = FALSE])
          mean(rowMeans(mat[available_genes, , drop = FALSE]/hkmean))
        } else {
          mean(rowMeans(mat[available_genes, , drop = FALSE]))
        }
      } else {
        NA
      }
    })
    
    # Store result using the ORIGINAL requested name 'drug' (not target_drug)
    # This keeps your output keys consistent with your input list
    all_metrics[[drug]] <- metrics
  }
  
  # Consolidate results
  if (length(all_metrics) == 0) {
    stop("No valid metrics calculated for any drug.")
  }
  
  all_metrics_df <- do.call(rbind, all_metrics)
  rownames(all_metrics_df) <- names(all_metrics)
  
  if (is.null(uni_scale)) {
    min_values <- apply(all_metrics_df, 2, min, na.rm = TRUE)
    max_values <- apply(all_metrics_df, 2, max, na.rm = TRUE)
  } else {
    min_values <- uni_scale$min_values
    max_values <- uni_scale$max_values
  }
  
  min_vals <- as.numeric(min_values)
  max_vals <- as.numeric(max_values)
  
  norm1 <- sweep(all_metrics_df, 2, min_vals, FUN = "-")
  # Avoid division by zero if max == min
  diff_vals <- max_vals - min_vals
  diff_vals[diff_vals == 0] <- 1 
  
  normalized_metrics_df <- sweep(norm1, 2, diff_vals, FUN = "/")
  
  return(normalized_metrics_df)
}


plot_radar_chart_A_drug_pdf <- function(normalized_metrics_df, selected_drug, 
                                    reference1 = 'control', reference2 = 'inflammatory', 
                                    pdf_file = "radar_chart.pdf", width = 4, height = 4) {
  # Get reference points for control and inflammatory
  ref1_values <- normalized_metrics_df[reference1, ]
  ref2_values <- normalized_metrics_df[reference2, ]
  
  # Select the specified drug condition
  if (!(selected_drug %in% rownames(normalized_metrics_df))) {
    stop("Selected drug condition not found in the dataset.")
  }
  selected_values <- normalized_metrics_df[selected_drug, ]

  # Combine reference points with the selected drug's metrics
  final_metrics_df <- rbind(rep(1, ncol(normalized_metrics_df)), 
                            rep(0, ncol(normalized_metrics_df)), 
                            ref1_values, 
                            ref2_values, 
                            selected_values)

  # Prepare data frame for radar chart
  plot_df <- as.data.frame(final_metrics_df[c(1:5), , drop = FALSE])
  colnames(plot_df) <- gsub("_", "\n", colnames(plot_df))  # Format column names
  
  # Save radar chart as a PDF
  pdf(pdf_file, width = width, height = height)
  par(mar = c(1, 1, 1, 1))
  radarchart(plot_df, axistype = 1,
             pcol = c("green", "red", "blue"),  # Colors for the reference and the drug condition
             pfcol = c(scales::alpha("green", 0.1), scales::alpha("red", 0.1), scales::alpha("blue", 0.5)),
             plwd = c(1, 1, 2), plty = c(2, 2, 1),  # Line styles and widths
             cglcol = "grey", cglty = 1, axislabcol = "grey", caxislabels = seq(0, 1, 0.5),
             cglwd = 0.8, vlcex = 0.8,
             title = paste("Radar Chart for", selected_drug))
  dev.off()
  
  message("Radar chart saved to ", pdf_file)
}

plot_radar_chart_A_drug <- function(normalized_metrics_df, selected_drug, 
                                    reference1 = 'control', reference2 = 'inflammatory') {
  # Get reference points for control and inflammatory
  ref1_values <- normalized_metrics_df[reference1, ]
  ref2_values <- normalized_metrics_df[reference2, ]
  
  # Select the specified drug condition
  if (!(selected_drug %in% rownames(normalized_metrics_df))) {
    stop("Selected drug condition not found in the dataset.")
  }
  selected_values <- normalized_metrics_df[selected_drug, ]

  # Combine reference points with the selected drug's metrics
  final_metrics_df <- rbind(rep(1, ncol(normalized_metrics_df)), 
                            rep(0, ncol(normalized_metrics_df)), 
                            ref1_values, 
                            ref2_values, 
                            selected_values)

  # Prepare data frame for radar chart
  plot_df <- as.data.frame(final_metrics_df[c(1:5), , drop = FALSE])
  colnames(plot_df) <- gsub("_", "\n", colnames(plot_df))  # Format column names
  
  par(mar = c(1, 1, 1, 1))
  radarchart(plot_df, axistype = 1,
             pcol = c("green", "red", "blue"),  # Colors for the reference and the drug condition
             pfcol = c(scales::alpha("green", 0.1), scales::alpha("red", 0.1), scales::alpha("blue", 0.5)),
             plwd = c(1, 1, 2), plty = c(2, 2, 1),  # Line styles and widths
             cglcol = "grey", cglty = 1, axislabcol = "grey", caxislabels = seq(0, 1, 0.5),
             cglwd = 0.8, vlcex = 0.8,
             title = paste("Radar Chart for", selected_drug))
}


### Figure 4f (2d GO map)
# GO sum for 2d coordinate display
# GO terms divided into two groups: anabolism and catabolism
library(ggplot2)

coordinates_GO_df_GSVA <- function(norm_cmb, norm_sng) {
    library(dplyr)

    # 2) assemble into one data frame
    df_cmb <- as_tibble(norm_cmb,    rownames = "drug") %>% mutate(rsc = "combi")
    df_sng <- as_tibble(norm_sng,    rownames = "drug") %>% mutate(rsc = "single")
    df_all <- bind_rows(df_cmb, df_sng)

    # 3) compute the raw 2 axes
    ana_terms <- c("Chondrocyte_Development",
                "Cartilage_Condensation",
                "Tissue_homeostasis")
    cat_terms <- c("NO_Synthase_Upregulation",
                "ECM_Disassembly",
                "Inflammatory_response")

    df_all <- df_all %>%
    rowwise() %>%
    mutate(
        anabolism = mean(c_across(all_of(ana_terms)),  na.rm = TRUE),
        catabolism = mean(c_across(all_of(cat_terms)), na.rm = TRUE)
    ) %>%
    ungroup()

    # 4) pull out the inflammatory reference coords
    ref <- df_all %>% filter(drug == "inflammatory") %>% 
    select(catabolism, anabolism)

    df_all <- df_all %>%
    # for combos, split on "&" and keep only those where the two parts differ
    filter(
        !(rsc == "combi" &
        sapply(strsplit(drug, "&"), function(x) length(x)==2 && x[1] == x[2]))
    )

    # 5) subtract to recenter on inflammatory = (0,0)
    df_all <- df_all %>%
    filter(!drug %in% c("control","inflammatory")) %>%
    mutate(
        cat_rel = catabolism - ref$catabolism[1],
        ana_rel = anabolism - ref$anabolism[1]
    )
  return(df_all)
}

df_to_pdf_GSVA <- function(GOsum_df, highlight1 = NULL, highlight2 = NULL, 
                                  x_scale = NULL, y_scale = NULL, 
                                  pdf_file = "coordinates_GO_plot.pdf", 
                                  width = 7, height = 7){  
  library(ggplot2)

  # Step 4: Set x and y scale based on provided arguments or data range if NULL
  if (is.null(x_scale)) {
    x_scale <- range(GOsum_df$cat_rel, na.rm = TRUE)
  }
  
  if (is.null(y_scale)) {
    y_scale <- range(GOsum_df$ana_rel, na.rm = TRUE)
  }
  
  # Ensure the scales include 0 to make axes cross at (0,0)
  x_scale <- range(c(0, x_scale))
  y_scale <- range(c(0, y_scale))
  
  # Step 5: Annotate the data points for highlighting specific drugs
  GOsum_df$highlight <- ifelse(GOsum_df$drug %in% highlight1, "Highlight1",
                         ifelse(GOsum_df$drug %in% highlight2, "Highlight2", "Other"))
  
  # Step 6: Create the plot
  p <- ggplot(GOsum_df, aes(x = cat_rel, y = ana_rel, color = highlight, shape = rsc)) +
    geom_point(size = 3) +   # Scatter plot with shapes determined by 'rsc' and colors by 'highlight'
    scale_x_continuous(limits = x_scale, expand = c(0, 0)) +  # x-axis crosses at 0
    scale_y_continuous(limits = y_scale, expand = c(0, 0)) +  # y-axis crosses at 0
    theme_minimal() +  # Use a clean theme
    labs(x = "Catabolism (relative to reference)", y = "Anabolism (relative to reference)", 
         title = "GO Term Comparison Plot") +  # Label axes and add title
    geom_hline(yintercept = 0, color = "black") +  # Add x-axis line at y=0
    geom_vline(xintercept = 0, color = "black") +  # Add y-axis line at x=0
    scale_color_manual(values = c("Highlight1" = "red", "Highlight2" = "blue", "Other" = "grey")) +  # Colors for highlighted drugs
    scale_shape_manual(values = c("single" = 16, "combi" = 17))  # Shapes for "single" and "combi" types
  
  # Step 7: Add drug names for highlighted points
  p <- p + geom_text(
    data = GOsum_df,
    aes(label = ifelse(highlight != "Other", as.character(drug), "")),  # Show labels only for highlighted points
    vjust = -1, hjust = 0.5, size = 3, check_overlap = FALSE  # Adjust the position and size of the labels
  )
  
  # Step 8: Save the plot to a PDF
  pdf(pdf_file, width = width, height = height)
  print(p)
  dev.off()
  
  message("Plot saved to ", pdf_file)
}

# Example Usage
# prepare datasets
df_all <- coordinates_GO_df_GSVA(
  norm_cmb           = radar_df_mac_cmb0.1,
  norm_sng           = radar_df_mac_sin0.1
)


# highlight a combinations and its composers.
df_to_pdf_GSVA(
  GOsum_df           = df_all,
  highlight1         = c("SB431542&SANT-1"), 
  highlight2         = c("SB431542", "SANT-1"), 
  pdf_file           = "./figures/figure4/GO_GSVA_2d_test26.pdf",
  width              = 5,
  height             = 4
)

# 2) Identify the Q2 (x < 0, y > 0) combinations
library(dplyr)
q2_combos <- df_all %>%
  filter(rsc == "combi", cat_rel < 0, ana_rel > 0) %>%
  pull(drug)

# 3) Extract their single-drug parts
comp_singles <- q2_combos %>%
  strsplit("&") %>%
  unlist() %>%
  unique()

# 4) Zoomed Q2 only plot
df_to_pdf_GSVA(
  GOsum_df           = df_all,
  highlight1         = q2_combos,     # label only these combos
  highlight2         = NULL,          # no single labels here
  x_scale            = c(min(df_all$cat_rel[df_all$cat_rel < 0]) * 1.1, 0),
  y_scale            = c(0, max(df_all$ana_rel[df_all$ana_rel > 0]) * 1.1),
  pdf_file           = "./figures/figure4/GO_GSVA_2d_Q2.pdf",
  width              = 5,
  height             = 4
)

# 5) Whole plot with component singles colored
df_to_pdf_GSVA(
  GOsum_df           = df_all,  
  highlight1         = q2_combos,     # combos in red
  highlight2         = comp_singles,  # their singles in blue
  pdf_file           = "./figures/figure4/GO_GSVA_2d_Q2_components.pdf",
  width              = 5,
  height             = 4
)

