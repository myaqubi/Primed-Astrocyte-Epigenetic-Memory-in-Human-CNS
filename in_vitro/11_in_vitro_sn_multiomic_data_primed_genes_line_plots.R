# ============================================================================== 
# Donor-level multiome visualization of selected primed genes
# ============================================================================== 
#
# Purpose:
#   This script generates donor-level RNA expression and GeneActivity
#   accessibility plots for selected primed-memory genes identified in the
#   previous analysis step. 
#
# Input:
#   1. Integrated Seurat object: In_vitro_prenatal_samples_QC_filtered_integrated_annotated_astrocytes
#   2. Primed-gene results table: primed_genes_all_cells.csv
#
# Output:
#   - One RNA plot per gene
#   - One GeneActivity plot per gene
#   - A legend plot for the donor styling
#
# Notes:
#   - Replace the generic placeholder paths below with your own project paths.
# ============================================================================== 

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(grid)
  library(readr)
})

# ------------------------------------------------------------------------------
# 0) Configuration
# ------------------------------------------------------------------------------
project_dir <- "/path/to/your/project"
output_dir <- file.path(project_dir, "multiome_qPCR_plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_object_path <- file.path(project_dir, "In_vitro_prenatal_samples_QC_filtered_integrated_annotated_astrocytes.rds")
input_results_path <- file.path(project_dir, "finding_memory_genes", "TRUE_genes", "primed_genes_all_cells.csv")

if (!file.exists(input_object_path)) {
  stop("Input object not found at: ", input_object_path)
}

if (!file.exists(input_results_path)) {
  stop("Primed-gene results file not found at: ", input_results_path)
}

obj <- readRDS(input_object_path)
res_B <- read_csv(input_results_path, show_col_types = FALSE)

# Remove the first column if it is an index column
if (ncol(res_B) > 0 && names(res_B)[1] == "...1") {
  res_B <- res_B[, -1]
}

# ------------------------------------------------------------------------------
# 1) Sample and donor metadata setup
# ------------------------------------------------------------------------------
condition_col <- "condition"
donor_col <- "Donor"

if (!donor_col %in% colnames(obj@meta.data)) {
  if ("donor_global_pretty" %in% colnames(obj@meta.data)) {
    obj@meta.data[[donor_col]] <- obj@meta.data[["donor_global_pretty"]]
  } else {
    stop("No donor column found. Expected 'Donor' or 'donor_global_pretty'.")
  }
}

# Harmonize conditions to Ctrl / Stim / Prime
obj@meta.data[[condition_col]] <- as.character(obj@meta.data[[condition_col]])
obj@meta.data[[condition_col]][obj@meta.data[[condition_col]] %in% c("CTRL", "Control_no_hit")] <- "Ctrl"
obj@meta.data[[condition_col]][obj@meta.data[[condition_col]] %in% c("Single_hit_D7", "Stim")] <- "Stim"
obj@meta.data[[condition_col]][obj@meta.data[[condition_col]] %in% c("Single_hit_D1", "Prime")] <- "Prime"

cond_levels <- c("Ctrl", "Stim", "Prime")
obj@meta.data[[condition_col]] <- factor(obj@meta.data[[condition_col]], levels = cond_levels)

cond_labels <- c(
  "Ctrl"  = "Ctrl",
  "Stim"  = "Stim",
  "Prime" = "Prime"
)

# Donor colors
all_donors <- sort(unique(obj@meta.data[[donor_col]][!is.na(obj@meta.data[[donor_col]])]))
donor_colors <- c(
  "Donor1" = "red",
  "Donor2" = "green",
  "Donor3" = "blue"
)

missing_donors <- setdiff(all_donors, names(donor_colors))
if (length(missing_donors) > 0) {
  extra_colors <- grDevices::hcl.colors(length(missing_donors), palette = "Set2")
  names(extra_colors) <- missing_donors
  donor_colors <- c(donor_colors, extra_colors)
}

# ------------------------------------------------------------------------------
# 2) Gene selection

gene_list <- c("IL6", "NOS2", "BACH1", "BIRC3", "IL32", "EGR2")
genes_to_plot <- gene_list

# ------------------------------------------------------------------------------
# 3) Helper functions
# ------------------------------------------------------------------------------
clean_filename <- function(x) {
  gsub("[^A-Za-z0-9_\\-]", "_", x)
}

star_label <- function(x) {
  if (is.na(x)) return("ns")
  if (isTRUE(x)) return("*")
  return("ns")
}

get_donor_means <- function(obj, assay_name, gene, condition_col, donor_col, cond_levels) {
  DefaultAssay(obj) <- assay_name
  vals <- FetchData(obj, vars = gene, slot = "data")
  
  df <- data.frame(
    condition = obj@meta.data[[condition_col]],
    donor     = obj@meta.data[[donor_col]],
    value     = vals[[gene]],
    stringsAsFactors = FALSE
  )
  
  df <- df[!is.na(df$condition) & !is.na(df$donor), , drop = FALSE]
  df <- df[df$condition %in% cond_levels, , drop = FALSE]
  df$condition <- factor(df$condition, levels = cond_levels)
  
  aggregate(value ~ condition + donor, data = df, FUN = mean)
}

make_ann_df <- function(res_row, assay_name, cond_levels, y_max, y_min) {
  if (assay_name == "RNA") {
    ann <- data.frame(
      group1 = c(cond_levels[1], cond_levels[2]),
      group2 = c(cond_levels[2], cond_levels[3]),
      sig = c(
        res_row$RNA_sig_Stim_vs_Ctrl,
        res_row$RNA_sig_Stim_vs_Prime
      ),
      stringsAsFactors = FALSE
    )
  } else if (assay_name == "GeneActivity") {
    ann <- data.frame(
      group1 = c(cond_levels[1], cond_levels[1]),
      group2 = c(cond_levels[2], cond_levels[3]),
      sig = c(
        res_row$GA_sig_Stim_vs_Ctrl,
        res_row$GA_sig_Prime_vs_Ctrl
      ),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Unknown assay: ", assay_name)
  }
  
  ann$label <- vapply(ann$sig, star_label, character(1))
  yr <- max(y_max - y_min, 1e-6)
  ann$y.position <- y_max + c(0.12, 0.30)[seq_len(nrow(ann))] * yr
  ann
}

make_one_multiome_plot <- function(obj, gene, assay_name, assay_title, y_lab, res_row,
                                   condition_col, donor_col, cond_levels, cond_labels) {
  df <- get_donor_means(
    obj = obj,
    assay_name = assay_name,
    gene = gene,
    condition_col = condition_col,
    donor_col = donor_col,
    cond_levels = cond_levels
  )
  
  if (nrow(df) == 0) return(NULL)
  
  colnames(df)[colnames(df) == "value"] <- "donor_mean"
  df$condition <- factor(df$condition, levels = cond_levels)
  df$donor <- factor(df$donor)
  
  means_df <- aggregate(donor_mean ~ condition, data = df, FUN = mean)
  means_df$condition <- factor(means_df$condition, levels = cond_levels)
  
  ann <- make_ann_df(
    res_row = res_row,
    assay_name = assay_name,
    cond_levels = cond_levels,
    y_max = max(df$donor_mean, na.rm = TRUE),
    y_min = min(df$donor_mean, na.rm = TRUE)
  )
  
  yr <- max(max(df$donor_mean, na.rm = TRUE) - min(df$donor_mean, na.rm = TRUE), 1e-6)
  
  ggplot(df, aes(x = condition, y = donor_mean, color = donor)) +
    geom_boxplot(
      aes(group = condition),
      width = 0.35,
      outlier.shape = NA,
      color = "black",
      fill = NA,
      linewidth = 1.2
    ) +
    geom_line(aes(group = donor), linewidth = 1.4, alpha = 0.9) +
    geom_point(size = 4.5) +
    geom_point(
      data = means_df,
      aes(x = condition, y = donor_mean),
      inherit.aes = FALSE,
      shape = 18,
      size = 5.5,
      color = "black"
    ) +
    stat_pvalue_manual(
      ann,
      label = "label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 1,
      size = 9
    ) +
    scale_x_discrete(labels = cond_labels[cond_levels]) +
    expand_limits(y = max(ann$y.position, na.rm = TRUE) + 0.15 * yr) +
    scale_color_manual(values = donor_colors) +
    labs(title = assay_title, x = NULL, y = y_lab, color = "Donor") +
    theme_classic(base_size = 28) +
    theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5, color = "black"),
      axis.text.x = element_text(size = 24, face = "bold", angle = 30, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 24, face = "bold", color = "black"),
      axis.title.y = element_text(size = 26, face = "bold", color = "black"),
      legend.position = "none",
      axis.line = element_line(color = "black", linewidth = 1.2),
      axis.ticks = element_line(color = "black", linewidth = 1.2),
      axis.ticks.length = unit(0.25, "cm")
    )
}

# ------------------------------------------------------------------------------
# 4) Create a legend plot
# ------------------------------------------------------------------------------
legend_gene <- "IL6"
res_row <- res_B[res_B$gene == legend_gene, , drop = FALSE]

p_legend <- make_one_multiome_plot(
  obj = obj,
  gene = legend_gene,
  assay_name = "RNA",
  assay_title = "RNA",
  y_lab = "Expression",
  res_row = res_row,
  condition_col = condition_col,
  donor_col = donor_col,
  cond_levels = cond_levels,
  cond_labels = cond_labels
) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 22),
    legend.key.size = unit(1.2, "cm")
  )

ggsave(file.path(output_dir, "Legend_plot_IL6_RNA.png"), p_legend, width = 6, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 5) Make and save RNA and GeneActivity plots for each selected gene
# ------------------------------------------------------------------------------
rna_genes <- rownames(obj[["RNA"]])
ga_genes  <- rownames(obj[["GeneActivity"]])

saved_files <- character(0)

for (g in genes_to_plot) {
  message("Plotting ", g)
  res_row <- res_B[res_B$gene == g, , drop = FALSE]
  
  if (nrow(res_row) == 0) {
    warning("Skipping ", g, ": not found in res_B")
    next
  }
  
  res_row <- res_row[1, , drop = FALSE]
  
  if (g %in% rna_genes) {
    p_rna <- make_one_multiome_plot(
      obj = obj,
      gene = g,
      assay_name = "RNA",
      assay_title = paste0(g, " — RNA"),
      y_lab = "Expression",
      res_row = res_row,
      condition_col = condition_col,
      donor_col = donor_col,
      cond_levels = cond_levels,
      cond_labels = cond_labels
    )
    
    out_rna <- file.path(output_dir, paste0(clean_filename(g), "_RNA.png"))
    ggsave(filename = out_rna, plot = p_rna, width = 6, height = 6, dpi = 300)
    saved_files <- c(saved_files, out_rna)
  }
  
  if (g %in% ga_genes) {
    p_ga <- make_one_multiome_plot(
      obj = obj,
      assay_name = "GeneActivity",
      gene = g,
      assay_title = paste0(g, " — Accessibility"),
      y_lab = "Accessibility",
      res_row = res_row,
      condition_col = condition_col,
      donor_col = donor_col,
      cond_levels = cond_levels,
      cond_labels = cond_labels
    )
    
    out_ga <- file.path(output_dir, paste0(clean_filename(g), "_GeneActivity.png"))
    ggsave(filename = out_ga, plot = p_ga, width = 6, height = 6, dpi = 300)
    saved_files <- c(saved_files, out_ga)
  }
}

saved_files
