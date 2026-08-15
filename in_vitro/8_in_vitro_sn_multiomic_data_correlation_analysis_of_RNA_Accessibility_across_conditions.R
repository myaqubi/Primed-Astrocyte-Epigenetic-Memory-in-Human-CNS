# ============================================================================== 
# Correlation analysis of RNA and accessibility across conditions
# ============================================================================== 
#
# Purpose:
#   This script compares Ctrl, Prime, and Stim conditions at both the
#   transcriptomic (RNA) and chromatin accessibility (ATAC) levels. It uses:
#   - top 3000 highly variable genes (RNA), and
#   - variable peaks from the ATAC assay.
#
# The analysis computes pseudobulk profiles by donor-condition and generates:
#   1) condition-level scatter plots comparing Ctrl vs Prime and Ctrl vs Stim,
#   2) Pearson correlations for RNA and ATAC,
#   3) donor-level correlation summaries.
#
# Notes:
#   - Replace the generic placeholder paths below with your own project paths.
# ============================================================================== 

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# 0) Configuration
# ------------------------------------------------------------------------------
project_dir <- "/path/to/your/project"
output_dir <- file.path(project_dir)
output_file <- file.path(output_dir, "Ctrl_vs_Stim_Prime_pseudobulk_correlation_RNA_ATAC.png")

input_object_path <- file.path(project_dir, "In_vitro_prenatal_samples_QC_filtered_integrated_annotated_astrocytes.rds")

if (!file.exists(input_object_path)) {
  stop("Input object not found at: ", input_object_path)
}

obj <- readRDS(input_object_path)

# Standardize donor metadata to Donor
if (!"Donor" %in% colnames(obj@meta.data)) {
  if ("donor_global_pretty" %in% colnames(obj@meta.data)) {
    obj@meta.data[["Donor"]] <- obj@meta.data[["donor_global_pretty"]]
  } else {
    stop("No donor column found. Expected 'Donor' or 'donor_global_pretty'.")
  }
}

# Standardize condition labels
if ("condition" %in% colnames(obj@meta.data)) {
  obj@meta.data[["condition"]] <- as.character(obj@meta.data[["condition"]])
  obj@meta.data[["condition"]][obj@meta.data[["condition"]] %in% c("CTRL", "Control_no_hit")] <- "Ctrl"
  obj@meta.data[["condition"]][obj@meta.data[["condition"]] %in% c("Single_hit_D7", "Stim")] <- "Stim"
  obj@meta.data[["condition"]][obj@meta.data[["condition"]] %in% c("Single_hit_D1", "Prime")] <- "Prime"
}

# ------------------------------------------------------------------------------
# 1) Helper functions
# ------------------------------------------------------------------------------
make_pseudobulk <- function(obj, assay_use, features_use, slot_use = "data") {
  mat <- GetAssayData(object = obj, assay = assay_use, slot = slot_use)
  features_use <- intersect(features_use, rownames(mat))
  mat <- mat[features_use, , drop = FALSE]
  meta <- obj@meta.data
  
  group_id <- paste(meta[["Donor"]], meta[["condition"]], sep = "_")
  pb <- sapply(
    unique(group_id),
    function(g) {
      cells_use <- rownames(meta)[group_id == g]
      Matrix::rowMeans(mat[, cells_use, drop = FALSE])
    }
  )
  
  pb <- as.data.frame(pb)
  return(pb)
}

average_condition <- function(pb_mat) {
  condition_names <- sub(".*_", "", colnames(pb_mat))
  avg <- sapply(
    unique(condition_names),
    function(cond) {
      rowMeans(pb_mat[, condition_names == cond, drop = FALSE])
    }
  )
  avg <- as.data.frame(avg)
  return(avg)
}

calculate_donor_correlations <- function(pb_mat, assay_name) {
  sample_names <- colnames(pb_mat)
  donor_names <- sub("_(Ctrl|Prime|Stim)$", "", sample_names)
  condition_names <- sub(".*_", "", sample_names)
  donors <- unique(donor_names)
  out <- data.frame()
  
  for (d in donors) {
    donor_cols <- sample_names[donor_names == d]
    donor_conditions <- condition_names[donor_names == d]
    
    if (!all(c("Ctrl", "Prime", "Stim") %in% donor_conditions)) next
    
    ctrl_col  <- donor_cols[donor_conditions == "Ctrl"]
    prime_col <- donor_cols[donor_conditions == "Prime"]
    stim_col  <- donor_cols[donor_conditions == "Stim"]
    
    r_ctrl_prime <- cor(log1p(pb_mat[, ctrl_col]), log1p(pb_mat[, prime_col]), method = "pearson")
    r_ctrl_stim  <- cor(log1p(pb_mat[, ctrl_col]), log1p(pb_mat[, stim_col]), method = "pearson")
    
    out <- rbind(
      out,
      data.frame(
        Donor = d,
        Assay = assay_name,
        Ctrl_vs_Prime = r_ctrl_prime,
        Ctrl_vs_Stim = r_ctrl_stim
      )
    )
  }
  out
}

# ------------------------------------------------------------------------------
# 2) RNA pseudobulk: top 3000 HVGs
# ------------------------------------------------------------------------------
DefaultAssay(obj) <- "RNA"
obj <- FindVariableFeatures(obj, assay = "RNA", selection.method = "vst", nfeatures = 3000)
rna_hvg <- VariableFeatures(obj[["RNA"]])

pb_rna <- make_pseudobulk(obj = obj, assay_use = "RNA", features_use = rna_hvg, slot_use = "data")
avg_rna <- average_condition(pb_rna)

# ------------------------------------------------------------------------------
# 3) ATAC pseudobulk: variable peaks
# ------------------------------------------------------------------------------
DefaultAssay(obj) <- "ATAC"
atac_peaks <- VariableFeatures(obj[["ATAC"]])

pb_atac <- make_pseudobulk(obj = obj, assay_use = "ATAC", features_use = atac_peaks, slot_use = "data")
avg_atac <- average_condition(pb_atac)

# ------------------------------------------------------------------------------
# 4) Create plotting data
# ------------------------------------------------------------------------------
rna_prime_df <- data.frame(Ctrl = log1p(avg_rna$Ctrl), Prime = log1p(avg_rna$Prime))
rna_stim_df  <- data.frame(Ctrl = log1p(avg_rna$Ctrl), Stim = log1p(avg_rna$Stim))
atac_prime_df <- data.frame(Ctrl = log1p(avg_atac$Ctrl), Prime = log1p(avg_atac$Prime))
atac_stim_df  <- data.frame(Ctrl = log1p(avg_atac$Ctrl), Stim = log1p(avg_atac$Stim))

# ------------------------------------------------------------------------------
# 5) Correlations
# ------------------------------------------------------------------------------
r_rna_prime <- cor(rna_prime_df$Ctrl, rna_prime_df$Prime)
r_rna_stim  <- cor(rna_stim_df$Ctrl, rna_stim_df$Stim)
r_atac_prime <- cor(atac_prime_df$Ctrl, atac_prime_df$Prime)
r_atac_stim  <- cor(atac_stim_df$Ctrl, atac_stim_df$Stim)

cor_df <- data.frame(
  Assay = c("RNA", "RNA", "ATAC", "ATAC"),
  Comparison = c("Ctrl vs Prime", "Ctrl vs Stim", "Ctrl vs Prime", "Ctrl vs Stim"),
  Correlation = c(r_rna_prime, r_rna_stim, r_atac_prime, r_atac_stim)
)

print(cor_df)

# ------------------------------------------------------------------------------
# 6) Plotting theme
# ------------------------------------------------------------------------------
my_theme <- theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, face = "bold", color = "black"),
    axis.text.y = element_text(size = 20, face = "bold", color = "black"),
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    aspect.ratio = 1
  )

# ------------------------------------------------------------------------------
# 7) Create and save the 2x2 figure
# ------------------------------------------------------------------------------
p_rna_prime <- ggplot(rna_prime_df, aes(Ctrl, Prime)) +
  geom_point(size = 0.5, alpha = 0.5, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  coord_equal() +
  labs(title = paste0("RNA (r = ", round(r_rna_prime, 3), ")"), x = "Ctrl", y = "Prime") +
  my_theme

p_rna_stim <- ggplot(rna_stim_df, aes(Ctrl, Stim)) +
  geom_point(size = 0.5, alpha = 0.5, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  coord_equal() +
  labs(title = paste0("RNA (r = ", round(r_rna_stim, 3), ")"), x = "Ctrl", y = "Stim") +
  my_theme

p_atac_prime <- ggplot(atac_prime_df, aes(Ctrl, Prime)) +
  geom_point(size = 0.3, alpha = 0.3, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  coord_equal() +
  labs(title = paste0("ATAC (r = ", round(r_atac_prime, 3), ")"), x = "Ctrl", y = "Prime") +
  my_theme

p_atac_stim <- ggplot(atac_stim_df, aes(Ctrl, Stim)) +
  geom_point(size = 0.3, alpha = 0.3, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  coord_equal() +
  labs(title = paste0("ATAC (r = ", round(r_atac_stim, 3), ")"), x = "Ctrl", y = "Stim") +
  my_theme

final_plot <- (p_rna_prime | p_rna_stim) / (p_atac_prime | p_atac_stim) + plot_layout(widths = c(1, 1), heights = c(1, 1))
print(final_plot)

ggsave(output_file, final_plot, width = 8, height = 6)

# ------------------------------------------------------------------------------
# 8) Donor-level correlation analysis
# ------------------------------------------------------------------------------
rna_donor_corr <- calculate_donor_correlations(pb_rna, assay_name = "RNA")
atac_donor_corr <- calculate_donor_correlations(pb_atac, assay_name = "ATAC")
donor_corr <- rbind(rna_donor_corr, atac_donor_corr)

print(donor_corr)
