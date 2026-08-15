# ============================================================================== 
# Quality control filtering and doublet removal
# ============================================================================== 
#
# Purpose:
#   quality-control filtering of the merged multiome object and doublet removal
#   using scDblFinder.
#
# Notes:
#   - The workflow is written to be portable and reusable.
#   - Replace the generic placeholder path below with your own output directory.
#   - This script assumes that the merged object from the previous merge step is
#     available in the workspace as `combined`.
# ============================================================================== 

# ------------------------------------------------------------------------------
# 0) Packages
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(scater)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(dplyr)
  library(Matrix)
})

# ------------------------------------------------------------------------------
# 1) Project settings
# ------------------------------------------------------------------------------
output_dir <- "/path/to/your/output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("combined") || !inherits(combined, "Seurat")) {
  stop("The merged Seurat object 'combined' was not found. Please run the merge step first.")
}

# ------------------------------------------------------------------------------
# 2) Sanity check required metadata columns
# ------------------------------------------------------------------------------
required_cols <- c(
  "nFeature_RNA", "nCount_RNA", "percent.mt",
  "nCount_ATAC", "nFeature_ATAC",
  "TSS.enrichment", "nucleosome_signal", "reads_in_peaks_frac",
  "orig.ident", "droplet_qc"
)

missing_cols <- setdiff(required_cols, colnames(combined@meta.data))
if (length(missing_cols) > 0) {
  stop("These required metadata columns are missing from 'combined': ",
       paste(missing_cols, collapse = ", "))
}

# ------------------------------------------------------------------------------
# 3) Define batch for sample-aware MAD filtering
# ------------------------------------------------------------------------------
batch <- combined$orig.ident

# ------------------------------------------------------------------------------
# 4) MAD-based outlier flags (within-sample) for count/feature/mito metrics
# ------------------------------------------------------------------------------
feature_flag <- scater::isOutlier(
  combined$nFeature_RNA,
  nmads = 5,
  type  = "both",
  batch = batch
)

umi_flag <- scater::isOutlier(
  combined$nCount_RNA,
  nmads = 5,
  type  = "both",
  batch = batch
)

mito_flag <- scater::isOutlier(
  combined$percent.mt,
  nmads = 5,
  type  = "higher",
  batch = batch
)

count_ATAC_flag <- scater::isOutlier(
  combined$nCount_ATAC,
  nmads = 5,
  type  = "both",
  batch = batch
)

feature_ATAC_flag <- scater::isOutlier(
  combined$nFeature_ATAC,
  nmads = 5,
  type  = "both",
  batch = batch
)

# ------------------------------------------------------------------------------
# 5) Fixed-threshold ATAC QC flags (global thresholds)
# ------------------------------------------------------------------------------
tss_low_flag <- !is.na(combined$TSS.enrichment) & (combined$TSS.enrichment <= 3)
nucleosome_high_flag <- !is.na(combined$nucleosome_signal) & (combined$nucleosome_signal >= 2)
frip_low_flag <- !is.na(combined$reads_in_peaks_frac) & (combined$reads_in_peaks_frac <= 0.3)

tss_na_flag <- is.na(combined$TSS.enrichment)
nucleosome_na_flag <- is.na(combined$nucleosome_signal)
frip_na_flag <- is.na(combined$reads_in_peaks_frac)

# ------------------------------------------------------------------------------
# 6) Combine all QC flags
# ------------------------------------------------------------------------------
is_outlier <- feature_flag |
  umi_flag |
  mito_flag |
  count_ATAC_flag |
  feature_ATAC_flag |
  tss_low_flag |
  nucleosome_high_flag |
  frip_low_flag |
  tss_na_flag |
  nucleosome_na_flag |
  frip_na_flag

message("Total cells in merged object: ", length(is_outlier))
message("QC outliers flagged (MAD + fixed thresholds): ", sum(is_outlier))

message("\nBreakdown of QC flags:")
print(c(
  nFeature_RNA_MAD5_both  = sum(feature_flag, na.rm = TRUE),
  nCount_RNA_MAD5_both    = sum(umi_flag, na.rm = TRUE),
  percent_mt_MAD5_high    = sum(mito_flag, na.rm = TRUE),
  nCount_ATAC_MAD5_both   = sum(count_ATAC_flag, na.rm = TRUE),
  nFeature_ATAC_MAD5_both = sum(feature_ATAC_flag, na.rm = TRUE),
  TSS_le_3                = sum(tss_low_flag, na.rm = TRUE),
  nucleosome_ge_2         = sum(nucleosome_high_flag, na.rm = TRUE),
  FRiP_le_0.3             = sum(frip_low_flag, na.rm = TRUE),
  TSS_NA                  = sum(tss_na_flag, na.rm = TRUE),
  nucleosome_NA           = sum(nucleosome_na_flag, na.rm = TRUE),
  FRiP_NA                 = sum(frip_na_flag, na.rm = TRUE)
))

# ------------------------------------------------------------------------------
# 7) Retention summary by sample before filtering
# ------------------------------------------------------------------------------
qc_df <- data.frame(
  sample = combined$orig.ident,
  is_outlier = is_outlier,
  stringsAsFactors = FALSE
)

qc_summary_by_sample <- qc_df %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    total_cells = dplyr::n(),
    lost_cells = sum(is_outlier, na.rm = TRUE),
    kept_cells = sum(!is_outlier, na.rm = TRUE),
    lost_pct = round(100 * lost_cells / total_cells, 2),
    kept_pct = round(100 * kept_cells / total_cells, 2),
    .groups = "drop"
  )

message("\nQC retention summary by sample:")
print(qc_summary_by_sample)

# ------------------------------------------------------------------------------
# 8) Subset non-outliers
# ------------------------------------------------------------------------------
combined_filtered <- subset(
  combined,
  cells = Cells(combined)[!is_outlier]
)

message("Cells after QC filtering: ", ncol(combined_filtered))

# ------------------------------------------------------------------------------
# 9) Apply DropletQC classification: keep only true cells
# ------------------------------------------------------------------------------
cleaned_object <- subset(
  combined_filtered,
  subset = droplet_qc == "cell"
)

message("Cells after DropletQC filter (droplet_qc == 'cell'): ",
        ncol(cleaned_object))

# ------------------------------------------------------------------------------
# 10) Save the QC-filtered object
# ------------------------------------------------------------------------------
saveRDS(
  cleaned_object,
  file = file.path(output_dir, "cleaned_object.rds")
)

# ============================================================================== 
# PART E — Doublet detection with scDblFinder, then keep singlets
# ============================================================================== 

# ------------------------------------------------------------------------------
# E1) Convert Seurat -> SingleCellExperiment (RNA counts used by default)
# ------------------------------------------------------------------------------
sce <- as.SingleCellExperiment(cleaned_object)

# ------------------------------------------------------------------------------
# E2) Remove zero-library cells before scDblFinder
# ------------------------------------------------------------------------------
libsizes <- colSums(counts(sce))
message("Cells with zero RNA library size: ", sum(libsizes == 0))

sce <- sce[, libsizes > 0]
message("Cells kept for scDblFinder: ", ncol(sce))

# ------------------------------------------------------------------------------
# E3) Run scDblFinder (sample-aware via orig.ident)
# ------------------------------------------------------------------------------
sce <- scDblFinder(
  sce,
  samples = sce$orig.ident,
  verbose = TRUE
)

# ------------------------------------------------------------------------------
# E4) Extract doublet results and map them back to the Seurat object
# ------------------------------------------------------------------------------
dbl_meta <- as.data.frame(colData(sce)[, c("scDblFinder.class", "scDblFinder.score")])
dbl_meta$cell <- rownames(dbl_meta)

cells <- Cells(cleaned_object)

cleaned_object$doublet_class <- NA_character_
cleaned_object$doublet_score <- NA_real_

common_cells <- intersect(cells, dbl_meta$cell)

cleaned_object$doublet_class[common_cells] <-
  dbl_meta[match(common_cells, dbl_meta$cell), "scDblFinder.class"]

cleaned_object$doublet_score[common_cells] <-
  dbl_meta[match(common_cells, dbl_meta$cell), "scDblFinder.score"]

message("\nRaw mapped scDblFinder classes:")
print(table(cleaned_object$doublet_class, useNA = "ifany"))

# ------------------------------------------------------------------------------
# E5) Relabel classes if they appear as "1"/"2" instead of text
# ------------------------------------------------------------------------------
dc <- as.character(cleaned_object$doublet_class)
dc[dc == "1"] <- "singlet"
dc[dc == "2"] <- "doublet"

cleaned_object$doublet_class <- factor(
  dc,
  levels = c("singlet", "doublet")
)

message("\nRelabeled scDblFinder classes:")
print(table(cleaned_object$doublet_class, useNA = "ifany"))

# ------------------------------------------------------------------------------
# E6) Keep singlets only
# ------------------------------------------------------------------------------
cleaned_object <- subset(
  cleaned_object,
  subset = doublet_class == "singlet"
)

message("Cells after scDblFinder singlet filter: ", ncol(cleaned_object))

# ------------------------------------------------------------------------------
# E7) Save post-doublet-filter object
# ------------------------------------------------------------------------------
saveRDS(
  cleaned_object,
  file = file.path(output_dir, "cleaned_object_noDoublets.rds")
)

message("\nQC filtering and doublet removal completed successfully.")
