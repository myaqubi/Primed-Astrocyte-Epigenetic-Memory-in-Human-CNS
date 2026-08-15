# ============================================================================== 
# Merge samples and build a consensus ATAC peak set
# ============================================================================== 
#
# Purpose:
#   This script merges the per-sample preprocessing objects for Ctrl,
#   Stim, and Prime, rebuilds a consensus ATAC peak set, computes QC metrics,
#   and performs initial filtering before downstream integration.
#
# Notes:
#   - This script is written to be portable and easy to reuse.
#   - Replace the generic placeholder path below with your project directory.
#   - The script assumes that the per-sample preprocessing step has already
#     created an `objs` list containing Seurat objects for Ctrl, Stim, and Prime.
# ============================================================================== 

# ------------------------------------------------------------------------------
# 0) Packages
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(scater)
  library(dplyr)
  library(Matrix)
})

# ------------------------------------------------------------------------------
# 1) Project settings
# ------------------------------------------------------------------------------
project_dir <- "/path/to/your/project"
selected_samples <- c("Ctrl", "Stim", "Prime")

if (!exists("objs") || !is.list(objs)) {
  stop("The 'objs' list was not found. Please run the per-sample preprocessing script first.")
}

if (!all(selected_samples %in% names(objs))) {
  stop("The expected sample objects were not found in 'objs'.")
}

# ------------------------------------------------------------------------------
# 2) Merge selected samples
# ------------------------------------------------------------------------------
merge_name_map <- c(Ctrl = "Ctrl", Stim = "Stim", Prime = "Prime")
merge_ids <- merge_name_map[names(objs)]
merge_ids[is.na(merge_ids)] <- names(objs)[is.na(merge_ids)]

objs_for_merge <- objs
names(objs_for_merge) <- unname(merge_ids)

print(names(objs_for_merge))

if (length(objs_for_merge) < 2) {
  stop("Need at least 2 samples to merge.")
}

combined <- merge(
  x = objs_for_merge[[1]],
  y = objs_for_merge[-1],
  add.cell.ids = names(objs_for_merge),
  project = "Multiome_Combined_NoDoubleHit"
)

# Ensure sample column exists
if (!"sample" %in% colnames(combined@meta.data)) {
  combined$sample <- combined$orig.ident
}

message("\nCells by sample x donor_global_pretty:")
print(table(combined$sample, combined$donor_global_pretty, useNA = "ifany"))

# ------------------------------------------------------------------------------
# 3) Build a consensus ATAC peak set across selected samples
# ------------------------------------------------------------------------------
gr_list <- lapply(objs, function(o) o@assays$ATAC@ranges)
peaks <- reduce(unlist(as(gr_list, "GRangesList")))

# Width filter
peaks <- peaks[width(peaks) > 20 & width(peaks) < 10000]

# Keep only standard chromosomes
standard_chroms <- GenomeInfoDb::standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
peaks <- keepSeqlevels(
  peaks,
  intersect(seqlevels(peaks), standard_chroms),
  pruning.mode = "coarse"
)

message("Consensus peak count after filtering: ", length(peaks))

# ------------------------------------------------------------------------------
# 4) Re-quantify ATAC fragments over consensus peaks and replace ATAC assay
# ------------------------------------------------------------------------------
DefaultAssay(combined) <- "ATAC"

frags <- Fragments(combined[["ATAC"]])
message("Number of Fragment objects in merged ATAC assay: ", length(frags))

counts_atac_merged <- FeatureMatrix(
  fragments = frags,
  features  = peaks,
  cells     = Cells(combined)
)

combined[["ATAC"]] <- CreateChromatinAssay(
  counts     = counts_atac_merged,
  fragments  = frags,
  annotation = combined@assays$ATAC@annotation,
  sep        = c(":", "-"),
  genome     = "hg38"
)

# Optional seqlevels cleanup
seqlevels(combined[["ATAC"]]@ranges) <- intersect(
  seqlevels(granges(combined[["ATAC"]])),
  unique(as.character(seqnames(granges(combined[["ATAC"]]))))
)

# ------------------------------------------------------------------------------
# 5) Basic RNA + ATAC QC metrics (pre-filter)
# ------------------------------------------------------------------------------
DefaultAssay(combined) <- "RNA"
combined <- JoinLayers(combined)
combined <- PercentageFeatureSet(combined, pattern = "^MT-", col.name = "percent.mt", assay = "RNA")

DefaultAssay(combined) <- "ATAC"
combined <- NucleosomeSignal(combined, assay = "ATAC")
combined <- TSSEnrichment(combined, assay = "ATAC")

# Quick QC plots
VlnPlot(
  combined,
  features = c("nFeature_RNA", "percent.mt", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 5,
  pt.size = 0,
  group.by = "sample"
)

# ------------------------------------------------------------------------------
# 6) FRiP-like metric for the rebuilt consensus ATAC assay
# ------------------------------------------------------------------------------
DefaultAssay(combined) <- "ATAC"

frags <- Fragments(combined[["ATAC"]])

frag_counts_list <- lapply(frags, function(frag_obj) {
  frag_path <- frag_obj@path
  message("Counting fragments in: ", frag_path)
  x <- Signac::CountFragments(fragments = frag_path)

  bc_col <- intersect(c("CB", "cell", "barcode", "barcodes"), colnames(x))
  if (length(bc_col) == 0) {
    stop("Could not find barcode column in CountFragments output. Columns were: ", paste(colnames(x), collapse = ", "))
  }
  bc_col <- bc_col[1]

  total_col_candidates <- c("frequency_count", "reads_count", "total", "total_fragments", "count")
  total_col <- intersect(total_col_candidates, colnames(x))
  if (length(total_col) == 0) {
    stop("Could not find total fragment count column in CountFragments output. Columns were: ", paste(colnames(x), collapse = ", "))
  }
  total_col <- total_col[1]

  out <- data.frame(
    cell = x[[bc_col]],
    total_atac_fragments = as.numeric(x[[total_col]]),
    stringsAsFactors = FALSE
  )

  out
})

frag_counts_df <- dplyr::bind_rows(frag_counts_list)
frag_counts_df <- frag_counts_df %>%
  dplyr::group_by(cell) %>%
  dplyr::summarise(total_atac_fragments = sum(total_atac_fragments, na.rm = TRUE), .groups = "drop")

combined_cells <- Cells(combined)
strip_merge_prefix <- function(x) sub("^[^_]+_", "", x)
combined_cells_unpref <- strip_merge_prefix(combined_cells)

lut_total <- setNames(frag_counts_df$total_atac_fragments, frag_counts_df$cell)
total_match_direct <- unname(lut_total[combined_cells])
total_match_unpref <- unname(lut_total[combined_cells_unpref])

total_atac_fragments <- ifelse(!is.na(total_match_direct), total_match_direct, total_match_unpref)

combined$total_atac_fragments <- total_atac_fragments
combined$reads_in_peaks_frac <- ifelse(
  is.na(combined$total_atac_fragments) | combined$total_atac_fragments <= 0,
  NA_real_,
  combined$nCount_ATAC / combined$total_atac_fragments
)

message("Computed 'reads_in_peaks_frac' = nCount_ATAC / total_atac_fragments")
print(summary(combined$reads_in_peaks_frac))

# ------------------------------------------------------------------------------
# 7) Save the merged object for the next step
# ------------------------------------------------------------------------------
message("Merged object created successfully.")
message("You can now use this object for QC filtering and doublet removal.")
