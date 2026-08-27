# 02_merge_consensus_peaks.R
# Merge preprocessed samples, construct a unified ATAC peak set,
# re-quantify fragments over consensus peaks, retain standard chromosomes,
# and calculate post-merge QC metrics.
#
# Run after 01_sample_preprocessing_QC.R in the same R session.

library(Seurat)
library(Signac)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

# Merge preprocessed sample objects
 ─────────────────────────────────
combined <- merge(
  x = sample1_CTRL,
  y = list(
    sample2_MS,
    sample3_CTRL,
    sample4_MS,
    sample7_CTRL,
    sample8_MS,
    sample9_MS,
    sample11_CTRL,
    sample12_MS
  ),
  add.cell.ids = names(samples),    # length 9: Sample1…Sample12
  project      = "Multiome_Combined"
)


# Check
combined  
# Should show all cells, both RNA & ATAC assays, and an “orig.ident” metadata column

peaks <- reduce(
  unlist(
    as(
      c(
        sample1_CTRL@assays$ATAC@ranges,
        sample2_MS@assays$ATAC@ranges,
        sample3_CTRL@assays$ATAC@ranges,
        sample4_MS@assays$ATAC@ranges,
        sample7_CTRL@assays$ATAC@ranges,
        sample8_MS@assays$ATAC@ranges,
        sample9_MS@assays$ATAC@ranges,
        sample11_CTRL@assays$ATAC@ranges,
        sample12_MS@assays$ATAC@ranges
      ),
      "GRangesList"
    )
  )
)


peakwidths <- width(peaks)
peaks <- peaks[peakwidths > 20 & peakwidths < 10000]




# 1) Quantify all fragments in `combined` over your consensus peaks
counts_atac_merged <- FeatureMatrix(
  fragments = combined@assays$ATAC@fragments,
  features  = peaks,
  cells     = Cells(combined)
)

# 2) Replace the ATAC assay in `combined` with this merged‐peak assay
combined[["ATAC"]] <- CreateChromatinAssay(
  counts     = counts_atac_merged,
  fragments  = combined@assays$ATAC@fragments,
  annotation = combined@assays$ATAC@annotation,
  sep        = c(":", "-"),
  genome     = "hg38"
)



library(BSgenome.Hsapiens.UCSC.hg38)

# 1) Grab the “standard” chromosomes for hg38
standard_chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)

# 2) Find which peaks lie on those chromosomes
#    (Here, combined@assays$ATAC@ranges is a GRanges of all peaks.)
keep_idx <- which(
  as.character(seqnames(granges(combined[["ATAC"]]))) %in% standard_chroms
)

# 3) Subset the ATAC assay to only those “standard‐chr” features
combined[["ATAC"]] <- subset(
  combined[["ATAC"]],
  features = rownames(combined[["ATAC"]])[keep_idx]
)

# 4) Finally, drop any unused seqlevels from the GRanges object
seqlevels(combined[["ATAC"]]@ranges) <- intersect(
  seqlevels(granges(combined[["ATAC"]])),
  unique(seqnames(granges(combined[["ATAC"]])))
