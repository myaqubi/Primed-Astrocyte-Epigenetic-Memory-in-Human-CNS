# 03_post_merge_QC.R
# Post-merge QC and per-sample MAD-based outlier filtering.
#
#### Quality Control
DefaultAssay(combined) <- "RNA"
combined <- JoinLayers(object    = combined)
combined <- PercentageFeatureSet(combined, pattern = "^MT-", col.name = "percent.mt", assay = "RNA")
DefaultAssay(combined) <- "ATAC"
combined <- NucleosomeSignal(combined, assay = "ATAC")
combined <- TSSEnrichment(combined, assay = "ATAC")
VlnPlot(combined,features = c("nFeature_RNA", "percent.mt", "nFeature_ATAC", "TSS.enrichment","nucleosome_signal"),
        ncol = 5, pt.size = 0, group.by = "orig.ident")

#### Temp save
saveRDS(combined, file = "combined_object_temp_save.rds")


# At this stage, define QC parameters

library(scater)
library(SingleCellExperiment)
library(Seurat)

# Load the merged object if this script is run in a fresh R session.
if (!exists("combined")) {
  combined <- readRDS("combined_object_temp_save.rds")
}

# 1) Convert Seurat → SingleCellExperiment
sce <- as.SingleCellExperiment(combined, assay = "RNA")

# 2) Flag outliers **per sample** for each QC metric: in the following order:
# a) nFeature_RNA, b) nCount_RNA, c) percent.mt, d) nCount_ATAC, e) nFeature_ATAC f) TSS.enrichment g) nucleosome_signal

feature_flag <- isOutlier(sce$nFeature_RNA, nmads = 5, type  = "both", batch = sce$orig.ident)

umi_flag <- isOutlier(sce$nCount_RNA,nmads = 5, type  = "both", batch = sce$orig.ident)

mito_flag <- isOutlier(sce$percent.mt, nmads = 5, type  = "higher", batch = sce$orig.ident)

count_ATAC_flag <- isOutlier(sce$nCount_ATAC, nmads = 5, type  = "both", batch = sce$orig.ident)

feature_ATAC_flag <- isOutlier(sce$nFeature_ATAC, nmads = 5, type  = "both", batch = sce$orig.ident)

TSS_flag <- isOutlier(sce$TSS.enrichment, nmads = 5, type  = "higher", batch = sce$orig.ident)

nucleosome_flag <- isOutlier(sce$nucleosome_signal, nmads = 5, type  = "higher", batch = sce$orig.ident)

# 3) Combine & subset
is_outlier    <- feature_flag | umi_flag | mito_flag | count_ATAC_flag | feature_ATAC_flag | TSS_flag | nucleosome_flag
sce_filtered  <- sce[, !is_outlier]

# 4) Back to Seurat
seurat_filtered <- as.Seurat(sce_filtered, counts="counts")

# 5) (Optional) Summary
message("Before: ", ncol(sce), " cells;\n",
        "After:  ", ncol(seurat_filtered), " cells;\n",
        "Dropped: ", sum(is_outlier), " outliers.")


#### save merged combined object in which I have done QC steps but not integration
saveRDS(seurat_filtered, file = "combined_object_seurat_QC_filtered.rds")
