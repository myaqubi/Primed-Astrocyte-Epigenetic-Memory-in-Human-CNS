# 05_ATAC_integration.R
# ATAC preprocessing, LSI, and ATAC-based CSS integration.

library(Seurat)
library(Signac)
library(simspec)

if (!exists("seurat_filtered")) {
  seurat_filtered <- readRDS("seurat_filtered_RNA_integrated.rds")
}

# ATAC preprocessing and initial LSI embedding
DefaultAssay(seurat_filtered) <- "ATAC"
seurat_filtered <- FindTopFeatures(seurat_filtered, min.cutoff = 50)
seurat_filtered <- RunTFIDF(seurat_filtered, method = 1)
seurat_filtered <- RunSVD(seurat_filtered, n = 50)
seurat_filtered <- RunUMAP(seurat_filtered,
                           reduction = "lsi",
                           dims = 2:30,
                           reduction.name = "umap_atac",
                           reduction.key = "UMAPATAC_")
saveRDS(seurat_filtered, file = "seurat_filtered_integrated.rds")

# ATAC integration using cluster similarity spectrum (CSS)
seurat_filtered <- cluster_sim_spectrum(
  seurat_filtered,
  label_tag        = "orig.ident",
  use_dr           = "lsi",
  dims_use         = 2:30,
  cluster_resolution = 0.6,
  reduction.name   = "css_atac",
  reduction.key    = "CSSATAC_"
)

seurat_filtered <- RunUMAP(
  seurat_filtered,
  reduction      = "css_atac",
  dims           = 1:ncol(Embeddings(seurat_filtered, "css_atac")),
  reduction.name = "umap_css_atac",
  reduction.key  = "UMAPCSSATAC_"

# Save the ATAC-integrated object for handoff to WNN integration.
saveRDS(seurat_filtered, file = "seurat_filtered_ATAC_integrated.rds")
