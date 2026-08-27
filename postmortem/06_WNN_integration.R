# 06_WNN_integration.R
# Weighted nearest-neighbor (WNN) integration of RNA and ATAC CSS reductions.

library(Seurat)

if (!exists("seurat_filtered")) {
  seurat_filtered <- readRDS("seurat_filtered_ATAC_integrated.rds")
}

# Weighted nearest-neighbor integration of RNA and ATAC modalities
seurat_filtered <- FindMultiModalNeighbors(seurat_filtered,
                                           reduction.list = list("css_rna", "css_atac"),
                                           dims.list = list(1:ncol(Embeddings(seurat_filtered,"css_rna")),
                                                            1:ncol(Embeddings(seurat_filtered,"css_atac"))),
                                           modality.weight.name = c("RNA.weight","ATAC.weight"),
                                           verbose = TRUE)

seurat_filtered <- RunUMAP(seurat_filtered, nn.name = "weighted.nn", assay = "RNA")
seurat_filtered <- FindClusters(seurat_filtered, graph.name = "wsnn", resolution = 0.5)
DimPlot(seurat_filtered, label = TRUE)
saveRDS(seurat_filtered, file = "seurat_filtered_integrated_final.rds")
