# 04_RNA_integration.R
# RNA preprocessing and RNA-based CSS integration.


library(Seurat)
library(simspec)

if (!exists("seurat_filtered")) {
  seurat_filtered <- readRDS("combined_object_seurat_QC_filtered.rds")
}

# RNA preprocessing and RNA-based integration.
# Normalization, variable-feature selection, scaling, PCA, and RNA UMAP are
# performed before cluster similarity spectrum (CSS) integration.

DefaultAssay(seurat_filtered) <- "RNA"

seurat_filtered <- NormalizeData(seurat_filtered) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  CellCycleScoring(s.features = cc.genes.updated.2019$s.genes,
                   g2m.features = cc.genes.updated.2019$g2m.genes) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20, reduction.name = "umap_rna", reduction.key = "UMAPRNA_")
DimPlot(seurat_filtered, group.by = "orig.ident", reduction = "umap_rna")
FeaturePlot(seurat_filtered, c("PLP1"), reduction = "umap_rna", min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(seurat_filtered, c("C1QA"), reduction = "umap_rna", min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(seurat_filtered, c("GFAP"), reduction = "umap_rna", min.cutoff = "q10", max.cutoff = "q90")


# RNA integration using cluster similarity spectrum (CSS)
library(simspec)
seurat_filtered <- cluster_sim_spectrum(seurat_filtered,
                                        label_tag = "orig.ident",
                                        cluster_resolution = 0.6,
                                        reduction.name = "css_rna",
                                        reduction.key = "CSSRNA_")
seurat_filtered <- RunUMAP(seurat_filtered,
                           reduction = "css_rna",
                           dims = 1:ncol(Embeddings(seurat_filtered,"css_rna")),
                           reduction.name = "umap_css_rna",
                           reduction.key = "UMAPCSSRNA_")

# Identify and drop cells with NA in CSS reduction
css_matrix <- Embeddings(seurat_filtered, "css_rna")
valid_cells <- rownames(css_matrix)[complete.cases(css_matrix)]
seurat_filtered <- subset(seurat_filtered, cells = valid_cells)

# Re-run UMAP
seurat_filtered <- RunUMAP(
  seurat_filtered,
  reduction = "css_rna",
  dims = 1:ncol(Embeddings(seurat_filtered, "css_rna")),
  reduction.name = "umap_css_rna",
  reduction.key = "UMAPCSSRNA_")
DimPlot(seurat_filtered, group.by = "orig.ident", reduction = "umap_css_rna")

# Save the RNA-integrated object for handoff to the ATAC integration script.
saveRDS(seurat_filtered, file = "seurat_filtered_RNA_integrated.rds")
