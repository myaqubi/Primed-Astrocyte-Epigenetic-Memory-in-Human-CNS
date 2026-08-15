# ============================================================================== 
# Batch correction and multimodal integration
# ============================================================================== 
#
# Purpose:
#   This script contains the batch-correction and multimodal integration 
#   It starts from the cleaned object after QC and
#   doublet removal, performs RNA/ATAC batch correction using CSS and WNN, and
#   calculates the GeneActivity assay.
# ============================================================================== 

# ------------------------------------------------------------------------------
# 0) Packages
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(simspec)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 1) Input object
# ------------------------------------------------------------------------------
# Replace this with the path to your cleaned object if needed.
# The expected input is the object produced after QC filtering and doublet removal.
if (!exists("In_vitro_combined_object_filtered") || !inherits(In_vitro_combined_object_filtered, "Seurat")) {
  stop("The object 'In_vitro_combined_object_filtered' was not found. Please load your cleaned object first.")
}

# Use the cleaned object as the starting point.
obj <- In_vitro_combined_object_filtered

# ------------------------------------------------------------------------------
# 2) RNA preprocessing before batch correction
# ------------------------------------------------------------------------------
DefaultAssay(obj) <- "RNA"

obj <- obj |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 3000) |>
  ScaleData() |>
  RunPCA(npcs = 50)

obj <- RunUMAP(
  obj,
  reduction      = "pca",
  dims           = 1:20,
  reduction.name = "umap_rna_pre",
  reduction.key  = "UMAPRNApre_"
)

p_rna_pre <- DimPlot(
  obj,
  reduction = "umap_rna_pre",
  group.by  = "orig.ident"
) + ggtitle("RNA • BEFORE (PCA→UMAP) ")

# ------------------------------------------------------------------------------
# 3) RNA batch correction using CSS
# ------------------------------------------------------------------------------
obj <- cluster_sim_spectrum(
  obj,
  label_tag          = "orig.ident",
  use_dr             = "pca",
  dims_use           = 1:30,
  cluster_resolution = 0.4,
  reduction.name     = "css_rna",
  reduction.key      = "CSSRNA_"
)

obj <- RunUMAP(
  obj,
  reduction      = "css_rna",
  dims           = 1:ncol(Embeddings(obj, "css_rna")),
  reduction.name = "umap_rna_post",
  reduction.key  = "UMAPRNApost_"
)

p_rna_post <- DimPlot(
  obj,
  reduction = "umap_rna_post",
  group.by  = "orig.ident"
) + ggtitle("RNA • AFTER (CSS→UMAP) ")

# Optional side-by-side view
p_rna_pre + p_rna_post

# ------------------------------------------------------------------------------
# 4) ATAC preprocessing before batch correction
# ------------------------------------------------------------------------------
DefaultAssay(obj) <- "ATAC"

obj <- obj |>
  FindTopFeatures(min.cutoff = 50) |>
  RunTFIDF(method = 1) |>
  RunSVD(n = 50)

obj <- RunUMAP(
  obj,
  reduction      = "lsi",
  dims           = 2:30,
  reduction.name = "umap_atac_pre",
  reduction.key  = "UMAPATACpre_"
)

p_atac_pre <- DimPlot(
  obj,
  reduction = "umap_atac_pre",
  group.by  = "orig.ident"
) + ggtitle("ATAC • BEFORE (LSI→UMAP) ")

# ------------------------------------------------------------------------------
# 5) ATAC batch correction using CSS
# ------------------------------------------------------------------------------
obj <- cluster_sim_spectrum(
  obj,
  label_tag          = "orig.ident",
  use_dr             = "lsi",
  dims_use           = 2:30,
  cluster_resolution = 0.4,
  reduction.name     = "css_atac",
  reduction.key      = "CSSATAC_"
)

obj <- RunUMAP(
  obj,
  reduction      = "css_atac",
  dims           = 1:ncol(Embeddings(obj, "css_atac")),
  reduction.name = "umap_atac_post",
  reduction.key  = "UMAPATACpost_"
)

p_atac_post <- DimPlot(
  obj,
  reduction = "umap_atac_post",
  group.by  = "orig.ident"
) + ggtitle("ATAC • AFTER (CSS→UMAP) ")

# Optional side-by-side view
p_atac_pre + p_atac_post

# ------------------------------------------------------------------------------
# 6) Multimodal integration using WNN over CSS reductions
# ------------------------------------------------------------------------------
stopifnot(all(c("css_rna", "css_atac") %in% Reductions(obj)))

rna_dims  <- seq_len(ncol(Embeddings(obj, "css_rna")))
atac_dims <- seq_len(ncol(Embeddings(obj, "css_atac")))

obj <- FindMultiModalNeighbors(
  object               = obj,
  reduction.list       = list("css_rna", "css_atac"),
  dims.list            = list(rna_dims, atac_dims),
  modality.weight.name = c("RNA.weight", "ATAC.weight"),
  verbose              = TRUE
)

obj <- RunUMAP(
  obj,
  nn.name        = "weighted.nn",
  reduction.name = "umap_css_rna_atac_wnn",
  reduction.key  = "UMAPCSSWNN_"
)

obj <- FindClusters(
  obj,
  graph.name = "wsnn",
  resolution = 0.4
)

p_wnn <- DimPlot(
  obj,
  reduction = "umap_css_rna_atac_wnn",
  group.by  = "orig.ident"
) + ggtitle("WNN (CSS-RNA + CSS-ATAC) ")

print(p_wnn)

# ------------------------------------------------------------------------------
# 7) Calculate the GeneActivity assay and add it to the integrated object
# ------------------------------------------------------------------------------
DefaultAssay(obj) <- "ATAC"
gene.activities <- GeneActivity(obj)

obj[["GeneActivity"]] <- CreateAssayObject(counts = gene.activities)

obj <- NormalizeData(
  object = obj,
  assay = "GeneActivity",
  normalization.method = "LogNormalize"
)

# ------------------------------------------------------------------------------
# 8) Save the integrated object
# ------------------------------------------------------------------------------
cleaned_object_noDoublets_integrated <- obj

saveRDS(
  cleaned_object_noDoublets_integrated,
  file = "cleaned_object_noDoublets_integrated.rds"
)

message("Integration completed and object saved as 'cleaned_object_noDoublets_integrated.rds'.")
