# ============================================================================== 
# Identification of Stim-upregulated genes and GO enrichment
# ============================================================================== 
#
# Purpose:
#   This script identifies genes significantly upregulated in Stim cells relative to the 
#   combined Ctrl and Prime populations using the MAST differential- expression framework. 
#   Donor identity is included as a latent covariate to adjust for donor-associated variation 
#   in gene expression and reduce the likelihood that differences between donors are interpreted 
#   as condition- specific effects. The resulting significant genes are then tested for Gene Ontology Biological Process enrichment.
#
# Notes:
#   - Replace the generic placeholder paths below with your own project paths.
# ============================================================================== 

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# ------------------------------------------------------------------------------
# 0) Configuration
# ------------------------------------------------------------------------------
project_dir <- "/path/to/your/project"
output_dir <- file.path(project_dir, "DEG_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_object_path <- file.path(project_dir, "In_vitro_prenatal_samples_QC_filtered_integrated_annotated_astrocytes.rds")
output_file <- file.path(output_dir, "Upreg_genes_stim_vs_CTRL_PRIME_all_cells_sig.csv")

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
# 1) Differential expression: Stim vs Ctrl + Prime
# ------------------------------------------------------------------------------
Idents(obj) <- "condition"
DefaultAssay(obj) <- "RNA"

DEG_stim_vs_CTRL_PRIME <- FindMarkers(
  object = obj,
  ident.1 = "Stim",
  ident.2 = c("Ctrl", "Prime"),
  assay = "RNA",
  slot = "data",
  test.use = "MAST",
  latent.vars = "Donor",
  only.pos = TRUE,
  logfc.threshold = 0.1,
  min.pct = 0.1
)

# ------------------------------------------------------------------------------
# 2) Format results
# ------------------------------------------------------------------------------
DEG_stim_vs_CTRL_PRIME <- DEG_stim_vs_CTRL_PRIME %>%
  rownames_to_column("gene") %>%
  arrange(p_val_adj, desc(avg_log2FC))

DEG_stim_vs_CTRL_PRIME_sig <- DEG_stim_vs_CTRL_PRIME %>%
  filter(p_val_adj < 0.05)

cat("Significant Stim vs Ctrl+Prime genes:", nrow(DEG_stim_vs_CTRL_PRIME_sig), "\n")

write.csv(DEG_stim_vs_CTRL_PRIME_sig, file = output_file, row.names = FALSE)

# ------------------------------------------------------------------------------
# 3) Gene list for GO analysis
# ------------------------------------------------------------------------------
gene_symbols <- unique(DEG_stim_vs_CTRL_PRIME_sig$gene)
message("Number of genes for GO analysis: ", length(gene_symbols))

# ------------------------------------------------------------------------------
# 4) Convert gene symbols to Entrez IDs
# ------------------------------------------------------------------------------
gene_conversion <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

head(gene_conversion)

# ------------------------------------------------------------------------------
# 5) GO Biological Process enrichment
# ------------------------------------------------------------------------------
go_bp <- enrichGO(
  gene = gene_conversion$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

go_bp_results <- as.data.frame(go_bp)

head(go_bp_results)
message("Number of GO BP terms: ", nrow(go_bp_results))

# ------------------------------------------------------------------------------
# 6) Basic plots
# ------------------------------------------------------------------------------
barplot(go_bp, showCategory = 10)
