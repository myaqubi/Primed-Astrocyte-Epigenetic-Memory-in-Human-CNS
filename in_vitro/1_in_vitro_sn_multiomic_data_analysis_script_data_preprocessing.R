# ============================================================================== 
# Part 1: Per-sample preprocessing for the in vitro multiome analysis
# ============================================================================== 
#
# Purpose:
#   per-sample preprocessing for the multiome RNA+ATAC data, SoupX ambient RNA
#   correction, DropletQC-based empty-droplet filtering, and Vireo donor labeling.
# ============================================================================== 

# ------------------------------------------------------------------------------
# 0) Packages
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(DropletQC)
  library(SoupX)
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(GenomicRanges)
  library(Matrix)
  library(readr)
  library(dplyr)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(scater)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(simspec)
  library(patchwork)
})

set.seed(1234)

# ==============================================================================
# PART A — Per-sample preprocessing (RNA+ATAC), SoupX, DropletQC, donor labeling
# ==============================================================================

# ------------------------------------------------------------------------------
# A1) Project root, sample metadata, and paths
#     Replace the placeholder paths with your actual local directories.
# ------------------------------------------------------------------------------
project_dir <- "/path/to/your/project"

# Names of the sample folders beneath project_dir.

sample_dir_names <- c("Ctrl", "Stim", "Prime")

# Names of Vireo output folders for each sample.

vireo_dir_names <- c("Ctrl_vireo_out", "Stim_vireo_out", "Prime_vireo_out")

sample_info <- data.frame(
  sample_id = c("Ctrl", "Stim", "Prime"),
  condition = c("Ctrl", "Stim", "Prime"),
  outs_dir  = file.path(project_dir, sample_dir_names, "multiome", "outs"),
  donor_tsv = file.path(project_dir, sample_dir_names, vireo_dir_names, "donor_ids.tsv"),
  stringsAsFactors = FALSE
)

print(sample_info[, c("sample_id", "condition")])

# ------------------------------------------------------------------------------
# A2) Annotation (load once; hg38 human)
# ------------------------------------------------------------------------------
gene_annot <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
seqlevelsStyle(gene_annot) <- "UCSC"
genome(gene_annot) <- "hg38"

# ------------------------------------------------------------------------------
# A3) Helper functions
# ------------------------------------------------------------------------------

# Detect how peak names are formatted (chr:start-end vs chr_start_end)
guess_sep <- function(peak_rownames) {
  if (any(grepl(":", peak_rownames))) return(c(":", "-"))
  if (any(grepl("_", peak_rownames))) return(c("_", "_"))
  c("-", "-")
}

# Find fragments file in common ARC locations
pick_frag <- function(outs_dir) {
  cands <- c(
    file.path(outs_dir, "atac_fragments.tsv.gz"),
    file.path(outs_dir, "fragments.tsv.gz"),
    file.path(outs_dir, "atac", "fragments.tsv.gz")
  )
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit)) stop("Fragments file not found in: ", outs_dir)
  hit
}

# Find ARC GEX BAM for DropletQC nuclear fraction
pick_bam <- function(outs_dir) {
  bam <- file.path(outs_dir, "gex_possorted_bam.bam")
  if (!file.exists(bam)) stop("Expected ARC GEX BAM not found: ", bam)
  bam
}

# Robustly attach Vireo best_singlet donor labels despite barcode formatting differences
attach_best_singlet <- function(obj, donor_tsv) {
  d <- read_tsv(donor_tsv, show_col_types = FALSE)
  stopifnot(all(c("cell", "best_singlet") %in% names(d)))

  obj_bcs <- colnames(obj)

  strip_prefix <- function(x) sub("^[^_]+_", "", x)
  add_minus1   <- function(x) ifelse(grepl("-\\d+$", x), x, paste0(x, "-1"))
  strip_minus1 <- function(x) sub("-\\d+$", "", x)

  variants_obj <- list(
    obj             = obj_bcs,
    obj_noprefix    = strip_prefix(obj_bcs),
    obj_nosuffix    = strip_minus1(obj_bcs),
    obj_nopre_nosuf = strip_minus1(strip_prefix(obj_bcs))
  )

  variants_don <- list(
    don          = d$cell,
    don_plus1    = add_minus1(d$cell),
    don_nosuffix = strip_minus1(d$cell)
  )

  best_hits <- -1
  best_ov <- NULL
  best_dv <- NULL

  for (ov in variants_obj) {
    for (dv in variants_don) {
      h <- sum(dv %in% ov)
      if (h > best_hits) {
        best_hits <- h
        best_ov <- ov
        best_dv <- dv
      }
    }
  }

  lut <- setNames(d$best_singlet, best_dv)

  donor_assign_variant <- unname(lut[best_ov])
  names(donor_assign_variant) <- best_ov

  donor_final <- unname(donor_assign_variant[variants_obj$obj])

  obj$donor <- donor_final
  obj
}

# ------------------------------------------------------------------------------
# A4) Per-sample processing function
#     - Load RNA+ATAC from 10x ARC H5
#     - Build Seurat + ChromatinAssay
#     - Do quick RNA clustering/UMAP for SoupX
#     - SoupX RNA decontamination
#     - DropletQC nuclear fraction + empty droplet labeling
# ------------------------------------------------------------------------------
process_one_multiome <- function(sample_id, condition, outs_dir, gene_annot) {
  message("\n==================== Processing ", sample_id, " (", condition, ") ====================")

  if (!dir.exists(outs_dir)) stop("outs_dir missing: ", outs_dir)

  filt_h5 <- file.path(outs_dir, "filtered_feature_bc_matrix.h5")
  raw_h5  <- file.path(outs_dir, "raw_feature_bc_matrix.h5")

  if (!file.exists(filt_h5)) stop("Missing filtered_feature_bc_matrix.h5 in ", outs_dir)
  if (!file.exists(raw_h5))  stop("Missing raw_feature_bc_matrix.h5 in ", outs_dir)

  # 1) Load 10x ARC matrices
  filt_list <- Read10X_h5(filt_h5)
  raw_list  <- Read10X_h5(raw_h5)

  stopifnot(
    "Gene Expression" %in% names(filt_list),
    "Peaks" %in% names(filt_list),
    "Gene Expression" %in% names(raw_list)
  )

  rna_filt  <- filt_list[["Gene Expression"]]
  atac_filt <- filt_list[["Peaks"]]
  rna_raw   <- raw_list[["Gene Expression"]]

  # 2) Ensure raw/filt RNA have same gene rows for SoupX
  common_genes <- intersect(rownames(rna_raw), rownames(rna_filt))
  rna_raw  <- rna_raw[common_genes, , drop = FALSE]
  rna_filt <- rna_filt[common_genes, , drop = FALSE]

  # 3) Create Seurat object (RNA) + add condition metadata
  so <- CreateSeuratObject(counts = rna_filt, assay = "RNA", project = sample_id)
  so$condition <- condition

  # 4) Attach ATAC ChromatinAssay
  sep_use <- guess_sep(rownames(atac_filt))
  frag    <- pick_frag(outs_dir)

  so[["ATAC"]] <- CreateChromatinAssay(
    counts     = atac_filt,
    sep        = sep_use,
    genome     = "hg38",
    fragments  = frag,
    annotation = gene_annot
  )

  # 5) Quick RNA processing ONLY to support SoupX contamination estimation
  DefaultAssay(so) <- "RNA"
  so <- NormalizeData(so)
  so <- FindVariableFeatures(so)
  so <- ScaleData(so, features = VariableFeatures(so), verbose = FALSE)
  so <- RunPCA(so, npcs = 30, verbose = FALSE)
  so <- FindNeighbors(so, dims = 1:30)
  so <- FindClusters(so, resolution = 0.3)
  so <- RunUMAP(so, dims = 1:30, verbose = FALSE)

  # 6) SoupX ambient RNA correction
  sc <- SoupChannel(tod = rna_raw, toc = rna_filt)
  sc <- setClusters(sc, setNames(as.character(so$seurat_clusters), Cells(so)))
  sc <- setDR(sc, Embeddings(so, "umap"))
  sc <- autoEstCont(sc)

  clean_rna <- adjustCounts(sc, roundToInt = TRUE)

  # Replace RNA counts with decontaminated counts
  so <- SetAssayData(so, assay = "RNA", slot = "counts", new.data = clean_rna)

  # 7) DropletQC nuclear fraction + empty droplet classification
  bam <- pick_bam(outs_dir)
  bcs <- colnames(so)

  nf <- nuclear_fraction_tags(
    bam = bam,
    barcodes = bcs,
    tiles = 100,
    cores = 4,
    verbose = TRUE
  )

  so$nuclear_fraction <- nf[bcs, "nuclear_fraction"]

  qc_df <- data.frame(
    nuclear_fraction = so$nuclear_fraction,
    umi = so$nCount_RNA
  )

  empty_qc <- identify_empty_drops(
    nf_umi       = qc_df,
    nf_rescue    = 0.05,
    umi_rescue   = 1000,
    include_plot = FALSE
  )

  so$droplet_qc <- empty_qc$cell_status

  DefaultAssay(so) <- "RNA"
  so
}

# ------------------------------------------------------------------------------
# A5) Run per-sample preprocessing for the three selected samples
# ------------------------------------------------------------------------------
objs <- setNames(vector("list", nrow(sample_info)), sample_info$sample_id)

for (i in seq_len(nrow(sample_info))) {
  objs[[sample_info$sample_id[i]]] <- process_one_multiome(
    sample_id  = sample_info$sample_id[i],
    condition  = sample_info$condition[i],
    outs_dir   = sample_info$outs_dir[i],
    gene_annot = gene_annot
  )
}

# ------------------------------------------------------------------------------
# A6) Attach donor labels (Vireo best_singlet) for each selected sample
# ------------------------------------------------------------------------------
stopifnot(all(file.exists(sample_info$donor_tsv)))

for (i in seq_len(nrow(sample_info))) {
  sid <- sample_info$sample_id[i]
  message("Attaching Vireo donor labels for ", sid)
  objs[[sid]] <- attach_best_singlet(objs[[sid]], sample_info$donor_tsv[i])
}

# Optional convenience aliases for the three samples
if ("Ctrl" %in% names(objs)) ctrl_obj <- objs[["Ctrl"]]
if ("Stim" %in% names(objs)) stim_obj <- objs[["Stim"]]
if ("Prime" %in% names(objs)) prime_obj <- objs[["Prime"]]

# ==============================================================================
# End of Part 1: per-sample preprocessing
# ==============================================================================
# The next section in the original workflow begins donor harmonization across samples

