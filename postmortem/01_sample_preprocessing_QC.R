# 01_sample_preprocessing_QC.R
# Per-sample preprocessing of multiome data.
# Steps: object construction, preliminary RNA clustering for SoupX,
# ambient-RNA correction, nuclear-fraction calculation, and DropletQC.
# Update `base_dir` before running.

library(DropletQC)
library(SoupX)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)

# ─── Parameters ─────────────────────────────────────────────────────────────
base_dir <- "PATH/TO/MULTIOME_DATA"

# map folder names → object names
samples <- list(
  Sample1  = "sample1_CTRL",
  Sample2  = "sample2_MS",
  Sample3  = "sample3_CTRL",
  Sample4  = "sample4_MS",
  Sample7  = "sample7_CTRL",
  Sample8  = "sample8_MS",
  Sample9  = "sample9_MS",
  Sample11 = "sample11_CTRL",
  Sample12 = "sample12_MS"
)

# build gene annotation once
gene_annot <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
seqlevelsStyle(gene_annot) <- "UCSC"
genome(gene_annot)       <- "hg38"

# ─── Loop over each sample ─────────────────────────────────────────────────
for (folder in names(samples)) {
  
  obj_name <- samples[[folder]]
  message("Processing ", folder, " → ", obj_name)
  
  # set sample‐specific paths
  samp_dir    <- file.path(base_dir, folder)
  filt_h5     <- file.path(samp_dir, "filtered_feature_bc_matrix.h5")
  raw_h5      <- file.path(samp_dir, "raw_feature_bc_matrix.h5")
  frag_path   <- file.path(samp_dir, "atac_fragments.tsv.gz")
  bam_path    <- file.path(samp_dir, "possorted_genome_bam.bam")
  
  # 1) Read in filtered RNA & ATAC
  all_counts  <- Read10X_h5(filt_h5)
  rna_counts  <- all_counts$`Gene Expression`
  peak_counts <- all_counts$Peaks
  
  # 2) Read raw RNA for SoupX
  raw_all <- Read10X_h5(raw_h5)
  raw_rna <- raw_all$`Gene Expression`
  
  # 3) Create Seurat & add assays
  so <- CreateSeuratObject(counts = rna_counts, assay = "RNA", project = obj_name)
  so[["ATAC"]] <- CreateChromatinAssay(
    counts     = peak_counts,
    sep        = c(":", "-"),
    genome     = 'hg38',
    fragments  = frag_path,
    annotation = gene_annot
  )
  
  # 4) RNA preprocessing & clustering
  so <- NormalizeData(so)
  so <- FindVariableFeatures(so)
  so <- ScaleData(so)
  so <- RunPCA(so, dims = 1:30, verbose = FALSE)
  so <- RunUMAP(so, dims = 1:30, verbose = FALSE)
  so <- FindNeighbors(so, dims = 1:30)
  so <- FindClusters(so, resolution = 0.5)
  
  # 5) SoupX ambient‐RNA correction
  sc <- SoupChannel(tod = raw_rna, toc = rna_counts)
  sc <- setClusters(sc, setNames(so$seurat_clusters, Cells(so)))
  sc <- setDR(sc, so@reductions$umap@cell.embeddings)
  sc <- autoEstCont(sc)
  clean_rna <- adjustCounts(sc, roundToInt = TRUE)
  so <- SetAssayData(so, assay = "RNA", layer = "counts", new.data = clean_rna)
  
  # 6) Compute nuclear fraction
  bc       <- colnames(so)
  nf       <- nuclear_fraction_tags(bam = bam_path, barcodes = bc, tiles = 100, cores = 4)
  so$nuclear_fraction <- nf[bc, "nuclear_fraction"]
  if (any(is.na(so$nuclear_fraction))) {
    warning("Some barcodes NA for nuclear_fraction in ", obj_name)
  }
  
  # 7) Identify & remove empty droplets
  qc_df   <- data.frame(nuclear_fraction = so$nuclear_fraction, umi = so$nCount_RNA)
  empty_qc <- identify_empty_drops(nf_umi = qc_df,
                                   nf_rescue = 0.05,
                                   umi_rescue = 1000,
                                   include_plot = TRUE)
  so$droplet_qc <- empty_qc$cell_status
  so <- subset(so, subset = droplet_qc == "cell")
  
  # 8) Save into global env
  assign(x = obj_name, value = so, envir = .GlobalEnv)
  
  message("→ Done: ", obj_name, "\n")
}
