#=======================================================================
# Shared logic for labelling Seurat clusters by a confidence-weighted
# majority vote of MapMyCells per-cell calls.
#
# Sourced by label_pv_integrated.R and label_geo_object.R - edit the method
# here once and both objects stay consistent.
#=======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
})

# The nine samples in the integrated object.
MAPMYCELLS_SAMPLES <- c("sample1", "sample2", "sample3", "sample4",
                        "sample7", "sample8", "sample9", "sample11", "sample12")

LABEL_COL  <- "supercluster_name"                      # the per-cell call
WEIGHT_COL <- "subcluster_bootstrapping_probability"   # the per-cell confidence

# ----------------------------------------------------------------------
# read the MapMyCells calls
# ----------------------------------------------------------------------
# One CSV per sample, keyed by raw 10x barcode. Rebuilds the integrated
# object's barcode namespace: sample1 + AAACAGCC...-1 -> Sample1_AAACAGCC...-1
read_mapmycells <- function(dir, samples = MAPMYCELLS_SAMPLES) {
  calls <- bind_rows(lapply(samples, function(s) {
    path <- file.path(dir, paste0(s, ".csv"))
    if (!file.exists(path)) stop("Missing MapMyCells file: ", path)
    df <- read.csv(path, comment.char = "#", header = TRUE,
                   stringsAsFactors = FALSE, check.names = FALSE)
    missing <- setdiff(c("cell_id", LABEL_COL, WEIGHT_COL), colnames(df))
    if (length(missing)) stop(s, " is missing column(s): ",
                              paste(missing, collapse = ", "))
    tibble(barcode = paste0(str_to_title(s), "_", df$cell_id),
           !!LABEL_COL  := df[[LABEL_COL]],
           !!WEIGHT_COL := as.numeric(df[[WEIGHT_COL]]))
  }))
  stopifnot(!any(duplicated(calls$barcode)))
  message("Read ", nrow(calls), " MapMyCells calls across ",
          length(samples), " samples")
  calls
}

# ----------------------------------------------------------------------
# collapse MapMyCells superclusters into broad categories
# ----------------------------------------------------------------------
ASTRO_SUPERCLUSTERS <- c("Astrocyte", "Bergmann glia")
NEURON_SUPERCLUSTERS <- c(
  "Cerebellar inhibitory", "CGE interneuron",
  "Deep-layer corticothalamic and 6b", "Deep-layer intratelencephalic",
  "Deep-layer near-projecting", "Eccentric medium spiny neuron",
  "Hippocampal CA1-3", "Hippocampal CA4",
  "LAMP5-LHX6 and Chandelier", "Medium spiny neuron",
  "MGE interneuron", "Midbrain-derived inhibitory",
  "Upper-layer intratelencephalic", "Upper rhombic lip"
)
STROMAL_SUPERCLUSTERS <- c("Fibroblast", "Vascular")

# ----------------------------------------------------------------------
# cluster ordering
# ----------------------------------------------------------------------
# Orders labels like "Astrocyte_1", "Astrocyte_2", ..., "Oligodendrocyte_10"
# by prefix then by NUMERIC suffix, so _10 sorts after _2 rather than after _1
# as plain alphabetical sorting would have it.
#
# Only the order *within* a category matters, because the _1.._n suffixes are
# handed out by row_number() inside each voted category.
natural_cluster_order <- function(labels) {
  u    <- unique(as.character(labels))
  base <- sub("_[0-9]+$", "", u)
  num  <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", u)))
  num[is.na(num)] <- 0L          # labels with no numeric suffix
  u[order(base, num)]
}

collapse_supercluster <- function(x) {
  case_when(
    x %in% ASTRO_SUPERCLUSTERS   ~ "Astrocyte",
    x %in% NEURON_SUPERCLUSTERS  ~ "Neuron",
    x %in% STROMAL_SUPERCLUSTERS ~ "Stromal",
    TRUE                         ~ x           # keep glia etc. as they are
  )
}

# ----------------------------------------------------------------------
# the vote
# ----------------------------------------------------------------------
# Transfers the calls onto `obj`, votes per cluster, and writes
# category / assigned_label / vote_score / norm_vote_score into the metadata.
#
# `cluster_col` must already exist in obj@meta.data. Its ordering matters: the
# _1.._n suffixes are handed out in cluster order, so a factor with numeric
# levels (seurat_clusters) keeps the historical numbering.
label_by_vote <- function(obj, calls, cluster_col = "seurat_clusters") {

  stopifnot(cluster_col %in% colnames(obj@meta.data))

  # Drop any labelling columns already present, so the vote is genuinely
  # re-derived rather than inherited (and left_join cannot make .x/.y columns).
  stale <- c(LABEL_COL, WEIGHT_COL, "category",
             "assigned_label", "vote_score", "norm_vote_score")
  for (cc in intersect(stale, colnames(obj@meta.data))) obj[[cc]] <- NULL

  meta <- obj@meta.data %>%
    rownames_to_column("barcode") %>%
    as_tibble() %>%
    left_join(calls %>% select(barcode, all_of(c(LABEL_COL, WEIGHT_COL))),
              by = "barcode")

  matched <- sum(!is.na(meta[[LABEL_COL]]))
  message("Matched ", matched, "/", nrow(meta), " cells (",
          round(100 * matched / nrow(meta), 2), "%)")
  if (matched == 0) stop("No barcodes matched - check the sample name prefixes.")

  meta <- meta %>% mutate(category = collapse_supercluster(.data[[LABEL_COL]]))

  # Each cell votes for its category with weight = its mapping confidence, so a
  # confident call counts for more than a marginal one. Unmapped cells (NA)
  # contribute nothing via na.rm = TRUE.
  votes <- meta %>%
    group_by(.data[[cluster_col]], category) %>%
    summarise(vote_score = sum(.data[[WEIGHT_COL]], na.rm = TRUE),
              .groups = "drop")

  # Winning category per cluster.
  assigned <- votes %>%
    group_by(.data[[cluster_col]]) %>%
    slice_max(vote_score, n = 1, with_ties = FALSE) %>%
    ungroup()

  # Where several clusters win the same category, suffix them _1.._n in
  # cluster order so every cluster keeps a distinct label.
  assigned <- assigned %>%
    group_by(category) %>%
    mutate(assigned_label = if (n() == 1) category
                            else paste0(category, "_", row_number())) %>%
    ungroup()

  # norm_vote_score = mean per-cell confidence behind the winning call, i.e.
  # how clean the cluster's identity is.
  assigned <- assigned %>%
    left_join(meta %>% count(.data[[cluster_col]], name = "n_cells"),
              by = cluster_col) %>%
    mutate(norm_vote_score = vote_score / n_cells)

  meta <- meta %>%
    left_join(assigned %>% select(all_of(cluster_col), assigned_label,
                                  vote_score, norm_vote_score),
              by = cluster_col)

  obj@meta.data <- meta %>% column_to_rownames("barcode") %>% as.data.frame()
  stopifnot(identical(rownames(obj@meta.data), colnames(obj)))

  assignment_table <- assigned %>%
    select(all_of(cluster_col), assigned_label, category,
           vote_score, n_cells, norm_vote_score)

  list(object = obj, assignments = assignment_table)
}

# ----------------------------------------------------------------------
# plots
# ----------------------------------------------------------------------
save_label_plots <- function(obj, reduction, out_dir, suffix = "") {
  if (!reduction %in% names(obj@reductions))
    stop("Reduction not found: ", reduction)
  f <- function(name) file.path(out_dir, paste0(name, suffix, ".png"))

  ggsave(f("umap_assigned_label"),
    DimPlot(obj, reduction = reduction, group.by = "assigned_label",
            label = TRUE, repel = TRUE) +
      ggtitle("Majority-vote cell type"),
    width = 12, height = 8, dpi = 300)

  ggsave(f("umap_supercluster_name"),
    DimPlot(obj, reduction = reduction, group.by = LABEL_COL,
            label = TRUE, repel = TRUE) +
      ggtitle("MapMyCells supercluster (per cell)"),
    width = 14, height = 8, dpi = 300)

  ggsave(f("norm_vote_score"),
    FeaturePlot(obj, features = "norm_vote_score",
                reduction = reduction, pt.size = 0.5) +
      scale_colour_viridis_c(option = "viridis", name = "Norm.\nvote score") +
      theme_minimal() +
      ggtitle("Per-cluster vote confidence"),
    width = 10, height = 8, dpi = 300)

  # Per-cell MapMyCells bootstrapping probability - the raw confidence each
  # cell contributes to the vote, before any per-cluster aggregation. Cells
  # with no MapMyCells call are grey.
  n_unmapped <- sum(is.na(obj@meta.data[[WEIGHT_COL]]))
  ggsave(f("bootstrap_probability"),
    FeaturePlot(obj, features = WEIGHT_COL,
                reduction = reduction, pt.size = 0.5) +
      scale_colour_viridis_c(option = "viridis", na.value = "lightgrey",
                             name = "Bootstrap\nprobability") +
      theme_minimal() +
      ggtitle("Per-cell MapMyCells bootstrapping probability",
              subtitle = if (n_unmapped > 0)
                paste0(n_unmapped, " unmapped cells in grey") else NULL),
    width = 10, height = 8, dpi = 300)

  message("Wrote plots on reduction '", reduction, "' to ", out_dir)
}

# Directory of the running script, so drivers can resolve paths relative to it.
this_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
