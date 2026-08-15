# ============================================================================== 
# Cell type annotation and astrocyte refinement
# ============================================================================== 
#
# Purpose:
#   1) exporting the integrated object for MapMyCells annotation,
#   2) importing the annotation output back into the Seurat object,
#   3) assigning cluster-level labels based on MapMyCells plus marker-gene review,
#   4) removing the postnatal donor and non-astrocyte clusters to generate an
#      astrocyte-focused object.
#
# Notes:
#   - The required annotation file for this step is generated from the Seurat
#     object and analyzed using the online MapMyCells workflow (standard
#     approach). The resulting annotation matrix is then imported back into the
#     object.
#   - Replace the generic placeholder paths below with your own project paths.
# ============================================================================== 

# ------------------------------------------------------------------------------
# 0) Packages
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(Matrix)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 1) Project paths and input object
# ------------------------------------------------------------------------------
project_dir <- "/path/to/your/project"
output_dir  <- "/path/to/your/output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_object_path <- file.path(output_dir, "cleaned_object_noDoublets_integrated.rds")
mapmycells_output_file <- file.path(project_dir, "In_vitro_three_condition_for_MapMyCells.csv")

if (!file.exists(input_object_path)) {
  stop("Input object not found at: ", input_object_path)
}

# Load the integrated object from the previous step
integrated_obj <- readRDS(input_object_path)

# Keep the sample labels consistent with the workflow
sample_names <- c("Ctrl", "Stim", "Prime")

# ------------------------------------------------------------------------------
# 2) Export the object for MapMyCells annotation
# ------------------------------------------------------------------------------
# The required file for this step is generated from the integrated object and
# analyzed using the online MapMyCells workflow (standard approach). The output
# annotation file is then read back into the Seurat object.

DefaultAssay(integrated_obj) <- "RNA"
counts_mat <- GetAssayData(
  object = integrated_obj,
  assay  = "RNA",
  slot   = "counts"
)

if (nrow(counts_mat) == 0 || ncol(counts_mat) == 0) {
  stop("Counts matrix is empty. Check that the RNA assay has counts.")
}

sce <- SingleCellExperiment(assays = list(counts = counts_mat))
colnames(sce) <- colnames(counts_mat)
rownames(sce) <- rownames(counts_mat)

writeH5AD(
  sce,
  file = file.path(project_dir, "In_vitro_three_condition_for_MapMyCells_raw_counts.h5ad"),
  X_name = "counts",
  compression = "gzip"
)

message("Written raw counts file for MapMyCells: ",
        file.path(project_dir, "In_vitro_three_condition_for_MapMyCells_raw_counts.h5ad"))

# ------------------------------------------------------------------------------
# 3) Read MapMyCells annotation results back into the Seurat object
# ------------------------------------------------------------------------------
if (!file.exists(mapmycells_output_file)) {
  stop("MapMyCells output file not found: ", mapmycells_output_file)
}

mapmycells_df <- read.csv(
  mapmycells_output_file,
  header = TRUE,
  check.names = FALSE,
  comment.char = "#"
)

rownames(mapmycells_df) <- mapmycells_df$cell_id
mapmycells_df <- mapmycells_df[colnames(integrated_obj), ]

integrated_obj <- AddMetaData(integrated_obj, metadata = mapmycells_df)

message("MapMyCells annotations have been added to the object.")

# ------------------------------------------------------------------------------
# 4) Assign cluster-level labels based on MapMyCells and manual marker review
# ------------------------------------------------------------------------------
# The main goal here is to create a consistent cluster-level annotation for the
# downstream analysis and to distinguish true astrocytes from other populations.

Idents(integrated_obj) <- "seurat_clusters"
clusters_to_drop <- c("15", "18") # these clusters composed of post-natal cells
all_clusters <- as.character(levels(Idents(integrated_obj)))
clusters_to_annotate <- setdiff(all_clusters, clusters_to_drop)

# Summary table for cluster review
cluster_df <- integrated_obj@meta.data %>%
  mutate(
    cluster = as.character(Idents(integrated_obj)),
    supercluster_name = as.character(.data$supercluster_name),
    supercluster_boot = suppressWarnings(as.numeric(.data$supercluster_bootstrapping_probability))
  ) %>%
  filter(cluster %in% clusters_to_annotate) %>%
  filter(!is.na(supercluster_name) & supercluster_name != "")

breakdown <- cluster_df %>%
  count(cluster, supercluster_name, name = "n_cells") %>%
  group_by(cluster) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  arrange(as.integer(cluster), desc(prop), desc(n_cells))

boot_median_by_sc <- cluster_df %>%
  group_by(cluster, supercluster_name) %>%
  summarise(
    n_cells = n(),
    median_boot = stats::median(supercluster_boot, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(as.integer(cluster), desc(median_boot), desc(n_cells))

summary_table <- breakdown %>%
  left_join(
    boot_median_by_sc %>% select(cluster, supercluster_name, median_boot),
    by = c("cluster", "supercluster_name")
  ) %>%
  arrange(as.integer(cluster), desc(prop), desc(n_cells))

print(summary_table, n = 200)

# ------------------------------------------------------------------------------
# 5) Cluster-to-label mapping
# ------------------------------------------------------------------------------
# Astrocyte clusters are assigned names based on cluster order.
astro_clusters <- c("0", "1", "3", "7", "8", "9", "10", "11", "13")
astro_clusters_sorted <- astro_clusters[order(as.integer(astro_clusters))]
astro_labels <- paste0("Astrocyte_", seq_along(astro_clusters_sorted))
astro_map <- setNames(astro_labels, astro_clusters_sorted)

# Other labels
other_map <- c(
  "2"  = "Oligodendrocyte",
  "5"  = "Microglia_1",
  "12" = "Immune_like_astrocyte",
  "14" = "Cycling",
  "4"  = "mixed_glial_stromal",
  "6"  = "mixed_ast-vas-like",
  "16" = "mixed_OPC-astro"
)

cluster_to_final <- c(astro_map, other_map)

missing_clusters <- setdiff(clusters_to_annotate, names(cluster_to_final))
if (length(missing_clusters) > 0) {
  stop("No final_annotation assigned for these clusters: ",
       paste(missing_clusters, collapse = ", "))
}

# ------------------------------------------------------------------------------
# 6) Remove the postnatal donor and associated non-astrocyte populations
# ------------------------------------------------------------------------------
# There is one postnatal sample in the object, and it is almost entirely composed
# of mature oligodendrocytes and some microglia cells. To focus the downstream
# analysis on prenatal-derived astrocytes, we remove that donor and exclude
# oligodendrocyte and microglial populations from the astrocyte-focused object.

# We remove the postnatal donor using the donor label from the harmonized donor
# metadata. If your donor mapping uses a different label, update this line.
if ("donor_global_pretty" %in% colnames(integrated_obj@meta.data)) {
  obj_prenatal <- subset(integrated_obj, subset = donor_global_pretty != "donor2")
} else if ("Donor" %in% colnames(integrated_obj@meta.data)) {
  obj_prenatal <- subset(integrated_obj, subset = Donor != "Astrocyte3")
} else {
  obj_prenatal <- integrated_obj
}

# Remove the clusters that are not considered bona fide astrocytes
clusters_to_remove <- c(
  "Oligodendrocyte",
  "Microglia_1",
  "Cycling",
  "mixed_ast-vas-like",
  "mixed_glial_stromal",
  "mixed_OPC-astro",
  "Immune_like_astrocyte"
)

# Keep only the astrocyte-focused populations
obj_astrocytes <- subset(
  obj_prenatal,
  subset = final_annotation %in% grep("^Astrocyte_", unique(obj_prenatal$final_annotation), value = TRUE)
)

# ------------------------------------------------------------------------------
# 8) Final object and a brief summary
# ------------------------------------------------------------------------------
message("Final astrocyte-focused object created.")
message("Cells retained: ", ncol(obj_astrocytes))

print(table(obj_astrocytes$final_annotation, useNA = "ifany"))

# ------------------------------------------------------------------------------
# 9) Optional final object save
# ------------------------------------------------------------------------------
# If you want to save the refined object for downstream analyses, uncomment the
# line below and adjust the output path.
# saveRDS(obj_astrocytes, file = file.path(output_dir, "In_vitro_prenatal_samples_QC_filtered_integrated_annotated_astrocytes.rds"))
