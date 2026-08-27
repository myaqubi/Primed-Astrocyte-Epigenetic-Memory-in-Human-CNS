#=======================================================================
# Label the MS/CTRL GEO object by confidence-weighted majority vote of
# MapMyCells calls per cluster, and plot on the umap_MP (metaprogram) embedding.
#
# STANDALONE: the GEO object plus the MapMyCells CSVs, nothing else.
#
# Inputs   data/MS_CTRL_final_object_for_GEO.rds   (52,080 cells)
#          data/mapmycells/sample*.csv             MapMyCells calls
#
# Outputs  output/MS_CTRL_final_object_for_GEO_labelled.rds
#          output/geo_cluster_label_assignments.csv
#          output/geo_comparison_vs_Cluster_label.csv
#          output/umap_assigned_label_MP.png
#          output/umap_supercluster_name_MP.png
#          output/norm_vote_score_MP.png
#          output/bootstrap_probability_MP.png
#
# Usage:   Rscript label_geo_object.R
#
# On clustering: this object has no seurat_clusters column (it was dropped when
# the object was prepared for GEO) and Idents() holds Cell_label, the 13 cell
# types - not the 22 clusters. What it does ship is Cluster_label, a 1:1
# relabelling of those 22 clusters, which carries the cluster memberships the
# vote needs.
#
# Cluster_label is stored as plain character, so it carries no ordering, and
# ordering matters here: the _1.._n suffixes are handed out in cluster order,
# and R sorts character labels alphabetically, putting "Oligodendrocyte_10"
# before "Oligodendrocyte_2". natural_cluster_order() sorts by prefix then by
# numeric suffix instead, which restores the original numbering from the
# object's own labels - no external ordering file needed.
#
# The consequence is worth being explicit about: the *suffix numbers* are
# inherited from Cluster_label. The cluster memberships, the winning category
# for each cluster, and every score are derived independently from the
# MapMyCells calls.
#=======================================================================

script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
})
source(file.path(script_dir, "vote_labelling.R"))

set.seed(1)

in_object   <- file.path(script_dir, "data", "MS_CTRL_final_object_for_GEO.rds")
mapmycells  <- file.path(script_dir, "data", "mapmycells")
out_dir     <- file.path(script_dir, "output")
cluster_col <- "cluster"
reduction   <- "umap_MP"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- load ------------------------------------------------------------
message("Loading ", in_object)
geo <- readRDS(in_object)
message("  ", ncol(geo), " cells x ", nrow(geo), " features")
message("  reductions: ", paste(names(geo@reductions), collapse = ", "))

for (needed in c("Cluster_label")) {
  if (!needed %in% colnames(geo@meta.data))
    stop("Expected metadata column missing from the object: ", needed)
}
if (!reduction %in% names(geo@reductions))
  stop("Expected reduction missing from the object: ", reduction)

# ---- cluster grouping, in cluster order ------------------------------
# Same grouping as Cluster_label, but as an ordered factor with _10 after _2.
geo$cluster <- factor(as.character(geo$Cluster_label),
                      levels = natural_cluster_order(geo$Cluster_label))
stopifnot(!any(is.na(geo$cluster)))
message("  ", nlevels(geo$cluster), " clusters, ordered: ",
        paste(head(levels(geo$cluster), 3), collapse = ", "), ", ...")

previous_labels <- geo@meta.data[, c(cluster_col, "Cluster_label"), drop = FALSE]

# ---- vote ------------------------------------------------------------
calls <- read_mapmycells(mapmycells)
res   <- label_by_vote(geo, calls, cluster_col = cluster_col)
geo   <- res$object

assignment_table <- res$assignments %>%
  arrange(cluster) %>%
  select(cluster, assigned_label, category, vote_score, n_cells, norm_vote_score)

message("\nCluster assignments:")
print(as.data.frame(assignment_table), row.names = FALSE)
write.csv(assignment_table,
          file.path(out_dir, "geo_cluster_label_assignments.csv"),
          row.names = FALSE)

# ---- compare against the Cluster_label shipped in the object ---------
cmp <- previous_labels %>%
  as_tibble() %>%
  rename(stored = Cluster_label) %>%
  mutate(recomputed = as.character(geo$assigned_label)) %>%
  distinct() %>%
  mutate(stored_category     = sub("_[0-9]+$", "", as.character(stored)),
         recomputed_category = sub("_[0-9]+$", "", recomputed),
         category_agrees     = stored_category == recomputed_category,
         label_agrees        = as.character(stored) == recomputed) %>%
  arrange(cluster)

message("\nComparison vs. Cluster_label shipped in the GEO object:")
message("  broad category agrees for ", sum(cmp$category_agrees), "/",
        nrow(cmp), " clusters")
message("  exact label agrees for ", sum(cmp$label_agrees), "/", nrow(cmp),
        " clusters")
if (any(!cmp$label_agrees)) {
  message("  clusters that differ:")
  print(as.data.frame(cmp %>% filter(!label_agrees) %>%
                        select(cluster, stored, recomputed)),
        row.names = FALSE)
}
write.csv(as.data.frame(cmp),
          file.path(out_dir, "geo_comparison_vs_Cluster_label.csv"),
          row.names = FALSE)

# ---- plots and save --------------------------------------------------
save_label_plots(geo, reduction, out_dir, suffix = "_MP")

out_object <- file.path(out_dir, "MS_CTRL_final_object_for_GEO_labelled.rds")
message("\nSaving ", out_object)
saveRDS(geo, file = out_object)
message("Done.")
