# ============================================================================== 
# Donor label harmonization across samples
# ============================================================================== 
#
# Purpose:
#   This script performs donor-label harmonization for the in vitro multiome
#   analysis after the per-sample preprocessing step.
#
# Notes:
#   - The original `donor` column contains sample-specific donor assignments
#     generated independently during demultiplexing. As a result, donor labels
#     (for example `donor0`–`donor3`) are only valid within an individual
#     sample and do not necessarily correspond to the same biological donor
#     across different samples.
#   - To enable donor-matched analyses across conditions, these local labels
#     are harmonized into a unified identifier stored in the `donor_global`
#     column.
#   - Using `Ctrl` as the reference sample, each donor is assigned a unique
#     global identifier (`D_A`–`D_D`), and donor mapping information is used to
#     translate the local donor labels in the remaining samples to the
#     corresponding global donor identity.
#   - Finally, the `donor_global_pretty` column provides a display-friendly
#     version of the harmonized donor labels, and the `Donor` column is added
#     to convert these labels into more interpretable names.
#
# Important:
#   - Replace the generic placeholder path below with your actual project path.
#   - This script assumes that the per-sample preprocessing step has already
#     produced the `objs` list containing the Seurat objects for Ctrl, Stim,
#     and Prime.
# ============================================================================== 

# ------------------------------------------------------------------------------
# 0) Packages
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# 1) Project root and sample information
# ------------------------------------------------------------------------------
project_dir <- "/path/to/your/project"
anchor_sample <- "Ctrl"
selected_samples <- c("Ctrl", "Stim", "Prime")

if (!exists("objs") || !is.list(objs)) {
  stop("The 'objs' list was not found. Please run the per-sample preprocessing script first.")
}

if (!all(selected_samples %in% names(objs))) {
  stop("The expected sample objects were not found in 'objs'.")
}

map_path <- file.path(project_dir, "Demultiplexing", "donor_mapping.csv")

if (!file.exists(map_path)) {
  stop("Donor mapping file not found at: ", map_path)
}

# ------------------------------------------------------------------------------
# 2) Helper functions
# ------------------------------------------------------------------------------
read_mapping <- function(path) {
  tries <- list(
    function(p) readr::read_csv(p, trim_ws = TRUE, show_col_types = FALSE),
    function(p) readr::read_delim(p, delim = "\t", trim_ws = TRUE, show_col_types = FALSE),
    function(p) readr::read_delim(p, delim = ";", trim_ws = TRUE, show_col_types = FALSE),
    function(p) readr::read_table(p, trim_ws = TRUE, show_col_types = FALSE)
  )

  for (fn in tries) {
    df <- tryCatch(fn(path), error = function(e) NULL)
    if (!is.null(df) && ncol(df) >= 4) return(df)
  }

  lines <- readLines(path, warn = FALSE)
  stopifnot(length(lines) >= 2)

  header_tokens <- strsplit(lines[1], "\\s+")[[1]]
  body_tokens   <- strsplit(lines[-1], "\\s+")
  ok <- vapply(body_tokens, length, integer(1)) == 4

  if (!all(ok)) stop("Fallback parse failed: not all rows split into 4 tokens.")

  mat <- do.call(rbind, body_tokens)
  df  <- as.data.frame(mat, stringsAsFactors = FALSE)

  colnames(df) <- if (length(header_tokens) == 4) header_tokens else
    c("ReferenceCondition", "ReferenceDonor", "OtherCondition", "MatchedDonor")

  df
}

norm_names <- function(x) {
  x <- gsub("[\uFEFF\u200B]", "", x)
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
}

pick_col <- function(df, aliases) {
  nm <- intersect(aliases, names(df))
  if (length(nm) == 0) return(NULL)
  nm[1]
}

# ------------------------------------------------------------------------------
# 3) Load and normalize the donor mapping table
# ------------------------------------------------------------------------------
don_map_raw <- read_mapping(map_path)
names(don_map_raw) <- norm_names(names(don_map_raw))

ref_cond_col    <- pick_col(don_map_raw, c("referencecondition", "reference_condition", "reference", "sample1", "anchor_condition"))
ref_donor_col   <- pick_col(don_map_raw, c("referencedonor", "reference_donor", "ref_donor", "donor_ref", "reference_donor_id"))
other_cond_col  <- pick_col(don_map_raw, c("othercondition", "other_condition", "other", "sample2", "target_condition"))
match_donor_col <- pick_col(don_map_raw, c("matcheddonor", "matched_donor", "otherdonor", "donor_other", "matched_donor_id"))

if (any(sapply(list(ref_cond_col, ref_donor_col, other_cond_col, match_donor_col), is.null))) {
  stop("Could not detect all donor mapping columns. Found: ", paste(names(don_map_raw), collapse = ", "))
}

don_map <- don_map_raw %>%
  transmute(
    ReferenceCondition = .data[[ref_cond_col]],
    ReferenceDonor     = .data[[ref_donor_col]],
    OtherCondition     = .data[[other_cond_col]],
    MatchedDonor       = .data[[match_donor_col]]
  )

# ------------------------------------------------------------------------------
# 4) Build global donor IDs anchored on Ctrl
# ------------------------------------------------------------------------------
if (!(anchor_sample %in% names(objs))) {
  stop("Anchor sample '", anchor_sample, "' is not present among selected samples.")
}

# Global donor IDs are anchored to the reference sample.
global_levels <- c("D_A", "D_B", "D_C", "D_D")
ref_to_global <- setNames(global_levels, c("donor0", "donor1", "donor2", "donor3"))

mk_map_for <- function(target_sample, don_map, ref_to_global, anchor_sample = "Ctrl") {
  if (target_sample == anchor_sample) {
    return(ref_to_global)
  }

  x <- don_map %>%
    dplyr::filter(
      ReferenceCondition == anchor_sample,
      OtherCondition == target_sample
    ) %>%
    dplyr::mutate(
      global = ref_to_global[ReferenceDonor]
    ) %>%
    dplyr::select(MatchedDonor, global)

  if (nrow(x) == 0) {
    stop("No rows found in donor mapping for ", anchor_sample, " -> ", target_sample)
  }

  setNames(x$global, x$MatchedDonor)
}

maps <- lapply(selected_samples, function(sid) {
  mk_map_for(
    target_sample = sid,
    don_map = don_map,
    ref_to_global = ref_to_global,
    anchor_sample = anchor_sample
  )
})

names(maps) <- selected_samples

# ------------------------------------------------------------------------------
# 5) Apply donor harmonization to each Seurat object
# ------------------------------------------------------------------------------
harmonize_one <- function(obj, sample_name, maps) {
  d_local <- as.character(obj$donor)
  obj$donor_global <- unname(maps[[sample_name]][d_local])
  if (!"sample" %in% colnames(obj@meta.data)) obj$sample <- sample_name
  obj
}

for (sid in selected_samples) {
  objs[[sid]] <- harmonize_one(objs[[sid]], sid, maps)
}

# Pretty labels for donor_global_pretty
pretty_map <- c(D_A = "donor0", D_B = "donor1", D_C = "donor2", D_D = "donor3")

add_pretty <- function(obj, map, new_col = "donor_global_pretty") {
  x <- unname(map[as.character(obj$donor_global)])
  obj[[new_col]] <- factor(x, levels = unname(map[c("D_A", "D_B", "D_C", "D_D")]))
  obj
}

for (sid in selected_samples) {
  objs[[sid]] <- add_pretty(objs[[sid]], pretty_map)
}

# ------------------------------------------------------------------------------
# 6) Add a simplified donor column for downstream visualization
# ------------------------------------------------------------------------------
# The requested display names are:
#   donor0 -> Astrocyte1
#   donor1 -> Astrocyte2
#   donor2 -> Astrocyte3
#   donor3 -> Astrocyte3
#
# If you want a different mapping for donor2/donor3, update the vector below.
donor_display_map <- c(
  "donor0" = "Astrocyte1",
  "donor1" = "Astrocyte2",
  "donor2" = "Astrocyte3",
  "donor3" = "Astrocyte3"
)

for (sid in selected_samples) {
  obj <- objs[[sid]]
  obj$Donor <- factor(
    donor_display_map[as.character(obj$donor_global_pretty)],
    levels = c("Astrocyte1", "Astrocyte2", "Astrocyte3")
  )
  objs[[sid]] <- obj
}

# ------------------------------------------------------------------------------
# 7) Sanity checks
# ------------------------------------------------------------------------------
for (sid in selected_samples) {
  message("\nDonor sanity table for ", sid)
  print(table(objs[[sid]]$donor, objs[[sid]]$donor_global, useNA = "ifany"))
}

message("\nDonor harmonization completed successfully.")
message("The 'Donor' column now contains simplified labels for visualization.")
