# 07_donor_demultiplexing_explicit.R
# Donor demultiplexing and donor identity assignment for multiplexed multiome samples
#
# This script:
#   1. Loads the post-QC multiplexed Seurat object.
#   2. Adds demultiplexing calls from matched barcode files.
#   3. Examines sex-associated RNA markers separately for donor0 and donor1.
#   4. Assigns each donor0/donor1 group to the corresponding biological donor.
#   5. Adds Group, Donor, Sample, and Sex metadata back to the master object.
#
# Unassigned cells are retained in the master object but do not receive
# donor-specific metadata.

library(Seurat)

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------

base_dir <- "PATH/TO/MULTIOME_DATA"

input_rds <- file.path(
  base_dir,
  "multiome_object_postQC.rds"
)

output_rds <- file.path(
  base_dir,
  "multiome_object_demultiplexed.rds"
)

# -------------------------------------------------------------------------
# Load object
# -------------------------------------------------------------------------

multiome_obj <- readRDS(input_rds)

# Demultiplexing is performed within each multiplexed sample.
Idents(multiome_obj) <- "orig.ident"

# Sex-associated genes used to distinguish male/female donor assignments.
sex_genes <- c("XIST", "DDX3Y", "SRY", "ZFY", "USP9Y")

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

# Read the matched demultiplexing results and assign the best singlet call.
add_demux_calls <- function(sample_obj, matched_file) {

  donor_data <- read.csv(
    matched_file,
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!"best_singlet" %in% colnames(donor_data)) {
    stop("Column 'best_singlet' not found in: ", matched_file)
  }

  # Use exact cell-name matching where possible.
  if (all(Cells(sample_obj) %in% rownames(donor_data))) {
    sample_obj$donor <- donor_data[Cells(sample_obj), "best_singlet"]
  } else {
    # Retain the row-order behavior of the original analysis if exact matching
    # is unavailable, but require equal numbers of cells/rows.
    if (nrow(donor_data) != ncol(sample_obj)) {
      stop(
        "Could not align donor calls for ",
        basename(matched_file),
        ": ",
        nrow(donor_data),
        " rows for ",
        ncol(sample_obj),
        " cells."
      )
    }

    warning(
      "Exact barcode matching unavailable for ",
      basename(matched_file),
      "; assigning by row order as in the original analysis."
    )

    sample_obj$donor <- donor_data$best_singlet
  }

  sample_obj
}

# Print the number of nuclei expressing each sex-associated marker.
check_sex_markers <- function(donor_obj, donor_label, sample_label) {

  available_genes <- intersect(
    sex_genes,
    rownames(donor_obj[["RNA"]])
  )

  if (length(available_genes) == 0) {
    warning("No sex-associated genes found for ", sample_label, " ", donor_label)
    return(invisible(NULL))
  }

  sex_genes_expression <- FetchData(
    donor_obj,
    vars = available_genes
  )

  num_cells_expressing_sex_genes <- colSums(
    sex_genes_expression > 0
  )

  cat(
    "\n",
    sample_label,
    " - ",
    donor_label,
    " (n = ",
    ncol(donor_obj),
    " nuclei)\n",
    sep = ""
  )

  print(num_cells_expressing_sex_genes)

  invisible(num_cells_expressing_sex_genes)
}

# -------------------------------------------------------------------------
# Sample 1
# -------------------------------------------------------------------------
# Multiplexed pool: sample1_CTRL
# Final donor assignment:
#   donor0 -> C1670, Female, CTRL
#   donor1 -> C1143, Male, CTRL

Sample1 <- subset(
  multiome_obj,
  ident = "sample1_CTRL"
)

Sample1 <- add_demux_calls(
  Sample1,
  file.path(
    base_dir,
    "S1_demultiplexing",
    "matched_barcodes_Sample1.csv"
  )
)

Idents(Sample1) <- "donor"

Donor1_S1 <- subset(Sample1, idents = "donor0")
Donor2_S1 <- subset(Sample1, idents = "donor1")

check_sex_markers(Donor1_S1, "donor0", "S1")
check_sex_markers(Donor2_S1, "donor1", "S1")

# Based on the sex-marker profile and known sample information:
# donor0 was assigned to female CTRL donor C1670.
# donor1 was assigned to male CTRL donor C1143.

Donor1_S1[["Group"]] <- "CTRL"
Donor1_S1[["Donor"]] <- "C1670"
Donor1_S1[["Sample"]] <- "S1"
Donor1_S1[["Sex"]] <- "Female"

Donor2_S1[["Group"]] <- "CTRL"
Donor2_S1[["Donor"]] <- "C1143"
Donor2_S1[["Sample"]] <- "S1"
Donor2_S1[["Sex"]] <- "Male"


# -------------------------------------------------------------------------
# Sample 2
# -------------------------------------------------------------------------
# Multiplexed pool: sample2_MS
# Final donor assignment:
#   donor0 -> MS1002, Male, MS
#   donor1 -> MS1736, Female, MS

Sample2 <- subset(
  multiome_obj,
  ident = "sample2_MS"
)

Sample2 <- add_demux_calls(
  Sample2,
  file.path(
    base_dir,
    "S2_demultiplexing",
    "matched_barcodes_Sample2.csv"
  )
)

Idents(Sample2) <- "donor"

Donor1_S2 <- subset(Sample2, idents = "donor0")
Donor2_S2 <- subset(Sample2, idents = "donor1")

check_sex_markers(Donor1_S2, "donor0", "S2")
check_sex_markers(Donor2_S2, "donor1", "S2")

# Based on the sex-marker profile and known sample information:
# donor0 was assigned to male MS donor MS1002.
# donor1 was assigned to female MS donor MS1736.

Donor1_S2[["Group"]] <- "MS"
Donor1_S2[["Donor"]] <- "MS1002"
Donor1_S2[["Sample"]] <- "S2"
Donor1_S2[["Sex"]] <- "Male"

Donor2_S2[["Group"]] <- "MS"
Donor2_S2[["Donor"]] <- "MS1736"
Donor2_S2[["Sample"]] <- "S2"
Donor2_S2[["Sex"]] <- "Female"


# -------------------------------------------------------------------------
# Sample 3
# -------------------------------------------------------------------------
# Multiplexed pool: sample3_CTRL
# Final donor assignment:
#   donor0 -> C1670, Female, CTRL
#   donor1 -> C1143, Male, CTRL

Sample3 <- subset(
  multiome_obj,
  ident = "sample3_CTRL"
)

Sample3 <- add_demux_calls(
  Sample3,
  file.path(
    base_dir,
    "S3_demultiplexing",
    "matched_barcodes_Sample3.csv"
  )
)

Idents(Sample3) <- "donor"

Donor1_S3 <- subset(Sample3, idents = "donor0")
Donor2_S3 <- subset(Sample3, idents = "donor1")

check_sex_markers(Donor1_S3, "donor0", "S3")
check_sex_markers(Donor2_S3, "donor1", "S3")

# Based on the sex-marker profile and known sample information:
# donor0 was assigned to female CTRL donor C1670.
# donor1 was assigned to male CTRL donor C1143.

Donor1_S3[["Group"]] <- "CTRL"
Donor1_S3[["Donor"]] <- "C1670"
Donor1_S3[["Sample"]] <- "S3"
Donor1_S3[["Sex"]] <- "Female"

Donor2_S3[["Group"]] <- "CTRL"
Donor2_S3[["Donor"]] <- "C1143"
Donor2_S3[["Sample"]] <- "S3"
Donor2_S3[["Sex"]] <- "Male"


# -------------------------------------------------------------------------
# Sample 4
# -------------------------------------------------------------------------
# Multiplexed pool: sample4_MS
# Final donor assignment:
#   donor0 -> MS1002, Male, MS
#   donor1 -> MS1736, Female, MS

Sample4 <- subset(
  multiome_obj,
  ident = "sample4_MS"
)

Sample4 <- add_demux_calls(
  Sample4,
  file.path(
    base_dir,
    "S4_demultiplexing",
    "matched_barcodes_Sample4.csv"
  )
)

Idents(Sample4) <- "donor"

Donor1_S4 <- subset(Sample4, idents = "donor0")
Donor2_S4 <- subset(Sample4, idents = "donor1")

check_sex_markers(Donor1_S4, "donor0", "S4")
check_sex_markers(Donor2_S4, "donor1", "S4")

# Based on the sex-marker profile and known sample information:
# donor0 was assigned to male MS donor MS1002.
# donor1 was assigned to female MS donor MS1736.

Donor1_S4[["Group"]] <- "MS"
Donor1_S4[["Donor"]] <- "MS1002"
Donor1_S4[["Sample"]] <- "S4"
Donor1_S4[["Sex"]] <- "Male"

Donor2_S4[["Group"]] <- "MS"
Donor2_S4[["Donor"]] <- "MS1736"
Donor2_S4[["Sample"]] <- "S4"
Donor2_S4[["Sex"]] <- "Female"


# -------------------------------------------------------------------------
# Sample 7
# -------------------------------------------------------------------------
# Multiplexed pool: sample7_CTRL
# Final donor assignment:
#   donor0 -> C1722, Male, CTRL
#   donor1 -> C1888, Female, CTRL

Sample7 <- subset(
  multiome_obj,
  ident = "sample7_CTRL"
)

Sample7 <- add_demux_calls(
  Sample7,
  file.path(
    base_dir,
    "S7_demultiplexing",
    "matched_barcodes_Sample7.csv"
  )
)

Idents(Sample7) <- "donor"

Donor1_S7 <- subset(Sample7, idents = "donor0")
Donor2_S7 <- subset(Sample7, idents = "donor1")

check_sex_markers(Donor1_S7, "donor0", "S7")
check_sex_markers(Donor2_S7, "donor1", "S7")

# Based on the sex-marker profile and known sample information:
# donor0 was assigned to male CTRL donor C1722.
# donor1 was assigned to female CTRL donor C1888.

Donor1_S7[["Group"]] <- "CTRL"
Donor1_S7[["Donor"]] <- "C1722"
Donor1_S7[["Sample"]] <- "S7"
Donor1_S7[["Sex"]] <- "Male"

Donor2_S7[["Group"]] <- "CTRL"
Donor2_S7[["Donor"]] <- "C1888"
Donor2_S7[["Sample"]] <- "S7"
Donor2_S7[["Sex"]] <- "Female"


# -------------------------------------------------------------------------
# Sample 8
# -------------------------------------------------------------------------
# Multiplexed pool: sample8_MS
# Final donor assignment:
#   donor0 -> MS1421, Female, MS
#   donor1 -> MS1011, Male, MS
#
# Note: the original narrative comment in the historical script was inconsistent
# for this sample. The assignments below reproduce the metadata that were
# actually written to Donor1_S8 and Donor2_S8.

Sample8 <- subset(
  multiome_obj,
  ident = "sample8_MS"
)

Sample8 <- add_demux_calls(
  Sample8,
  file.path(
    base_dir,
    "S8_demultiplexing",
    "matched_barcodes_Sample8.csv"
  )
)

Idents(Sample8) <- "donor"

Donor1_S8 <- subset(Sample8, idents = "donor0")
Donor2_S8 <- subset(Sample8, idents = "donor1")

check_sex_markers(Donor1_S8, "donor0", "S8")
check_sex_markers(Donor2_S8, "donor1", "S8")

# Based on the sex-marker profile and known sample information:
# donor0 was assigned to female MS donor MS1421.
# donor1 was assigned to male MS donor MS1011.

Donor1_S8[["Group"]] <- "MS"
Donor1_S8[["Donor"]] <- "MS1421"
Donor1_S8[["Sample"]] <- "S8"
Donor1_S8[["Sex"]] <- "Female"

Donor2_S8[["Group"]] <- "MS"
Donor2_S8[["Donor"]] <- "MS1011"
Donor2_S8[["Sample"]] <- "S8"
Donor2_S8[["Sex"]] <- "Male"


# -------------------------------------------------------------------------
# Sample 9
# -------------------------------------------------------------------------
# Multiplexed pool: sample9_MS
# Final donor assignment:
#   donor0 -> MS1736, Female, MS
#   donor1 -> MS1090, Male, MS

Sample9 <- subset(
  multiome_obj,
  ident = "sample9_MS"
)

Sample9 <- add_demux_calls(
  Sample9,
  file.path(
    base_dir,
    "S9_demultiplexing",
    "matched_barcodes_Sample9.csv"
  )
)

Idents(Sample9) <- "donor"

Donor1_S9 <- subset(Sample9, idents = "donor0")
Donor2_S9 <- subset(Sample9, idents = "donor1")

check_sex_markers(Donor1_S9, "donor0", "S9")
check_sex_markers(Donor2_S9, "donor1", "S9")

# Based on the sex-marker profile and known sample information:
# donor0 was assigned to female MS donor MS1736.
# donor1 was assigned to male MS donor MS1090.

Donor1_S9[["Group"]] <- "MS"
Donor1_S9[["Donor"]] <- "MS1736"
Donor1_S9[["Sample"]] <- "S9"
Donor1_S9[["Sex"]] <- "Female"

Donor2_S9[["Group"]] <- "MS"
Donor2_S9[["Donor"]] <- "MS1090"
Donor2_S9[["Sample"]] <- "S9"
Donor2_S9[["Sex"]] <- "Male"


# -------------------------------------------------------------------------
# Sample 11
# -------------------------------------------------------------------------
# Multiplexed pool: sample11_CTRL
# Final donor assignment:
#   donor0 -> C1670, Female, CTRL
#   donor1 -> C1143, Male, CTRL

Sample11 <- subset(
  multiome_obj,
  ident = "sample11_CTRL"
)

Sample11 <- add_demux_calls(
  Sample11,
  file.path(
    base_dir,
    "S11_demultiplexing",
    "matched_barcodes_Sample11.csv"
  )
)

Idents(Sample11) <- "donor"

Donor1_S11 <- subset(Sample11, idents = "donor0")
Donor2_S11 <- subset(Sample11, idents = "donor1")

check_sex_markers(Donor1_S11, "donor0", "S11")
check_sex_markers(Donor2_S11, "donor1", "S11")

# Based on the sex-marker profile and known sample information:
# donor0 was assigned to female CTRL donor C1670.
# donor1 was assigned to male CTRL donor C1143.

Donor1_S11[["Group"]] <- "CTRL"
Donor1_S11[["Donor"]] <- "C1670"
Donor1_S11[["Sample"]] <- "S11"
Donor1_S11[["Sex"]] <- "Female"

Donor2_S11[["Group"]] <- "CTRL"
Donor2_S11[["Donor"]] <- "C1143"
Donor2_S11[["Sample"]] <- "S11"
Donor2_S11[["Sex"]] <- "Male"


# -------------------------------------------------------------------------
# Sample 12
# -------------------------------------------------------------------------
# Multiplexed pool: sample12_MS
# Final donor assignment:
#   donor0 -> MS1736, Female, MS
#   donor1 -> MS1090, Male, MS

Sample12 <- subset(
  multiome_obj,
  ident = "sample12_MS"
)

Sample12 <- add_demux_calls(
  Sample12,
  file.path(
    base_dir,
    "S12_demultiplexing",
    "matched_barcodes_Sample12.csv"
  )
)

Idents(Sample12) <- "donor"

Donor1_S12 <- subset(Sample12, idents = "donor0")
Donor2_S12 <- subset(Sample12, idents = "donor1")

check_sex_markers(Donor1_S12, "donor0", "S12")
check_sex_markers(Donor2_S12, "donor1", "S12")

# Based on the sex-marker profile and known sample information:
# donor0 was assigned to female MS donor MS1736.
# donor1 was assigned to male MS donor MS1090.

Donor1_S12[["Group"]] <- "MS"
Donor1_S12[["Donor"]] <- "MS1736"
Donor1_S12[["Sample"]] <- "S12"
Donor1_S12[["Sex"]] <- "Female"

Donor2_S12[["Group"]] <- "MS"
Donor2_S12[["Donor"]] <- "MS1090"
Donor2_S12[["Sample"]] <- "S12"
Donor2_S12[["Sex"]] <- "Male"


# -------------------------------------------------------------------------
# Add demultiplexing calls to the master object
# -------------------------------------------------------------------------

sample_objects <- list(
  S1  = Sample1,
  S2  = Sample2,
  S3  = Sample3,
  S4  = Sample4,
  S7  = Sample7,
  S8  = Sample8,
  S9  = Sample9,
  S11 = Sample11,
  S12 = Sample12
)

if (!"demux_call" %in% colnames(multiome_obj@meta.data)) {
  multiome_obj$demux_call <- NA_character_
}

for (nm in names(sample_objects)) {

  obj <- sample_objects[[nm]]
  common <- intersect(
    Cells(multiome_obj),
    Cells(obj)
  )

  multiome_obj@meta.data[
    common,
    "demux_call"
  ] <- obj$donor[common]
}


# -------------------------------------------------------------------------
# Add donor metadata back to the master object
# -------------------------------------------------------------------------

get_four_cols <- function(
    obj,
    cols = c("Group", "Donor", "Sample", "Sex")
) {

  meta <- obj@meta.data

  miss <- setdiff(cols, colnames(meta))

  if (length(miss) > 0) {
    meta[, miss] <- NA
  }

  meta[, cols, drop = FALSE]
}


update_from_subset <- function(master, subset_obj) {

  subm <- get_four_cols(subset_obj)

  common <- intersect(
    rownames(master@meta.data),
    rownames(subm)
  )

  if (length(common) > 0) {
    master@meta.data[
      common,
      colnames(subm)
    ] <- subm[
      common,
      ,
      drop = FALSE
    ]
  }

  master
}


donor_objects <- list(
  Donor1_S1,
  Donor2_S1,
  Donor1_S2,
  Donor2_S2,
  Donor1_S3,
  Donor2_S3,
  Donor1_S4,
  Donor2_S4,
  Donor1_S7,
  Donor2_S7,
  Donor1_S8,
  Donor2_S8,
  Donor1_S9,
  Donor2_S9,
  Donor1_S11,
  Donor2_S11,
  Donor1_S12,
  Donor2_S12
)


for (donor_obj in donor_objects) {
  multiome_obj <- update_from_subset(
    multiome_obj,
    donor_obj
  )
}


# -------------------------------------------------------------------------
# Final donor labels used for downstream single-cell analyses
# -------------------------------------------------------------------------
#
# Biological donor IDs are converted to the simplified donor labels used
# throughout downstream single-cell analyses and figures:
#
#   C1143  -> CTRL1
#   C1670  -> CTRL2
#   C1722  -> CTRL3
#   C1888  -> CTRL4
#
#   MS1002 -> MS1
#   MS1736 -> MS2
#   MS1011 -> MS3
#   MS1421 -> MS4
#   MS1090 -> MS5
#
# Cells without a confident donor assignment remain NA.

multiome_obj$Donor_labels_for_singlocell <- dplyr::recode(
  as.character(multiome_obj$Donor),

  "C1143"  = "CTRL1",
  "C1670"  = "CTRL2",
  "C1722"  = "CTRL3",
  "C1888"  = "CTRL4",

  "MS1002" = "MS1",
  "MS1736" = "MS2",
  "MS1011" = "MS3",
  "MS1421" = "MS4",
  "MS1090" = "MS5",

  .default = NA_character_
)

cat("\nFinal donor ID -> downstream label mapping:\n")

print(
  unique(
    na.omit(
      multiome_obj@meta.data[
        ,
        c("Donor", "Donor_labels_for_singlocell")
      ]
    )
  )
)


# -------------------------------------------------------------------------
# QC summaries
# -------------------------------------------------------------------------

cat("\nDemultiplexing calls by multiplexed sample:\n")

print(
  table(
    multiome_obj$orig.ident,
    multiome_obj$demux_call,
    useNA = "ifany"
  )
)


cat("\nFinal donor annotation:\n")

print(
  table(
    multiome_obj$Sample,
    multiome_obj$Donor,
    useNA = "ifany"
  )
)


cat("\nDonor metadata combinations:\n")

print(
  unique(
    na.omit(
      multiome_obj@meta.data[
        ,
        c("Group", "Donor", "Sample", "Sex")
      ]
    )
  )
)


# -------------------------------------------------------------------------
# Save final demultiplexed object
# -------------------------------------------------------------------------

saveRDS(
  multiome_obj,
  file = output_rds
)

message(
  "\nSaved demultiplexed object to:\n",
  output_rds
)
