# ============================================================================== 
# Primed gene identification and bootstrap stability analysis
# ============================================================================== 
#
# Purpose:
#   This script contains identification of in-vitro primed genes workflow from the integrated
#   object using RNA and GeneActivity data. The script includes:
#   1) differential testing with limma using condition + donor as covariates,
#   2) bootstrap stability analysis across 1000 iterations,
#   3) export of summary tables and iteration-level results.
#
#   Primed gene criteria identification used in this script:
#   RNA: Stim > Ctrl and Stim > Prime.
#   Gene Activity: Stim > Ctrl and Prime > Ctrl.
#
# Notes:
#   - Replace the generic placeholder paths below with your own project paths.
# ============================================================================== 

suppressPackageStartupMessages({
  library(Seurat)
  library(limma)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
  library(openxlsx)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# 0) Configuration
# ------------------------------------------------------------------------------
project_dir <- "/path/to/your/project"
output_dir  <- file.path(project_dir, "finding_primed_genes", "TRUE_genes")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_object_path <- file.path(project_dir, "In_vitro_prenatal_samples_QC_filtered_integrated_annotated_astrocytes.rds")

if (!file.exists(input_object_path)) {
  stop("Input object not found at: ", input_object_path)
}

obj <- readRDS(input_object_path)

condition_col <- "condition"
donor_col     <- "Donor"
cond_levels   <- c("Ctrl", "Stim", "Prime")

q_cut    <- 0.20
alpha_ns <- 0.05

excel_file <- file.path(output_dir, "Astrocytes_primed_score_all_genes.xlsx")
true_txt   <- file.path(output_dir, "TRUE_genes.txt")

# Ensure condition order is consistent with the analysis
obj@meta.data[[condition_col]] <- factor(obj@meta.data[[condition_col]], levels = cond_levels)

# If the object uses donor_global_pretty instead of Donor, create a portable alias
if (!donor_col %in% colnames(obj@meta.data)) {
  if ("donor_global_pretty" %in% colnames(obj@meta.data)) {
    obj@meta.data[[donor_col]] <- obj@meta.data[["donor_global_pretty"]]
  } else {
    stop("No donor column found. Expected 'Donor'.")
  }
}

# ------------------------------------------------------------------------------
# 1) Helper functions
# ------------------------------------------------------------------------------
clean_filename <- function(x) {
  x <- gsub("[/\\?<>\\:*|\"']", "_", x)
  x <- gsub("\\s+", "_", x)
  x
}

sig_q_to_label <- function(q, est, q_cut = 0.20) {
  if (is.na(q) || is.na(est)) return("NA")
  if (!(est > 0) || q >= q_cut) return("ns")
  if (q < 1e-4) return("****")
  if (q < 1e-3) return("***")
  if (q < 1e-2) return("**")
  if (q < 5e-2) return("*")
  return("*")
}

ns_p_to_label <- function(p, alpha_ns = 0.05) {
  if (is.na(p)) return("NA")
  if (p >= alpha_ns) return("ns")
  if (p < 1e-4) return("****")
  if (p < 1e-3) return("***")
  if (p < 1e-2) return("**")
  if (p < 5e-2) return("*")
  return("*")
}

get_donor_means_one_gene <- function(obj, assay, gene, condition_col, donor_col, cond_levels, slot = "data") {
  DefaultAssay(obj) <- assay
  expr <- FetchData(obj, vars = gene, slot = slot)
  
  df <- data.frame(
    condition = obj@meta.data[[condition_col]],
    donor     = obj@meta.data[[donor_col]],
    value     = expr[[gene]],
    stringsAsFactors = FALSE
  )
  
  df <- df[
    !is.na(df$condition) &
      !is.na(df$donor) &
      df$condition %in% cond_levels,
    , drop = FALSE
  ]
  
  df$condition <- factor(df$condition, levels = cond_levels)
  df$donor     <- factor(df$donor)
  
  dm <- aggregate(value ~ condition + donor, data = df, FUN = mean)
  dm$gene  <- gene
  dm$assay <- assay
  dm
}

# ------------------------------------------------------------------------------
# 2) limma function
#    expression ~ condition + donor
# ------------------------------------------------------------------------------
run_limma_contrasts <- function(expr_mat, meta_df) {
  design <- model.matrix(~ 0 + condition + donor, data = meta_df)
  colnames(design) <- sub("^condition", "", colnames(design))
  colnames(design) <- sub("^donor", "donor_", colnames(design))
  
  contrast_mat <- makeContrasts(
    Stim_vs_Ctrl = Stim - Ctrl,
    Stim_vs_Prime = Stim - Prime,
    Prime_vs_Ctrl = Prime - Ctrl,
    levels = design
  )
  
  fit <- lmFit(expr_mat, design)
  fit <- contrasts.fit(fit, contrast_mat)
  fit <- eBayes(fit)
  
  extract <- function(coef_name) {
    tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
    data.frame(
      gene = rownames(tt),
      est  = tt$logFC,
      p    = tt$P.Value,
      q    = tt$adj.P.Val,
      stringsAsFactors = FALSE
    )
  }
  
  t1 <- extract("Stim_vs_Ctrl")
  names(t1)[2:4] <- c("est_Stim_vs_Ctrl", "p_Stim_vs_Ctrl", "q_Stim_vs_Ctrl")
  
  t2 <- extract("Stim_vs_Prime")
  names(t2)[2:4] <- c("est_Stim_vs_Prime", "p_Stim_vs_Prime", "q_Stim_vs_Prime")
  
  t3 <- extract("Prime_vs_Ctrl")
  names(t3)[2:4] <- c("est_Prime_vs_Ctrl", "p_Prime_vs_Ctrl", "q_Prime_vs_Ctrl")
  
  Reduce(function(x, y) merge(x, y, by = "gene"), list(t1, t2, t3))
}

# ------------------------------------------------------------------------------
# 3) Pass rules
# ------------------------------------------------------------------------------
apply_pass_rule_original <- function(df, q_cut, alpha_ns) {
  sig_part <- (
    (df$est_Stim_vs_Ctrl__RNA > 0 & df$q_Stim_vs_Ctrl__RNA < q_cut) &
      (df$est_Stim_vs_Prime__RNA > 0 & df$q_Stim_vs_Prime__RNA < q_cut) &
      (df$est_Stim_vs_Ctrl__GA > 0 & df$q_Stim_vs_Ctrl__GA < q_cut) &
      (df$est_Prime_vs_Ctrl__GA > 0 & df$q_Prime_vs_Ctrl__GA < q_cut)
  )
  
  ns_part <- (
    df$p_Stim_vs_Prime__GA >= alpha_ns
  )
  
  df$PASS_primed_ORIGINAL <- sig_part & ns_part
  df$RNA_sig_Stim_vs_Ctrl <- df$est_Stim_vs_Ctrl__RNA > 0 & df$q_Stim_vs_Ctrl__RNA < q_cut
  df$RNA_sig_Stim_vs_Prime <- df$est_Stim_vs_Prime__RNA > 0 & df$q_Stim_vs_Prime__RNA < q_cut
  df$GA_sig_Stim_vs_Ctrl <- df$est_Stim_vs_Ctrl__GA > 0 & df$q_Stim_vs_Ctrl__GA < q_cut
  df$GA_sig_Prime_vs_Ctrl <- df$est_Prime_vs_Ctrl__GA > 0 & df$q_Prime_vs_Ctrl__GA < q_cut
  df$GA_ns_Stim_vs_Prime <- df$p_Stim_vs_Prime__GA >= alpha_ns
  df
}

apply_pass_rule_modified <- function(df, q_cut) {
  df$PASS_RNA_MODIFIED <- (
    df$est_Stim_vs_Ctrl__RNA > 0 & df$q_Stim_vs_Ctrl__RNA < q_cut &
      df$est_Stim_vs_Prime__RNA > 0 & df$q_Stim_vs_Prime__RNA < q_cut
  )
  
  df$PASS_GA_MODIFIED <- (
    df$est_Stim_vs_Ctrl__GA > 0 & df$q_Stim_vs_Ctrl__GA < q_cut &
      df$est_Prime_vs_Ctrl__GA > 0 & df$q_Prime_vs_Ctrl__GA < q_cut
  )
  
  df$PASS_primed_MODIFIED <- df$PASS_RNA_MODIFIED & df$PASS_GA_MODIFIED
  df
}

# ------------------------------------------------------------------------------
# 4) Prepare metadata and matrices
# ------------------------------------------------------------------------------
genes_rna <- rownames(GetAssayData(obj, assay = "RNA", slot = "data"))
genes_ga  <- rownames(GetAssayData(obj, assay = "GeneActivity", slot = "data"))

valid_genes <- intersect(genes_rna, genes_ga)
cat("Genes present in both RNA and GeneActivity assays:", length(valid_genes), "\n")

meta <- obj@meta.data

keep_cells <- rownames(meta)[
  !is.na(meta[[condition_col]]) &
    !is.na(meta[[donor_col]]) &
    meta[[condition_col]] %in% cond_levels
]

meta_sub <- meta[keep_cells, , drop = FALSE]
meta_sub[[condition_col]] <- factor(meta_sub[[condition_col]], levels = cond_levels)
meta_sub[[donor_col]]     <- factor(meta_sub[[donor_col]])

cat("Cells retained for analysis:", nrow(meta_sub), "\n")
cat("Condition counts:\n")
print(table(meta_sub[[condition_col]]))
cat("Donor counts:\n")
print(table(meta_sub[[donor_col]]))

rna_mat <- as.matrix(GetAssayData(obj, assay = "RNA", slot = "data")[valid_genes, keep_cells])
ga_mat  <- as.matrix(GetAssayData(obj, assay = "GeneActivity", slot = "data")[valid_genes, keep_cells])

# ------------------------------------------------------------------------------
# 5) Run limma once
# ------------------------------------------------------------------------------
cat("\nRunning limma: expression ~ condition + donor\n")

rna_B <- run_limma_contrasts(rna_mat, meta_sub)
ga_B  <- run_limma_contrasts(ga_mat, meta_sub)

res_B <- merge(rna_B, ga_B, by = "gene", suffixes = c("__RNA", "__GA"))

# ------------------------------------------------------------------------------
# 6) Modified rule
# ------------------------------------------------------------------------------
res_B <- apply_pass_rule_modified(res_B, q_cut = q_cut)

n_rna_pass <- sum(res_B$PASS_RNA_MODIFIED, na.rm = TRUE)
n_ga_pass_among_rna <- sum(res_B$PASS_RNA_MODIFIED & res_B$PASS_GA_MODIFIED, na.rm = TRUE)

cat("\n========================================\n")
cat("MODIFIED RULE\n")
cat("========================================\n")
cat("Genes passing RNA criteria only:", n_rna_pass, "\n")
cat("Among RNA-pass genes, genes also passing GA criteria:", n_ga_pass_among_rna, "\n")
cat("Final TRUE primed genes:", sum(res_B$PASS_primed_MODIFIED, na.rm = TRUE), "\n")
cat("\nModified TRUE/FALSE table:\n")
print(table(res_B$PASS_primed_MODIFIED, useNA = "ifany"))

# Save the full result table
write.csv(res_B, file = file.path(output_dir, "primed_genes_all_cells.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 7) Bootstrap stability of modified primed genes
# ------------------------------------------------------------------------------
n_iter <- 1000
set.seed(123)

primed_genes_modified <- res_B$gene[res_B$PASS_primed_MODIFIED]
cat("\nNumber of modified primed genes:", length(primed_genes_modified), "\n")

sampling_table <- meta_sub %>%
  dplyr::count(
    condition = .data[[condition_col]],
    donor = .data[[donor_col]],
    name = "n_cells"
  )

print(sampling_table)

make_bootstrap_cells <- function(meta, sampling_table, condition_col, donor_col) {
  sampled_cells <- c()
  for (i in seq_len(nrow(sampling_table))) {
    this_condition <- as.character(sampling_table$condition[i])
    this_donor <- as.character(sampling_table$donor[i])
    this_n <- sampling_table$n_cells[i]
    candidate_cells <- rownames(meta)[
      as.character(meta[[condition_col]]) == this_condition &
        as.character(meta[[donor_col]]) == this_donor
    ]
    sampled_cells <- c(
      sampled_cells,
      sample(candidate_cells, size = this_n, replace = TRUE)
    )
  }
  sampled_cells
}

all_pass_lists <- list()
iter_summary <- list()

for (iter in seq_len(n_iter)) {
  if (iter %% 50 == 0) {
    cat("\nBootstrap iteration:", iter, "of", n_iter, "\n")
  }
  
  sampled_cells <- make_bootstrap_cells(
    meta = meta_sub,
    sampling_table = sampling_table,
    condition_col = condition_col,
    donor_col = donor_col
  )
  
  boot_meta <- meta_sub[sampled_cells, , drop = FALSE]
  boot_meta[[condition_col]] <- factor(boot_meta[[condition_col]], levels = cond_levels)
  boot_meta[[donor_col]] <- factor(boot_meta[[donor_col]])
  
  boot_rna <- rna_mat[, sampled_cells, drop = FALSE]
  boot_ga  <- ga_mat[, sampled_cells, drop = FALSE]
  
  rna_res <- run_limma_contrasts(boot_rna, boot_meta)
  ga_res  <- run_limma_contrasts(boot_ga, boot_meta)
  
  boot_res <- merge(rna_res, ga_res, by = "gene", suffixes = c("__RNA", "__GA"))
  boot_res <- apply_pass_rule_modified(boot_res, q_cut = q_cut)
  
  pass_genes <- boot_res$gene[boot_res$PASS_primed_MODIFIED]
  all_pass_lists[[iter]] <- data.frame(iteration = iter, gene = pass_genes)
  iter_summary[[iter]] <- data.frame(iteration = iter, n_primed_genes = length(pass_genes))
}

pass_df <- dplyr::bind_rows(all_pass_lists)
iter_summary_df <- dplyr::bind_rows(iter_summary)

# ------------------------------------------------------------------------------
# 8) Bootstrap summary
# ------------------------------------------------------------------------------
bootstrap_summary <- data.frame(gene = primed_genes_modified, stringsAsFactors = FALSE)

freq_tbl <- pass_df %>% dplyr::count(gene, name = "n_pass")
bootstrap_summary <- bootstrap_summary %>% dplyr::left_join(freq_tbl, by = "gene")
bootstrap_summary$n_pass[is.na(bootstrap_summary$n_pass)] <- 0
bootstrap_summary$pass_frequency <- bootstrap_summary$n_pass / n_iter
bootstrap_summary$pct_pass <- 100 * bootstrap_summary$pass_frequency

write.csv(bootstrap_summary, file = file.path(output_dir, "primed_genes_all_cells_bootstrap_iterations_resul.csv"), row.names = FALSE)

bootstrap_summary$stability_class <- dplyr::case_when(
  bootstrap_summary$pct_pass >= 75 ~ "Stable",
  bootstrap_summary$pct_pass >= 50 ~ "Moderately stable",
  bootstrap_summary$pct_pass >= 25 ~ "Relatively unstable",
  TRUE ~ "Unstable"
)

cat("\n========================================\n")
cat("BOOTSTRAP STABILITY SUMMARY\n")
cat("========================================\n")
print(table(bootstrap_summary$stability_class))
cat("\nPercentages:\n")
print(round(100 * prop.table(table(bootstrap_summary$stability_class)), 2))

write.csv(iter_summary_df, file = file.path(output_dir, "primed_genes_all_cells_bootstrap_iterations_result.csv"), row.names = FALSE)

summary(iter_summary_df$n_primed_genes)

