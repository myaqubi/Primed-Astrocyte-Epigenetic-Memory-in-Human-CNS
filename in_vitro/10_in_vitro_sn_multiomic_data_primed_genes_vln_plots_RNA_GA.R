# ============================================================================== 
# 405-gene prime program UCell violin plots for RNA and GeneActivity
# ============================================================================== 
#
# Purpose:
#   This script adds a final 405-gene primed-profile score to the astrocyte
#   object and uses it to generate two violin plots for RNA
#   and GeneActivity. .
#
# Analysis:
#   - The 405 genes passing the primed screen are read from the
#     output table generated in the primed-gene identification script.
#   - UCell module scores are computed independently for RNA and GeneActivity.
#   - Mixed-effects models are fit for condition and donor effects.
#   - Violin plots display the distribution of single-cell scores with donor-
#     level means overlaid.
#
# Notes:
#   - Replace the generic placeholder paths below with your own project paths.
# ============================================================================== 

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(dplyr)
  library(readr)
  library(lme4)
  library(emmeans)
  library(ggplot2)
  library(ggpubr)
})

# ------------------------------------------------------------------------------
# 0) Configuration
# ------------------------------------------------------------------------------
project_dir <- "/path/to/your/project"
output_dir  <- file.path(project_dir, "Epigenetic_primed_anchor_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_object_path <- file.path(project_dir, "In_vitro_prenatal_samples_QC_filtered_integrated_annotated_astrocytes.rds")
primed_file <- file.path(project_dir, "finding_primed_genes", "TRUE_genes", "primed_genes_all_cells.csv")

if (!file.exists(input_object_path)) {
  stop("Input object not found at: ", input_object_path)
}

if (!file.exists(primed_file)) {
  stop("primed gene table not found at: ", primed_file)
}

obj <- readRDS(input_object_path)

# Standardize condition labels to Ctrl / Stim / Prime
condition_col <- "condition"
donor_col <- "Donor"
condition_levels <- c("Ctrl", "Stim", "Prime")

if (condition_col %in% colnames(obj@meta.data)) {
  obj@meta.data[[condition_col]] <- as.character(obj@meta.data[[condition_col]])
  obj@meta.data[[condition_col]][obj@meta.data[[condition_col]] %in% c("CTRL", "Control_no_hit")] <- "Ctrl"
  obj@meta.data[[condition_col]][obj@meta.data[[condition_col]] %in% c("Single_hit_D7", "Stim")] <- "Stim"
  obj@meta.data[[condition_col]][obj@meta.data[[condition_col]] %in% c("Single_hit_D1", "Prime")] <- "Prime"
  obj@meta.data[[condition_col]] <- factor(obj@meta.data[[condition_col]], levels = condition_levels)
}

# Standardize donor metadata to Donor
if (!donor_col %in% colnames(obj@meta.data)) {
  if ("donor_global_pretty" %in% colnames(obj@meta.data)) {
    obj@meta.data[[donor_col]] <- obj@meta.data[["donor_global_pretty"]]
  } else {
    stop("No donor column found. Expected 'Donor' or 'donor_global_pretty'.")
  }
}

obj@meta.data[[donor_col]] <- as.character(obj@meta.data[[donor_col]])
obj@meta.data$condition_pretty <- obj@meta.data[[condition_col]]
obj@meta.data$converted_donor_ID <- obj@meta.data[[donor_col]]

# ------------------------------------------------------------------------------
# 1) Read the primed gene list
# ------------------------------------------------------------------------------
primed_df <- read_csv(primed_file)

primed_genes <- primed_df %>%
  filter(PASS_primed_MODIFIED == TRUE) %>%
  pull(gene)

message("Number of primed genes loaded: ", length(primed_genes))

# ------------------------------------------------------------------------------
# 2) Add UCell module scores for RNA and GeneActivity
# ------------------------------------------------------------------------------
DefaultAssay(obj) <- "RNA"
obj <- AddModuleScore_UCell(
  obj,
  features = list(Final_405_primed_genes_RNA = primed_genes)
)

DefaultAssay(obj) <- "GeneActivity"
 ga_genes <- rownames(obj[["GeneActivity"]])
primed_genes_GA <- intersect(primed_genes, ga_genes)

message("Number of intersecting primed genes available in GeneActivity: ", length(primed_genes_GA))

obj <- AddModuleScore_UCell(
  obj,
  features = list(Final_405_primed_genes_GA = primed_genes_GA)
)

# ------------------------------------------------------------------------------
# 3) Plotting settings
# ------------------------------------------------------------------------------
condition_colors <- c(
  "Ctrl"  = "orange",
  "Prime" = "black",
  "Stim"  = "purple"
)

condition_colors <- condition_colors[condition_levels]

# Use donor colors that are robust to the donor labels present in the object
all_donors <- sort(unique(obj@meta.data[[donor_col]][!is.na(obj@meta.data[[donor_col]])]))
donor_colors <- c(
  "Donor1" = "red",
  "Donor2" = "green",
  "Donor3" = "blue"
)

missing_donors <- setdiff(all_donors, names(donor_colors))
if (length(missing_donors) > 0) {
  extra_colors <- grDevices::hcl.colors(length(missing_donors), palette = "Set2")
  names(extra_colors) <- missing_donors
  donor_colors <- c(donor_colors, extra_colors)
}

extra_theme <- if (exists("my_theme")) my_theme else theme()

# ------------------------------------------------------------------------------
# 4) Helper function for violin plots
# ------------------------------------------------------------------------------
make_lmer_vln_plot <- function(feature, plot_title, y_label, comparison_mode = c("GA", "RNA")) {
  comparison_mode <- match.arg(comparison_mode)
  
  meta <- obj@meta.data %>%
    mutate(
      condition_pretty = factor(condition_pretty, levels = condition_levels),
      converted_donor_ID = factor(converted_donor_ID, levels = sort(unique(converted_donor_ID)))
    ) %>%
    filter(
      !is.na(.data[[feature]]),
      !is.na(condition_pretty),
      !is.na(.data[[donor_col]]),
      !is.na(converted_donor_ID)
    )
  
  fit <- lmer(
    as.formula(paste0("`", feature, "` ~ condition_pretty + (1 | converted_donor_ID)")),
    data = meta
  )
  
  print(summary(fit))
  print(anova(fit))
  
  donor_means <- meta %>%
    group_by(condition_pretty, converted_donor_ID) %>%
    summarise(
      donor_mean = mean(.data[[feature]], na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    )
  
  emm <- emmeans(fit, ~ condition_pretty)
  
  if (comparison_mode == "GA") {
    ctrl_contrasts <- contrast(
      emm,
      method = list(
        "Stim_vs_Ctrl"  = c(-1, 1, 0),
        "Prime_vs_Ctrl" = c(-1, 0, 1)
      )
    )
  } else {
    ctrl_contrasts <- contrast(
      emm,
      method = list(
        "Stim_vs_Ctrl"  = c(-1, 1, 0),
        "Stim_vs_Prime" = c(0, 1, -1)
      )
    )
  }
  
  contrast_df <- as.data.frame(summary(ctrl_contrasts, adjust = "holm"))
  
  y_range <- range(meta[[feature]], na.rm = TRUE)
  y_step <- 0.08 * diff(y_range)
  if (y_step == 0) y_step <- 0.1
  
  p_annot <- contrast_df %>%
    mutate(
      group1 = case_when(
        contrast == "Stim_vs_Ctrl"  ~ "Ctrl",
        contrast == "Prime_vs_Ctrl" ~ "Ctrl",
        contrast == "Stim_vs_Prime" ~ "Stim"
      ),
      group2 = case_when(
        contrast == "Stim_vs_Ctrl"  ~ "Stim",
        contrast == "Prime_vs_Ctrl" ~ "Prime",
        contrast == "Stim_vs_Prime" ~ "Prime"
      ),
      p.signif = case_when(
        p.value <= 0.0001 ~ "****",
        p.value <= 0.001  ~ "***",
        p.value <= 0.01   ~ "**",
        p.value <= 0.05   ~ "*",
        TRUE            ~ "ns"
      ),
      y.position = max(meta[[feature]], na.rm = TRUE) + row_number() * y_step
    )
  
  p <- ggplot(
    meta,
    aes(
      x = condition_pretty,
      y = .data[[feature]],
      fill = condition_pretty
    )
  ) +
    geom_violin(scale = "width", trim = FALSE, color = "grey30", alpha = 0.8) +
    geom_boxplot(
      data = donor_means,
      aes(x = condition_pretty, y = donor_mean),
      inherit.aes = FALSE,
      width = 0.16,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      alpha = 0.85
    ) +
    geom_point(
      data = donor_means,
      aes(x = condition_pretty, y = donor_mean, color = converted_donor_ID),
      inherit.aes = FALSE,
      position = position_jitter(width = 0.08, height = 0, seed = 1),
      size = 2,
      shape = 16
    ) +
    stat_pvalue_manual(
      p_annot,
      label = "p.signif",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.5,
      size = 6
    ) +
    scale_fill_manual(values = condition_colors) +
    scale_color_manual(values = donor_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    ggtitle(plot_title) +
    xlab(NULL) +
    ylab(y_label) +
    theme_classic(base_size = 14) +
    extra_theme +
    theme(legend.position = "none")
  
  return(list(
    plot = p,
    fit = fit,
    anova = anova(fit),
    contrasts = contrast_df,
    p_annot = p_annot,
    donor_means = donor_means
  ))
}

# ------------------------------------------------------------------------------
# 5) Generate the RNA and GeneActivity plots
# ------------------------------------------------------------------------------
GA_results <- make_lmer_vln_plot(
  feature = "Final_405_primed_genes_GA_UCell",
  plot_title = "Accessibility_GA",
  y_label = "Final 405 primed Genes GA UCell",
  comparison_mode = "GA"
)

p_GA_405 <- GA_results$plot
print(p_GA_405)

RNA_results <- make_lmer_vln_plot(
  feature = "Final_405_primed_genes_RNA_UCell",
  plot_title = "Final 405 primed Genes (RNA)",
  y_label = "Final 405 primed Genes RNA UCell",
  comparison_mode = "RNA"
)

p_RNA_405 <- RNA_results$plot
print(p_RNA_405)

# ------------------------------------------------------------------------------
# 6) Save plots
# ------------------------------------------------------------------------------
ggsave(file.path(output_dir, "Final_405_primed_genes_GA_UCell_vln_plot.png"), p_GA_405, width = 3.7, height = 3)
ggsave(file.path(output_dir, "Final_405_primed_genes_RNA_UCell_vln_plot.png"), p_RNA_405, width = 3.7, height = 3)
