library(readxl)
library(ggplot2)
library(viridis)

input_file <- "Astrocyte_memory_genes_qPCR_data.xlsx"
output_dir <- "qPCR_memory_heatmap"
output_file <- file.path(
  output_dir,
  "qPCR_donor_panel_log_cividis_gene_donor_scaled_memory_stars.png"
)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

raw <- as.data.frame(read_excel(input_file, col_names = FALSE), stringsAsFactors = FALSE)

source_donor_labels <- paste0("Donor", 1:4)
source_sample_labels <- c(
  Donor1 = "HFA598",
  Donor2 = "HFA596",
  Donor3 = "HFA593",
  Donor4 = "HFA579"
)

display_order <- c("Donor3", "Donor4", "Donor2", "Donor1")
display_labels <- c(
  Donor3 = "Donor1",
  Donor4 = "Donor2",
  Donor2 = "Donor3",
  Donor1 = "Donor4"
)

condition_order <- c("Ctrl", "Stim", "Prime", "Dhit")

scale_01 <- function(z) {
  rng <- range(z, na.rm = TRUE)
  if (!is.finite(rng[1]) || diff(rng) == 0) return(rep(0.5, length(z)))
  (z - rng[1]) / diff(rng)
}

gene_blocks <- list()
idx <- 1
for (i in seq_len(nrow(raw))) {
  gene <- raw[i, 2]
  if (is.na(gene) || !nzchar(as.character(gene))) next

  block_rows <- i + 2:5
  if (max(block_rows) > nrow(raw)) next

  condition_names <- as.character(raw[block_rows, 1])
  if (!all(condition_order %in% condition_names)) next

  values <- matrix(NA_real_, nrow = length(condition_order), ncol = length(source_donor_labels))
  rownames(values) <- condition_order
  colnames(values) <- source_donor_labels

  for (d in seq_along(source_donor_labels)) {
    vals <- setNames(as.numeric(raw[block_rows, d + 1]), condition_names)
    values[, d] <- vals[condition_order]
  }

  gene_blocks[[idx]] <- list(gene = as.character(gene), values = values)
  idx <- idx + 1
}

gene_blocks <- gene_blocks[
  vapply(gene_blocks, function(x) toupper(x$gene) != "CCL2", logical(1))
]

count_crosses_for_gene <- function(vals) {
  donor_crosses <- numeric(ncol(vals))
  for (d in seq_len(ncol(vals))) {
    donor_vals <- vals[, d]
    donor_crosses[d] <- sum(c(
      donor_vals["Stim"] < donor_vals["Ctrl"],
      donor_vals["Prime"] > donor_vals["Stim"],
      donor_vals["Dhit"] < donor_vals["Stim"]
    ), na.rm = TRUE)
  }
  donor_crosses
}

order_stats <- data.frame(
  gene = vapply(gene_blocks, `[[`, character(1), "gene"),
  mean_crosses = NA_real_,
  total_crosses = NA_integer_,
  mean_log_memory = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(gene_blocks)) {
  vals <- gene_blocks[[i]]$values
  donor_crosses <- count_crosses_for_gene(vals)
  log_vals <- log2(vals)
  memory_diff <- log_vals["Dhit", ] - log_vals["Stim", ]

  order_stats$mean_crosses[i] <- mean(donor_crosses, na.rm = TRUE)
  order_stats$total_crosses[i] <- sum(donor_crosses, na.rm = TRUE)
  order_stats$mean_log_memory[i] <- mean(memory_diff, na.rm = TRUE)
}

order_stats <- order_stats[order(
  order_stats$mean_crosses,
  order_stats$total_crosses,
  -order_stats$mean_log_memory,
  order_stats$gene
), ]
gene_order <- order_stats$gene

records <- list()
idx <- 1
for (block in gene_blocks) {
  raw_vals <- block$values
  log_vals <- log2(raw_vals)

  for (source_donor in source_donor_labels) {
    raw_donor <- raw_vals[, source_donor]
    crossed <- c(
      Ctrl = FALSE,
      Stim = unname(raw_donor["Stim"] < raw_donor["Ctrl"]),
      Prime = unname(raw_donor["Prime"] > raw_donor["Stim"]),
      Dhit = unname(raw_donor["Dhit"] < raw_donor["Stim"])
    )
    row_has_no_crosses <- !any(crossed)

    for (x in seq_along(condition_order)) {
      condition <- condition_order[x]
      records[[idx]] <- data.frame(
        gene = block$gene,
        source_donor = source_donor,
        donor = unname(display_labels[source_donor]),
        x = x,
        condition = condition,
        value = as.numeric(log_vals[condition, source_donor]),
        star = condition == "Dhit" && row_has_no_crosses,
        slash = crossed[condition],
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
}

plot_data <- do.call(rbind, records)
plot_data$gene <- factor(plot_data$gene, levels = rev(gene_order))
plot_data$y <- as.numeric(plot_data$gene)
plot_data$source_donor <- factor(plot_data$source_donor, levels = display_order)
plot_data$donor <- factor(plot_data$donor, levels = unname(display_labels[display_order]))
plot_data$scaled <- ave(
  plot_data$value,
  plot_data$source_donor,
  plot_data$gene,
  FUN = scale_01
)

slash_data <- plot_data[plot_data$slash, ]
slash_data$x_start <- slash_data$x - 0.36
slash_data$x_end <- slash_data$x + 0.36
slash_data$y_start <- slash_data$y - 0.36
slash_data$y_end <- slash_data$y + 0.36

star_data <- plot_data[plot_data$star, ]
n_genes <- length(gene_order)
donor_headers <- data.frame(
  donor = factor(unname(display_labels[display_order]), levels = unname(display_labels[display_order])),
  x = 2.5,
  y = n_genes + 0.82,
  label = unname(display_labels[display_order])
)

p <- ggplot(plot_data, aes(x = x, y = y, fill = scaled)) +
  geom_tile(width = 1, height = 1) +
  geom_segment(
    data = slash_data,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    color = "#2b2b2b",
    linewidth = 0.9
  ) +
  geom_text(
    data = star_data,
    aes(x = x, y = y, label = "*"),
    inherit.aes = FALSE,
    color = "black",
    fontface = "bold",
    size = 8,
    vjust = 0.63
  ) +
  geom_text(
    data = donor_headers,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 6.2
  ) +
  facet_grid(. ~ donor, switch = "x") +
  scale_x_continuous(
    breaks = seq_along(condition_order),
    labels = parse(text = paste0("bold(log[2](", condition_order, "))")),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(plot_data$gene)),
    labels = levels(plot_data$gene),
    expand = c(0, 0)
  ) +
  scale_fill_viridis_c(
    option = "cividis",
    limits = c(0, 1),
    name = "scaled log value",
    guide = guide_colorbar(
      title.position = "right",
      title.theme = element_text(
        angle = 90,
        hjust = 0.5,
        face = "bold",
        size = 19
      ),
      label.theme = element_text(size = 15),
      barheight = unit(2.35, "in"),
      barwidth = unit(0.3, "in")
    )
  ) +
  coord_cartesian(
    xlim = c(0.5, 4.5),
    ylim = c(0.5, n_genes + 1.08),
    clip = "off"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(
      angle = 38,
      hjust = 1,
      vjust = 1,
      face = "bold",
      size = 15,
      color = "black"
    ),
    axis.text.y = element_text(face = "bold", size = 15, color = "black"),
    panel.grid = element_blank(),
    panel.spacing.x = unit(0.26, "in"),
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.title = element_text(face = "bold", size = 19),
    legend.text = element_text(size = 15, color = "black"),
    plot.margin = margin(t = 18, r = 18, b = 16, l = 16)
  )

ggsave(
  output_file,
  p,
  width = 14,
  height = 8,
  dpi = 128,
  bg = "white"
)

cat("Wrote final heatmap:", output_file, "\n")
cat("Displayed donor mapping:\n")
for (source_donor in display_order) {
  cat(
    unname(display_labels[source_donor]),
    "=",
    unname(source_sample_labels[source_donor]),
    "(original",
    source_donor,
    ")\n"
  )
}
