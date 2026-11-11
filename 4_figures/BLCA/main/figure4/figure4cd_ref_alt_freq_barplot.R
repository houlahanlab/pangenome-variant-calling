# Figure 4c-d
# Tumour and normal REF / ALT allele frequency of the HIGH and MODERATE impact Strelka2 FN (MC3-exclusive) and FP
# (alignment-exclusive) SNVs in bladder cancer driver genes (BLCA list, including genes also on the PANCAN list),
# one bar per variant, gene:position labels coloured by ancestry. Every plot is plot_width_in wide; the bar width
# is solved per plot to reach that width.
#   Figure 4c / 4d: VG pangenome, tumour / normal. The BWA GRCh38 plots are saved too, but are not in a figure.
#
# Input:  <data_dir>/False_Positive_Analysis/Allele_Freq/*_high_impact_plot_data.txt (build_ref_alt_freq.py)
#         driver gene list from config.yaml
# Output: <plots_dir>/Figure4/*.{pdf,png}; Source Data for the Figure 4 panels
# Prints the variants per gene category and the plot sizes.

# Setup
library(ggplot2)
library(dplyr)
library(readxl)
library(ggtext)

source("4_figures/figure_helpers.R")

ancestry_text_colours <- c(
  "AA_Patients" = "#2b6ca3",
  "EA_Patients" = "#b7950b",
  "EAA_Patients" = "#c2185b"
)

part_colours <- c(
  "REF" = "white",
  "ALT" = "steelblue"
)

category_levels <- c("BLCA", "BOTH", "PANCAN", "other")

# Gene categories plotted: bladder cancer driver genes only
plot_categories <- c("BLCA", "BOTH")

# Full plot width in inches
plot_width_in <- 20

# Segments shorter than this (allele frequency) are too small to hold their read count
min_label_freq <- 0.06

# Read in files
data_dir <- file.path(config$data_dir, "False_Positive_Analysis", "Allele_Freq")

genes <- read_xlsx(config$resources$driver_genes)

# Cancer gene categories
pancan <- genes$Gene[genes$Cancer == "PANCAN"]
blca <- genes$Gene[genes$Cancer == "BLCA"]

both <- intersect(pancan, blca)
pancan_only <- setdiff(pancan, blca)
blca_only <- setdiff(blca, pancan)

# Useful functions
read_plot_data <- function(file) {
  read.table(file.path(data_dir, file), header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
}

# FN and FP tables combined, ordered by type, gene category, gene, position
# Rows come in REF / ALT pairs, so each pair is one bar; the same gene:position in two donors gets two bars
build_plot_df <- function(fn_file, fp_file, fp_label) {
  fn <- read_plot_data(fn_file) %>% mutate(type = "MC3-exclusive")
  fp <- read_plot_data(fp_file) %>% mutate(type = fp_label)

  bind_rows(fn, fp) %>%
    mutate(
      pair = paste(type, ceiling(row_number() / 2)),
      gene = sub(":.*", "", gene_pos),
      position = as.numeric(sub(".*:", "", gene_pos)),
      category = case_when(
        gene %in% blca_only ~ "BLCA",
        gene %in% both ~ "BOTH",
        gene %in% pancan_only ~ "PANCAN",
        TRUE ~ "other"
      ),
      category = factor(category, levels = category_levels),
      type = factor(type, levels = c("MC3-exclusive", fp_label)),
      part = factor(part, levels = c("REF", "ALT"))
    ) %>%
    filter(category %in% plot_categories) %>%
    arrange(type, category, gene, position, pair) %>%
    mutate(bar_id = factor(pair, levels = unique(pair)))
}

# Variant counts and ALT frequency by ancestry, printed with each plot
print_counts <- function(df, title) {
  cat("\n========================================\n")
  cat(title, "\n")
  cat("========================================\n")

  print(
    df %>%
      group_by(type, ancestry, bar_id) %>%
      summarise(alt_freq = freq[part == "ALT"], depth = sum(count), .groups = "drop") %>%
      group_by(type, ancestry) %>%
      summarise(
        Variants = n(),
        Median_ALT_freq = median(alt_freq),
        Median_depth = median(depth),
        .groups = "drop"
      ) %>%
      as.data.frame()
  )
}

# Shared legend design, same as the Figure 2 / Figure 4 plots
legend_theme <- theme(
  text = element_text(family = "Helvetica"),
  legend.background = element_rect(fill = "white", colour = NA),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key.size = unit(12, "pt"),
  legend.key.spacing.y = unit(4, "pt")
)

# Width / height in inches of a size 6 text label or of the legend
text_size_in <- function(label) {
  g <- grid::textGrob(label, gp = grid::gpar(fontsize = 6 * .pt, fontfamily = "Helvetica"))
  c(
    width = grid::convertWidth(grid::grobWidth(g), "in", valueOnly = TRUE),
    height = grid::convertHeight(grid::grobHeight(g), "in", valueOnly = TRUE)
  )
}

legend_size <- function() {
  dummy <- data.frame(x = "a", y = 1, part = factor("REF", levels = names(part_colours)))

  g <- ggplotGrob(
    ggplot(dummy, aes(x, y, fill = part)) +
      geom_col(colour = "black", show.legend = TRUE) +
      scale_fill_manual(values = part_colours, name = "Group", limits = names(part_colours), drop = FALSE) +
      legend_theme +
      theme(legend.position = "right", legend.box.spacing = unit(0, "pt"))
  )

  box <- g$grobs[[which(g$layout$name == "guide-box-right")]]
  c(
    width = grid::convertWidth(sum(box$widths), "in", valueOnly = TRUE),
    height = grid::convertHeight(sum(box$heights), "in", valueOnly = TRUE)
  )
}

legend_in <- legend_size()

# Shared barplot design; df from build_plot_df()
# Every bar gets gene_width inches of panel; bars span 0-1 over bar_height inches
# Type brackets sit just above the bars; sample label (top left) and legend (top right) sit above the brackets
freq_barplot <- function(df, sample_label, gene_width, bar_height = 2.85) {

  n_bars <- nlevels(df$bar_id)

  type_brackets <- df %>%
    distinct(bar_id, type) %>%
    mutate(bar_index = as.numeric(bar_id)) %>%
    group_by(type) %>%
    summarise(
      xmin = min(bar_index) - 0.45,
      xmax = max(bar_index) + 0.45,
      .groups = "drop"
    ) %>%
    mutate(xmid = (xmin + xmax) / 2)

  # Widen bars if neighbouring type labels would overlap
  label_in <- sapply(as.character(type_brackets$type), function(label) text_size_in(label)[["width"]])
  n_types <- nrow(type_brackets)
  if (n_types > 1) {
    label_gap <- (label_in[-1] + label_in[-n_types]) / 2 + 0.15
    gene_width <- max(gene_width, label_gap / diff(type_brackets$xmid))
  }

  # Extra space at either end so the outer type labels stay inside the plot
  x_left <- max(0.6, (label_in[1] / 2 + 0.15) / gene_width - (type_brackets$xmid[1] - 1))
  x_right <- max(0.6, (label_in[n_types] / 2 + 0.15) / gene_width - (n_bars - type_brackets$xmid[n_types]))
  x_lo <- 1 - x_left
  x_hi <- n_bars + x_right

  panel_width <- gene_width * (x_hi - x_lo)

  # Heights in data units (1 unit = bar_height inches)
  text_h <- text_size_in("Tumour")[["height"]] / bar_height
  y_bracket_down <- 1.02
  y_bracket <- 1.05
  y_text <- 1.08
  y_band <- y_text + text_h + 0.12 / bar_height
  y_top <- y_band + (max(legend_in[["height"]], text_h * bar_height) + 0.05) / bar_height

  bar_labels <- df %>%
    distinct(bar_id, gene_pos, ancestry) %>%
    mutate(label = paste0("<span style='color:", ancestry_text_colours[ancestry], "'>", gene_pos, "</span>"))
  label_map <- setNames(bar_labels$label, bar_labels$bar_id)

  ggplot(df, aes(bar_id, freq, fill = part)) +
    geom_col(
      width = 0.85,
      colour = "black",
      linewidth = 0.4,
      show.legend = TRUE
    ) +
    geom_text(
      aes(label = ifelse(count != 0 & freq >= min_label_freq, count, "")),
      position = position_stack(vjust = 0.5),
      size = 4,
      colour = "black",
      family = "Helvetica"
    ) +
    # ALT segments too short to hold their count get it just above the segment
    geom_text(
      data = df %>% filter(part == "ALT", count != 0, freq < min_label_freq),
      aes(bar_id, freq, label = count),
      inherit.aes = FALSE,
      vjust = -0.4,
      size = 4,
      colour = "black",
      family = "Helvetica"
    ) +
    geom_segment(
      data = type_brackets,
      aes(x = xmin, xend = xmax, y = y_bracket, yend = y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_segment(
      data = type_brackets,
      aes(x = xmin, xend = xmin, y = y_bracket, yend = y_bracket_down),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_segment(
      data = type_brackets,
      aes(x = xmax, xend = xmax, y = y_bracket, yend = y_bracket_down),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_text(
      data = type_brackets,
      aes(x = xmid, y = y_text, label = type),
      inherit.aes = FALSE,
      size = 6,
      vjust = 0,
      family = "Helvetica"
    ) +
    annotate(
      "text",
      x = x_lo + 0.15 / gene_width, y = y_top,
      label = sample_label,
      size = 6,
      colour = "black",
      hjust = 0,
      vjust = 1,
      family = "Helvetica"
    ) +
    labs(
      x = "Gene:Position",
      y = "Allele Frequency"
    ) +
    scale_fill_manual(
      values = part_colours,
      name = "Group",
      limits = names(part_colours),
      drop = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, y_top),
      breaks = seq(0, 1, 0.2),
      expand = c(0, 0)
    ) +
    scale_x_discrete(
      labels = label_map,
      expand = expansion(add = c(x_left, x_right))
    ) +
    coord_cartesian(clip = "off") +
    legend_theme +
    theme(
      legend.position = "inside",
      legend.position.inside = c(1, 1),
      legend.justification = c(1, 1),
      panel.background = element_rect(fill = "white"),
      panel.widths = unit(panel_width, "in"),
      panel.heights = unit(y_top * bar_height, "in"),
      plot.title = element_text(hjust = 0.5, size = 15),
      axis.text.x = element_markdown(size = 14, angle = 60, hjust = 1),
      axis.text.y = element_text(size = 15, colour = "black"),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 16),
      axis.ticks.x = element_blank(),
      panel.border = element_blank(),
      axis.line.x.bottom = element_line(colour = "black"),
      axis.line.y.left = element_line(colour = "black"),
      plot.margin = margin(10, 15, 5, 5)
    )
}

# Full plot size in inches once the panel size is fixed
plot_size <- function(p) {
  g <- ggplotGrob(p)
  c(
    width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
    height = grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  )
}

# Plots
# One row per plot: FN / FP files, FP label, sample, output file name
plot_specs <- data.frame(
  fn_file = character(),
  fp_file = character(),
  fp_label = character(),
  sample = character(),
  file = character()
)

add_spec <- function(fn_file, fp_file, fp_label, sample, file) {
  plot_specs <<- rbind(plot_specs, data.frame(fn_file = fn_file, fp_file = fp_file, fp_label = fp_label,
                                              sample = sample, file = file))
}

add_spec("tumour_fn_high_impact_plot_data.txt", "tumour_fp_high_impact_plot_data.txt",
         "Pangenome-exclusive", "Tumour", "figure4c_tumour_allele_freq_high_impact")
add_spec("normal_fn_high_impact_plot_data.txt", "normal_fp_high_impact_plot_data.txt",
         "Pangenome-exclusive", "Normal", "figure4d_normal_allele_freq_high_impact")
add_spec("tumour_fn_bwa_linear_high_impact_plot_data.txt", "tumour_fp_bwa_linear_high_impact_plot_data.txt",
         "BWA-exclusive", "Tumour", "tumour_allele_freq_bwa_linear_high_impact")
add_spec("normal_fn_bwa_linear_high_impact_plot_data.txt", "normal_fp_bwa_linear_high_impact_plot_data.txt",
         "BWA-exclusive", "Normal", "normal_allele_freq_bwa_linear_high_impact")

plot_dfs <- list()
for (i in seq_len(nrow(plot_specs))) {
  spec <- plot_specs[i, ]
  plot_dfs[[spec$file]] <- build_plot_df(spec$fn_file, spec$fp_file, spec$fp_label)
  print_counts(plot_dfs[[spec$file]], spec$file)
}

# Bar width that makes a plot plot_width_in wide
# Plot width is close to linear in bar width (end padding shifts it a little), so refine with a few secant steps
fit_gene_width <- function(df, sample) {
  width_off <- function(w) plot_size(freq_barplot(df, sample, w))[["width"]] - plot_width_in
  w0 <- 0.3
  w1 <- 0.6
  for (step in 1:6) {
    f0 <- width_off(w0)
    f1 <- width_off(w1)
    if (abs(f1) < 0.005) break
    w2 <- w1 - f1 * (w1 - w0) / (f1 - f0)
    w0 <- w1
    w1 <- w2
  }
  w1
}

plots <- list()
for (i in seq_len(nrow(plot_specs))) {
  spec <- plot_specs[i, ]
  gene_width <- fit_gene_width(plot_dfs[[spec$file]], spec$sample)
  cat(spec$file, "bar width:", round(gene_width, 3), "in\n")
  plots[[spec$file]] <- freq_barplot(plot_dfs[[spec$file]], spec$sample, gene_width)
}

for (p in plots) print(p)

# Save plots and Source Data
source_data_panels <- c(
  "figure4c_tumour_allele_freq_high_impact" = "Figure 4c",
  "figure4d_normal_allele_freq_high_impact" = "Figure 4d"
)
for (i in seq_len(nrow(plot_specs))) {
  spec <- plot_specs[i, ]
  size <- plot_size(plots[[spec$file]])
  cat(spec$file, ":", round(size[["width"]], 2), "x", round(size[["height"]], 2), "in\n")

  out_dir <- file.path(config$plots_dir, "Figure4")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (ext in c("pdf", "png")) {
    ggsave(file.path(out_dir, paste0(spec$file, ".", ext)),
           plots[[spec$file]], width = size[["width"]], height = size[["height"]], dpi = 300, limitsize = FALSE)
  }

  if (spec$file %in% names(source_data_panels)) {
    write_source_data(select(plot_dfs[[spec$file]], type, gene_pos, gene, ancestry, impact, part, freq, count),
                      source_data_panels[[spec$file]])
  }
}
