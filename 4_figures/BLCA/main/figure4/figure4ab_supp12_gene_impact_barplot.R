# Figure 4a-b, Supplementary Figure 12a-b
# Strelka2 FN (MC3-exclusive) and FP (alignment-exclusive) SNVs in cancer driver genes per gene, stacked by ancestry.
# Genes are grouped into BLCA only, both BLCA and PANCAN, and PANCAN only. Uses the first SnpEff annotation of each SNV.
#   Figure 4a / 4b: VG pangenome MC3-exclusive (FN) / VG pangenome-exclusive (FP) SNVs, HIGH and MODERATE impact
#   Supplementary Figure 12a / 12b: BWA GRCh38 MC3-exclusive (FN) / BWA GRCh38-exclusive (FP) SNVs, HIGH and MODERATE impact
#   The plots with all impacts are saved too, but are not in a figure
#
# Input:  <data_dir>/Variant_Caller_Overlap/<alignment>/Shared_Variants_Annotations/{False_Negatives,False_Positives}/
#             *_filtered_shared_by_1_caller_in_capture_annotated.vcf (annotate_variant_sets.py)
#         driver gene list from config.yaml
# Output: <plots_dir>/Figure4, <plots_dir>/Supp12 (*.pdf, *.png); Source Data for each plot
# Prints the variants per gene category for each plot and the plot widths.

# Setup
library(ggplot2)
library(dplyr)
library(stringr)
library(readxl)

source("4_figures/figure_helpers.R")

ancestry_colours <- c(
  "AA_Patients" = "#add8e6",
  "EA_Patients" = "#FEF5C3",
  "EAA_Patients" = "#fdc0cb"
)

ancestry_labels <- c(
  "AA_Patients" = "African",
  "EA_Patients" = "European",
  "EAA_Patients" = "East Asian"
)

ancestry_levels <- c("AA_Patients", "EA_Patients", "EAA_Patients")

category_levels <- c("BLCA", "BOTH", "PANCAN", "other")

# Read in files
overlap_dir <- file.path(config$data_dir, "Variant_Caller_Overlap")

vcf_files <- list(
  "T2T_pangenome" = list(
    "FN" = file.path(overlap_dir, "T2T_pangenome/Shared_Variants_Annotations/False_Negatives/FN_filtered_shared_by_1_caller_in_capture_annotated.vcf"),
    "FP" = file.path(overlap_dir, "T2T_pangenome/Shared_Variants_Annotations/False_Positives/FP_filtered_shared_by_1_caller_in_capture_annotated.vcf")
  ),
  "GRCh38_bwa" = list(
    "FN" = file.path(overlap_dir, "GRCh38_bwa/Shared_Variants_Annotations/False_Negatives/FN_filtered_shared_by_1_caller_in_capture_annotated.vcf"),
    "FP" = file.path(overlap_dir, "GRCh38_bwa/Shared_Variants_Annotations/False_Positives/FP_filtered_shared_by_1_caller_in_capture_annotated.vcf")
  )
)

genes <- read_xlsx(config$resources$driver_genes)

samples <- read.table(
  "metadata/BLCA_Patients.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)

sample_ancestry <- samples %>%
  transmute(SAMPLE = trimws(V1), ancestry = trimws(V4)) %>%
  distinct()

if (anyDuplicated(sample_ancestry$SAMPLE) > 0) {
  stop("Donor IDs with more than one ancestry in BLCA_Patients.txt")
}

cat("Patients per ancestry:\n")
print(table(sample_ancestry$ancestry))

# Cancer gene categories
pancan <- genes$Gene[genes$Cancer == "PANCAN"]
blca <- genes$Gene[genes$Cancer == "BLCA"]

all_genes <- union(pancan, blca)
both <- intersect(pancan, blca)
pancan_only <- setdiff(pancan, blca)
blca_only <- setdiff(blca, pancan)

# Useful functions
# Strelka records of an annotated VCF with SAMPLE, ancestry and the first SnpEff annotation
read_annotated_vcf <- function(vcf_file) {
  vcf_data <- readLines(vcf_file)
  vcf_data <- vcf_data[!grepl("^#", vcf_data)]

  if (length(vcf_data) == 0) {
    stop("The VCF contains no variant records: ", vcf_file)
  }

  df <- read.table(
    text = vcf_data,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  if (ncol(df) < 8) {
    stop("The VCF does not contain the expected eight VCF columns: ", vcf_file)
  }

  colnames(df)[1:8] <- c("chr", "pos", "ID", "ref", "alt", "QUAL", "FILTER", "INFO")

  df$SAMPLE <- trimws(str_extract(df$INFO, "(?<=SAMPLE=)[^;]+"))
  df$TOOL <- str_extract(df$INFO, "(?<=TOOL=)[^;]+")

  # Keep only the first SnpEff annotation for each variant
  ann <- str_split_fixed(str_extract(df$INFO, "(?<=ANN=)[^;]+"), ",", 2)[, 1]
  ann_fields <- str_split_fixed(ann, "\\|", 16)

  df <- df %>%
    mutate(
      annotation = ann_fields[, 2],
      impact = ann_fields[, 3],
      gene = ann_fields[, 4]
    ) %>%
    filter(TOOL == "Strelka") %>%
    left_join(sample_ancestry, by = "SAMPLE")

  unmatched_samples <- unique(df$SAMPLE[is.na(df$ancestry)])
  if (length(unmatched_samples) > 0) {
    stop(length(unmatched_samples), " donor IDs could not be matched to ancestry information: ",
         paste(unmatched_samples, collapse = ", "))
  }

  df$ancestry <- factor(df$ancestry, levels = ancestry_levels)
  df
}

# Variants per gene and ancestry in driver genes; impact_only keeps HIGH + MODERATE
count_gene_variants <- function(df, impact_only = FALSE) {
  df <- df %>% filter(gene %in% all_genes)
  if (impact_only) df <- df %>% filter(impact %in% c("HIGH", "MODERATE"))

  df %>%
    mutate(
      category = case_when(
        gene %in% blca_only ~ "BLCA",
        gene %in% both ~ "BOTH",
        gene %in% pancan_only ~ "PANCAN",
        TRUE ~ "other"
      ),
      category = factor(category, levels = category_levels)
    ) %>%
    dplyr::count(gene, ancestry, category, name = "variants") %>%
    arrange(category, gene) %>%
    mutate(gene = factor(gene, levels = unique(gene)))
}

# Variant counts by ancestry, printed with each plot
print_counts <- function(df, title) {
  cat("\n========================================\n")
  cat(title, "\n")
  cat("========================================\n")

  cat("Genes:", n_distinct(df$gene), "  Variants:", sum(df$variants), "\n")
  print(
    df %>%
      group_by(ancestry, .drop = FALSE) %>%
      summarise(Genes = n_distinct(gene), Variants = sum(variants), .groups = "drop") %>%
      as.data.frame()
  )
}

# Shared legend design, same as the Figure 2 ancestry boxplots
legend_theme <- theme(
  text = element_text(family = "Helvetica"),
  legend.background = element_rect(fill = "white", colour = NA),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key.size = unit(12, "pt"),
  legend.key.spacing.y = unit(4, "pt")
)

# Legend width / height in inches, used to find a free spot for it inside the panel
legend_size <- function() {
  dummy <- data.frame(gene = "a", variants = 1, ancestry = factor(ancestry_levels[1], levels = ancestry_levels))

  g <- ggplotGrob(
    ggplot(dummy, aes(gene, variants, fill = ancestry)) +
      geom_col(colour = "black", show.legend = TRUE) +
      scale_fill_manual(values = ancestry_colours, name = "Ancestry", limits = ancestry_levels,
                        labels = ancestry_labels, drop = FALSE) +
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

# Shared barplot design; count_df needs gene, ancestry, category, variants
# Every gene gets gene_width inches of panel so bars are the same width across plots
# Category brackets sit above the bars; the legend goes in the emptiest spot below the brackets
# show.legend = TRUE keeps every ancestry in the legend even if it has no bars
# legend_over_bars = TRUE keeps the y scale and puts the legend just under the brackets, even over bars
gene_barplot <- function(count_df, y_label, gene_width, panel_height = 4, y_min_max = 1,
                         legend_over_bars = FALSE) {

  gene_totals <- count_df %>%
    group_by(gene) %>%
    summarise(total = sum(variants), .groups = "drop") %>%
    mutate(gene_index = as.numeric(gene))

  n_genes <- nrow(gene_totals)
  y_max <- max(max(gene_totals$total), y_min_max)
  y_step <- if (y_max > 20) 5 else if (y_max > 10) 2 else 1

  category_brackets <- count_df %>%
    distinct(gene, category) %>%
    mutate(gene_index = as.numeric(gene)) %>%
    group_by(category) %>%
    summarise(
      xmin = min(gene_index) - 0.45,
      xmax = max(gene_index) + 0.45,
      .groups = "drop"
    ) %>%
    mutate(xmid = (xmin + xmax) / 2)

  # Extra left space when the first category label sits over a single gene
  x_left <- if (category_brackets$xmax[1] - category_brackets$xmin[1] < 1) 1.2 else 0.6
  x_right <- 0.6
  x_lo <- 1 - x_left
  x_hi <- n_genes + x_right

  # Panel at least wide enough for the legend
  panel_width <- max(gene_width * (x_hi - x_lo), legend_in[["width"]] + 0.4)
  x_per_in <- (x_hi - x_lo) / panel_width

  # Legend window in gene units, slide it along and keep the lowest bars underneath (rightmost on ties)
  legend_genes <- (legend_in[["width"]] + 0.1) * x_per_in
  starts <- seq(x_lo, max(x_lo, x_hi - legend_genes), by = 0.1)
  under_max <- sapply(starts, function(s) {
    covered <- gene_totals$total[gene_totals$gene_index + 0.425 > s &
                                   gene_totals$gene_index - 0.425 < s + legend_genes]
    if (length(covered) == 0) 0 else max(covered)
  })
  best <- max(which(under_max == min(under_max)))
  legend_x <- starts[best] + legend_genes / 2
  under_max <- under_max[best]

  # Raise the brackets if the legend does not fit between the bars and the brackets
  # (bracket_base * 1.02 - under_max) * panel_height / (bracket_base * 1.28) >= legend height + padding
  legend_need <- (legend_in[["height"]] + 0.15) * 1.28 / panel_height
  bracket_base <- if (legend_over_bars) y_max else max(y_max, under_max / (1.02 - legend_need))

  y_bracket <- bracket_base * 1.08
  y_bracket_down <- bracket_base * 1.02
  y_text <- bracket_base * 1.15
  y_top <- bracket_base * 1.28

  legend_y <- if (legend_over_bars) {
    y_bracket_down - (legend_in[["height"]] / 2 + 0.05) * y_top / panel_height
  } else {
    (under_max + y_bracket_down) / 2
  }

  ggplot(count_df, aes(gene, variants, fill = ancestry)) +
    geom_col(
      width = 0.85,
      colour = "black",
      linewidth = 0.4,
      show.legend = TRUE
    ) +
    geom_segment(
      data = category_brackets,
      aes(x = xmin, xend = xmax, y = y_bracket, yend = y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_segment(
      data = category_brackets,
      aes(x = xmin, xend = xmin, y = y_bracket, yend = y_bracket_down),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_segment(
      data = category_brackets,
      aes(x = xmax, xend = xmax, y = y_bracket, yend = y_bracket_down),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_text(
      data = category_brackets,
      aes(x = xmid, y = y_text, label = category),
      inherit.aes = FALSE,
      size = 5,
      vjust = 0,
      family = "Helvetica"
    ) +
    labs(
      x = "Gene",
      y = y_label
    ) +
    scale_fill_manual(
      values = ancestry_colours,
      name = "Ancestry",
      limits = ancestry_levels,
      labels = ancestry_labels,
      drop = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, y_top),
      breaks = seq(0, y_max, y_step),
      expand = c(0, 0)
    ) +
    scale_x_discrete(
      expand = expansion(add = c(x_left, x_right))
    ) +
    # Category labels over a single gene can run past the panel
    coord_cartesian(clip = "off") +
    legend_theme +
    theme(
      legend.position = "inside",
      legend.position.inside = c((legend_x - x_lo) / (x_hi - x_lo), legend_y / y_top),
      legend.justification = c(0.5, 0.5),
      panel.background = element_rect(fill = "white"),
      panel.widths = unit(panel_width, "in"),
      panel.heights = unit(panel_height, "in"),
      plot.title = element_text(hjust = 0.5, size = 15),
      axis.text.x = element_text(size = 14, angle = 60, hjust = 1, colour = "black"),
      axis.text.y = element_text(size = 15, colour = "black"),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
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
# One row per plot: alignment, FP or FN, high impact only or not, output folder and file name
# Bar width (inches per gene) is shared within a figure, unless gene_width sets it for one plot,
# width sets the full plot width in inches, match_width names another plot whose full width
# this plot should take, or match_bars names another plot whose bar width this plot should use
# legend_over_bars lets the legend cover bars instead of squashing the y scale to make room for it
plot_specs <- data.frame(
  alignment = character(),
  status = character(),
  impact_only = logical(),
  folder = character(),
  file = character(),
  match_width = character(),
  gene_width = numeric(),
  width = numeric(),
  match_bars = character(),
  legend_over_bars = logical()
)

add_spec <- function(alignment, status, impact_only, folder, file, match_width = NA_character_,
                     gene_width = NA_real_, width = NA_real_, match_bars = NA_character_,
                     legend_over_bars = FALSE) {
  plot_specs <<- rbind(plot_specs, data.frame(alignment = alignment, status = status,
                                              impact_only = impact_only, folder = folder, file = file,
                                              match_width = match_width, gene_width = gene_width,
                                              width = width, match_bars = match_bars,
                                              legend_over_bars = legend_over_bars))
}

add_spec("T2T_pangenome", "FN", TRUE, "Figure4", "figure4a_fn_high_impact_variants_pangenome", width = 10)
add_spec("T2T_pangenome", "FP", TRUE, "Figure4", "figure4b_fp_high_impact_variants_pangenome",
         match_bars = "figure4a_fn_high_impact_variants_pangenome")

add_spec("T2T_pangenome", "FN", FALSE, "Supp12", "fn_all_variants_pangenome")
add_spec("T2T_pangenome", "FP", FALSE, "Supp12", "fp_all_variants_pangenome")
add_spec("GRCh38_bwa", "FN", FALSE, "Supp12", "fn_all_variants_bwa")
add_spec("GRCh38_bwa", "FN", TRUE, "Supp12", "supp12a_fn_high_impact_variants_bwa", legend_over_bars = TRUE)
add_spec("GRCh38_bwa", "FP", FALSE, "Supp12", "fp_all_variants_bwa")
add_spec("GRCh38_bwa", "FP", TRUE, "Supp12", "supp12b_fp_high_impact_variants_bwa", gene_width = 0.2)

# FN are MC3-exclusive; FP are exclusive to the alignment
y_labels <- c(
  "FN" = "Number of MC3-exclusive SNVs",
  "T2T_pangenome_FP" = "Number of pangenome-exclusive SNVs",
  "GRCh38_bwa_FP" = "Number of BWA-exclusive SNVs"
)

vcfs <- list()
for (alignment in names(vcf_files)) {
  for (status in names(vcf_files[[alignment]])) {
    vcfs[[paste(alignment, status, sep = "_")]] <- read_annotated_vcf(vcf_files[[alignment]][[status]])
  }
}

gene_widths <- c(
  "Figure4" = 0.2,
  "Supp12" = 0.3
)

plots <- list()
plot_builders <- list()
n_genes <- list()
used_gene_width <- list()

for (i in seq_len(nrow(plot_specs))) {
  spec <- plot_specs[i, ]
  vcf_name <- paste(spec$alignment, spec$status, sep = "_")

  count_df <- count_gene_variants(vcfs[[vcf_name]], spec$impact_only)
  print_counts(count_df, spec$file)

  y_label <- if (spec$status == "FN") y_labels[["FN"]] else y_labels[[vcf_name]]

  plot_builders[[spec$file]] <- local({
    count_df <- count_df
    y_label <- y_label
    y_min_max <- if (spec$impact_only) 2 else 3
    legend_over_bars <- spec$legend_over_bars
    function(gene_width) gene_barplot(count_df, y_label, gene_width, y_min_max = y_min_max,
                                      legend_over_bars = legend_over_bars)
  })

  n_genes[[spec$file]] <- n_distinct(count_df$gene)
  gene_width <- if (is.na(spec$gene_width)) gene_widths[[spec$folder]] else spec$gene_width
  used_gene_width[[spec$file]] <- gene_width
  plots[[spec$file]] <- plot_builders[[spec$file]](gene_width)
}

# Gene width that hits the target width: plot width is close to linear in gene width once the bars,
# not the legend, set the width, so start near target / genes and refine with a few secant steps
for (i in which(!is.na(plot_specs$match_width) | !is.na(plot_specs$width))) {
  spec <- plot_specs[i, ]
  build <- plot_builders[[spec$file]]
  target <- if (!is.na(spec$width)) spec$width else plot_size(plots[[spec$match_width]])[["width"]]
  width_off <- function(g) plot_size(build(g))[["width"]] - target

  g1 <- target / n_genes[[spec$file]]
  g0 <- g1 * 0.7
  for (step in 1:6) {
    f0 <- width_off(g0)
    f1 <- width_off(g1)
    if (abs(f1) < 0.005) break
    g2 <- g1 - f1 * (g1 - g0) / (f1 - f0)
    g0 <- g1
    g1 <- g2
  }
  gene_width <- g1

  cat("\n", spec$file, "width:", round(target, 2), "in, gene width", round(gene_width, 3), "in\n")

  used_gene_width[[spec$file]] <- gene_width
  plots[[spec$file]] <- build(gene_width)
}

for (i in which(!is.na(plot_specs$match_bars))) {
  spec <- plot_specs[i, ]
  gene_width <- used_gene_width[[spec$match_bars]]
  cat("\n", spec$file, "uses the bar width of", spec$match_bars, ":", round(gene_width, 3), "in\n")
  plots[[spec$file]] <- plot_builders[[spec$file]](gene_width)
}

for (p in plots) print(p)

# Save plots and Source Data
source_data_panels <- c(
  "figure4a_fn_high_impact_variants_pangenome" = "Figure 4a",
  "figure4b_fp_high_impact_variants_pangenome" = "Figure 4b",
  "supp12a_fn_high_impact_variants_bwa" = "Supplementary Figure 12a",
  "supp12b_fp_high_impact_variants_bwa" = "Supplementary Figure 12b"
)
for (i in seq_len(nrow(plot_specs))) {
  spec <- plot_specs[i, ]
  out_dir <- file.path(config$plots_dir, spec$folder)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  size <- plot_size(plots[[spec$file]])

  for (ext in c("pdf", "png")) {
    ggsave(file.path(out_dir, paste0(spec$file, ".", ext)),
           plots[[spec$file]], width = size[["width"]], height = size[["height"]], dpi = 300)
  }

  if (spec$file %in% names(source_data_panels)) {
    write_source_data(plots[[spec$file]]$data, source_data_panels[[spec$file]])
  }
}
