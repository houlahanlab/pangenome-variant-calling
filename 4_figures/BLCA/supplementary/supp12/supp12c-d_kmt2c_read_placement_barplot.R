# Supplementary Figure 12c-d
# Where the reads carrying the ALT allele of the four KMT2C FN SNVs align: one stacked barplot for BWA GRCh38 and one for
# VG pangenome (the same reads), coloured by chromosome (12c VG pangenome, 12d BWA GRCh38), with the % of reads with MAPQ 0 above each bar.
#
# Input:  <data_dir>/False_Negatives_Analysis/kmt2c_read_placement/kmt2c_alt_read_placement.txt (build_kmt2c_read_placement.py)
# Output: <plots_dir>/Supp12/supp12{c,d}_*.{pdf,png}; Source Data

# Setup
library(ggplot2)
library(dplyr)

source("4_figures/figure_helpers.R")

chrom_colours <- c(
  "chr7 (KMT2C)" = "#5b7db1",
  "chr13" = "#fdc0cb",
  "chr14" = "#FEF5C3",
  "chr15" = "#b5e7a0",
  "chr21" = "#add8e6",
  "chr22" = "#d7b5d8",
  "unmapped" = "#d9d9d9",
  "not found" = "#ffffff"
)

alignment_titles <- c(
  "GRCh38_bwa" = "BWA[GRCh38]",
  "T2T_pangenome" = "VG[pangenome]"
)

# Read in files
data_dir <- file.path(config$data_dir, "False_Negatives_Analysis", "kmt2c_read_placement")
out_dir <- file.path(config$plots_dir, "Supp12")

placement <- read.table(file.path(data_dir, "kmt2c_alt_read_placement.txt"), header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)

# Build plot data
site_levels <- unique(placement$site)
chrom_levels <- intersect(names(chrom_colours), placement$chrom)

plot_df <- placement %>%
  group_by(site, alignment, chrom) %>%
  summarise(reads = sum(reads), .groups = "drop") %>%
  mutate(
    site = factor(site, levels = site_levels),
    chrom = factor(chrom, levels = chrom_levels)
  )

bar_labels <- placement %>%
  group_by(site, alignment) %>%
  summarise(
    total = sum(reads),
    mapq0 = sum(reads[mapq == "MAPQ 0"]),
    .groups = "drop"
  ) %>%
  mutate(
    site = factor(site, levels = site_levels),
    label = paste0(round(100 * mapq0 / total), "%\nMAPQ 0")
  )

print(as.data.frame(bar_labels))

# Same y axis for both plots
y_max <- max(bar_labels$total)
y_top <- y_max * 1.35

# Bar spacing in inches per site
site_width <- 0.56

# Useful functions
site_barplot <- function(alignment_name) {
  ggplot(plot_df %>% filter(alignment == alignment_name), aes(site, reads, fill = chrom)) +
    geom_col(
      width = 0.7,
      colour = "black",
      linewidth = 0.4,
      show.legend = TRUE
    ) +
    geom_text(
      data = bar_labels %>% filter(alignment == alignment_name),
      aes(site, total, label = label),
      inherit.aes = FALSE,
      vjust = -0.3,
      size = 4,
      lineheight = 0.9,
      family = "Helvetica"
    ) +
    annotate(
      "text",
      x = 0.3, y = Inf,
      label = alignment_titles[[alignment_name]],
      parse = TRUE,
      size = 6,
      colour = "black",
      hjust = 0,
      vjust = 1.3,
      family = "Helvetica"
    ) +
    labs(
      x = "KMT2C position",
      y = "ALT supporting reads",
      fill = "Aligned to"
    ) +
    scale_fill_manual(
      values = chrom_colours,
      limits = chrom_levels,
      drop = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, y_top),
      breaks = seq(0, ceiling(y_max / 10) * 10, 10),
      expand = c(0, 0)
    ) +
    scale_x_discrete(
      expand = expansion(add = c(0.8, 0.6))
    ) +
    theme(
      text = element_text(family = "Helvetica"),
      panel.background = element_rect(fill = "white"),
      panel.widths = unit(site_width * (length(site_levels) + 1.4), "in"),
      panel.heights = unit(4, "in"),
      axis.text.x = element_text(size = 14, angle = 60, hjust = 1, colour = "black"),
      axis.text.y = element_text(size = 15, colour = "black"),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      axis.ticks.x = element_blank(),
      panel.border = element_blank(),
      axis.line.x.bottom = element_line(colour = "black"),
      axis.line.y.left = element_line(colour = "black"),
      legend.position = "inside",
      legend.position.inside = c(1, 0.88),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.key.size = unit(9, "pt"),
      legend.key.spacing.y = unit(2, "pt"),
      legend.margin = margin(2, 2, 2, 2),
      plot.margin = margin(10, 15, 5, 5)
    )
}

plots <- list(
  "supp12c_kmt2c_read_placement_vg_pangenome" = site_barplot("T2T_pangenome"),
  "supp12d_kmt2c_read_placement_bwa_grch38" = site_barplot("GRCh38_bwa")
)

for (p in plots) print(p)

# Save plots
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (name in names(plots)) {
  g <- ggplotGrob(plots[[name]])
  width <- grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE)
  height <- grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE)

  for (ext in c("pdf", "png")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), plots[[name]], width = width, height = height, dpi = 300)
  }
}

write_source_data(placement, "Supplementary Figure 12c-d")
