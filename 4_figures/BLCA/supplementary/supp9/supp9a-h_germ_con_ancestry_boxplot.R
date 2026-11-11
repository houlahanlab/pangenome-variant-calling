# Supplementary Figure 9a-h
# Per donor change in somatic calls (TP + FP) matching a gnomAD germline variant, by ancestry, for Strelka2 (9a-d) and
# Mutect2 (9e-h): VG pangenome - BWA GRCh38 and VG GRCh38 - BWA GRCh38, for exact gnomAD matches and exact matches
# with MAF >= 0.01.
#
# Input:  <data_dir>/False_Positive_Analysis/Germline_Contamination/
#             Megatable_all_alignments_all_callers_germline_annotated.tsv (germline_contamination.py)
# Output: <plots_dir>/Supp9/supp9{a-h}_*.{pdf,png}; Source Data
# Prints the median change and the tests for each plot.

# Setup
library(ggplot2)
library(dplyr)

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

reference_names <- c(
  "GRCh38_bwa" = "BWA_GRCh38",
  "GRCh38_linear" = "VG_GRCh38",
  "T2T_pangenome" = "VG_Pangenome"
)

# Caller name in the megatable = caller name used in the output file names
callers <- c("Strelka", "Mutect2")

# Space on each side of the boxes (x axis units). European stays centred; the extra room
# left of the African box is for arrows added in Adobe. Box, whisker and jitter widths are
# scaled so the boxes keep the same size on the page as with the original 0.30 padding
x_padding <- 0.75
width_scale <- (2 + 2 * x_padding) / (2 + 2 * 0.30)

# Read in files
patients <- read_patients()

germline <- read.delim(
  file.path(config$data_dir, "False_Positive_Analysis", "Germline_Contamination",
            "Megatable_all_alignments_all_callers_germline_annotated.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

germline <- germline %>%
  mutate(
    STATUS = toupper(trimws(STATUS)),
    Reference = reference_names[ALIGNMENT],
    GNOMAD_EXACT_MATCH = toupper(trimws(GNOMAD_EXACT_MATCH)) == "TRUE",
    COMMON_GERMLINE_EXACT = toupper(trimws(COMMON_GERMLINE_EXACT)) == "TRUE"
  ) %>%
  filter(
    STATUS %in% c("TP", "FP"),
    TOOL %in% callers,
    !is.na(Reference),
    DONOR_ID %in% patients
  )

cat("Patients per caller and ancestry:\n")
print(table(distinct(germline, TOOL, DONOR_ID, ANCESTRY)[, c("TOOL", "ANCESTRY")]))

# Useful functions
# VG minus BWA germline overlap count per patient for one caller, alignment and metric
make_delta_table <- function(caller, vg_reference, metric) {
  counts <- germline %>%
    filter(TOOL == caller) %>%
    group_by(DONOR_ID, ANCESTRY, Reference) %>%
    summarise(
      Germline_Overlap = sum(GNOMAD_EXACT_MATCH, na.rm = TRUE),
      Germline_Overlap_MAF_001 = sum(COMMON_GERMLINE_EXACT, na.rm = TRUE),
      .groups = "drop"
    )

  bwa <- counts %>%
    filter(Reference == "BWA_GRCh38") %>%
    select(DONOR_ID, Ancestry = ANCESTRY, BWA = all_of(metric))

  vg <- counts %>%
    filter(Reference == vg_reference) %>%
    select(DONOR_ID, VG = all_of(metric))

  inner_join(bwa, vg, by = "DONOR_ID") %>%
    mutate(
      Delta = VG - BWA,
      Ancestry = factor(Ancestry, levels = ancestry_levels)
    )
}

print_delta_summary <- function(delta_df, title) {
  cat("\n========================================\n")
  cat(title, "\n")
  cat("========================================\n")
  print(
    delta_df %>%
      group_by(Ancestry) %>%
      summarise(
        N = n(),
        Mean_Delta = mean(Delta),
        Median_Delta = median(Delta),
        Minimum = min(Delta),
        Maximum = max(Delta),
        .groups = "drop"
      )
  )
}

# Tick spacing that gives at most 4 steps below 0, with one labelled tick above 0
# and half a step more for the title (room for arrows / labels added in Adobe)
y_axis_settings <- function(values) {
  steps <- c(1, 2, 5, 10, 20, 50, 100, 200, 500)
  step <- steps[which(-min(values, 0) / steps <= 4)[1]]

  y_lower <- floor(min(values, 0) / step) * step
  y_top_tick <- max(step, ceiling(max(values) / step) * step)

  list(
    breaks = seq(y_lower, y_top_tick, by = step),
    limits = c(y_lower, y_top_tick + step / 2)
  )
}

make_delta_plot <- function(delta_df, plot_title, y_title, output_name) {

  # Invisible points so the legend shows all three ancestries as filled squares
  legend_data <- data.frame(
    Ancestry = factor(ancestry_levels, levels = ancestry_levels),
    x = 1,
    y = 0
  )

  y_axis <- y_axis_settings(delta_df$Delta)

  set.seed(1)

  p <- ggplot(delta_df, aes(Ancestry, Delta, fill = Ancestry)) +
    geom_hline(
      yintercept = 0,
      colour = "black",
      linewidth = 0.65
    ) +
    stat_boxplot(
      geom = "errorbar",
      width = 0.18 * width_scale,
      colour = "black",
      linewidth = 0.55,
      show.legend = FALSE
    ) +
    geom_boxplot(
      width = 0.42 * width_scale,
      outlier.shape = NA,
      colour = "black",
      linewidth = 0.55,
      show.legend = FALSE
    ) +
    geom_jitter(
      width = 0.07 * width_scale,
      height = 0,
      shape = 21,
      size = 1.8,
      stroke = 0.45,
      fill = NA,
      colour = "black",
      alpha = 0.8,
      show.legend = FALSE
    ) +
    geom_point(
      data = legend_data,
      aes(x = x, y = y, fill = Ancestry),
      inherit.aes = FALSE,
      shape = 22,
      size = 0,
      alpha = 0,
      show.legend = TRUE
    ) +
    # Title centred over the middle box
    annotate(
      "text",
      x = 2, y = Inf,
      label = plot_title,
      parse = TRUE,
      size = 5.5,
      colour = "black",
      hjust = 0.5,
      vjust = 1.5,
      family = "Helvetica"
    ) +
    labs(
      x = NULL,
      y = parse(text = y_title)[[1]],
      fill = "Ancestry"
    ) +
    scale_x_discrete(
      labels = ancestry_labels,
      expand = expansion(add = c(x_padding, x_padding))
    ) +
    scale_y_continuous(
      breaks = y_axis$breaks,
      expand = expansion(mult = c(0.03, 0))
    ) +
    coord_cartesian(ylim = y_axis$limits, clip = "on") +
    scale_fill_manual(
      values = ancestry_colours,
      breaks = ancestry_levels,
      labels = ancestry_labels,
      drop = FALSE
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(shape = 22, size = 5, alpha = 1, colour = "black", linewidth = 0.7),
        keyheight = grid::unit(0.75, "lines"),
        keywidth = grid::unit(0.75, "lines")
      )
    ) +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "inside",
      legend.position.inside = c(0.95, 0.08),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      legend.spacing.y = grid::unit(0.10, "cm"),
      legend.title = element_text(size = 12.5),
      legend.text = element_text(size = 12.5),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      # hjust below 0.5 nudges "East Asian" right so it does not touch "European"
      axis.text.x = element_text(size = 11, colour = "black", margin = margin(t = 6), hjust = c(0.5, 0.5, 0.3)),
      axis.text.y = element_text(size = 11.5, colour = "black"),
      axis.title.y = element_text(size = 12.5, margin = margin(r = 8), hjust = 0.5),
      panel.border = element_blank(),
      axis.line.x.bottom = element_line(colour = "black", linewidth = 0.6),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.6),
      plot.margin = margin(t = 10, r = 8, b = 8, l = 8)
    )

  print(p)

  file_name <- paste0(panel_letters[[output_name]], "_", output_name)
  ggsave(file.path(out_dir, paste0(file_name, ".pdf")),
         p, width = 3.5, height = 4.5, units = "in")
  ggsave(file.path(out_dir, paste0(file_name, ".png")),
         p, width = 3.5, height = 4.5, units = "in", dpi = 600, bg = "white")

  p
}

# Panel letter of each plot in Supplementary Figure 9
panel_letters <- c(
  "Strelka_germline_overlap_delta_VG_pangenome_vs_BWA" = "supp9a",
  "Strelka_germline_overlap_delta_VG_GRCh38_vs_BWA" = "supp9b",
  "Strelka_germline_overlap_MAF001_delta_VG_pangenome_vs_BWA" = "supp9c",
  "Strelka_germline_overlap_MAF001_delta_VG_GRCh38_vs_BWA" = "supp9d",
  "Mutect2_germline_overlap_delta_VG_pangenome_vs_BWA" = "supp9e",
  "Mutect2_germline_overlap_delta_VG_GRCh38_vs_BWA" = "supp9f",
  "Mutect2_germline_overlap_MAF001_delta_VG_pangenome_vs_BWA" = "supp9g",
  "Mutect2_germline_overlap_MAF001_delta_VG_GRCh38_vs_BWA" = "supp9h"
)

# Plots
out_dir <- file.path(config$plots_dir, "Supp9")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# One row per plot for each caller: VG alignment, metric, title, file name ending
comparisons <- data.frame(
  vg_reference = c("VG_Pangenome", "VG_GRCh38", "VG_Pangenome", "VG_GRCh38"),
  metric = c("Germline_Overlap", "Germline_Overlap", "Germline_Overlap_MAF_001", "Germline_Overlap_MAF_001"),
  title = c("VG[pangenome] - BWA[GRCh38]", "VG[GRCh38] - BWA[GRCh38]",
            "VG[pangenome] - BWA[GRCh38]", "VG[GRCh38] - BWA[GRCh38]"),
  y_title = c(rep('Delta*" Number of Somatic-Germline Overlaps"', 2),
              rep('Delta*" Number of Somatic-Germline Overlaps"~"MAF" >= 0.01', 2)),
  file = c("germline_overlap_delta_VG_pangenome_vs_BWA", "germline_overlap_delta_VG_GRCh38_vs_BWA",
           "germline_overlap_MAF001_delta_VG_pangenome_vs_BWA", "germline_overlap_MAF001_delta_VG_GRCh38_vs_BWA")
)

plots <- list()

for (caller in callers) {
  for (i in seq_len(nrow(comparisons))) {
    comparison <- comparisons[i, ]
    output_name <- paste0(caller, "_", comparison$file)

    delta_df <- make_delta_table(caller, comparison$vg_reference, comparison$metric)

    if (nrow(delta_df) != length(patients)) {
      stop(output_name, ": ", nrow(delta_df), " of the ", length(patients), " patients")
    }

    print_delta_summary(delta_df, output_name)
    write_source_data(delta_df, paste0("Supplementary Figure ", sub("supp", "", panel_letters[[output_name]])))
    plots[[output_name]] <- make_delta_plot(delta_df, comparison$title, comparison$y_title, output_name)
  }
}
