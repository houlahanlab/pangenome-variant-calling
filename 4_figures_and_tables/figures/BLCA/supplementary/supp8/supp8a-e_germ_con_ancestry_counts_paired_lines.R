# Supplementary Figure 8a-e
# Somatic calls (TP + FP) matching a gnomAD germline variant per donor, by alignment and ancestry:
#   8a / 8b: Strelka2 number of overlaps (exact gnomAD match / exact match with MAF >= 0.01)
#   8c / 8d: the same for Mutect2
#   8e: Mutect2 somatic-germline overlap (%), one line per donor across the three alignments (as Figure 3b for Strelka2)
# Brackets: one-sided paired Wilcoxon tests (BWA GRCh38 > VG) within each ancestry.
#
# Input:  <data_dir>/False_Positive_Analysis/Germline_Contamination/
#             Megatable_all_alignments_all_callers_germline_annotated.tsv (germline_contamination.py)
# Output: <plots_dir>/Supp8/supp8{a-e}_*.{pdf,png}; Source Data
# Prints the Wilcoxon tests for each plot.

# Setup
library(ggplot2)
library(dplyr)
library(tidyr)

source("4_figures/figure_helpers.R")

reference_colours <- c(
  "BWA_GRCh38" = "#C0DEA4",
  "VG_GRCh38" = "#228B22",
  "VG_Pangenome" = "#72A643"
)

reference_levels <- c("BWA_GRCh38", "VG_GRCh38", "VG_Pangenome")
reference_labels <- expression(BWA[GRCh38], VG[GRCh38], VG[Pangenome])

reference_names <- c(
  "GRCh38_bwa" = "BWA_GRCh38",
  "GRCh38_linear" = "VG_GRCh38",
  "T2T_pangenome" = "VG_Pangenome"
)

ancestry_levels <- c("AA_Patients", "EA_Patients", "EAA_Patients")

ancestry_labels <- c(
  "AA_Patients" = "African",
  "EA_Patients" = "European",
  "EAA_Patients" = "East Asian"
)

# Caller name in the megatable = caller name in the plot title
callers <- c("Strelka" = "Strelka2", "Mutect2" = "Mutect2")

# Box positions: one ancestry per whole number, the three alignments side by side
ancestry_x <- c("AA_Patients" = 1, "EA_Patients" = 2, "EAA_Patients" = 3)
reference_offsets <- c("BWA_GRCh38" = -0.26, "VG_GRCh38" = 0, "VG_Pangenome" = 0.26)

# Paired lines: grey = lower with VG pangenome than BWA, red = higher
direction_colours <- c("Decrease" = "grey55", "Increase" = "#C0392B")

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
    TOOL %in% names(callers),
    !is.na(Reference),
    DONOR_ID %in% patients
  )

# Per patient counts for each caller and alignment
counts <- germline %>%
  group_by(TOOL, DONOR_ID, ANCESTRY, Reference) %>%
  summarise(
    Germline_Overlap = sum(GNOMAD_EXACT_MATCH, na.rm = TRUE),
    Germline_Overlap_MAF_001 = sum(COMMON_GERMLINE_EXACT, na.rm = TRUE),
    Somatic_Calls = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Germline_Overlap_Percent = 100 * Germline_Overlap / Somatic_Calls,
    Ancestry = factor(ANCESTRY, levels = ancestry_levels),
    Reference = factor(Reference, levels = reference_levels),
    X_Position = ancestry_x[as.character(Ancestry)] + reference_offsets[as.character(Reference)]
  )

cat("Patients per caller, ancestry and alignment:\n")
print(table(counts$TOOL, counts$Ancestry, counts$Reference))

for (caller in names(callers)) {
  for (reference in reference_levels) {
    n_patients <- sum(counts$TOOL == caller & counts$Reference == reference)
    if (n_patients != length(patients)) {
      stop(caller, " ", reference, ": ", n_patients, " of the ", length(patients), " patients")
    }
  }
}

# Useful functions
# One-sided paired Wilcoxon (BWA GRCh38 > VG) for each ancestry and VG alignment
paired_pvalues <- function(caller_df, metric) {
  wide <- caller_df %>%
    select(DONOR_ID, Ancestry, Reference, all_of(metric)) %>%
    pivot_wider(names_from = Reference, values_from = all_of(metric))

  expand.grid(
    Ancestry = ancestry_levels,
    Second_Reference = c("VG_GRCh38", "VG_Pangenome"),
    stringsAsFactors = FALSE
  ) %>%
    rowwise() %>%
    mutate(
      P_Value = {
        ancestry_df <- wide[wide$Ancestry == Ancestry, ]
        differences <- ancestry_df$BWA_GRCh38 - ancestry_df[[Second_Reference]]
        if (all(differences == 0)) 1 else
          wilcox.test(ancestry_df$BWA_GRCh38, ancestry_df[[Second_Reference]],
                      paired = TRUE, alternative = "greater", exact = FALSE)$p.value
      },
      Decreased = sum(wide$BWA_GRCh38[wide$Ancestry == Ancestry] > wide[[Second_Reference]][wide$Ancestry == Ancestry]),
      Increased = sum(wide$BWA_GRCh38[wide$Ancestry == Ancestry] < wide[[Second_Reference]][wide$Ancestry == Ancestry])
    ) %>%
    ungroup()
}

# p = 8.28 x 10^-4 style labels
format_p <- function(p) {
  if (p < 0.001) {
    parts <- strsplit(formatC(p, format = "e", digits = 2), "e")[[1]]
    paste0("italic(p) == ", parts[1], " %*% 10^", as.integer(parts[2]))
  } else {
    paste0("italic(p) == \"", signif(p, 3), "\"")
  }
}

# Tick step giving about 5 ticks up to the top of the data
nice_step <- function(top) {
  steps <- c(1, 2, 5, 10, 20, 50, 100, 200, 500)
  steps[which(top / steps <= 6)[1]]
}

# Shared layout: grey shading behind the VG alignments, dashed lines between ancestries,
# p-value brackets, caller title, reference legend
make_reference_plot <- function(plot_df, metric, pvalue_df, caller_title, y_title, y_upper, y_step,
                                bracket_text_size, p_parse, whisker_caps = TRUE) {

  vg_background <- data.frame(
    xmin = c(0.87, 1.87, 2.87),
    xmax = c(1.50, 2.50, 3.50)
  )

  bracket_data <- pvalue_df %>%
    mutate(
      x = ancestry_x[Ancestry] + reference_offsets[["BWA_GRCh38"]],
      xend = ancestry_x[Ancestry] + reference_offsets[Second_Reference],
      y = ifelse(Second_Reference == "VG_GRCh38", y_upper * 0.76, y_upper * 0.86),
      label = if (p_parse) sapply(P_Value, format_p) else P_Value
    )

  legend_data <- data.frame(
    Reference = factor(reference_levels, levels = reference_levels),
    x = 1,
    y = 0
  )

  ggplot(plot_df, aes(Ancestry, .data[[metric]])) +
    # Sets up the ancestry axis; every layer below is placed by X_Position
    geom_blank() +
    geom_rect(
      data = vg_background,
      aes(xmin = xmin, xmax = xmax, ymin = -0.03 * y_upper, ymax = y_upper),
      inherit.aes = FALSE,
      fill = "grey80",
      alpha = 0.45,
      colour = NA
    ) +
    geom_vline(
      xintercept = c(1.5, 2.5),
      linetype = "dashed",
      colour = "grey65",
      linewidth = 0.55
    ) +
    stat_boxplot(
      aes(group = interaction(Ancestry, Reference), x = X_Position),
      geom = "errorbar",
      width = if (whisker_caps) 0.10 else 0,
      colour = "black",
      linewidth = 0.5
    ) +
    geom_boxplot(
      aes(group = interaction(Ancestry, Reference), x = X_Position, fill = Reference),
      width = 0.20,
      outlier.shape = NA,
      colour = "black",
      linewidth = 0.5,
      show.legend = FALSE
    ) +
    # p-value brackets
    geom_segment(
      data = bracket_data,
      aes(x = x, xend = xend, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_segment(
      data = bracket_data,
      aes(x = x, xend = x, y = y - y_upper * 0.015, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_segment(
      data = bracket_data,
      aes(x = xend, xend = xend, y = y - y_upper * 0.015, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_text(
      data = bracket_data,
      aes(x = (x + xend) / 2, y = y + y_upper * 0.025, label = label),
      inherit.aes = FALSE,
      parse = p_parse,
      size = bracket_text_size,
      family = "Helvetica"
    ) +
    # Invisible points so the legend shows the three alignments as filled squares
    geom_point(
      data = legend_data,
      aes(x = x, y = y, fill = Reference),
      inherit.aes = FALSE,
      shape = 22,
      size = 0,
      alpha = 0,
      show.legend = TRUE
    ) +
    annotate(
      "text",
      x = 0.68, y = Inf,
      label = caller_title,
      size = 7,
      colour = "black",
      hjust = 0,
      vjust = 1.7,
      family = "Helvetica"
    ) +
    labs(
      x = NULL,
      y = y_title,
      fill = "Reference"
    ) +
    scale_x_discrete(
      labels = ancestry_labels,
      expand = expansion(add = c(0.45, 0.45))
    ) +
    scale_y_continuous(
      breaks = seq(0, y_upper, by = y_step),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(ylim = c(-0.03 * y_upper, y_upper), clip = "on") +
    scale_fill_manual(
      values = reference_colours,
      breaks = reference_levels,
      labels = reference_labels,
      drop = FALSE
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(shape = 22, size = 6, alpha = 1, colour = "black", linewidth = 0.7),
        keyheight = grid::unit(0.85, "lines"),
        keywidth = grid::unit(0.85, "lines")
      )
    ) +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "inside",
      legend.position.inside = c(0.995, 0.045),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      legend.spacing.y = grid::unit(0.15, "cm"),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 16, colour = "black", margin = margin(t = 7)),
      axis.text.y = element_text(size = 13, colour = "black"),
      axis.title.y = element_text(size = 15, margin = margin(r = 10)),
      panel.border = element_blank(),
      axis.line.x.bottom = element_line(colour = "black", linewidth = 0.6),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.6),
      plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
    )
}

# Panel letter of each plot in Supplementary Figure 8
panel_letters <- c(
  "Strelka_germline_overlap_counts" = "supp8a",
  "Strelka_germline_overlap_MAF001_counts" = "supp8b",
  "Mutect2_germline_overlap_counts" = "supp8c",
  "Mutect2_germline_overlap_MAF001_counts" = "supp8d",
  "Mutect2_germline_overlap_paired_lines" = "supp8e"
)

save_supp8_plot <- function(p, output_name, width = 10, height = 6) {
  print(p)
  file_name <- paste0(panel_letters[[output_name]], "_", output_name)
  ggsave(file.path(out_dir, paste0(file_name, ".pdf")),
         p, width = width, height = height, units = "in")
  ggsave(file.path(out_dir, paste0(file_name, ".png")),
         p, width = width, height = height, units = "in", dpi = 300, bg = "white")
}

print_tests <- function(pvalue_df, title) {
  cat("\n========================================\n")
  cat(title, "\n")
  cat("========================================\n")
  print(as.data.frame(pvalue_df), row.names = FALSE)
}

# Plots
out_dir <- file.path(config$plots_dir, "Supp8")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Supp 8a-d: number of somatic-germline overlaps
count_plots <- data.frame(
  metric = c("Germline_Overlap", "Germline_Overlap_MAF_001"),
  file = c("germline_overlap_counts", "germline_overlap_MAF001_counts")
)

count_y_titles <- list(
  "Germline_Overlap" = "Number of Somatic-Germline Overlaps",
  "Germline_Overlap_MAF_001" = expression(atop("Number of Somatic-Germline Overlaps", "(MAF" >= "0.01)"))
)

plots <- list()

for (caller in names(callers)) {
  caller_df <- counts %>% filter(TOOL == caller)

  for (i in seq_len(nrow(count_plots))) {
    metric <- count_plots$metric[i]
    output_name <- paste0(caller, "_", count_plots$file[i])

    pvalue_df <- paired_pvalues(caller_df, metric)
    print_tests(pvalue_df, output_name)

    # Brackets sit at 76% and 86% of the axis, above the highest point
    y_top <- max(caller_df[[metric]]) / 0.72
    y_step <- nice_step(y_top)
    y_upper <- ceiling(y_top / y_step) * y_step

    set.seed(1)
    plot_df <- caller_df %>% mutate(Jitter_X = X_Position + runif(n(), -0.05, 0.05))

    p <- make_reference_plot(plot_df, metric, pvalue_df, callers[[caller]],
                             count_y_titles[[metric]], y_upper, y_step,
                             bracket_text_size = 4.2, p_parse = TRUE) +
      geom_point(
        aes(x = Jitter_X),
        shape = 21,
        size = 1.8,
        stroke = 0.45,
        fill = NA,
        colour = "black",
        alpha = 0.8
      )

    save_supp8_plot(p, output_name)
    plots[[output_name]] <- p
  }
}

# Supp 8e: Mutect2 overlap (%) with one line per donor
mutect_df <- counts %>% filter(TOOL == "Mutect2")

pvalue_df <- paired_pvalues(mutect_df, "Germline_Overlap_Percent")
print_tests(pvalue_df, "Mutect2_germline_overlap_paired_lines")

# Line colour from the overall BWA GRCh38 -> VG pangenome change
donor_direction <- mutect_df %>%
  select(DONOR_ID, Reference, Germline_Overlap_Percent) %>%
  pivot_wider(names_from = Reference, values_from = Germline_Overlap_Percent) %>%
  mutate(Direction = ifelse(BWA_GRCh38 - VG_Pangenome > 0, "Decrease", "Increase")) %>%
  select(DONOR_ID, Direction)

line_df <- mutect_df %>% left_join(donor_direction, by = "DONOR_ID")

y_upper <- ceiling(max(line_df$Germline_Overlap_Percent) / 5) * 5

# Same p = 0.074 / p = 8.3 x 10^-4 text as the original paired lines script
pvalue_df$P_Value <- sapply(pvalue_df$P_Value, function(p) {
  if (p < 0.001) paste0("p = ", sub("e-0*", " × 10^-", formatC(p, format = "e", digits = 2)))
  else paste0("p = ", signif(p, 2))
})

p <- make_reference_plot(line_df, "Germline_Overlap_Percent", pvalue_df, "Mutect2",
                         "Somatic-Germline Overlap (%)", y_upper, 10,
                         bracket_text_size = 4.5, p_parse = FALSE, whisker_caps = FALSE) +
  geom_line(
    aes(x = X_Position, group = DONOR_ID, colour = Direction),
    linewidth = 0.55,
    alpha = 0.75,
    show.legend = FALSE
  ) +
  geom_point(
    aes(x = X_Position, colour = Direction),
    size = 1.5,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = direction_colours, guide = "none")

save_supp8_plot(p, "Mutect2_germline_overlap_paired_lines")
plots[["Mutect2_germline_overlap_paired_lines"]] <- p

# Source Data: overlap counts and % per donor, alignment and caller
source_cols <- function(df) select(df, -any_of(c("ANCESTRY", "X_Position", "Jitter_X", "Ancestry_X", "Direction")))
write_source_data(source_cols(filter(counts, TOOL == "Strelka")), "Supplementary Figure 8a-b")
write_source_data(source_cols(filter(counts, TOOL == "Mutect2")), "Supplementary Figure 8c-e")
