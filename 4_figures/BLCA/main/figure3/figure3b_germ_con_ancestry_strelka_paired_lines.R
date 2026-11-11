# Figure 3b
# Strelka2 somatic-germline overlap per donor: % of TP + FP calls with an exact gnomAD match, by ancestry and
# alignment, one line per donor (BWA GRCh38 -> VG GRCh38 -> VG pangenome). Brackets: one-sided paired Wilcoxon
# tests (BWA GRCh38 > VG) within each ancestry. The test is paired, so a small decrease that is consistent across
# donors can be significant while the medians barely move.
#
# Input:  <data_dir>/False_Positive_Analysis/Germline_Contamination/
#             Megatable_all_alignments_all_callers_germline_annotated.tsv (germline_contamination.py)
# Output: <plots_dir>/Figure3/figure3b_Strelka_germline_contamination_paired_lines_full.{pdf,png}; Source Data
# Prints the donors per ancestry and alignment, the median per group and the Wilcoxon tests.

# Setup
library(ggplot2)
library(dplyr)
library(tidyr)

source("4_figures/figure_helpers.R")

caller_name <- "Strelka2"

reference_colours <- c(
  "BWA_GRCh38" = "#C0DEA4",
  "VG_GRCh38" = "#228B22",
  "VG_Pangenome" = "#72A643"
)

reference_levels <- names(reference_colours)
reference_labels <- expression(BWA[GRCh38], VG[GRCh38], VG[Pangenome])

ancestry_labels <- c(
  "AA" = "African",
  "EA" = "European",
  "EAA" = "East Asian"
)

ancestry_levels <- names(ancestry_labels)

# Per-donor change from BWA GRCh38 to VG pangenome: grey = lower with VG (expected), red = higher with VG
direction_colours <- c(
  "Decrease" = "grey55",
  "Increase" = "#C0392B"
)

# Each donor's three points sit at the ancestry position plus a per-alignment offset
ancestry_x_positions <- c("AA" = 1, "EA" = 2, "EAA" = 3)
reference_offsets <- c("BWA_GRCh38" = -0.26, "VG_GRCh38" = 0, "VG_Pangenome" = 0.26)

# Read in files
patients <- read_patients()

# Same filtering as Supplementary Figure 9, so both describe the same numbers
germline <- read.delim(
  file.path(config$data_dir, "False_Positive_Analysis", "Germline_Contamination",
            "Megatable_all_alignments_all_callers_germline_annotated.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  mutate(
    DONOR_ID = trimws(DONOR_ID),
    STATUS = toupper(trimws(STATUS)),
    TOOL = trimws(TOOL),
    GNOMAD_EXACT_MATCH = toupper(trimws(as.character(GNOMAD_EXACT_MATCH))) == "TRUE",
    Ancestry = sub("_Patients$", "", trimws(ANCESTRY)),
    Reference = recode(trimws(ALIGNMENT), GRCh38_bwa = "BWA_GRCh38", GRCh38_linear = "VG_GRCh38",
                       T2T_pangenome = "VG_Pangenome", .default = NA_character_)
  ) %>%
  filter(STATUS %in% c("TP", "FP"), Ancestry %in% ancestry_levels, !is.na(Reference), DONOR_ID %in% patients)

# Strelka2: % of TP + FP calls with an exact gnomAD match, per donor and alignment
strelka_proportions <- germline %>%
  filter(TOOL == "Strelka") %>%
  group_by(DONOR_ID, Ancestry, Reference) %>%
  summarise(
    Germline_Overlap_Count = sum(GNOMAD_EXACT_MATCH, na.rm = TRUE),
    Total_Somatic_SNVs = n(),
    Germline_Contamination_Percent = 100 * Germline_Overlap_Count / Total_Somatic_SNVs,
    .groups = "drop"
  ) %>%
  mutate(
    Ancestry = factor(Ancestry, levels = ancestry_levels),
    Reference = factor(Reference, levels = reference_levels)
  )

# Only donors with all three alignments can be drawn as a paired line and enter the paired test
complete_donors <- strelka_proportions %>%
  count(DONOR_ID) %>%
  filter(n == length(reference_levels)) %>%
  pull(DONOR_ID)

dropped_donors <- setdiff(unique(strelka_proportions$DONOR_ID), complete_donors)
if (length(dropped_donors) > 0) {
  cat("\nDonors without all three references, dropped from the paired figure:\n")
  print(dropped_donors)
}

strelka_proportions <- filter(strelka_proportions, DONOR_ID %in% complete_donors)

cat("\nDonors retained per ancestry and reference:\n")
print(count(strelka_proportions, Ancestry, Reference))

wide_proportions <- strelka_proportions %>%
  select(DONOR_ID, Ancestry, Reference, Germline_Contamination_Percent) %>%
  pivot_wider(names_from = Reference, values_from = Germline_Contamination_Percent)

# Useful functions
# Same p-value format as Figure 3d
format_p <- function(p) {
  if (length(p) == 0 || is.na(p)) return("p = NA")
  if (p < 0.001) {
    sci <- formatC(p, format = "e", digits = 2)
    sci <- sub("e\\+0*", " × 10^", sub("e-0*", " × 10^-", sci))
    return(paste0("p = ", sci))
  }
  paste0("p = ", signif(p, 2))
}

# One-sided paired Wilcoxon, BWA GRCh38 > VG
paired_wilcox <- function(df, reference_1, reference_2) {
  paired_df <- filter(df, !is.na(.data[[reference_1]]), !is.na(.data[[reference_2]]))
  if (nrow(paired_df) < 2) return(NA_real_)
  if (all(paired_df[[reference_1]] - paired_df[[reference_2]] == 0)) return(1)
  wilcox.test(paired_df[[reference_1]], paired_df[[reference_2]], paired = TRUE, alternative = "greater",
              exact = FALSE)$p.value
}

# Tests and donor direction counts per ancestry
comparison_second_reference <- c(
  "BWA_vs_VG_GRCh38" = "VG_GRCh38",
  "BWA_vs_VG_Pangenome" = "VG_Pangenome"
)

pvalue_df <- bind_rows(lapply(ancestry_levels, function(ancestry_name) {
  ancestry_df <- filter(wide_proportions, Ancestry == ancestry_name)
  bind_rows(lapply(names(comparison_second_reference), function(comparison_name) {
    second_reference <- comparison_second_reference[[comparison_name]]
    differences <- ancestry_df$BWA_GRCh38 - ancestry_df[[second_reference]]
    differences <- differences[!is.na(differences)]
    data.frame(
      Ancestry = ancestry_name,
      Comparison = comparison_name,
      Second_Reference = second_reference,
      P_Value = paired_wilcox(ancestry_df, "BWA_GRCh38", second_reference),
      Donors_Total = length(differences),
      Donors_Decreased = sum(differences > 0),
      Donors_Increased = sum(differences < 0),
      Donors_Unchanged = sum(differences == 0),
      Median_Difference = median(differences),
      stringsAsFactors = FALSE
    )
  }))
})) %>%
  mutate(
    Ancestry = factor(Ancestry, levels = ancestry_levels),
    Comparison = factor(Comparison, levels = names(comparison_second_reference))
  )

cat("\n========================================\n", caller_name, " paired germline-overlap tests\n",
    "========================================\n", sep = "")
print(pvalue_df, row.names = FALSE)

# The medians barely move while most donors move in the same direction
cat("\nMedian percent per ancestry and reference:\n")
print(
  strelka_proportions %>%
    group_by(Ancestry, Reference) %>%
    summarise(Median_Percent = round(median(Germline_Contamination_Percent), 2), .groups = "drop") %>%
    pivot_wider(names_from = Reference, values_from = Median_Percent),
  row.names = FALSE
)

# Plot
# One colour per donor line, from the net BWA GRCh38 -> VG pangenome change
donor_direction <- wide_proportions %>%
  mutate(Direction = ifelse(BWA_GRCh38 - VG_Pangenome > 0, "Decrease", "Increase")) %>%
  select(DONOR_ID, Direction)

line_df <- strelka_proportions %>%
  left_join(donor_direction, by = "DONOR_ID") %>%
  mutate(
    X_Position = unname(ancestry_x_positions[as.character(Ancestry)]) +
      unname(reference_offsets[as.character(Reference)]),
    Direction = factor(Direction, levels = names(direction_colours))
  )

y_upper <- ceiling(max(line_df$Germline_Contamination_Percent, na.rm = TRUE) / 5) * 5
cat("Full y-axis maximum:", y_upper, "%\n")

# Same dodge width and x expansion as Supplementary Figure 9
dodge <- position_dodge(width = 0.78)

# Grey background behind the VG alignments
vg_background <- data.frame(
  xmin = c(0.87, 1.87, 2.87),
  xmax = c(1.50, 2.50, 3.50),
  ymin = -0.03 * y_upper,
  ymax = y_upper
)

# P-value brackets
bracket_data <- pvalue_df %>%
  mutate(
    Ancestry_X = unname(ancestry_x_positions[as.character(Ancestry)]),
    x = Ancestry_X + reference_offsets[["BWA_GRCh38"]],
    xend = Ancestry_X + unname(reference_offsets[Second_Reference]),
    y = ifelse(Comparison == "BWA_vs_VG_GRCh38", y_upper * 0.76, y_upper * 0.86),
    label = vapply(P_Value, format_p, character(1))
  )

# Invisible points that give the fill legend square swatches instead of boxplot keys
legend_data <- data.frame(Reference = factor(reference_levels, levels = reference_levels), x = 1, y = 0)

plot_lines <- ggplot(line_df, aes(Ancestry, Germline_Contamination_Percent)) +
  geom_rect(
    data = vg_background,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey80",
    alpha = 0.45,
    colour = NA
  ) +
  geom_segment(
    data = data.frame(x = c(1.5, 2.5)),
    aes(x = x, xend = x, y = -0.03 * y_upper, yend = y_upper),
    inherit.aes = FALSE,
    linetype = "dashed",
    colour = "grey65",
    linewidth = 0.55
  ) +
  geom_boxplot(
    aes(fill = Reference),
    width = 0.60,
    outlier.shape = NA,
    colour = "black",
    linewidth = 0.5,
    alpha = 1,
    position = dodge,
    show.legend = FALSE
  ) +
  # One line per donor: the pairing the Wilcoxon test uses
  geom_line(
    aes(x = X_Position, group = DONOR_ID, colour = Direction),
    linewidth = 0.55,
    alpha = 0.75,
    show.legend = FALSE
  ) +
  geom_point(
    aes(x = X_Position, group = DONOR_ID, colour = Direction),
    size = 1.5,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  geom_segment(
    data = bracket_data,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_segment(
    data = bracket_data,
    aes(x = x, xend = x, y = y - y_upper * 0.015, yend = y),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_segment(
    data = bracket_data,
    aes(x = xend, xend = xend, y = y - y_upper * 0.015, yend = y),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_text(
    data = bracket_data,
    aes(x = (x + xend) / 2, y = y + y_upper * 0.020, label = label),
    inherit.aes = FALSE,
    size = 4.5,
    family = "Helvetica"
  ) +
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
    label = caller_name,
    size = 7,
    colour = "black",
    hjust = 0,
    vjust = 1.7,
    family = "Helvetica"
  ) +
  scale_x_discrete(labels = ancestry_labels[ancestry_levels], expand = expansion(add = c(0.45, 0.45))) +
  scale_y_continuous(breaks = seq(0, y_upper, by = 10), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(-0.03 * y_upper, y_upper), clip = "on") +
  scale_fill_manual(values = reference_colours, breaks = reference_levels, labels = reference_labels, drop = FALSE) +
  scale_colour_manual(values = direction_colours, breaks = names(direction_colours), drop = FALSE, guide = "none") +
  labs(x = NULL, y = "Somatic-Germline Overlap (%)", fill = "Reference") +
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

plot_lines

# Save plots and Source Data
save_plot(plot_lines, "figure3b_Strelka_germline_contamination_paired_lines_full", "Figure3", width = 10, height = 6)
write_source_data(select(line_df, DONOR_ID, Ancestry, Reference, Germline_Overlap_Count, Total_Somatic_SNVs,
                         Germline_Contamination_Percent), "Figure 3b")
