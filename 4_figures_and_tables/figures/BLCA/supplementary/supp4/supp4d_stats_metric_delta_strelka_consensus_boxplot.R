# Supplementary Figure 4d
# Strelka2 - >= 3 caller consensus difference in recall, precision and F1 per donor, one box per alignment
# (VG pangenome, VG GRCh38, BWA GRCh38), paired by donor and alignment.
#
# Input:  Supplementary Table 1 (BLCA; Strelka_* and Consensus_* columns)
# Output: <plots_dir>/Supp4/supp4d_*.{pdf,png}; Source Data
# Prints the matched donors and the median difference per alignment.

# Setup
library(ggplot2)
library(dplyr)

source("4_figures/figure_helpers.R")

reference_colours <- c(
  "T2T_pangenome" = "#72a643",
  "GRCh38_linear" = "#228B22",
  "GRCh38_bwa" = "#C0DEA4"
)

black_colours <- c(
  "T2T_pangenome" = "black",
  "GRCh38_linear" = "black",
  "GRCh38_bwa" = "black"
)

reference_labels <- c(
  "T2T_pangenome" = "VG pangenome",
  "GRCh38_linear" = "VG GRCh38",
  "GRCh38_bwa" = "BWA GRCh38"
)

reference_levels <- names(reference_colours)

# Read in files
metrics <- supp_table_1_metrics("BLCA", c("Strelka", "Consensus")) %>%
  filter(Reference %in% reference_levels)
strelka <- filter(metrics, Variant_Caller == "Strelka")
consensus <- filter(metrics, Variant_Caller == "Consensus") %>% mutate(Variant_Caller = "Consensus_ge3")

# Strelka2 minus Consensus >=3, matched by donor and alignment
keep_cols <- c("Donor_ID", "Reference", "F1_Score", "Precision", "Recall")

merged <- inner_join(
  strelka[, keep_cols],
  consensus[, keep_cols],
  by = c("Donor_ID", "Reference"),
  suffix = c("_strelka", "_consensus")
)

cat("Matched donors per alignment:\n")
print(table(merged$Reference))

delta <- data.frame(
  Donor_ID = merged$Donor_ID,
  Reference = merged$Reference,
  F1_Score = merged$F1_Score_strelka - merged$F1_Score_consensus,
  Precision = merged$Precision_strelka - merged$Precision_consensus,
  Recall = merged$Recall_strelka - merged$Recall_consensus
)

delta_long <- data.frame(
  Reference = rep(delta$Reference, times = 3),
  Metric = rep(c("F1_Score", "Precision", "Recall"), each = nrow(delta)),
  Value = c(delta$F1_Score, delta$Precision, delta$Recall)
)

delta_long$Metric <- factor(delta_long$Metric, levels = c("F1_Score", "Precision", "Recall"))
delta_long$Reference <- factor(delta_long$Reference, levels = reference_levels)

cat("\n========================================\n")
cat("Strelka2 - Consensus >=3 summary\n")
cat("========================================\n")
print(
  delta_long %>%
    group_by(Reference, Metric) %>%
    summarise(
      N = sum(!is.na(Value)),
      Mean = mean(Value, na.rm = TRUE),
      Median = median(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      .groups = "drop"
    ),
  n = Inf
)

# Plot
dodge <- position_dodge(width = 0.9)
jitterdodge <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9)

set.seed(1)

plot_delta <- ggplot(delta_long, aes(Metric, Value, fill = Reference)) +
  annotate(
    "rect",
    xmin = -Inf, xmax = 1.5,
    ymin = -Inf, ymax = Inf,
    alpha = 0.2,
    fill = "grey"
  ) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.36,
    colour = "black",
    linewidth = 0.6,
    position = dodge
  ) +
  geom_hline(
    yintercept = 0,
    colour = "black",
    linewidth = 1
  ) +
  geom_boxplot(
    width = 0.72,
    outlier.shape = NA,
    alpha = 1,
    colour = "black",
    linewidth = 0.6,
    position = dodge
  ) +
  geom_jitter(
    size = 2,
    shape = 21,
    aes(colour = Reference),
    fill = NA,
    alpha = 0.9,
    position = jitterdodge
  ) +
  labs(
    x = NULL,
    y = expression(Delta ~ "Metric Percentage (%)")
  ) +
  scale_x_discrete(labels = c(
    "F1_Score" = "F1 Score",
    "Precision" = "Precision",
    "Recall" = "Recall"
  )) +
  scale_y_continuous(
    limits = c(-80, 80),
    breaks = seq(-80, 80, 20),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  scale_fill_manual(
    values = reference_colours,
    name = "Alignment"
  ) +
  scale_colour_manual(
    values = black_colours,
    name = "Alignment",
    labels = reference_labels
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        shape = 22,
        size = 6,
        fill = reference_colours,
        colour = "black"
      )
    ),
    fill = "none"
  ) +
  theme(
    text = element_text(family = "Helvetica"),
    legend.position = "inside",
    legend.position.inside = c(0.88, 0.85),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 15, colour = "black"),
    axis.title.y = element_text(size = 15),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black")
  ) +
  annotate(
    "text",
    x = 0.5, y = Inf,
    label = "Strelka2 - Consensus ≥3",
    size = 6,
    colour = "black",
    hjust = 0,
    vjust = 2,
    family = "Helvetica"
  ) +
  geom_vline(
    xintercept = c(1.5, 2.5),
    linetype = "dashed",
    colour = "grey",
    linewidth = 0.6
  )

plot_delta

# Save plots and Source Data
save_plot(plot_delta, "supp4d_consensus_ge3_vs_strelka_delta_all_alignments", "Supp4")
write_source_data(delta, "Supplementary Figure 4d")
