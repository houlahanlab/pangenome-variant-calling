# Supplementary Figure 4a-c
# Recall, precision and F1 per donor of the >= 3 caller consensus against MC3:
# 4a VG pangenome, 4b BWA GRCh38, 4c VG pangenome - BWA GRCh38 (paired by donor).
#
# Input:  Supplementary Table 1 (BLCA; Consensus_* columns)
# Output: <plots_dir>/Supp4/supp4{a,b,c}_*.{pdf,png}; Source Data
# Prints the donors per alignment and the median and IQR of each metric.

# Setup
library(ggplot2)
library(dplyr)

source("4_figures/figure_helpers.R")

consensus_colour <- c("Consensus_ge3" = "#5592D9")

# Read in files
consensus <- supp_table_1_metrics("BLCA", "Consensus") %>%
  filter(Reference %in% c("T2T_pangenome", "GRCh38_bwa")) %>%
  mutate(Variant_Caller = "Consensus_ge3")

print(table(consensus$Reference))

# Useful functions
to_long <- function(df) {
  long_df <- data.frame(
    Donor_ID = rep(df$Donor_ID, times = 3),
    Variant_Caller = rep(df$Variant_Caller, times = 3),
    Metric = rep(c("F1_Score", "Precision", "Recall"), each = nrow(df)),
    Value = c(df$F1_Score, df$Precision, df$Recall)
  )

  long_df$Metric <- factor(
    long_df$Metric,
    levels = c("F1_Score", "Precision", "Recall")
  )

  long_df
}

print_metric_summary <- function(long_df, title) {
  cat("\n========================================\n")
  cat(title, "\n")
  cat("========================================\n")
  print(
    long_df %>%
      group_by(Metric) %>%
      summarise(
        N = sum(!is.na(Value)),
        Mean = mean(Value, na.rm = TRUE),
        Median = median(Value, na.rm = TRUE),
        SD = sd(Value, na.rm = TRUE),
        .groups = "drop"
      )
  )
}

# One Consensus >=3 box per metric
consensus_boxplot <- function(long_df, plot_label, y_label, y_limits, y_breaks,
                              legend_pos, zero_line = FALSE) {

  set.seed(1)

  ggplot(long_df, aes(Metric, Value, fill = Variant_Caller)) +
    annotate(
      "rect",
      xmin = -Inf, xmax = 1.5,
      ymin = -Inf, ymax = Inf,
      alpha = 0.2,
      fill = "grey"
    ) +
    stat_boxplot(
      geom = "errorbar",
      width = 0.3,
      colour = "black",
      linewidth = 0.6
    ) +
    {if (zero_line) geom_hline(
      yintercept = 0,
      colour = "black",
      linewidth = 1
    )} +
    geom_boxplot(
      width = 0.6,
      outlier.shape = NA,
      alpha = 1,
      colour = "black",
      linewidth = 0.6
    ) +
    geom_jitter(
      aes(colour = Variant_Caller),
      width = 0.1,
      height = 0,
      size = 2,
      shape = 21,
      fill = NA,
      alpha = 0.9
    ) +
    labs(
      x = NULL,
      y = y_label
    ) +
    scale_x_discrete(
      labels = c(
        "F1_Score" = "F1 Score",
        "Precision" = "Precision",
        "Recall" = "Recall"
      ),
      expand = expansion(mult = c(0.4, 0.3))
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks
    ) +
    scale_fill_manual(values = consensus_colour) +
    scale_colour_manual(
      values = c("Consensus_ge3" = "black"),
      name = NULL,
      labels = c("Consensus_ge3" = "Consensus ≥3")
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(
          shape = 22,
          size = 6,
          fill = consensus_colour,
          colour = "black"
        )
      ),
      fill = "none"
    ) +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "inside",
      legend.position.inside = legend_pos,
      legend.background = element_rect(fill = NA, colour = NA),
      legend.key = element_rect(fill = NA, colour = NA),
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
      label = plot_label,
      parse = TRUE,
      size = 6,
      colour = "black",
      hjust = 0,
      vjust = 2,
      family = "Helvetica"
    )
}

# Supp 4a / 4b: metrics per reference
pan_long <- to_long(consensus %>% filter(Reference == "T2T_pangenome"))
bwa_long <- to_long(consensus %>% filter(Reference == "GRCh38_bwa"))

print_metric_summary(pan_long, "Consensus >=3 - VG pangenome")
print_metric_summary(bwa_long, "Consensus >=3 - BWA GRCh38")

plot_pan <- consensus_boxplot(
  pan_long,
  "VG[pangenome]",
  y_label = "Metric Percentage (%)",
  y_limits = c(-1, 100),
  y_breaks = seq(0, 100, 20),
  legend_pos = c(0.75, 0.08)
)

plot_bwa <- consensus_boxplot(
  bwa_long,
  "BWA[GRCh38]",
  y_label = "Metric Percentage (%)",
  y_limits = c(-1, 100),
  y_breaks = seq(0, 100, 20),
  legend_pos = c(0.75, 0.08)
)

# Supp 4c: VG pangenome minus BWA GRCh38, matched by donor
pan <- consensus %>%
  filter(Reference == "T2T_pangenome") %>%
  select(Donor_ID, Variant_Caller, F1_Score, Precision, Recall)

bwa <- consensus %>%
  filter(Reference == "GRCh38_bwa") %>%
  select(Donor_ID, Variant_Caller, F1_Score, Precision, Recall)

if (anyDuplicated(pan$Donor_ID) > 0 || anyDuplicated(bwa$Donor_ID) > 0) {
  stop("Duplicate donor rows found in the Consensus >=3 data")
}

merged <- inner_join(pan, bwa, by = c("Donor_ID", "Variant_Caller"), suffix = c("_pan", "_bwa"))
cat("\nMatched donors for the delta:", nrow(merged), "\n")

delta <- data.frame(
  Donor_ID = merged$Donor_ID,
  Variant_Caller = merged$Variant_Caller,
  F1_Score = merged$F1_Score_pan - merged$F1_Score_bwa,
  Precision = merged$Precision_pan - merged$Precision_bwa,
  Recall = merged$Recall_pan - merged$Recall_bwa
)

delta_long <- to_long(delta)
print_metric_summary(delta_long, "Consensus >=3 - VG pangenome minus BWA GRCh38")

plot_delta <- consensus_boxplot(
  delta_long,
  "VG[pangenome] - BWA[GRCh38]",
  y_label = expression(Delta ~ "Metric Percentage (%)"),
  y_limits = c(-90, 50),
  y_breaks = seq(-80, 40, 20),
  legend_pos = c(0.75, 0.05),
  zero_line = TRUE
)

plot_pan
plot_bwa
plot_delta

# Save plots and Source Data
plots_to_save <- list(
  "supp4a_consensus_ge3_boxplot_vg_pangenome" = plot_pan,
  "supp4b_consensus_ge3_boxplot_bwa_grch38" = plot_bwa,
  "supp4c_consensus_ge3_delta_pangenome_minus_bwa" = plot_delta
)
for (plot_name in names(plots_to_save)) {
  save_plot(plots_to_save[[plot_name]], plot_name, "Supp4", width = 5, height = 6)
}

write_source_data(pan_long, "Supplementary Figure 4a")
write_source_data(bwa_long, "Supplementary Figure 4b")
write_source_data(delta_long, "Supplementary Figure 4c")
