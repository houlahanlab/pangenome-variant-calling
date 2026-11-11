# Figure 5b-d (LUAD)
# Recall, precision and F1 of Strelka2 and Mutect2 against MC3, per LUAD donor (as Figure 1b-d for BLCA).
#   Figure 5b: VG pangenome; Figure 5c: VG pangenome - BWA GRCh38; Figure 5d: VG pangenome - VG GRCh38
# The BWA GRCh38 and VG GRCh38 plots are saved to Extras.
#
# Input:  Supplementary Table 1 (LUAD)
# Output: <plots_dir>/LUAD/Figure5, <plots_dir>/LUAD/Extras (*.pdf, *.png); Source Data for Figure 5b-d
# Prints the donors per alignment and caller, and the median and mean differences for Figure 5c-d.

# Setup
library(ggplot2)
library(dplyr)
library(patchwork)

source("4_figures/figure_helpers.R")

caller_colours <- c(
  "Strelka" = "#f1a1ff",
  "Mutect2" = "#9c8fff"
)

black_colours <- c(
  "Strelka" = "black",
  "Mutect2" = "black"
)

# Read in files
combined <- supp_table_1_metrics("LUAD", c("Strelka", "Mutect2"))

combined$Variant_Caller <- factor(
  combined$Variant_Caller,
  levels = c("Strelka", "Mutect2")
)

combined$Reference <- factor(
  combined$Reference,
  levels = c("T2T_pangenome", "GRCh38_bwa", "GRCh38_linear")
)

print(table(combined$Reference, combined$Variant_Caller))

# Useful functions
# Shared boxplot design; long_df needs Variant_Caller, Metric, Value
plot_metric_boxplot <- function(long_df, plot_label, y_label, y_limits, y_breaks,
                                legend_pos, zero_line = FALSE,
                                y_expand = expansion(mult = 0.05)) {

  dodge <- position_dodge(width = 0.65)
  jitterdodge <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.65)

  set.seed(1)

  p <- ggplot(long_df, aes(Metric, Value, fill = Variant_Caller)) +
    annotate(
      "rect",
      xmin = -Inf, xmax = 1.5,
      ymin = -Inf, ymax = Inf,
      alpha = 0.2,
      fill = "grey"
    ) +
    stat_boxplot(
      geom = "errorbar",
      width = 0.25,
      colour = "black",
      linewidth = 0.6,
      position = dodge
    ) +
    {if (zero_line) geom_hline(
      yintercept = 0,
      colour = "black",
      linewidth = 1
    )} +
    geom_boxplot(
      width = 0.5,
      outlier.shape = NA,
      alpha = 1,
      colour = "black",
      linewidth = 0.6,
      position = dodge
    ) +
    geom_jitter(
      size = 2,
      shape = 21,
      aes(colour = Variant_Caller),
      fill = NA,
      alpha = 0.9,
      position = jitterdodge
    ) +
    labs(
      x = NULL,
      y = y_label
    ) +
    scale_x_discrete(labels = c(
      "F1_Score" = "F1 Score",
      "Precision" = "Precision",
      "Recall" = "Recall"
    )) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      expand = y_expand
    ) +
    scale_fill_manual(
      values = caller_colours,
      name = "Variant Caller"
    ) +
    scale_colour_manual(
      values = black_colours,
      name = "Variant Caller",
      labels = c(
        "Strelka" = "Strelka2",
        "Mutect2" = "Mutect2"
      )
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(
          shape = 22,
          size = 6,
          fill = caller_colours,
          colour = "black"
        )
      ),
      fill = "none"
    ) +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "inside",
      legend.position.inside = legend_pos,
      legend.background = element_rect(fill = "white", colour = NA),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      panel.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5, size = 15),
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
    ) +
    geom_vline(
      xintercept = c(1.5, 2.5),
      linetype = "dashed",
      colour = "grey",
      linewidth = 0.6
    )

  return(p)
}

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

  long_df[!is.na(long_df$Value), ]
}

# Figure 5b: raw metrics for one reference
make_boxplot <- function(df, ref_name, plot_label) {

  df_ref <- df %>%
    filter(Reference == ref_name)

  plot_metric_boxplot(
    to_long(df_ref),
    plot_label,
    y_label = "Metric Percentage (%)",
    y_limits = c(-1, 100),
    y_breaks = seq(0, 100, 20),
    legend_pos = c(0.88, 0.19)
  )
}

# Pangenome minus another reference, matched by donor and caller
calc_delta <- function(df, ref_a, ref_b) {

  df_a <- df %>%
    filter(Reference == ref_a) %>%
    select(Donor_ID, Variant_Caller, Recall, Precision, F1_Score)

  df_b <- df %>%
    filter(Reference == ref_b) %>%
    select(Donor_ID, Variant_Caller, Recall, Precision, F1_Score)

  if (anyDuplicated(df_a[, c("Donor_ID", "Variant_Caller")]) > 0) {
    stop(paste("Duplicate donor rows found for", ref_a))
  }

  if (anyDuplicated(df_b[, c("Donor_ID", "Variant_Caller")]) > 0) {
    stop(paste("Duplicate donor rows found for", ref_b))
  }

  merged <- inner_join(df_a, df_b, by = c("Donor_ID", "Variant_Caller"), suffix = c("_a", "_b"))

  cat("\n", ref_a, "-", ref_b, "matched donors:\n")
  print(table(merged$Variant_Caller))

  data.frame(
    Donor_ID = merged$Donor_ID,
    Variant_Caller = merged$Variant_Caller,
    F1_Score = merged$F1_Score_a - merged$F1_Score_b,
    Precision = merged$Precision_a - merged$Precision_b,
    Recall = merged$Recall_a - merged$Recall_b
  )
}

summarise_delta <- function(delta_df, title) {

  delta_summary <- delta_df %>%
    group_by(Variant_Caller) %>%
    summarise(
      N = n(),
      Median_Delta_F1 = median(F1_Score, na.rm = TRUE),
      Median_Delta_Precision = median(Precision, na.rm = TRUE),
      Median_Delta_Recall = median(Recall, na.rm = TRUE),
      Mean_Delta_F1 = mean(F1_Score, na.rm = TRUE),
      Mean_Delta_Precision = mean(Precision, na.rm = TRUE),
      Mean_Delta_Recall = mean(Recall, na.rm = TRUE)
    )

  cat("\n========================================\n")
  cat(title, "summary\n")
  cat("========================================\n")
  print(delta_summary, width = Inf)
}

# Figure 5c / 5d: delta metrics
make_delta_boxplot <- function(delta_df, plot_label) {

  plot_metric_boxplot(
    to_long(delta_df),
    plot_label,
    y_label = expression(Delta ~ "Metric Percentage (%)"),
    y_limits = c(-80, 80),
    y_breaks = seq(-80, 80, 20),
    legend_pos = c(0.88, 0.85),
    zero_line = TRUE,
    y_expand = expansion(mult = c(0.02, 0.03))
  )
}

# Figure 5b: metrics per reference
plot_pangenome <- make_boxplot(
  combined,
  "T2T_pangenome",
  "VG[pangenome]"
)

plot_bwa <- make_boxplot(
  combined,
  "GRCh38_bwa",
  "BWA[GRCh38]"
)

plot_linear <- make_boxplot(
  combined,
  "GRCh38_linear",
  "VG[GRCh38]"
)

plot_pangenome
plot_bwa
plot_linear

# Figure 5c + 5d: pangenome deltas
delta_bwa <- calc_delta(combined, "T2T_pangenome", "GRCh38_bwa")
delta_linear <- calc_delta(combined, "T2T_pangenome", "GRCh38_linear")

summarise_delta(delta_bwa, "VG pangenome - BWA GRCh38")
summarise_delta(delta_linear, "VG pangenome - VG GRCh38")

plot_delta_bwa <- make_delta_boxplot(
  delta_bwa,
  "VG[pangenome] - BWA[GRCh38]"
)

plot_delta_linear <- make_delta_boxplot(
  delta_linear,
  "VG[pangenome] - VG[GRCh38]"
)

plot_delta_bwa
plot_delta_linear

# Save plots and Source Data
plots_to_save <- list(
  "LUAD_figure5b_variant_caller_boxplot_pangenome" = plot_pangenome,
  "LUAD_variant_caller_boxplot_grch38_bwa" = plot_bwa,
  "LUAD_variant_caller_boxplot_grch38_linear" = plot_linear,
  "LUAD_figure5c_variant_caller_boxplot_delta_pangenome_bwa" = plot_delta_bwa,
  "LUAD_figure5d_variant_caller_boxplot_delta_pangenome_linear" = plot_delta_linear
)
for (plot_name in names(plots_to_save)) {
  save_plot(plots_to_save[[plot_name]], plot_name, file.path("LUAD", if (grepl("figure5", plot_name)) "Figure5" else "Extras"))
}

write_source_data(to_long(filter(combined, Reference == "T2T_pangenome")), "Figure 5b")
write_source_data(to_long(delta_bwa), "Figure 5c")
write_source_data(to_long(delta_linear), "Figure 5d")
