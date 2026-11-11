# Supplementary Figure 3c
# Strelka2 and Mutect2 on the VG GRCh38 pangenome (HPRC v1.1 surjected onto GRCh38), as Figure 1b-d:
# recall, precision and F1 per donor; VG GRCh38 pangenome - BWA GRCh38; VG pangenome - VG GRCh38 pangenome (all
# metrics and F1 only). The VG GRCh38 pangenome is not in Supplementary Table 1, so the stats tables are read.
#
# Input:  <data_dir>/Stats_Summary/{Strelka,Mutect2}/mega_stats_table.txt (compare_to_mc3.py)
# Output: <plots_dir>/Supp3/supp3c_*.{pdf,png}; Source Data
# Prints the donors per alignment and caller, and the median and mean differences.

# Setup
library(ggplot2)
library(dplyr)

source("4_figures/figure_helpers.R")

caller_colours <- c(
  "Strelka" = "#f1a1ff",
  "Mutect2" = "#9c8fff"
)

black_colours <- c(
  "Strelka" = "black",
  "Mutect2" = "black"
)

rename_col <- function(df) {
  colnames(df)[colnames(df) == "Recall..."] <- "Recall"
  colnames(df)[colnames(df) == "Precision..."] <- "Precision"
  colnames(df)[colnames(df) == "F1_Score..."] <- "F1_Score"
  df
}

read_stats <- function(file) {
  df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df <- rename_col(df)
  df
}

# Read in files
patients <- read_patients()

strelka <- read_stats(file.path(config$data_dir, "Stats_Summary", "Strelka", "mega_stats_table.txt"))
mutect <- read_stats(file.path(config$data_dir, "Stats_Summary", "Mutect2", "mega_stats_table.txt"))

combined <- bind_rows(
  strelka,
  mutect
) %>%
  filter(Donor_ID %in% patients)

combined$F1_Score <- as.numeric(combined$F1_Score)
combined$Precision <- as.numeric(combined$Precision)
combined$Recall <- as.numeric(combined$Recall)

combined$Variant_Caller <- factor(
  combined$Variant_Caller,
  levels = c("Strelka", "Mutect2")
)

print(table(combined$Reference, combined$Variant_Caller))

# Useful functions
# Shared boxplot design; long_df needs Variant_Caller, Metric, Value
# A single metric (F1 only) gets the whole panel shaded and no dashed dividers
plot_metric_boxplot <- function(long_df, plot_label, y_label, y_limits, y_breaks,
                                legend_pos, zero_line = FALSE,
                                y_expand = expansion(mult = 0.05)) {

  single_metric <- length(unique(long_df$Metric)) == 1

  dodge <- position_dodge(width = 0.65)
  jitterdodge <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.65)

  set.seed(1)

  p <- ggplot(long_df, aes(Metric, Value, fill = Variant_Caller)) +
    annotate(
      "rect",
      xmin = -Inf, xmax = if (single_metric) Inf else 1.5,
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
      label = plot_label,
      parse = TRUE,
      size = 6,
      colour = "black",
      hjust = 0,
      vjust = 2,
      family = "Helvetica"
    )

  if (!single_metric) {
    p <- p +
      geom_vline(
        xintercept = c(1.5, 2.5),
        linetype = "dashed",
        colour = "grey",
        linewidth = 0.6
      )
  }

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

# ref_a minus ref_b, matched by donor and caller
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

make_delta_boxplot <- function(delta_df, plot_label, f1_only = FALSE) {
  long_df <- to_long(delta_df)
  if (f1_only) long_df <- long_df %>% filter(Metric == "F1_Score")

  plot_metric_boxplot(
    long_df,
    plot_label,
    y_label = expression(Delta ~ "Metric Percentage (%)"),
    y_limits = c(-80, 80),
    y_breaks = seq(-80, 80, 20),
    legend_pos = if (f1_only) c(0.82, 0.15) else c(0.88, 0.85),
    zero_line = TRUE,
    y_expand = expansion(mult = c(0.02, 0.03))
  )
}

# VG GRCh38 HPRC metrics
plot_hprc <- plot_metric_boxplot(
  to_long(combined %>% filter(Reference == "GRCh38_hprc_v1.1")),
  "VG[GRCh38~pangenome]",
  y_label = "Metric Percentage (%)",
  y_limits = c(-1, 100),
  y_breaks = seq(0, 100, 20),
  legend_pos = c(0.88, 0.19)
)

# Deltas
delta_hprc_bwa <- calc_delta(combined, "GRCh38_hprc_v1.1", "GRCh38_bwa")
delta_pan_hprc <- calc_delta(combined, "T2T_pangenome", "GRCh38_hprc_v1.1")

summarise_delta(delta_hprc_bwa, "VG GRCh38 HPRC - BWA GRCh38")
summarise_delta(delta_pan_hprc, "VG pangenome - VG GRCh38 HPRC")

plot_delta_hprc_bwa <- make_delta_boxplot(delta_hprc_bwa, "VG[GRCh38~pangenome] - BWA[GRCh38]")
plot_delta_pan_hprc <- make_delta_boxplot(delta_pan_hprc, "VG[pangenome] - VG[GRCh38~pangenome]")
plot_delta_pan_hprc_f1 <- make_delta_boxplot(delta_pan_hprc, "VG[pangenome] - VG[GRCh38~pangenome]", f1_only = TRUE)

plot_hprc
plot_delta_hprc_bwa
plot_delta_pan_hprc
plot_delta_pan_hprc_f1

# Save plots and Source Data
plots_to_save <- list(
  "supp3c_BLCA_GRCh38_HPRC_all_metrics" = plot_hprc,
  "supp3c_BLCA_GRCh38_HPRC_minus_BWA_GRCh38_all_metrics" = plot_delta_hprc_bwa,
  "supp3c_BLCA_VG_pangenome_minus_GRCh38_HPRC_all_metrics" = plot_delta_pan_hprc
)
for (plot_name in names(plots_to_save)) {
  save_plot(plots_to_save[[plot_name]], plot_name, "Supp3")
}
save_plot(plot_delta_pan_hprc_f1, "supp3c_F1_VG_pangenome_minus_GRCh38_HPRC", "Supp3", width = 5.5, height = 6)

write_source_data(to_long(combined %>% filter(Reference == "GRCh38_hprc_v1.1")), "Supplementary Figure 3c VG GRCh38 pangenome")
write_source_data(to_long(delta_hprc_bwa), "Supplementary Figure 3c VG GRCh38 pangenome minus BWA GRCh38")
write_source_data(to_long(delta_pan_hprc), "Supplementary Figure 3c VG pangenome minus VG GRCh38 pangenome")
