# Figure 2c-e, Supplementary Figures 5 and 6 (ancestry panels)
# Recall, precision and F1 per donor by ancestry for Strelka2, Mutect2, DeepSomatic and SomaticSniper.
#   Figure 2c: Strelka2 VG pangenome F1; Figure 2d / 2e: Strelka2 VG pangenome - BWA GRCh38 / - VG GRCh38 F1
#   Supplementary Figures 5 and 6: the other callers, BWA GRCh38 and the differences
# Plots without a panel prefix are saved too, but are not in a figure.
#
# Input:  Supplementary Table 1 (BLCA)
# Output: <plots_dir>/ancestry/*.{pdf,png}; Source Data for the figure2*, supp5* and supp6* plots
# Prints the median, mean and IQR per ancestry and the Kruskal-Wallis and pairwise Wilcoxon tests for each plot.

# Setup
library(ggplot2)
library(dplyr)

source("4_figures/figure_helpers.R")

ancestry_colours <- c(
  "AA_Patients" = "#add8e6",
  "EA_Patients" = "#FEF5C3",
  "EAA_Patients" = "#fdc0cb"
)

black_colours <- c(
  "AA_Patients" = "black",
  "EA_Patients" = "black",
  "EAA_Patients" = "black"
)

ancestry_labels <- c(
  "AA_Patients" = "African",
  "EA_Patients" = "European",
  "EAA_Patients" = "East Asian"
)

ancestry_levels <- c("AA_Patients", "EA_Patients", "EAA_Patients")

# Read in files
metrics <- supp_table_1_metrics("BLCA", c("Strelka", "Mutect2", "DeepSomatic", "SomaticSniper"))

callers <- list(
  "Strelka2" = filter(metrics, Variant_Caller == "Strelka"),
  "Mutect2" = filter(metrics, Variant_Caller == "Mutect2"),
  "DeepSomatic" = filter(metrics, Variant_Caller == "DeepSomatic"),
  "SomaticSniper" = filter(metrics, Variant_Caller == "SomaticSniper")
)

required_references <- c("T2T_pangenome", "GRCh38_bwa", "GRCh38_linear")

for (caller in names(callers)) {
  df <- callers[[caller]]

  missing_references <- setdiff(required_references, unique(df$Reference))
  if (length(missing_references) > 0) {
    stop("Missing references for ", caller, ": ", paste(missing_references, collapse = ", "),
         "\nAvailable references: ", paste(unique(df$Reference), collapse = ", "))
  }

  df <- df %>%
    filter(Reference %in% required_references)

  df$F1_Score <- as.numeric(df$F1_Score)
  df$Precision <- as.numeric(df$Precision)
  df$Recall <- as.numeric(df$Recall)
  df$Ancestry <- factor(df$Ancestry, levels = ancestry_levels)

  callers[[caller]] <- df
}

# Useful functions
# Shared boxplot design; long_df needs Ancestry, Metric, Value
# A single metric (F1 only) gets the whole panel shaded and no dashed dividers
# Legend sits below the plot; ancestry boxes fill the whole panel so an inside legend covers points
ancestry_boxplot <- function(long_df, plot_label, y_label, y_limits, y_breaks,
                             zero_line = FALSE,
                             y_expand = expansion(mult = 0.05)) {

  single_metric <- length(unique(long_df$Metric)) == 1

  dodge <- position_dodge(width = 0.9)
  jitterdodge <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9)

  set.seed(1)

  p <- ggplot(long_df, aes(Metric, Value, fill = Ancestry)) +
    annotate(
      "rect",
      xmin = -Inf, xmax = if (single_metric) Inf else 1.5,
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
    {if (zero_line) geom_hline(
      yintercept = 0,
      colour = "black",
      linewidth = 1
    )} +
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
      aes(colour = Ancestry),
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
      values = ancestry_colours,
      name = "Ancestry"
    ) +
    scale_colour_manual(
      values = black_colours,
      name = "Ancestry",
      labels = ancestry_labels
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(
          shape = 22,
          size = 6,
          fill = ancestry_colours,
          colour = "black"
        )
      ),
      fill = "none"
    ) +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "bottom",
      legend.margin = margin(0, 0, 0, 0),
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
    Ancestry = rep(df$Ancestry, times = 3),
    Metric = rep(c("F1_Score", "Precision", "Recall"), each = nrow(df)),
    Value = c(df$F1_Score, df$Precision, df$Recall)
  )

  long_df$Metric <- factor(
    long_df$Metric,
    levels = c("F1_Score", "Precision", "Recall")
  )

  long_df[!is.na(long_df$Value), ]
}

# ref_a minus ref_b, matched by donor
calc_delta <- function(df, ref_a, ref_b) {

  df_a <- df %>%
    filter(Reference == ref_a) %>%
    select(Donor_ID, Ancestry, Recall, Precision, F1_Score)

  df_b <- df %>%
    filter(Reference == ref_b) %>%
    select(Donor_ID, Recall, Precision, F1_Score)

  if (anyDuplicated(df_a$Donor_ID) > 0) {
    stop(paste("Duplicate donor rows found for", ref_a))
  }

  if (anyDuplicated(df_b$Donor_ID) > 0) {
    stop(paste("Duplicate donor rows found for", ref_b))
  }

  merged <- inner_join(df_a, df_b, by = "Donor_ID", suffix = c("_a", "_b"))

  if (nrow(merged) == 0) {
    stop(paste("No matching donors between", ref_a, "and", ref_b))
  }

  data.frame(
    Donor_ID = merged$Donor_ID,
    Ancestry = merged$Ancestry,
    F1_Score = merged$F1_Score_a - merged$F1_Score_b,
    Precision = merged$Precision_a - merged$Precision_b,
    Recall = merged$Recall_a - merged$Recall_b
  )
}

# Mean / median / SD / min / max of each metric, overall and by ancestry
print_stats <- function(df, title) {

  summarise_metrics <- function(x) {
    x %>%
      summarise(
        N = n(),
        across(
          c(F1_Score, Precision, Recall),
          list(
            Mean = ~ mean(.x, na.rm = TRUE),
            Median = ~ median(.x, na.rm = TRUE),
            SD = ~ sd(.x, na.rm = TRUE),
            Min = ~ min(.x, na.rm = TRUE),
            Max = ~ max(.x, na.rm = TRUE)
          ),
          .names = "{.fn}_{.col}"
        ),
        .groups = "drop"
      )
  }

  cat("\n========================================\n")
  cat(title, "\n")
  cat("========================================\n")

  cat("\nOverall:\n")
  print(as_tibble(summarise_metrics(df)), width = Inf)

  cat("\nBy ancestry:\n")
  print(summarise_metrics(df %>% group_by(Ancestry)), width = Inf)
}

make_boxplot <- function(df, plot_label, f1_only = FALSE) {
  long_df <- to_long(df)
  if (f1_only) long_df <- long_df %>% filter(Metric == "F1_Score")

  ancestry_boxplot(
    long_df,
    plot_label,
    y_label = if (f1_only) "F1 Score (%)" else "Metric Percentage (%)",
    y_limits = c(-1, 100),
    y_breaks = seq(0, 100, 20)
  )
}

make_delta_boxplot <- function(delta_df, plot_label, f1_only = FALSE) {
  long_df <- to_long(delta_df)
  if (f1_only) long_df <- long_df %>% filter(Metric == "F1_Score")

  ancestry_boxplot(
    long_df,
    plot_label,
    y_label = if (f1_only) expression(Delta ~ "F1 Score (%)") else expression(Delta ~ "Metric Percentage (%)"),
    y_limits = c(-80, 80),
    y_breaks = seq(-80, 80, 20),
    zero_line = TRUE,
    y_expand = expansion(mult = c(0.02, 0.03))
  )
}

# Deltas and stats per caller
deltas <- list()

for (caller in names(callers)) {
  df <- callers[[caller]]

  deltas[[caller]] <- list(
    bwa = calc_delta(df, "T2T_pangenome", "GRCh38_bwa"),
    vg = calc_delta(df, "T2T_pangenome", "GRCh38_linear"),
    vg_bwa = calc_delta(df, "GRCh38_linear", "GRCh38_bwa")
  )

  print_stats(df %>% filter(Reference == "T2T_pangenome"), paste(caller, "- T2T_pangenome"))
  print_stats(df %>% filter(Reference == "GRCh38_bwa"), paste(caller, "- GRCh38_bwa"))
  print_stats(deltas[[caller]]$bwa, paste(caller, "- VG pangenome minus BWA GRCh38"))
  print_stats(deltas[[caller]]$vg, paste(caller, "- VG pangenome minus VG GRCh38"))
}

# Plots
# One row per plot: caller, what to plot, F1 only or not, output file name
plot_specs <- data.frame(
  caller = character(),
  data = character(),
  f1_only = logical(),
  file = character()
)

add_spec <- function(caller, data, f1_only, file) {
  plot_specs <<- rbind(plot_specs, data.frame(caller = caller, data = data, f1_only = f1_only, file = file))
}

add_spec("Strelka2", "pangenome", FALSE, "BLCA_Strelka_pangenome")
add_spec("Strelka2", "pangenome", TRUE, "figure2c_BLCA_Strelka_pangenome_F1_only")
add_spec("Strelka2", "bwa_ref", FALSE, "supp5_BLCA_Strelka_BWA_GRCh38")
add_spec("Strelka2", "bwa", FALSE, "BLCA_Strelka_delta_pangenome_minus_bwa")
add_spec("Strelka2", "vg", FALSE, "BLCA_Strelka_delta_pangenome_minus_vg_linear")
add_spec("Strelka2", "bwa", TRUE, "figure2d_BLCA_Strelka_delta_F1_pangenome_minus_bwa")
add_spec("Strelka2", "vg", TRUE, "figure2e_BLCA_Strelka_delta_F1_pangenome_minus_vg_linear")

add_spec("Mutect2", "pangenome", FALSE, "supp5_BLCA_Mutect2_pangenome")
add_spec("Mutect2", "pangenome", TRUE, "BLCA_Mutect2_pangenome_F1_only")
add_spec("Mutect2", "bwa_ref", FALSE, "supp6_BLCA_Mutect2_BWA_GRCh38")
add_spec("Mutect2", "bwa", FALSE, "BLCA_Mutect2_delta_pangenome_minus_bwa")
add_spec("Mutect2", "vg", FALSE, "BLCA_Mutect2_delta_pangenome_minus_vg_linear")
add_spec("Mutect2", "bwa", TRUE, "supp6_BLCA_Mutect2_delta_F1_pangenome_minus_bwa")
add_spec("Mutect2", "vg", TRUE, "BLCA_Mutect2_delta_F1_pangenome_minus_vg_linear")
add_spec("Mutect2", "vg_bwa", TRUE, "supp6_BLCA_Mutect2_delta_F1_VG_GRCh38_minus_BWA_GRCh38")

add_spec("DeepSomatic", "pangenome", FALSE, "supp5_BLCA_DeepSomatic_pangenome")
add_spec("DeepSomatic", "pangenome", TRUE, "BLCA_DeepSomatic_pangenome_F1_only")
add_spec("DeepSomatic", "bwa_ref", FALSE, "supp6_BLCA_DeepSomatic_BWA_GRCh38")
add_spec("DeepSomatic", "bwa", FALSE, "BLCA_DeepSomatic_delta_pangenome_minus_bwa")
add_spec("DeepSomatic", "vg", FALSE, "BLCA_DeepSomatic_delta_pangenome_minus_vg_linear")
add_spec("DeepSomatic", "bwa", TRUE, "BLCA_DeepSomatic_delta_F1_pangenome_minus_bwa")

add_spec("SomaticSniper", "pangenome", FALSE, "supp5_BLCA_SomaticSniper_pangenome")
add_spec("SomaticSniper", "pangenome", TRUE, "BLCA_SomaticSniper_pangenome_F1_only")
add_spec("SomaticSniper", "bwa_ref", FALSE, "supp6_BLCA_SomaticSniper_BWA_GRCh38")
add_spec("SomaticSniper", "bwa", FALSE, "BLCA_SomaticSniper_delta_pangenome_minus_bwa")
add_spec("SomaticSniper", "vg", FALSE, "BLCA_SomaticSniper_delta_pangenome_minus_vg_linear")
add_spec("SomaticSniper", "bwa", TRUE, "supp6_BLCA_SomaticSniper_delta_F1_pangenome_minus_bwa")

plot_labels <- c(
  "pangenome" = "VG[pangenome]",
  "bwa_ref" = "BWA[GRCh38]",
  "bwa" = "VG[pangenome] - BWA[GRCh38]",
  "vg" = "VG[pangenome] - VG[GRCh38]",
  "vg_bwa" = "VG[GRCh38] - BWA[GRCh38]"
)

plots <- list()

for (i in seq_len(nrow(plot_specs))) {
  spec <- plot_specs[i, ]
  df <- callers[[spec$caller]]
  plot_label <- paste0(spec$caller, "~", plot_labels[[spec$data]])

  plots[[spec$file]] <- switch(
    spec$data,
    "pangenome" = make_boxplot(df %>% filter(Reference == "T2T_pangenome"), plot_label, spec$f1_only),
    "bwa_ref" = make_boxplot(df %>% filter(Reference == "GRCh38_bwa"), plot_label, spec$f1_only),
    make_delta_boxplot(deltas[[spec$caller]][[spec$data]], plot_label, spec$f1_only)
  )
}

for (p in plots) print(p)

# Save plots and Source Data
for (i in seq_len(nrow(plot_specs))) {
  spec <- plot_specs[i, ]
  plot_width <- if (spec$f1_only) 5.5 else 10
  plot_height <- if (spec$f1_only) 6 else 5
  save_plot(plots[[spec$file]], spec$file, "ancestry", width = plot_width, height = plot_height)

  if (grepl("^(figure2|supp5|supp6)", spec$file)) {
    panel <- sub("^figure2([a-z])_.*", "Figure 2\\1", spec$file)
    if (panel == spec$file) panel <- paste("Supplementary Figure", sub("^supp([0-9]+)_BLCA_", "\\1 ", spec$file))
    write_source_data(plots[[spec$file]]$data, panel)
  }
}
