# Supplementary Figure 6d
# Mutect2 VG GRCh38 pangenome (HPRC v1.1 surjected onto GRCh38) - BWA GRCh38 F1 per donor by ancestry (6d), as Figure 2.
# Also saves the other VG GRCh38 pangenome plots (Strelka2 and Mutect2: VG GRCh38 pangenome, - BWA GRCh38,
# - VG pangenome), which are not in a figure. The VG GRCh38 pangenome is not in Supplementary Table 1, so the stats
# tables are read.
#
# Input:  <data_dir>/Stats_Summary/{Strelka,Mutect2}/mega_stats_table.txt (compare_to_mc3.py)
# Output: <plots_dir>/Supp6/supp6d_* and the other plots (*.pdf, *.png); Source Data for 6d
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
stats_dir <- file.path(config$data_dir, "Stats_Summary")

patients <- read.table(
  "metadata/BLCA_Patients.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)$V1

callers <- list(
  "Strelka2" = read_stats(file.path(stats_dir, "Strelka/mega_stats_table.txt")),
  "Mutect2" = read_stats(file.path(stats_dir, "Mutect2/mega_stats_table.txt"))
)

required_references <- c("GRCh38_hprc_v1.1", "GRCh38_bwa", "T2T_pangenome")

for (caller in names(callers)) {
  df <- callers[[caller]] %>%
    filter(Reference %in% required_references, Donor_ID %in% patients)

  missing_references <- setdiff(required_references, unique(df$Reference))
  if (length(missing_references) > 0) {
    stop("Missing references for ", caller, ": ", paste(missing_references, collapse = ", "))
  }

  df$F1_Score <- as.numeric(df$F1_Score)
  df$Precision <- as.numeric(df$Precision)
  df$Recall <- as.numeric(df$Recall)
  df$Ancestry <- factor(df$Ancestry, levels = ancestry_levels)

  cat("\n", caller, "patients per reference:\n")
  print(table(df$Reference, df$Ancestry))

  callers[[caller]] <- df
}

# Useful functions
# Shared boxplot design; long_df needs Ancestry, Metric, Value
# A single metric (F1 only) gets the whole panel shaded and no dashed dividers
ancestry_boxplot <- function(long_df, plot_label, y_label, y_limits, y_breaks,
                             legend_pos, zero_line = FALSE,
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

  if (anyDuplicated(df_a$Donor_ID) > 0 || anyDuplicated(df_b$Donor_ID) > 0) {
    stop(paste("Duplicate donor rows found for", ref_a, "or", ref_b))
  }

  merged <- inner_join(df_a, df_b, by = "Donor_ID", suffix = c("_a", "_b"))

  if (nrow(merged) != length(patients)) {
    stop(ref_a, " - ", ref_b, ": only ", nrow(merged), " of ", length(patients), " patients matched")
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

make_boxplot <- function(df, plot_label) {
  ancestry_boxplot(
    to_long(df),
    plot_label,
    y_label = "Metric Percentage (%)",
    y_limits = c(-1, 100),
    y_breaks = seq(0, 100, 20),
    legend_pos = c(0.88, 0.15)
  )
}

make_delta_boxplot <- function(delta_df, plot_label, f1_only = FALSE) {
  long_df <- to_long(delta_df)
  if (f1_only) long_df <- long_df %>% filter(Metric == "F1_Score")

  ancestry_boxplot(
    long_df,
    plot_label,
    y_label = expression(Delta ~ "Metric Percentage (%)"),
    y_limits = c(-80, 80),
    y_breaks = seq(-80, 80, 20),
    legend_pos = if (f1_only) c(0.82, 0.15) else c(0.88, 0.15),
    zero_line = TRUE,
    y_expand = expansion(mult = c(0.02, 0.03))
  )
}

# Plots
plots <- list()
plot_widths <- c()

for (caller in names(callers)) {
  df <- callers[[caller]]
  file_caller <- if (caller == "Strelka2") "Strelka" else caller

  hprc <- df %>% filter(Reference == "GRCh38_hprc_v1.1")
  delta_bwa <- calc_delta(df, "GRCh38_hprc_v1.1", "GRCh38_bwa")
  delta_pan <- calc_delta(df, "GRCh38_hprc_v1.1", "T2T_pangenome")

  print_stats(hprc, paste(caller, "- VG GRCh38 HPRC"))
  print_stats(delta_bwa, paste(caller, "- VG GRCh38 HPRC minus BWA GRCh38"))
  print_stats(delta_pan, paste(caller, "- VG GRCh38 HPRC minus VG pangenome"))

  plots[[paste0("supp6_", file_caller, "_grch38_hprc")]] <-
    make_boxplot(hprc, paste0(caller, "~VG[GRCh38~pangenome]"))

  plots[[paste0("supp6_", file_caller, "_delta_grch38_hprc_minus_bwa")]] <-
    make_delta_boxplot(delta_bwa, paste0(caller, "~VG[GRCh38~pangenome] - BWA[GRCh38]"))

  plots[[paste0("supp6_", file_caller, "_delta_grch38_hprc_minus_vg_pangenome")]] <-
    make_delta_boxplot(delta_pan, paste0(caller, "~VG[GRCh38~pangenome] - VG[pangenome]"))

  if (caller == "Mutect2") {
    plots[["supp6d_Mutect2_delta_F1_grch38_hprc_minus_bwa"]] <-
      make_delta_boxplot(delta_bwa, "Mutect2~VG[GRCh38~pangenome] - BWA[GRCh38]", f1_only = TRUE)
  }
}

for (p in plots) print(p)

# Save plots and Source Data
for (plot_name in names(plots)) {
  f1_only <- grepl("_F1_", plot_name)
  save_plot(plots[[plot_name]], plot_name, "Supp6", width = if (f1_only) 5.5 else 10, height = if (f1_only) 6 else 5)
}

write_source_data(plots[["supp6d_Mutect2_delta_F1_grch38_hprc_minus_bwa"]]$data, "Supplementary Figure 6d")
