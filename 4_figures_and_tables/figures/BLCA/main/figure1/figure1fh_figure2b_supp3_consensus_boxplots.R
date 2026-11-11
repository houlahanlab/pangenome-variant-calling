# Figure 1f, 1h, 2b and Supplementary Figure 3 (VG GRCh38 consensus panels)
# F1 of Strelka2, Mutect2 and the >= 3 caller consensus per donor, with paired Wilcoxon tests (two-sided).
#   Figure 2b: VG pangenome Strelka2, Mutect2 and consensus against MC3
#   Figure 1f / 1h: VG pangenome - BWA GRCh38 F1 of Strelka2 and Mutect2 against the BWA GRCh38 / VG pangenome
#                   consensus >= 3; bracket = Strelka2 vs Mutect2 difference
#   Supplementary Figure 3: VG GRCh38 Strelka2 and Mutect2 against each consensus
# Also saves the VG GRCh38+ALT plots and the plots not in a figure to <plots_dir>/extras/consensus.
#
# Input:  <data_dir>/Stats_Summary/{Strelka,Mutect2}/stats_table_*.txt, MC3_Consensus/stats_table_T2T_pangenome.txt
#         (compare_to_mc3.py, build_consensus.py); <data_dir>/Stats_Summary/{VG_Pangenome,BWA_GRCh38}_Consensus/
#         stats_table_*_vs_consensus_ge3.txt (compare_to_consensus.py)
# Output: <plots_dir>/Figure1, Figure2, Supp3, extras/consensus; Source Data for each figure panel
# Prints the donors per group and every paired Wilcoxon test.

# Setup
library(ggplot2)
library(dplyr)
library(tidyr)

source("4_figures/figure_helpers.R")

caller_colours <- c(
  "Mutect2" = "#9c8fff",
  "Strelka2" = "#f1a1ff",
  "Consensus_ge3" = "#6F97D8"
)

black_colours <- c(
  "Mutect2" = "black",
  "Strelka2" = "black",
  "Consensus_ge3" = "black"
)

mc3_order <- c("Strelka2", "Mutect2", "Consensus_ge3")
groundtruth_order <- c("Strelka2", "Mutect2")

mc3_labels <- c(
  "Mutect2" = "Mutect2",
  "Strelka2" = "Strelka2",
  "Consensus_ge3" = "Consensus ≥3"
)

groundtruth_labels <- c(
  "Mutect2" = "Mutect2",
  "Strelka2" = "Strelka2"
)

title_mc3_consensus <- "plain(MC3)[plain(Consensus)]"
title_vg_consensus <- "plain(VG)[plain(Consensus)]"
title_bwa_consensus <- "plain(BWA)[plain(Consensus)]"

stats_dir <- file.path(config$data_dir, "Stats_Summary")

plots_dir <- config$plots_dir

# Panel folder for each plot; anything not in a figure goes to extras/consensus
plot_folders <- c(
  "vg_t2t_pangenome_mutect_strelka_pangenome_consensus_boxplot" = "Figure1",
  "vg_t2t_pangenome_mutect_strelka_bwa_consensus_boxplot" = "Figure1",
  "figure2b_vg_t2t_pangenome_mutect_strelka_mc3_consensus_boxplot" = "Figure2",
  "vg_grch38_linear_mutect_strelka_pangenome_consensus_boxplot" = "Supp3",
  "vg_grch38_linear_mutect_strelka_bwa_consensus_boxplot" = "Supp3",
  "figure1h_vg_t2t_pangenome_minus_bwa_mutect_strelka_pangenome_consensus_delta_boxplot" = "Figure1",
  "figure1f_vg_t2t_pangenome_minus_bwa_mutect_strelka_bwa_consensus_delta_boxplot" = "Figure1"
)

mc3_plot_width <- 5.5
groundtruth_plot_width <- 4.5
plot_height <- 6
plot_dpi <- 300
y_limit <- 120
# Same delta axis for the VG and BWA consensus plots, so the 0 line sits at the same height
delta_y_limits <- c(-30, 70)

# Useful functions
read_stats <- function(file) {
  read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = TRUE)
}

rename_metric_columns <- function(df) {
  recall_column <- grep("^(Recall|ecall)", colnames(df), value = TRUE, ignore.case = TRUE)
  precision_column <- grep("^Precision", colnames(df), value = TRUE, ignore.case = TRUE)
  f1_column <- grep("^F1_Score", colnames(df), value = TRUE, ignore.case = TRUE)

  if (length(recall_column) == 0) {
    stop("No Recall column found. Columns are: ", paste(colnames(df), collapse = ", "))
  }
  if (length(precision_column) == 0) {
    stop("No Precision column found. Columns are: ", paste(colnames(df), collapse = ", "))
  }
  if (length(f1_column) == 0) {
    stop("No F1_Score column found. Columns are: ", paste(colnames(df), collapse = ", "))
  }

  colnames(df)[colnames(df) == recall_column[1]] <- "Recall"
  colnames(df)[colnames(df) == precision_column[1]] <- "Precision"
  colnames(df)[colnames(df) == f1_column[1]] <- "F1_Score"
  df
}

check_columns <- function(df, required_columns, table_name) {
  missing_columns <- setdiff(required_columns, colnames(df))
  if (length(missing_columns) > 0) {
    stop("Missing columns in ", table_name, ": ", paste(missing_columns, collapse = ", "))
  }
}

standardize_caller_names <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("Strelka", "Strelka2")] <- "Strelka2"
  x[x %in% c("Mutect", "Mutect2")] <- "Mutect2"
  x[x %in% c("Consensus_ge3", "Consensus >=3", "Consensus ≥3")] <- "Consensus_ge3"
  x
}

metrics_to_numeric <- function(df) {
  df %>%
    mutate(
      Recall = suppressWarnings(as.numeric(Recall)),
      Precision = suppressWarnings(as.numeric(Precision)),
      F1_Score = suppressWarnings(as.numeric(F1_Score))
    )
}

check_numeric_metrics <- function(df, table_name) {
  bad_rows <- df %>%
    filter(is.na(Recall) | is.na(Precision) | is.na(F1_Score))

  if (nrow(bad_rows) > 0) {
    print(bad_rows %>%
            select(any_of(c("Donor_ID", "Variant_Caller", "Reference", "Recall", "Precision", "F1_Score"))))
    stop("Some metric values could not be converted to numeric in ", table_name, ".")
  }
  df
}

check_plot_groups <- function(df, expected_groups, plot_name) {
  group_counts <- df %>%
    count(Variant_Caller, name = "N")

  cat("\n", plot_name, " group counts:\n", sep = "")
  print(group_counts)

  missing_groups <- setdiff(expected_groups, as.character(group_counts$Variant_Caller))
  if (length(missing_groups) > 0) {
    stop(plot_name, " is missing: ", paste(missing_groups, collapse = ", "))
  }
}

# Single caller table (MC3 section); Variant_Caller and Reference set explicitly
prepare_mc3_caller_table <- function(df, caller_name, plot_reference) {
  df <- rename_metric_columns(df)
  check_columns(df, c("Donor_ID", "Ancestry", "Recall", "Precision", "F1_Score"), caller_name)

  df <- df %>%
    mutate(Variant_Caller = caller_name, Reference = plot_reference) %>%
    metrics_to_numeric()

  df <- check_numeric_metrics(df, paste(caller_name, plot_reference))

  df %>%
    select(Donor_ID, Ancestry, Variant_Caller, Reference, F1_Score, Precision, Recall)
}

# MC3 consensus >=3 rows; same rows are reused for every reference
prepare_mc3_consensus_table <- function(df) {
  df <- rename_metric_columns(df)
  check_columns(df, c("Donor_ID", "Ancestry", "Variant_Caller", "Recall", "Precision", "F1_Score"),
                "consensus table")

  df <- df %>%
    mutate(Variant_Caller = standardize_caller_names(Variant_Caller)) %>%
    metrics_to_numeric() %>%
    filter(Variant_Caller == "Consensus_ge3")

  df <- check_numeric_metrics(df, "MC3 consensus table")

  df %>%
    select(Donor_ID, Ancestry, Variant_Caller, F1_Score, Precision, Recall)
}

# Mutect2 + Strelka2 evaluated against a consensus >=3 ground truth
prepare_groundtruth_table <- function(df, expected_reference) {
  df <- rename_metric_columns(df)
  check_columns(df, c("Donor_ID", "Ancestry", "Variant_Caller", "Reference", "Recall", "Precision", "F1_Score"),
                paste(expected_reference, "ground-truth table"))

  df <- df %>%
    mutate(
      Variant_Caller = standardize_caller_names(Variant_Caller),
      Reference = trimws(as.character(Reference))
    ) %>%
    metrics_to_numeric()

  unexpected_references <- df %>%
    filter(Reference != expected_reference) %>%
    distinct(Reference)

  if (nrow(unexpected_references) > 0) {
    warning("Unexpected Reference values in ", expected_reference, " table: ",
            paste(unexpected_references$Reference, collapse = ", "))
  }

  df <- df %>%
    filter(Variant_Caller %in% groundtruth_order)

  df <- check_numeric_metrics(df, paste(expected_reference, "ground-truth table"))

  df %>%
    select(Donor_ID, Ancestry, Variant_Caller, Reference, F1_Score, Precision, Recall)
}

p_label <- function(p) {
  if (is.na(p)) {
    return("italic(p) == NA")
  }

  if (p < 0.001) {
    parts <- strsplit(formatC(p, format = "e", digits = 2), "e")[[1]]
    coefficient <- as.numeric(parts[1])
    exponent <- as.integer(parts[2])
    return(paste0("italic(p) == ", sprintf("%.2f", coefficient), " %*% 10^{", exponent, "}"))
  }

  paste0("italic(p) == ", sprintf("%.3f", p))
}

# Paired Wilcoxon test on F1 score for each comparison; y = bracket height
make_paired_pvals <- function(df, comparisons, caller_order, y_step = 8) {

  y_max <- max(df$F1_Score, na.rm = TRUE)
  pvals <- data.frame()

  for (i in seq_along(comparisons)) {
    group_1 <- comparisons[[i]][1]
    group_2 <- comparisons[[i]][2]

    paired_data <- df %>%
      filter(Variant_Caller %in% c(group_1, group_2)) %>%
      select(Donor_ID, Variant_Caller, F1_Score) %>%
      distinct() %>%
      pivot_wider(names_from = Variant_Caller, values_from = F1_Score) %>%
      filter(!is.na(.data[[group_1]]), !is.na(.data[[group_2]]))

    if (nrow(paired_data) < 2) {
      stop("Not enough paired donors for ", group_1, " versus ", group_2, ".")
    }

    p_value <- wilcox.test(
      paired_data[[group_1]],
      paired_data[[group_2]],
      paired = TRUE,
      exact = FALSE
    )$p.value

    pvals <- bind_rows(
      pvals,
      data.frame(
        comparison = paste(group_1, "vs", group_2),
        x1 = which(caller_order == group_1),
        x2 = which(caller_order == group_2),
        y = y_max + y_step * i,
        label = p_label(p_value),
        p_value = p_value,
        N_paired = nrow(paired_data)
      )
    )
  }

  pvals
}

make_f1_boxplot <- function(df, caller_order, caller_labels, comparisons, plot_title) {

  pvals <- make_paired_pvals(df, comparisons, caller_order)

  if (max(pvals$y + 5) > y_limit) {
    warning("A p-value annotation exceeds the y-axis limit of ", y_limit,
            ". Increase y_limit if the annotation is clipped.")
  }

  set.seed(1)

  p <- ggplot(df, aes(Variant_Caller, F1_Score, fill = Variant_Caller)) +
    stat_boxplot(
      geom = "errorbar",
      width = 0.3,
      colour = "black",
      linewidth = 0.6
    ) +
    geom_boxplot(
      width = 0.6,
      outlier.shape = NA,
      alpha = 1,
      colour = "black",
      linewidth = 0.6
    ) +
    geom_jitter(
      aes(colour = Variant_Caller),
      size = 2,
      # horizontal jitter only, fixed seed so every saved file matches and F1 = 0 points stay on the plot
      position = position_jitter(width = 0.08, height = 0, seed = 1),
      shape = 21,
      fill = NA,
      alpha = 0.9
    ) +
    scale_fill_manual(values = caller_colours, drop = FALSE) +
    scale_colour_manual(values = black_colours, drop = FALSE) +
    scale_x_discrete(drop = FALSE, labels = caller_labels) +
    scale_y_continuous(
      limits = c(0, y_limit),
      breaks = seq(0, y_limit, by = 20),
      expand = expansion(mult = c(0, 0))
    ) +
    # p-value brackets
    geom_segment(
      data = pvals,
      aes(x = x1, xend = x2, y = y, yend = y),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.6
    ) +
    geom_segment(
      data = pvals,
      aes(x = x1, xend = x1, y = y, yend = y - 2),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.6
    ) +
    geom_segment(
      data = pvals,
      aes(x = x2, xend = x2, y = y, yend = y - 2),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.6
    ) +
    geom_text(
      data = pvals,
      aes(x = (x1 + x2) / 2, y = y + 2.2, label = label),
      inherit.aes = FALSE,
      size = 4.5,
      family = "Helvetica",
      parse = TRUE
    ) +
    annotate(
      "text",
      x = 0.55, y = Inf,
      label = plot_title,
      parse = TRUE,
      size = 6,
      colour = "black",
      hjust = 0,
      vjust = 2,
      family = "Helvetica"
    ) +
    labs(
      x = NULL,
      y = "F1 Score (%)"
    ) +
    guides(
      colour = "none",
      fill = "none"
    ) +
    theme(
      text = element_text(family = "Helvetica"),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 14, colour = "black", angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 16, colour = "black"),
      axis.title.y = element_text(size = 15),
      panel.border = element_blank(),
      axis.line.x.bottom = element_line(colour = "black"),
      axis.line.y.left = element_line(colour = "black"),
      plot.margin = margin(10, 12, 10, 10)
    )

  return(list(plot = p, pvals = pvals))
}

# Delta F1 boxplot, same style as make_f1_boxplot; bracket = paired Wilcoxon of Strelka2 vs Mutect2 delta F1
make_delta_boxplot <- function(df, pvals, caller_order, caller_labels, plot_title) {

  if (min(df$Delta_F1) < delta_y_limits[1] || max(pvals$y + 5) > delta_y_limits[2]) {
    stop("Delta F1 (", paste(round(range(df$Delta_F1), 1), collapse = " to "),
         ") or a p-value annotation falls outside the y axis ", delta_y_limits[1], " to ", delta_y_limits[2])
  }

  set.seed(1)

  ggplot(df, aes(Variant_Caller, Delta_F1, fill = Variant_Caller)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 1) +
    stat_boxplot(
      geom = "errorbar",
      width = 0.3,
      colour = "black",
      linewidth = 0.6
    ) +
    geom_boxplot(
      width = 0.6,
      outlier.shape = NA,
      alpha = 1,
      colour = "black",
      linewidth = 0.6
    ) +
    geom_jitter(
      aes(colour = Variant_Caller),
      size = 2,
      position = position_jitter(width = 0.08, height = 0, seed = 1),
      shape = 21,
      fill = NA,
      alpha = 0.9
    ) +
    scale_fill_manual(values = caller_colours, drop = FALSE) +
    scale_colour_manual(values = black_colours, drop = FALSE) +
    scale_x_discrete(drop = FALSE, labels = caller_labels) +
    scale_y_continuous(
      limits = delta_y_limits,
      breaks = seq(-20, 60, by = 20),
      expand = expansion(mult = c(0, 0))
    ) +
    # p-value bracket
    geom_segment(
      data = pvals,
      aes(x = x1, xend = x2, y = y, yend = y),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.6
    ) +
    geom_segment(
      data = pvals,
      aes(x = x1, xend = x1, y = y, yend = y - 2),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.6
    ) +
    geom_segment(
      data = pvals,
      aes(x = x2, xend = x2, y = y, yend = y - 2),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.6
    ) +
    geom_text(
      data = pvals,
      aes(x = (x1 + x2) / 2, y = y + 2.2, label = label),
      inherit.aes = FALSE,
      size = 4.5,
      family = "Helvetica",
      parse = TRUE
    ) +
    annotate(
      "text",
      x = 0.55, y = Inf,
      label = plot_title,
      parse = TRUE,
      size = 6,
      colour = "black",
      hjust = 0,
      vjust = 2,
      family = "Helvetica"
    ) +
    labs(
      x = NULL,
      y = expression(Delta ~ "F1 Score (%)")
    ) +
    guides(
      colour = "none",
      fill = "none"
    ) +
    theme(
      text = element_text(family = "Helvetica"),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 14, colour = "black", angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 16, colour = "black"),
      axis.title.y = element_text(size = 15),
      panel.border = element_blank(),
      axis.line.x.bottom = element_line(colour = "black"),
      axis.line.y.left = element_line(colour = "black"),
      plot.margin = margin(10, 12, 10, 10)
    )
}

save_plot_both <- function(plot, filename, width) {
  for (ext in c("pdf", "png")) {
    out_dir <- file.path(plots_dir, ifelse(filename %in% names(plot_folders),
                                           plot_folders[filename], "extras/consensus"))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(out_dir, paste0(filename, ".", ext)),
           plot, width = width, height = plot_height, dpi = plot_dpi)
  }
}

# Builds, prints and saves the VG GRCh38, VG GRCh38 plus and VG T2T pangenome plots for one section
# data_list / file_names are named by reference: linear, linear_plus, pangenome
plot_section <- function(data_list, caller_order, caller_labels, comparisons, plot_title,
                         section_name, file_names, width) {

  reference_names <- c(
    "linear" = "VG GRCh38",
    "linear_plus" = "VG GRCh38 Plus",
    "pangenome" = "VG T2T Pangenome"
  )

  results <- lapply(names(reference_names), function(ref) {
    make_f1_boxplot(data_list[[ref]], caller_order, caller_labels, comparisons, plot_title)
  })
  names(results) <- names(reference_names)

  for (ref in names(results)) {
    print(results[[ref]]$plot)
  }

  for (ref in names(results)) {
    cat("\n", section_name, " ", reference_names[[ref]], " paired Wilcoxon ",
        ifelse(length(comparisons) > 1, "tests", "test"), ":\n", sep = "")
    print(results[[ref]]$pvals)
  }

  for (ref in names(results)) {
    save_plot_both(results[[ref]]$plot, file_names[[ref]], width)
  }

  results
}

# Mutect2 and Strelka2 vs a consensus >=3 ground truth (VG pangenome or BWA GRCh38 consensus)
run_groundtruth_section <- function(consensus_dir, plot_title, section_name, check_name, file_suffix) {

  groundtruth <- list(
    linear = prepare_groundtruth_table(
      read_stats(file.path(consensus_dir, "stats_table_GRCh38_linear_vs_consensus_ge3.txt")),
      "GRCh38_linear"
    ),
    linear_plus = prepare_groundtruth_table(
      read_stats(file.path(consensus_dir, "stats_table_GRCh38_linear_plus_vs_consensus_ge3.txt")),
      "GRCh38_linear_plus"
    ),
    pangenome = prepare_groundtruth_table(
      read_stats(file.path(consensus_dir, "stats_table_T2T_pangenome_vs_consensus_ge3.txt")),
      "T2T_pangenome"
    )
  )

  for (ref in names(groundtruth)) {
    groundtruth[[ref]]$Variant_Caller <- factor(groundtruth[[ref]]$Variant_Caller, levels = groundtruth_order)
  }

  check_plot_groups(groundtruth$linear, groundtruth_order, paste0(check_name, ": VG GRCh38"))
  check_plot_groups(groundtruth$linear_plus, groundtruth_order, paste0(check_name, ": VG GRCh38 Plus"))
  check_plot_groups(groundtruth$pangenome, groundtruth_order, paste0(check_name, ": VG T2T Pangenome"))

  plot_section(
    groundtruth,
    groundtruth_order,
    groundtruth_labels,
    list(c("Mutect2", "Strelka2")),
    plot_title,
    section_name,
    file_names = c(
      "linear" = paste0("vg_grch38_linear_mutect_strelka_", file_suffix, "_consensus_boxplot"),
      "linear_plus" = paste0("vg_grch38_linear_plus_mutect_strelka_", file_suffix, "_consensus_boxplot"),
      "pangenome" = paste0("vg_t2t_pangenome_mutect_strelka_", file_suffix, "_consensus_boxplot")
    ),
    width = groundtruth_plot_width
  )
}

# MC3 comparisons
mc3_consensus_data <- prepare_mc3_consensus_table(
  read_stats(file.path(stats_dir, "MC3_Consensus/stats_table_T2T_pangenome.txt"))
)

mc3_files <- list(
  linear = c(
    "Mutect2" = "Mutect2/stats_table_grch38_linear.txt",
    "Strelka2" = "Strelka/stats_table_grch38_linear.txt"
  ),
  linear_plus = c(
    "Mutect2" = "Mutect2/stats_table_grch38_linear_plus.txt",
    "Strelka2" = "Strelka/stats_table_grch38_linear_plus.txt"
  ),
  pangenome = c(
    "Mutect2" = "Mutect2/stats_table_pangenome.txt",
    "Strelka2" = "Strelka/stats_table_pangenome.txt"
  )
)

mc3_references <- c(
  "linear" = "GRCh38_linear",
  "linear_plus" = "GRCh38_linear_plus",
  "pangenome" = "T2T_pangenome"
)

mc3 <- list()

for (ref in names(mc3_files)) {
  callers <- bind_rows(lapply(names(mc3_files[[ref]]), function(caller) {
    prepare_mc3_caller_table(
      read_stats(file.path(stats_dir, mc3_files[[ref]][[caller]])),
      caller,
      mc3_references[[ref]]
    )
  }))

  mc3[[ref]] <- bind_rows(callers, mc3_consensus_data %>% mutate(Reference = mc3_references[[ref]]))
  mc3[[ref]]$Variant_Caller <- factor(mc3[[ref]]$Variant_Caller, levels = mc3_order)
}

check_plot_groups(mc3$linear, mc3_order, "MC3: VG GRCh38")
check_plot_groups(mc3$linear_plus, mc3_order, "MC3: VG GRCh38 Plus")
check_plot_groups(mc3$pangenome, mc3_order, "MC3: VG T2T Pangenome")

mc3_results <- plot_section(
  mc3,
  mc3_order,
  mc3_labels,
  list(c("Mutect2", "Consensus_ge3"), c("Strelka2", "Consensus_ge3")),
  title_mc3_consensus,
  "MC3",
  file_names = c(
    "linear" = "vg_grch38_linear_mutect_strelka_mc3_consensus_boxplot",
    "linear_plus" = "vg_grch38_linear_plus_mutect_strelka_mc3_consensus_boxplot",
    "pangenome" = "figure2b_vg_t2t_pangenome_mutect_strelka_mc3_consensus_boxplot"
  ),
  width = mc3_plot_width
)

# Mutect2 and Strelka2 evaluated against VG consensus >=3
vg_results <- run_groundtruth_section(
  file.path(stats_dir, "VG_Pangenome_Consensus"),
  title_vg_consensus,
  "Consensus >=3 ground truth",
  "Consensus >=3 ground truth",
  "pangenome"
)

# Mutect2 and Strelka2 evaluated against BWA GRCh38 consensus >=3
bwa_results <- run_groundtruth_section(
  file.path(stats_dir, "BWA_GRCh38_Consensus"),
  title_bwa_consensus,
  "BWA Consensus >=3 ground truth",
  "Consensus >=3 ground truth",
  "bwa"
)

# Delta F1 (VG pangenome - BWA GRCh38) for Mutect2 and Strelka2 vs a consensus >=3 ground truth
run_delta_section <- function(consensus_dir, plot_title, section_name, file_suffix) {

  vg <- prepare_groundtruth_table(
    read_stats(file.path(consensus_dir, "stats_table_T2T_pangenome_vs_consensus_ge3.txt")),
    "T2T_pangenome"
  )
  bwa <- prepare_groundtruth_table(
    read_stats(file.path(consensus_dir, "stats_table_GRCh38_bwa_vs_consensus_ge3.txt")),
    "GRCh38_bwa"
  )

  delta <- inner_join(
    vg %>% select(Donor_ID, Variant_Caller, F1_VG = F1_Score),
    bwa %>% select(Donor_ID, Variant_Caller, F1_BWA = F1_Score),
    by = c("Donor_ID", "Variant_Caller")
  ) %>%
    mutate(
      Delta_F1 = F1_VG - F1_BWA,
      Variant_Caller = factor(Variant_Caller, levels = groundtruth_order)
    )

  check_plot_groups(delta, groundtruth_order, paste0(section_name, ": VG pangenome - BWA GRCh38"))

  # Printed only: each caller's VG pangenome vs BWA GRCh38 paired Wilcoxon test
  caller_pvals <- delta %>%
    group_by(Variant_Caller) %>%
    summarise(
      N_paired = n(),
      Median_Delta = median(Delta_F1),
      p_value = wilcox.test(F1_VG, F1_BWA, paired = TRUE, exact = FALSE)$p.value,
      .groups = "drop"
    )

  cat("\n", section_name, " VG pangenome vs BWA GRCh38 paired Wilcoxon tests:\n", sep = "")
  print(as.data.frame(caller_pvals))

  # Plotted bracket: Strelka2 vs Mutect2 delta F1, paired by donor
  paired_delta <- delta %>%
    select(Donor_ID, Variant_Caller, Delta_F1) %>%
    pivot_wider(names_from = Variant_Caller, values_from = Delta_F1)

  p_value <- wilcox.test(paired_delta$Strelka2, paired_delta$Mutect2, paired = TRUE, exact = FALSE)$p.value

  pvals <- data.frame(
    comparison = "Strelka2 vs Mutect2 delta F1",
    x1 = which(groundtruth_order == "Strelka2"),
    x2 = which(groundtruth_order == "Mutect2"),
    y = max(delta$Delta_F1) + 8,
    label = p_label(p_value),
    p_value = p_value,
    N_paired = nrow(paired_delta)
  )

  cat("\n", section_name, " Strelka2 vs Mutect2 delta F1 paired Wilcoxon test:\n", sep = "")
  print(pvals)

  p <- make_delta_boxplot(delta, pvals, groundtruth_order, groundtruth_labels, plot_title)
  print(p)
  # BWA consensus = Figure 1f, VG pangenome consensus = Figure 1h
  panel <- c("bwa" = "figure1f_", "pangenome" = "figure1h_")[[file_suffix]]
  save_plot_both(p, paste0(panel, "vg_t2t_pangenome_minus_bwa_mutect_strelka_", file_suffix, "_consensus_delta_boxplot"),
                 groundtruth_plot_width)

  list(plot = p, pvals = pvals)
}

vg_delta_results <- run_delta_section(
  file.path(stats_dir, "VG_Pangenome_Consensus"),
  title_vg_consensus,
  "Consensus >=3 ground truth",
  "pangenome"
)

bwa_delta_results <- run_delta_section(
  file.path(stats_dir, "BWA_GRCh38_Consensus"),
  title_bwa_consensus,
  "BWA Consensus >=3 ground truth",
  "bwa"
)

# Source Data: the values each figure panel draws
write_source_data(mc3_results$pangenome$plot$data, "Figure 2b")
write_source_data(bwa_delta_results$plot$data, "Figure 1f")
write_source_data(vg_delta_results$plot$data, "Figure 1h")
write_source_data(vg_results$linear$plot$data, "Supplementary Figure 3 VG GRCh38 vs VG pangenome consensus")
write_source_data(bwa_results$linear$plot$data, "Supplementary Figure 3 VG GRCh38 vs BWA GRCh38 consensus")
