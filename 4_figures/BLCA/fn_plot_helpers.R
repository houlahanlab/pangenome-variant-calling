# Shared helpers for the false-negative boxplots (Figure 3e-f, Supplementary Figure 11): Strelka2 EVS, Mutect2 TLOD,
# VAF, coverage depth and mapping quality. Sourced after 4_figures/figure_helpers.R.
#
# Variant groups are defined by calling status:
#   BWA TP + VG TP = Found by both
#   BWA FN + VG FN = Missed by both
#   BWA TP + VG FN = Missed by VG pangenome
# Each point is one donor: the median of the metric over that donor's variants in the group, measured in each of the
# three alignments. Only variants with a measurement in all three are kept.

library(ggplot2)
library(dplyr)
library(tidyr)

vaf_file <- file.path(config$data_dir, "False_Negatives_Analysis", "VAF_results.tsv")

# Plots go to <plots_dir>/Supp11 unless the calling script sets out_dir
out_dir <- file.path(config$plots_dir, "Supp11")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

alignments <- c(
  "BWA_GRCh38_linear",
  "VG_GRCh38_linear",
  "VG_T2T_pangenome"
)

reference_colours <- c(
  "BWA_GRCh38_linear" = "#C0DEA4",
  "VG_GRCh38_linear" = "#228B22",
  "VG_T2T_pangenome" = "#72a643"
)

group_levels <- c(
  "Found by both",
  "Missed by both",
  "Missed by VG pangenome"
)

# ------------------------------------------------------------
# Read VAF_results.tsv (TP/FN status plus VAF and depth per BAM)
# ------------------------------------------------------------

read_vaf_results <- function() {
  read.delim(
    vaf_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      across(
        c(Variant_ID, DONOR_ID, TUMOUR_ID, TOOL, Defined_By, STATUS, Counted_In),
        ~ trimws(as.character(.x))
      ),
      TOOL = ifelse(TOOL == "Strelka2", "Strelka", TOOL),
      VAF = suppressWarnings(as.numeric(VAF)),
      MPILEUP_DEPTH = suppressWarnings(as.numeric(MPILEUP_DEPTH))
    ) %>%
    filter(
      TOOL %in% c("Strelka", "Mutect2"),
      STATUS %in% c("TP", "FN")
    )
}

# ------------------------------------------------------------
# Assign each variant to a group (BWA vs VG pangenome status)
# ------------------------------------------------------------

variant_groups <- function(vaf_data, caller) {

  status_data <- vaf_data %>%
    filter(
      TOOL == caller,
      Defined_By %in% c("BWA_GRCh38_linear", "VG_T2T_pangenome")
    ) %>%
    mutate(
      Defined_Set = ifelse(Defined_By == "BWA_GRCh38_linear", "BWA", "Pangenome")
    ) %>%
    distinct(Variant_ID, DONOR_ID, TUMOUR_ID, Defined_Set, STATUS)

  conflicts <- status_data %>%
    count(Variant_ID, DONOR_ID, TUMOUR_ID, Defined_Set) %>%
    filter(n > 1)

  if (nrow(conflicts) > 0) {
    stop(paste("Conflicting TP/FN statuses found for", caller))
  }

  status_data %>%
    pivot_wider(names_from = Defined_Set, values_from = STATUS) %>%
    mutate(
      Variant_Group = case_when(
        BWA == "TP" & Pangenome == "TP" ~ "Found by both",
        BWA == "FN" & Pangenome == "FN" ~ "Missed by both",
        BWA == "TP" & Pangenome == "FN" ~ "Missed by VG pangenome",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Variant_Group)) %>%
    select(Variant_ID, DONOR_ID, TUMOUR_ID, Variant_Group)
}

# ------------------------------------------------------------
# Join a metric onto the groups and take per-patient medians
#
# metric_long: one row per Variant_ID, DONOR_ID, Counted_In with
# a numeric Value column.
# ------------------------------------------------------------

patient_medians <- function(groups, metric_long, caller, metric_name) {

  metric_long <- metric_long %>%
    filter(Counted_In %in% alignments) %>%
    distinct(Variant_ID, DONOR_ID, Counted_In, .keep_all = TRUE)

  variant_metric <- groups %>%
    inner_join(metric_long, by = c("Variant_ID", "DONOR_ID")) %>%
    select(Variant_ID, DONOR_ID, TUMOUR_ID, Variant_Group, Counted_In, Value) %>%
    pivot_wider(names_from = Counted_In, values_from = Value)

  for (bam in alignments) {
    if (!bam %in% colnames(variant_metric)) {
      variant_metric[[bam]] <- NA_real_
    }
  }

  variant_metric <- variant_metric %>%
    filter(
      !is.na(BWA_GRCh38_linear),
      !is.na(VG_GRCh38_linear),
      !is.na(VG_T2T_pangenome)
    )

  medians <- variant_metric %>%
    pivot_longer(
      cols = all_of(alignments),
      names_to = "Counted_In",
      values_to = "Value"
    ) %>%
    group_by(DONOR_ID, TUMOUR_ID, Variant_Group, Counted_In) %>%
    summarise(
      Median_Value = median(Value, na.rm = TRUE),
      N_Variants = n_distinct(Variant_ID),
      .groups = "drop"
    ) %>%
    filter(is.finite(Median_Value)) %>%
    mutate(
      Variant_Group = factor(Variant_Group, levels = group_levels),
      Counted_In = factor(Counted_In, levels = alignments)
    )

  cat("\n====", caller, "-", metric_name, "====\n")
  cat("Variants with a measurement in all three BAMs:\n")
  print(as.data.frame(count(variant_metric, Variant_Group, name = "Variants")))

  cat("Patients per group:\n")
  print(as.data.frame(
    medians %>%
      group_by(Variant_Group) %>%
      summarise(Patients = n_distinct(DONOR_ID), .groups = "drop")
  ))

  medians
}

# ------------------------------------------------------------
# Paired Wilcoxon tests
# ------------------------------------------------------------

safe_paired_wilcox <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)

  if (sum(keep) < 2) {
    return(NA_real_)
  }

  wilcox.test(
    x[keep],
    y[keep],
    paired = TRUE,
    exact = FALSE
  )$p.value
}

paired_pvalues <- function(medians) {

  pvalues <- medians %>%
    select(
      DONOR_ID,
      TUMOUR_ID,
      Variant_Group,
      Counted_In,
      Median_Value
    ) %>%
    pivot_wider(
      names_from = Counted_In,
      values_from = Median_Value
    ) %>%
    filter(
      !is.na(BWA_GRCh38_linear),
      !is.na(VG_GRCh38_linear),
      !is.na(VG_T2T_pangenome)
    ) %>%
    group_by(Variant_Group) %>%
    summarise(
      N_Patients = n(),
      p_bwa_vs_vg38 = safe_paired_wilcox(
        BWA_GRCh38_linear,
        VG_GRCh38_linear
      ),
      p_bwa_vs_pangenome = safe_paired_wilcox(
        BWA_GRCh38_linear,
        VG_T2T_pangenome
      ),
      .groups = "drop"
    )

  cat("Paired Wilcoxon tests:\n")
  print(as.data.frame(pvalues), row.names = FALSE)

  pvalues
}

superscript_digits <- c(
  "0" = "\u2070",
  "1" = "\u00b9",
  "2" = "\u00b2",
  "3" = "\u00b3",
  "4" = "\u2074",
  "5" = "\u2075",
  "6" = "\u2076",
  "7" = "\u2077",
  "8" = "\u2078",
  "9" = "\u2079",
  "-" = "\u207b"
)

format_p <- function(p) {

  if (is.na(p)) {
    return("p = NA")
  }

  if (p < 0.001) {
    parts <- strsplit(
      formatC(p, format = "e", digits = 2),
      "e"
    )[[1]]

    exponent <- as.character(as.integer(parts[2]))

    return(
      paste0(
        "p = ",
        parts[1],
        " \u00d7 10",
        paste(
          superscript_digits[strsplit(exponent, "")[[1]]],
          collapse = ""
        )
      )
    )
  }

  paste0("p = ", signif(p, 2))
}

# ------------------------------------------------------------
# Plot and save as PNG and PDF
# ------------------------------------------------------------

make_plot <- function(data, caller_label, y_label, output_file,
                      y_max = NULL, y_breaks = NULL,
                      legend_position = c(0.82, 0.68),
                      pvalues = NULL) {

  max_value <- max(data$Median_Value, na.rm = TRUE)

  if (is.null(y_max)) {
    y_breaks <- pretty(c(0, max_value))
    y_breaks <- y_breaks[y_breaks <= max_value * 1.05]
    y_max <- max_value * 1.18
  }

  dodge <- position_dodge(width = 0.8)

  jitterdodge <- position_jitterdodge(
    jitter.width = 0.08,
    dodge.width = 0.8
  )

  set.seed(1)

  bracket_layer <- NULL

  if (!is.null(pvalues)) {

    centers <- setNames(
      seq_along(group_levels),
      group_levels
    )

    tick <- y_max * 0.012

    bracket_stats <- bind_rows(
      pvalues %>%
        transmute(
          x = centers[as.character(Variant_Group)] - 0.27,
          xend = centers[as.character(Variant_Group)],
          y = y_max * 0.80,
          label = sapply(p_bwa_vs_vg38, format_p)
        ),

      pvalues %>%
        transmute(
          x = centers[as.character(Variant_Group)] - 0.27,
          xend = centers[as.character(Variant_Group)] + 0.27,
          y = y_max * 0.85,
          label = sapply(p_bwa_vs_pangenome, format_p)
        )
    )

    bracket_layer <- list(
      geom_segment(
        data = bracket_stats,
        aes(
          x = x,
          xend = xend,
          y = y,
          yend = y
        ),
        inherit.aes = FALSE,
        colour = "black",
        linewidth = 0.6
      ),

      geom_segment(
        data = bracket_stats,
        aes(
          x = x,
          xend = x,
          y = y - tick,
          yend = y
        ),
        inherit.aes = FALSE,
        colour = "black",
        linewidth = 0.6
      ),

      geom_segment(
        data = bracket_stats,
        aes(
          x = xend,
          xend = xend,
          y = y - tick,
          yend = y
        ),
        inherit.aes = FALSE,
        colour = "black",
        linewidth = 0.6
      ),

      geom_text(
        data = bracket_stats,
        aes(
          x = (x + xend) / 2,
          y = y + y_max * 0.015,
          label = label
        ),
        inherit.aes = FALSE,
        size = 3.3
      )
    )
  }

  p <- ggplot(
    data,
    aes(
      x = Variant_Group,
      y = Median_Value,
      fill = Counted_In
    )
  ) +

    annotate(
      "rect",
      xmin = -Inf,
      xmax = 1.5,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.2,
      fill = "grey"
    ) +

    geom_vline(
      xintercept = 1.5,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.6
    ) +

    stat_boxplot(
      geom = "errorbar",
      width = 0.18,
      colour = "black",
      linewidth = 0.6,
      position = dodge,
      show.legend = FALSE
    ) +

    geom_boxplot(
      width = 0.65,
      outlier.shape = NA,
      colour = "black",
      linewidth = 0.6,
      position = dodge,
      show.legend = FALSE
    ) +

    geom_point(
      aes(group = Counted_In),
      size = 1.5,
      shape = 21,
      colour = "black",
      fill = NA,
      alpha = 0.8,
      position = jitterdodge,
      show.legend = FALSE
    ) +

    geom_point(
      aes(fill = Counted_In),
      shape = 22,
      size = 0,
      alpha = 0,
      show.legend = TRUE
    ) +

    bracket_layer +

    annotate(
      "text",
      x = 0.50,
      y = y_max * 0.986,
      label = caller_label,
      hjust = 0,
      vjust = 1,
      size = 6
    ) +

    annotate(
      "text",
      x = 1,
      y = y_max * 0.915,
      label = "Shared",
      size = 5
    ) +

    annotate(
      "text",
      x = 2.5,
      y = y_max * 0.915,
      label = "MC3-exclusive",
      size = 5
    ) +

    scale_x_discrete(
      drop = FALSE,
      labels = c(
        "Found by both" = "Found by both",
        "Missed by both" = "Missed by both",
        "Missed by VG pangenome" =
          expression(paste("Missed by ", VG[pangenome]))
      )
    ) +

    scale_fill_manual(
      name = "Alignment",
      values = reference_colours,
      breaks = alignments,
      labels = expression(
        BWA[GRCh38],
        VG[GRCh38],
        VG[pangenome]
      ),
      drop = FALSE
    ) +

    guides(
      fill = guide_legend(
        override.aes = list(
          shape = 22,
          size = 6,
          alpha = 1,
          colour = "black",
          linewidth = 0.6
        ),
        keyheight = grid::unit(0.65, "cm"),
        keywidth = grid::unit(0.65, "cm")
      )
    ) +

    scale_y_continuous(
      limits = c(0, y_max),
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0))
    ) +

    labs(
      x = NULL,
      y = y_label
    ) +

    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "inside",
      legend.position.inside = legend_position,
      legend.justification = c(0.5, 0.5),
      legend.background = element_rect(
        fill = "transparent",
        colour = NA
      ),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      legend.key.height = grid::unit(0.65, "cm"),
      legend.key.width = grid::unit(0.65, "cm"),
      legend.spacing.y = grid::unit(0.15, "cm"),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 14, colour = "black"),
      axis.text.y = element_text(size = 14, colour = "black"),
      axis.title.y = element_text(size = 20),
      axis.line = element_line(colour = "black"),
      panel.border = element_blank(),
      plot.margin = margin(
        t = 10,
        r = 10,
        b = 10,
        l = 10
      )
    )

  if (max_value > y_max) {
    warning(
      output_file,
      ": largest median (",
      max_value,
      ") is above y_max (",
      y_max,
      ")"
    )
  }

  ggsave(
    file.path(out_dir, paste0(output_file, ".png")),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )

  ggsave(
    file.path(out_dir, paste0(output_file, ".pdf")),
    p,
    width = 8,
    height = 6,
    # cairo_pdf: the basic pdf device cannot store the superscript p-value exponents (they printed as "10..")
    device = cairo_pdf,
    family = "Helvetica"
  )

  invisible(p)
}