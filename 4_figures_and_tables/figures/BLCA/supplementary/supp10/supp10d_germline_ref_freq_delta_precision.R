# Supplementary Figure 10d
# Per donor change in median germline REF allele frequency (VG pangenome - BWA GRCh38) against the change in precision
# (VG pangenome - BWA GRCh38), with Spearman correlations overall and within each ancestry (a null result).
# 10d is Strelka2; the same plot for Mutect2 is saved too but is not in a figure. Germline sites are the heterozygous
# SNVs within 1-74 bp of each caller's VG pangenome TP + FP calls, with normal depth > 8 in all three alignments and the
# same reference base in GRCh38 and T2T-CHM13 (as in Figure 3d). Donors with no usable germline sites are left out.
#
# Input:  <data_dir>/Reference_Bias/paired_sites[_mutect2]/germline_paired_counts.tsv (combine_paired_counts.py)
#         Supplementary Table 1 (BLCA; Strelka_Precision, Mutect_Precision)
# Output: <plots_dir>/Supp10/supp10d_*.{pdf,png} and the Mutect2 plot; Source Data for 10d
# Prints the Spearman correlations for each caller.

# Setup
library(ggplot2)
library(dplyr)

source("4_figures/figure_helpers.R")

ancestry_colours <- c(
  "AA_Patients" = "#add8e6",
  "EA_Patients" = "#FEF5C3",
  "EAA_Patients" = "#fdc0cb"
)

ancestry_labels <- c(
  "AA_Patients" = "African",
  "EA_Patients" = "European",
  "EAA_Patients" = "East Asian"
)

ancestry_levels <- c("AA_Patients", "EA_Patients", "EAA_Patients")

alignments <- c("bwa_grch38", "vg_grch38", "vg_pangenome")

# Top at 55 so the title clears the highest Strelka2 patients (about 42)
y_limits <- c(-45, 55)

# Read in files
reference_bias_dir <- file.path(config$data_dir, "Reference_Bias")

callers <- list(
  "Strelka2" = list(paired_dir = "paired_sites", st1_caller = "Strelka",
                    plot_name = "supp10d_germline_ref_freq_delta_precision"),
  "Mutect2" = list(paired_dir = "paired_sites_mutect2", st1_caller = "Mutect2",
                   plot_name = "supp10_mutect2_germline_ref_freq_delta_precision")
)

patients <- read_patients()

# Useful functions
# Per patient delta Precision and delta median germline REF allele frequency (VG pangenome - BWA GRCh38) for one caller
delta_data <- function(caller) {
  info <- callers[[caller]]

  stats <- supp_table_1_metrics("BLCA", info$st1_caller)

  pan <- stats %>% filter(Reference == "T2T_pangenome") %>% select(Donor_ID, Ancestry, Precision)
  bwa <- stats %>% filter(Reference == "GRCh38_bwa") %>% select(Donor_ID, Precision)

  delta_precision <- inner_join(pan, bwa, by = "Donor_ID", suffix = c("_pangenome", "_bwa")) %>%
    filter(Donor_ID %in% patients) %>%
    mutate(Delta_Precision = Precision_pangenome - Precision_bwa)

  germline <- read.delim(file.path(reference_bias_dir, info$paired_dir, "germline_paired_counts.tsv"), stringsAsFactors = FALSE) %>%
    filter(donor %in% patients)

  # Sites with normal depth > 8 in all three alignments and the same reference base in GRCh38 and T2T, as in Figure 3d
  covered <- Reduce(`&`, lapply(alignments, function(a) germline[[paste0(a, "_normal_depth")]] > 8)) &
    germline$grch38_ref == germline$t2t_ref

  delta_ref_freq <- germline[covered, ] %>%
    group_by(donor) %>%
    summarise(
      N_Sites = n(),
      Delta_Ref_Freq = median(vg_pangenome_ref_freq, na.rm = TRUE) - median(bwa_grch38_ref_freq, na.rm = TRUE),
      .groups = "drop"
    )

  plot_data <- inner_join(delta_precision, delta_ref_freq, by = c("Donor_ID" = "donor"))
  plot_data$Ancestry <- factor(plot_data$Ancestry, levels = ancestry_levels)

  cat("\n========================================\n")
  cat(caller, "\n")
  cat("========================================\n")
  cat("Germline sites kept:", sum(covered), "of", nrow(germline), "\n")
  cat("Patients plotted:", nrow(plot_data), "of", length(patients), "\n")
  cat("No usable germline sites (left out):", paste(setdiff(patients, plot_data$Donor_ID), collapse = ", "), "\n")
  print(table(plot_data$Ancestry))
  cat("Delta median REF allele frequency range:", round(range(plot_data$Delta_Ref_Freq), 3), "\n")
  cat("Delta Precision range:", round(range(plot_data$Delta_Precision), 2), "\n")

  if (any(plot_data$Delta_Precision < y_limits[1] | plot_data$Delta_Precision > y_limits[2])) {
    stop(caller, ": some delta Precision values fall outside the y axis ", y_limits[1], " to ", y_limits[2])
  }

  plot_data
}

# Spearman correlation overall and within each ancestry
print_correlations <- function(df) {
  cat("\nDelta median germline REF allele frequency vs delta Precision (Spearman):\n")

  groups <- c(list(Overall = df), split(df, df$Ancestry))

  for (group in names(groups)) {
    x <- groups[[group]]
    corr <- cor.test(x$Delta_Ref_Freq, x$Delta_Precision, method = "spearman", exact = FALSE)
    cat(sprintf("%-13s n = %2d  rho = %6.3f  p = %.3g\n", group, nrow(x), corr$estimate, corr$p.value))
  }
}

ref_freq_scatterplot <- function(df, caller) {
  ggplot(df, aes(Delta_Ref_Freq, Delta_Precision, fill = Ancestry)) +
    geom_hline(
      yintercept = 0,
      colour = "black",
      linewidth = 1
    ) +
    geom_point(
      size = 3,
      shape = 21,
      colour = "black"
    ) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste0(caller, "~VG[pangenome] - BWA[GRCh38]"),
      parse = TRUE,
      size = 6,
      colour = "black",
      hjust = -0.1,
      vjust = 2,
      family = "Helvetica"
    ) +
    labs(
      x = expression(Delta ~ "Median Germline REF Allele Frequency"),
      y = expression(Delta ~ "Precision (%)")
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = seq(-40, 40, 20)
    ) +
    scale_fill_manual(
      values = ancestry_colours,
      name = NULL,
      labels = ancestry_labels
    ) +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "inside",
      legend.position.inside = c(0.85, 0.15),
      legend.background = element_rect(fill = NA, colour = NA),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.text = element_text(size = 14),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 15, colour = "black"),
      axis.text.y = element_text(size = 15, colour = "black"),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      panel.border = element_blank(),
      axis.line.x.bottom = element_line(colour = "black"),
      axis.line.y.left = element_line(colour = "black")
    )
}

# Plots
plots <- list()
for (caller in names(callers)) {
  plot_data <- delta_data(caller)
  print_correlations(plot_data)
  plots[[callers[[caller]]$plot_name]] <- ref_freq_scatterplot(plot_data, caller)
}

for (p in plots) print(p)

# Save plots and Source Data
for (plot_name in names(plots)) {
  save_plot(plots[[plot_name]], plot_name, "Supp10", width = 6, height = 6)  # same size as the Supp 10 boxplots
}
write_source_data(plots[["supp10d_germline_ref_freq_delta_precision"]]$data, "Supplementary Figure 10d")
