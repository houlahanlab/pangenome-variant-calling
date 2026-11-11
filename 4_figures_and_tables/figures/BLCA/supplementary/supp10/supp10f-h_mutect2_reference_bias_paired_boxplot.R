# Supplementary Figure 10f-h
# Mutect2 version of the Figure 3d germline reference bias, measured at the same sites in every alignment: heterozygous
# germline SNVs within 1-74 bp of the Mutect2 VG pangenome TP + FP calls, normal depth > 8 in all three alignments, and
# the same reference base in GRCh38 and T2T-CHM13. Per donor, the median REF allele frequency per alignment, then:
#   10f: VG pangenome - VG GRCh38 (same aligner, only the reference differs), one-sided
#   10g: BWA GRCh38 - VG GRCh38 (same reference, only the aligner differs), two-sided
#   10h: VG pangenome - BWA GRCh38 (aligner and reference both differ), two-sided
# Donors with no usable germline sites are left out.
#
# Input:  <data_dir>/Reference_Bias/paired_sites_mutect2/germline_paired_counts.tsv (combine_paired_counts.py --caller mutect2)
# Output: <plots_dir>/Supp10/supp10{f,g,h}_*.{pdf,png}; Source Data
# Prints the donors per ancestry, the median REF allele frequencies and the tests for each difference.

# Setup
library(ggplot2)
library(dplyr)

source("4_figures/figure_helpers.R")

ancestry_colours <- c(
  "AA" = "#add8e6",
  "EA" = "#FEF5C3",
  "EAA" = "#fdc0cb"
)

ancestry_labels <- c(
  "AA" = "African",
  "EA" = "European",
  "EAA" = "East Asian"
)

alignments <- c("bwa_grch38", "vg_grch38", "vg_pangenome")

# Same axis for both plots
y_limits <- c(-0.1, 0.1)

# Read in files
paired_dir <- file.path(config$data_dir, "Reference_Bias", "paired_sites_mutect2")

patients <- read.table(
  "metadata/BLCA_Patients.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)$V1

# Useful functions
# Sites passing the depth filter in all three alignments, then the median REF allele frequency per patient and alignment
patient_ref_freq <- function(file) {
  sites <- read.delim(file, stringsAsFactors = FALSE) %>%
    filter(donor %in% patients)

  # Somatic sites carry tumour and normal depth, germline sites normal depth only
  passes_depth <- function(a) {
    if (paste0(a, "_tumour_depth") %in% names(sites)) {
      sites[[paste0(a, "_tumour_depth")]] >= 5 & sites[[paste0(a, "_normal_depth")]] > 8
    } else {
      sites[[paste0(a, "_normal_depth")]] > 8
    }
  }
  covered <- Reduce(`&`, lapply(alignments, passes_depth))
  cat(basename(file), ": ", sum(covered), " of ", nrow(sites), " sites pass the depth filter in all three alignments\n", sep = "")

  same_ref <- sites$grch38_ref == sites$t2t_ref
  covered <- covered & same_ref
  cat("  of those, ", sum(covered), " have the same reference base in GRCh38 and T2T (kept)\n", sep = "")

  sites[covered, ] %>%
    group_by(donor, ancestry) %>%
    summarise(
      N_Sites = n(),
      BWA_Ref_Freq = median(bwa_grch38_ref_freq, na.rm = TRUE),
      VG_GRCh38_Ref_Freq = median(vg_grch38_ref_freq, na.rm = TRUE),
      Pangenome_Ref_Freq = median(vg_pangenome_ref_freq, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Donor_ID = donor,
      Ancestry = factor(sub("_Patients", "", ancestry), levels = names(ancestry_labels)),
      Delta_Ref_Freq_Pangenome_VG_GRCh38 = Pangenome_Ref_Freq - VG_GRCh38_Ref_Freq,
      Delta_Ref_Freq_BWA_VG_GRCh38 = BWA_Ref_Freq - VG_GRCh38_Ref_Freq,
      Delta_Ref_Freq_Pangenome_BWA = Pangenome_Ref_Freq - BWA_Ref_Freq
    )
}

# Statistical tests, same as the original Figure 3 script
# first_col - other_col = delta_col; alternative = "less" tests whether first_col is lower (less reference bias)
run_delta_stats <- function(ref_freq, first_col, other_col, delta_col, label, alternative = "less") {
  cat("\n========================================\n")
  cat(label, "\n")
  cat("========================================\n")

  cat("\nOverall paired test (", alternative, "):\n", sep = "")
  print(wilcox.test(ref_freq[[first_col]], ref_freq[[other_col]], paired = TRUE, alternative = alternative, exact = FALSE))

  ancestry_tests <- ref_freq %>%
    group_by(Ancestry) %>%
    summarise(
      N = n(),
      Median_Delta = median(.data[[delta_col]]),
      P_Value = wilcox.test(.data[[delta_col]], mu = 0, alternative = alternative, exact = FALSE)$p.value,
      .groups = "drop"
    )
  cat("\nTests within each ancestry:\n")
  print(as.data.frame(ancestry_tests))

  cat("\nKruskal-Wallis test between ancestries:\n")
  print(kruskal.test(ref_freq[[delta_col]] ~ ref_freq$Ancestry))

  cat("\nPairwise ancestry comparisons (BH adjusted):\n")
  print(pairwise.wilcox.test(ref_freq[[delta_col]], ref_freq$Ancestry, p.adjust.method = "BH", exact = FALSE))

  cat("\nEast Asian vs European delta comparison (", alternative, "):\n", sep = "")
  print(wilcox.test(ref_freq[[delta_col]][ref_freq$Ancestry == "EAA"],
                    ref_freq[[delta_col]][ref_freq$Ancestry == "EA"],
                    alternative = alternative, exact = FALSE))

  cat("\nAfrican vs European delta comparison (", alternative, "):\n", sep = "")
  print(wilcox.test(ref_freq[[delta_col]][ref_freq$Ancestry == "AA"],
                    ref_freq[[delta_col]][ref_freq$Ancestry == "EA"],
                    alternative = alternative, exact = FALSE))
}

# legend_y: vertical position of the top of the legend (Supp 10 plots sit it lower, under the longer title)
# y_lims: y axis; the shared y_limits unless a plot has a patient outside it
make_delta_boxplot <- function(df, delta_col, plot_title, y_title, legend_y = 0.97, y_lims = y_limits) {
  if (any(df[[delta_col]] < y_lims[1] | df[[delta_col]] > y_lims[2])) {
    stop(delta_col, ": some deltas fall outside the y axis ", y_lims[1], " to ", y_lims[2])
  }

  set.seed(1)

  ggplot(df, aes(Ancestry, .data[[delta_col]], fill = Ancestry)) +
    stat_boxplot(
      geom = "errorbar",
      width = 0.3,
      colour = "black",
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    geom_hline(
      yintercept = 0,
      colour = "black",
      linewidth = 1
    ) +
    geom_boxplot(
      width = 0.6,
      outlier.shape = NA,
      colour = "black",
      linewidth = 0.6,
      key_glyph = "polygon"
    ) +
    geom_jitter(
      width = 0.1,
      height = 0,
      size = 2,
      shape = 21,
      fill = NA,
      alpha = 0.9,
      show.legend = FALSE
    ) +
    labs(
      x = NULL,
      y = y_title,
      fill = "Ancestry"
    ) +
    scale_x_discrete(labels = ancestry_labels) +
    scale_y_continuous(
      limits = y_lims,
      breaks = seq(-1, 1, 0.05)
    ) +
    scale_fill_manual(values = ancestry_colours, labels = ancestry_labels) +
    annotate(
      "text",
      x = 0.85, y = Inf,
      label = plot_title,
      parse = TRUE,
      size = 6,
      colour = "black",
      hjust = 0,
      vjust = 1.3,
      family = "Helvetica"
    ) +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "inside",
      legend.position.inside = c(0.97, legend_y),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
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
}

# Per patient REF allele frequency
germline <- patient_ref_freq(file.path(paired_dir, "germline_paired_counts.tsv"))

cat("\nGermline: ", nrow(germline), " of ", length(patients), " patients; left out: ",
    paste(setdiff(patients, germline$Donor_ID), collapse = ", "), "\n", sep = "")
print(table(germline$Ancestry))
cat("Median REF allele frequency: BWA GRCh38", round(median(germline$BWA_Ref_Freq), 3),
    " VG GRCh38", round(median(germline$VG_GRCh38_Ref_Freq), 3),
    " VG pangenome", round(median(germline$Pangenome_Ref_Freq), 3), "\n")

# Figure 3d: does the pangenome lower the REF allele frequency compared with the same aligner on GRCh38? (one-sided)
run_delta_stats(germline, "Pangenome_Ref_Freq", "VG_GRCh38_Ref_Freq", "Delta_Ref_Freq_Pangenome_VG_GRCh38",
                "Germline REF allele frequency: VG Pangenome vs VG GRCh38")
# Supp 10g: do the two aligners differ on the same GRCh38 reference? No expected direction, so two-sided
run_delta_stats(germline, "BWA_Ref_Freq", "VG_GRCh38_Ref_Freq", "Delta_Ref_Freq_BWA_VG_GRCh38",
                "Germline REF allele frequency: BWA GRCh38 vs VG GRCh38", alternative = "two.sided")
# Supp 10h: VG pangenome vs the standard BWA GRCh38 pipeline (aligner and reference both change), two-sided
run_delta_stats(germline, "Pangenome_Ref_Freq", "BWA_Ref_Freq", "Delta_Ref_Freq_Pangenome_BWA",
                "Germline REF allele frequency: VG Pangenome vs BWA GRCh38", alternative = "two.sided")

# Plots
germline_y_title <- expression(Delta ~ "Median Germline REF Allele Frequency")

plots <- list(
  "supp10f_mutect2_germline_ref_frequency_delta_vg_pangenome_vg_grch38" =
    make_delta_boxplot(germline, "Delta_Ref_Freq_Pangenome_VG_GRCh38", "Mutect2~VG[pangenome] - VG[GRCh38]", germline_y_title,
                       legend_y = 0.86),
  # Supp 10g
  "supp10g_mutect2_germline_ref_frequency_delta_bwa_grch38_vg_grch38" =
    make_delta_boxplot(germline, "Delta_Ref_Freq_BWA_VG_GRCh38", "Mutect2~BWA[GRCh38] - VG[GRCh38]", germline_y_title,
                       legend_y = 0.86),
  # Supp 10h
  "supp10h_mutect2_germline_ref_frequency_delta_vg_pangenome_bwa_grch38" =
    make_delta_boxplot(germline, "Delta_Ref_Freq_Pangenome_BWA", "Mutect2~VG[pangenome] - BWA[GRCh38]", germline_y_title,
                       legend_y = 0.86,
                       # TCGA-FD-A3B7 is at -0.117 (3 nearby sites with about twice the ALT reads in the pangenome)
                       y_lims = c(-0.15, 0.1))
)

for (p in plots) print(p)

# Save plots and Source Data
for (plot_name in names(plots)) {
  save_plot(plots[[plot_name]], plot_name, "Supp10", width = 6, height = 6)
  write_source_data(select(plots[[plot_name]]$data, -donor, -ancestry),
                    paste0("Supplementary Figure ", sub("^supp(10[fgh])_.*", "\\1", plot_name)))
}
