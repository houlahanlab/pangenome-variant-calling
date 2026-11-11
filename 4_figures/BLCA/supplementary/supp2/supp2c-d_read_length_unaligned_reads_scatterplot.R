# Supplementary Figure 2c-d
# Mean aligned read length against % unaligned reads in blood samples (45 BLCA donors), coloured by ancestry:
# 2c VG pangenome (projected BAM; the vg statistics have no read length), 2d BWA GRCh38, on the same axes.
# Spearman correlations for each alignment and for the VG pangenome - BWA GRCh38 differences are printed and saved.
#
# Input:  <alignments_dir>/T2T_pangenome/2026-10-08_T2Thprc1.1_gatk_alignment_statistics.tsv,
#         <alignments_dir>/GRCh38_bwa/2026-10-09_GRCh38bwa_gatk_alignment_statistics.tsv (from pangenome-pipeline)
# Output: <plots_dir>/Supp2/supp2{c,d}_*.{pdf,png}, read_length_vs_unaligned_reads_correlations.txt; Source Data

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

alignment_dir <- config$alignments_dir
out_dir <- file.path(config$plots_dir, "Supp2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Read in files
patients <- read.table(
  "metadata/BLCA_Patients.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)
patients <- data.frame(Donor_ID = patients$V1, Ancestry = patients$V4)

pangenome <- read.delim(
  file.path(alignment_dir, "T2T_pangenome/2026-10-08_T2Thprc1.1_gatk_alignment_statistics.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

bwa <- read.delim(
  file.path(alignment_dir, "GRCh38_bwa/2026-10-09_GRCh38bwa_gatk_alignment_statistics.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Useful functions
# Blood samples only (TCGA sample type 10/11), unaligned reads as a percentage of all reads
prepare_blood <- function(data) {
  data %>%
    mutate(
      Donor_ID = sub("^(TCGA-[^-]+-[^-]+)-.*", "\\1", SAMPLE),
      Sample_Type = substr(sub("^TCGA-[^-]+-[^-]+-", "", SAMPLE), 1, 2),
      Read_Length = MEAN_ALIGNED_READ_LENGTH,
      Unaligned_Percent = 100 * UNALIGNED_READS / TOTAL_READS
    ) %>%
    filter(Sample_Type %in% c("10", "11"), Donor_ID %in% patients$Donor_ID) %>%
    left_join(patients, by = "Donor_ID") %>%
    mutate(Ancestry = factor(Ancestry, levels = ancestry_levels)) %>%
    select(SAMPLE, Donor_ID, Ancestry, Read_Length, Unaligned_Percent)
}

pangenome_blood <- prepare_blood(pangenome)
bwa_blood <- prepare_blood(bwa)

for (df_name in c("pangenome_blood", "bwa_blood")) {
  df <- get(df_name)
  if (nrow(df) != nrow(patients) || anyDuplicated(df$Donor_ID) > 0) {
    stop(df_name, ": ", nrow(df), " blood samples for the ", nrow(patients), " patients")
  }
}

cat("Blood samples per ancestry:\n")
print(table(pangenome_blood$Ancestry))

# Spearman correlations
deltas <- inner_join(pangenome_blood, bwa_blood, by = c("SAMPLE", "Donor_ID", "Ancestry"),
                     suffix = c("_Pangenome", "_BWA")) %>%
  mutate(
    Delta_Read_Length = Read_Length_Pangenome - Read_Length_BWA,
    Delta_Unaligned = Unaligned_Percent_Pangenome - Unaligned_Percent_BWA
  )

spearman <- function(x, y) cor.test(x, y, method = "spearman", exact = FALSE)

correlations <- list(
  "VG pangenome" = spearman(pangenome_blood$Read_Length, pangenome_blood$Unaligned_Percent),
  "BWA GRCh38" = spearman(bwa_blood$Read_Length, bwa_blood$Unaligned_Percent),
  "VG pangenome - BWA GRCh38" = spearman(deltas$Delta_Read_Length, deltas$Delta_Unaligned)
)

correlation_table <- data.frame(
  Comparison = names(correlations),
  Correlation = sapply(correlations, function(x) unname(x$estimate)),
  P_value = sapply(correlations, function(x) x$p.value)
)

cat("\nSpearman correlations, mean aligned read length vs unaligned reads (blood):\n")
print(correlation_table, row.names = FALSE)

write.table(
  correlation_table,
  file.path(out_dir, "read_length_vs_unaligned_reads_correlations.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Shared axes for both panels
x_limits <- c(floor(min(pangenome_blood$Read_Length, bwa_blood$Read_Length)),
              ceiling(max(pangenome_blood$Read_Length, bwa_blood$Read_Length)))
y_limits <- c(0, ceiling(max(pangenome_blood$Unaligned_Percent, bwa_blood$Unaligned_Percent) / 2) * 2)
cat("\nx axis:", x_limits, " y axis:", y_limits, "\n")

scatterplot <- function(data, title) {
  ggplot(data, aes(Read_Length, Unaligned_Percent, fill = Ancestry)) +
    geom_point(
      size = 2.5,
      shape = 21,
      colour = "black",
      stroke = 0.5
    ) +
    # Title in the top left
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = title,
      parse = TRUE,
      size = 5.5,
      colour = "black",
      hjust = -0.1,
      vjust = 1.5,
      family = "Helvetica"
    ) +
    labs(
      x = "Mean Aligned Read Length (bp)",
      y = "Unaligned Reads (%)",
      fill = "Ancestry"
    ) +
    # Ticks on even read lengths, as in the original plot
    scale_x_continuous(limits = x_limits, breaks = seq(ceiling(x_limits[1] / 2) * 2, x_limits[2], by = 2)) +
    scale_y_continuous(limits = y_limits, breaks = seq(y_limits[1], y_limits[2], by = 2)) +
    scale_fill_manual(values = ancestry_colours, labels = ancestry_labels, drop = FALSE) +
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "inside",
      legend.position.inside = c(0.97, 0.6),
      legend.justification = c(1, 0.5),
      legend.background = element_rect(fill = NA, colour = NA),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.title = element_text(size = 12.5),
      legend.text = element_text(size = 12.5),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 11.5, colour = "black"),
      axis.text.y = element_text(size = 11.5, colour = "black"),
      axis.title.x = element_text(size = 12.5, margin = margin(t = 6)),
      axis.title.y = element_text(size = 12.5, margin = margin(r = 8)),
      panel.border = element_blank(),
      axis.line.x.bottom = element_line(colour = "black"),
      axis.line.y.left = element_line(colour = "black"),
      plot.margin = margin(t = 5.5, r = 8, b = 5.5, l = 5.5)
    )
}

# Plots

plots <- list(
  "supp2c_read_length_vs_unaligned_reads_blood_vg_pangenome" = scatterplot(pangenome_blood, "VG[pangenome]"),
  "supp2d_read_length_vs_unaligned_reads_blood_bwa_grch38" = scatterplot(bwa_blood, "BWA[GRCh38]")
)

for (plot_name in names(plots)) {
  print(plots[[plot_name]])
  for (ext in c("pdf", "png")) {
    ggsave(file.path(out_dir, paste0(plot_name, ".", ext)),
           plots[[plot_name]], width = 3.5, height = 4.7, dpi = 300, bg = "white")
  }
}

write_source_data(bind_rows(mutate(pangenome_blood, Alignment = "VG_Pangenome"), mutate(bwa_blood, Alignment = "BWA_GRCh38")),
                  "Supplementary Figure 2c-d")
