# Supplementary Figure 7g-i
# Strelka2 VG pangenome - BWA GRCh38 F1 against the proportion of properly paired reads in the VG pangenome tumour
# sample (vg statistics, PROPER_PAIRS / TOTAL_READS), properly paired reads by ancestry, and the F1 difference by
# ancestry in donors above the median proportion (Wilcoxon East Asian vs African + European).
#
# Input:  <alignments_dir>/T2T_pangenome/2026-10-08_T2Thprc1.1_vg_alignment_statistics.tsv (from pangenome-pipeline)
#         Supplementary Table 1 (BLCA; Strelka_F1_Score)
# Output: <plots_dir>/Supp7/supp7g_* (scatterplot), supp7h_* (boxplot), supp7i_* (subset boxplot),
#         as PDF and PNG; Source Data
# Prints the donors per ancestry, the correlation and the Wilcoxon test.

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

supp7_theme <- theme(
  text = element_text(family = "Helvetica"),
  legend.position = "inside",
  legend.position.inside = c(0.85, 0.15),
  legend.background = element_rect(fill = NA, colour = NA),
  legend.key = element_rect(fill = NA, colour = NA),
  legend.title = element_text(size = 14),
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

# Read in files
patients <- read_patients()

alignment <- read.delim(
  file.path(config$alignments_dir, "T2T_pangenome", "2026-10-08_T2Thprc1.1_vg_alignment_statistics.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

strelka <- read_supp_table_1("BLCA") %>% transmute(Donor_ID, Ancestry, Reference, F1_Score = Strelka_F1_Score)

# Proportion of properly paired reads in each tumour sample (TCGA sample type 01)
tumour_paired <- alignment %>%
  mutate(
    Donor_ID = sub("^(TCGA-[^-]+-[^-]+)-.*", "\\1", Sample),
    Sample_Type = substr(sub("^TCGA-[^-]+-[^-]+-", "", Sample), 1, 2),
    Proportion_Paired = PROPER_PAIRS / TOTAL_READS
  ) %>%
  filter(Sample_Type == "01", Donor_ID %in% patients) %>%
  select(Donor_ID, Proportion_Paired)

if (anyDuplicated(tumour_paired$Donor_ID) > 0) {
  stop("More than one tumour sample for: ",
       paste(unique(tumour_paired$Donor_ID[duplicated(tumour_paired$Donor_ID)]), collapse = ", "))
}

# Delta F1, VG pangenome minus BWA GRCh38
pan <- strelka %>% filter(Reference == "T2T_pangenome") %>% select(Donor_ID, Ancestry, F1_Score)
bwa <- strelka %>% filter(Reference == "GRCh38_bwa") %>% select(Donor_ID, F1_Score)

paired <- inner_join(pan, bwa, by = "Donor_ID", suffix = c("_pangenome", "_bwa")) %>%
  filter(Donor_ID %in% patients) %>%
  mutate(Delta_F1_Score = F1_Score_pangenome - F1_Score_bwa) %>%
  inner_join(tumour_paired, by = "Donor_ID") %>%
  mutate(Ancestry = factor(Ancestry, levels = ancestry_levels))

if (nrow(paired) != length(patients)) {
  stop("Properly paired reads and delta F1 found for ", nrow(paired), " of the ", length(patients), " patients")
}

cat("Patients per ancestry:\n")
print(table(paired$Ancestry))
cat("Proportion properly paired range:", signif(range(paired$Proportion_Paired), 4), "\n")

# Delta F1 vs properly paired reads
plot_scatter <- ggplot(paired, aes(Proportion_Paired, Delta_F1_Score, fill = Ancestry)) +
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
    label = "VG[pangenome] - BWA[GRCh38]",
    parse = TRUE,
    size = 6,
    colour = "black",
    hjust = -0.1,
    vjust = 2,
    family = "Helvetica"
  ) +
  labs(
    x = "Proportion of Properly Paired Reads",
    y = expression(Delta ~ "F1 Score (%)")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.09))) +
  scale_fill_manual(values = ancestry_colours, name = NULL, labels = ancestry_labels) +
  supp7_theme

# Properly paired reads by ancestry
set.seed(1)

plot_paired <- ggplot(paired, aes(Ancestry, Proportion_Paired, fill = Ancestry)) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.3,
    colour = "black",
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 1,
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
    colour = "black",
    alpha = 0.9,
    show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = "Proportion of Properly Paired Reads"
  ) +
  scale_x_discrete(labels = ancestry_labels) +
  scale_fill_manual(values = ancestry_colours, name = "Ancestry", labels = ancestry_labels) +
  supp7_theme

# Samples above the median proportion of properly paired reads, EAA vs AA + EA
median_paired <- median(paired$Proportion_Paired)

above_median <- paired %>%
  filter(Proportion_Paired > median_paired) %>%
  mutate(Group = ifelse(Ancestry == "EAA_Patients", "EAA", "AA_EA"))

wilcox_result <- wilcox.test(Delta_F1_Score ~ Group, data = above_median, exact = FALSE)

cat("\n========================================\n")
cat("Median proportion properly paired:", median_paired, "\n")
cat("Samples above the median:", nrow(above_median), "\n")
cat("========================================\n")
print(table(above_median$Ancestry))
cat("\nWilcoxon EAA vs AA + EA, samples > median properly paired reads\n")
print(wilcox_result)

set.seed(1)

plot_subset <- ggplot(above_median, aes(Ancestry, Delta_F1_Score, fill = Ancestry)) +
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
    alpha = 1,
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
    colour = "black",
    alpha = 0.9,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = "VG[pangenome] - BWA[GRCh38]",
    parse = TRUE,
    size = 6,
    colour = "black",
    hjust = -0.1,
    vjust = 2,
    family = "Helvetica"
  ) +
  labs(
    x = NULL,
    y = expression(Delta ~ "F1 (%) in Samples > Median Paired Reads")
  ) +
  scale_x_discrete(
    labels = ancestry_labels,
    expand = expansion(mult = c(0.4, 0.3))
  ) +
  scale_fill_manual(values = ancestry_colours, name = "Ancestry", labels = ancestry_labels) +
  supp7_theme

plot_scatter
plot_paired
plot_subset

# Save plots and Source Data
save_plot(plot_scatter, "supp7g_proportion_of_properly_paired_reads_scatterplot", "Supp7", width = 8, height = 6)
save_plot(plot_paired, "supp7h_proportion_of_properly_paired_reads_boxplot", "Supp7", width = 5, height = 6)
save_plot(plot_subset, "supp7i_proportion_of_properly_paired_reads_subset_boxplot", "Supp7", width = 5, height = 6)
write_source_data(select(paired, Donor_ID, Ancestry, F1_Score_pangenome, F1_Score_bwa, Delta_F1_Score, Proportion_Paired),
                  "Supplementary Figure 7g-i")
