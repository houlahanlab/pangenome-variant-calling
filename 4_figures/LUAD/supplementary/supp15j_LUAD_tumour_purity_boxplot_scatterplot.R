# Supplementary Figure 15j (LUAD)
# Strelka2 VG pangenome - BWA GRCh38 F1 against tumour purity (ABSOLUTE), tumour purity by ancestry, and the F1
# difference by ancestry in donors at or above the median purity (LUAD donors with a purity estimate; as
# Supplementary Figure 7d-f).
#
# Input:  Supplementary Table 1 (LUAD; Strelka_F1_Score, Tumour_Purity)
# Output: <plots_dir>/LUAD/Supp15/supp15j_LUAD_tumour_purity_boxplot (in the figure), LUAD_tumour_purity_scatterplot
#         and LUAD_tumour_purity_subset_boxplot (not in the figure), as PDF and PNG; Source Data
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

supp15_theme <- theme(
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
patients <- read_patients("metadata/LUAD_Patients.txt")
purity <- strelka_delta_f1_covariates("LUAD") %>% select(Donor_ID, Ancestry, Delta_F1_Score, Purity)

purity$Purity <- as.numeric(purity$Purity)
purity$Delta_F1_Score <- as.numeric(purity$Delta_F1_Score)
purity$Ancestry <- factor(purity$Ancestry, levels = ancestry_levels)
purity <- purity %>% filter(!is.na(Purity))

cat("Patients with tumour purity:", nrow(purity), "of", length(patients), "\n")
cat("Missing:", paste(setdiff(patients, purity$Donor_ID), collapse = ", "), "\n")
cat("Patients per ancestry:\n")
print(table(purity$Ancestry))

# Delta F1 vs tumour purity
plot_scatter <- ggplot(purity, aes(Purity, Delta_F1_Score, fill = Ancestry)) +
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
    x = "Tumour Purity",
    y = expression(Delta ~ "F1 Score (%)")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.09))) +
  scale_fill_manual(values = ancestry_colours, name = NULL, labels = ancestry_labels) +
  supp15_theme

# Tumour purity by ancestry
set.seed(1)

plot_purity <- ggplot(purity, aes(Ancestry, Purity, fill = Ancestry)) +
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
    y = "Tumour Purity"
  ) +
  scale_x_discrete(labels = ancestry_labels) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(values = ancestry_colours, name = "Ancestry", labels = ancestry_labels) +
  supp15_theme

# Samples at or above the median tumour purity, EAA vs AA + EA
median_purity <- median(purity$Purity)

above_median <- purity %>%
  filter(Purity >= median_purity) %>%
  mutate(Group = ifelse(Ancestry == "EAA_Patients", "EAA", "AA_EA"))

wilcox_result <- wilcox.test(Delta_F1_Score ~ Group, data = above_median, exact = FALSE)

cat("\n========================================\n")
cat("Median tumour purity:", median_purity, "\n")
cat("Samples at or above the median:", nrow(above_median), "\n")
cat("========================================\n")
print(table(above_median$Ancestry))
cat("\nWilcoxon EAA vs AA + EA, samples >= median tumour purity\n")
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
    y = expression(Delta ~ "F1 (%) in Samples" >= "Median Tumour Purity")
  ) +
  scale_x_discrete(
    labels = ancestry_labels,
    expand = expansion(mult = c(0.4, 0.3))
  ) +
  scale_fill_manual(values = ancestry_colours, name = "Ancestry", labels = ancestry_labels) +
  supp15_theme

plot_scatter
plot_purity
plot_subset

# Save plots and Source Data
save_plot(plot_scatter, "LUAD_tumour_purity_scatterplot", file.path("LUAD", "Supp15"), width = 8, height = 6)
save_plot(plot_purity, "supp15j_LUAD_tumour_purity_boxplot", file.path("LUAD", "Supp15"), width = 5, height = 6)
save_plot(plot_subset, "LUAD_tumour_purity_subset_boxplot", file.path("LUAD", "Supp15"), width = 6, height = 7)
write_source_data(purity, "Supplementary Figure 15j")
