# Supplementary Figure 15i (LUAD)
# Strelka2 VG pangenome - BWA GRCh38 F1 against somatic mutation burden (PanCancer Atlas), mutation burden by ancestry,
# and the F1 difference by ancestry in donors below the median burden (29 LUAD donors; as Supplementary Figure 7a-c).
#
# Input:  Supplementary Table 1 (LUAD; Strelka_F1_Score, Somatic_Mutation_Burden)
# Output: <plots_dir>/LUAD/Supp15/supp15i_LUAD_mutation_burden_boxplot (in the figure), LUAD_mutation_burden_scatterplot
#         and LUAD_mutation_burden_subset_boxplot (not in the figure), as PDF and PNG; Source Data
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
mutation <- strelka_delta_f1_covariates("LUAD") %>% select(Donor_ID, Ancestry, Delta_F1_Score, Somatic)

mutation$Somatic <- as.numeric(mutation$Somatic)
mutation$Delta_F1_Score <- as.numeric(mutation$Delta_F1_Score)
mutation$Ancestry <- factor(mutation$Ancestry, levels = ancestry_levels)

if (nrow(mutation) != length(patients) || any(is.na(mutation$Somatic))) {
  stop("Mutation burden found for ", sum(!is.na(mutation$Somatic)), " of the ", length(patients), " patients")
}

cat("Patients per ancestry:\n")
print(table(mutation$Ancestry))

# Delta F1 vs mutation burden
plot_scatter <- ggplot(mutation, aes(Somatic, Delta_F1_Score, fill = Ancestry)) +
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
    x = "Mutation Burden",
    y = expression(Delta ~ "F1 Score (%)")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  scale_fill_manual(values = ancestry_colours, name = NULL, labels = ancestry_labels) +
  supp15_theme

# Mutation burden by ancestry
set.seed(1)

plot_burden <- ggplot(mutation, aes(Ancestry, Somatic, fill = Ancestry)) +
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
    y = "Mutation Burden"
  ) +
  scale_x_discrete(labels = ancestry_labels) +
  scale_fill_manual(values = ancestry_colours, name = "Ancestry", labels = ancestry_labels) +
  supp15_theme

# Samples below the median mutation burden, EAA vs AA + EA
median_burden <- median(mutation$Somatic)

below_median <- mutation %>%
  filter(Somatic < median_burden) %>%
  mutate(Group = ifelse(Ancestry == "EAA_Patients", "EAA", "AA_EA"))

wilcox_result <- wilcox.test(Delta_F1_Score ~ Group, data = below_median, exact = FALSE)

cat("\n========================================\n")
cat("Median mutation burden:", median_burden, "\n")
cat("Samples below the median:", nrow(below_median), "\n")
cat("========================================\n")
print(table(below_median$Ancestry))
cat("\nWilcoxon EAA vs AA + EA, samples < median mutation burden\n")
print(wilcox_result)

set.seed(1)

plot_subset <- ggplot(below_median, aes(Ancestry, Delta_F1_Score, fill = Ancestry)) +
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
    y = expression(Delta ~ "F1 (%) in Samples < Median Mutation Burden")
  ) +
  scale_x_discrete(
    labels = ancestry_labels,
    expand = expansion(mult = c(0.4, 0.3))
  ) +
  scale_fill_manual(values = ancestry_colours, name = "Ancestry", labels = ancestry_labels) +
  supp15_theme

plot_scatter
plot_burden
plot_subset

# Save plots and Source Data
save_plot(plot_scatter, "LUAD_mutation_burden_scatterplot", file.path("LUAD", "Supp15"), width = 8, height = 6)
save_plot(plot_burden, "supp15i_LUAD_mutation_burden_boxplot", file.path("LUAD", "Supp15"), width = 5, height = 6)
save_plot(plot_subset, "LUAD_mutation_burden_subset_boxplot", file.path("LUAD", "Supp15"), width = 6, height = 7)
write_source_data(mutation, "Supplementary Figure 15i")
