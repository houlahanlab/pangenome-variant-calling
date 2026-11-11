# Supplementary Figure 10a
# Number of heterozygous germline SNVs within 1-74 bp of a Strelka2 VG pangenome TP / FP call per donor, by ancestry,
# counted in the VG pangenome normal at sites with normal depth > 8 and the same reference base in GRCh38 and
# T2T-CHM13 (as in Figure 3d). Donors with no nearby germline SNVs count as 0.
#
# Input:  <data_dir>/Reference_Bias/paired_sites/germline_paired_counts.tsv (combine_paired_counts.py)
# Output: <plots_dir>/Supp10/supp10a_germline_variant_count.{pdf,png}; Source Data

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

# Read in files
patients <- read.table(
  "metadata/BLCA_Patients.txt",
  header = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(Donor_ID = V1, Ancestry = V4)

germline <- read.delim(
  file.path(config$data_dir, "Reference_Bias", "paired_sites", "germline_paired_counts.tsv"),
  stringsAsFactors = FALSE
) %>%
  filter(vg_pangenome_normal_depth > 8, grch38_ref == t2t_ref)  # same reference base in GRCh38 and T2T, as in Figure 3d

# Germline variants per patient
germline_count <- patients %>%
  left_join(count(germline, donor, name = "Germline_Count"), by = c("Donor_ID" = "donor")) %>%
  mutate(
    Germline_Count = ifelse(is.na(Germline_Count), 0, Germline_Count),
    Ancestry = factor(Ancestry, levels = ancestry_levels)
  ) %>%
  arrange(Ancestry, Donor_ID)

germline_count$Donor_ID <- factor(germline_count$Donor_ID, levels = germline_count$Donor_ID)

cat("Patients per ancestry:\n")
print(table(germline_count$Ancestry))
cat("Germline count range:", range(germline_count$Germline_Count), "\n")
cat("Patients with 0 germline variants nearby:",
    paste(germline_count$Donor_ID[germline_count$Germline_Count == 0], collapse = ", "), "\n")

# Plot
plot_count <- ggplot(germline_count, aes(Donor_ID, Germline_Count, fill = Ancestry)) +
  geom_col(colour = "black") +
  # Caller name inside the top left, like the other Supp 10 plots (the first bars are short)
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = "Strelka2",
    size = 6,
    colour = "black",
    hjust = -0.1,
    vjust = 2,
    family = "Helvetica"
  ) +
  labs(
    x = "Patient",
    y = "Number of Germline Variants"
  ) +
  scale_fill_manual(
    values = ancestry_colours,
    name = "Ancestry",
    labels = ancestry_labels
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    text = element_text(family = "Helvetica"),
    legend.position = "inside",
    legend.position.inside = c(0.97, 0.97),  # top right, above the short East Asian bars
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = 11, colour = "black", angle = 60, hjust = 1),
    axis.text.y = element_text(size = 15, colour = "black"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.ticks.x = element_blank(),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black")
  )

plot_count

# Save plots and Source Data
save_plot(plot_count, "supp10a_germline_variant_count", "Supp10", width = 10, height = 6)
write_source_data(germline_count, "Supplementary Figure 10a")
