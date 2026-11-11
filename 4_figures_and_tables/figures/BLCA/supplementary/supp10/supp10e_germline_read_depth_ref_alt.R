# Supplementary Figure 10e
# Per donor mean read depth (REF + ALT) against the median germline REF - ALT read count in the VG pangenome normal,
# at heterozygous germline SNVs within 1-74 bp of the Strelka2 VG pangenome TP / FP calls (normal depth > 8, same
# reference base in GRCh38 and T2T-CHM13, as in Figure 3d). Depth is on the x axis because it drives the raw
# REF - ALT count. Donors with no usable germline sites are left out.
#
# Input:  <data_dir>/Reference_Bias/paired_sites/germline_paired_counts.tsv (combine_paired_counts.py)
# Output: <plots_dir>/Supp10/supp10e_germline_read_depth_ref_alt.{pdf,png}; Source Data

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
patients <- read_patients()

germline <- read.delim(
  file.path(config$data_dir, "Reference_Bias", "paired_sites", "germline_paired_counts.tsv"),
  stringsAsFactors = FALSE
) %>%
  filter(donor %in% patients, vg_pangenome_normal_depth > 8, grch38_ref == t2t_ref)  # same reference base in GRCh38 and T2T, as in Figure 3d

# Median REF - ALT and mean depth per patient
read_depth <- germline %>%
  group_by(donor, ancestry) %>%
  summarise(
    Median_Ref_Alt = median(vg_pangenome_ref_count - vg_pangenome_alt_count),
    Mean_Read_Depth = mean(vg_pangenome_ref_count + vg_pangenome_alt_count),
    .groups = "drop"
  ) %>%
  mutate(Ancestry = factor(ancestry, levels = ancestry_levels))

cat("Patients plotted:", nrow(read_depth), "of", length(patients), "\n")
cat("No usable germline sites (left out):", paste(setdiff(patients, read_depth$donor), collapse = ", "), "\n")
print(table(read_depth$Ancestry))

corr <- cor.test(read_depth$Mean_Read_Depth, read_depth$Median_Ref_Alt, method = "spearman", exact = FALSE)
cat("Mean read depth vs median REF - ALT: rho =", round(corr$estimate, 3), ", p =", signif(corr$p.value, 3), "\n")

# Plot
plot_depth <- ggplot(read_depth, aes(Mean_Read_Depth, Median_Ref_Alt, fill = Ancestry)) +
  geom_point(
    size = 3,
    shape = 21,
    colour = "black"
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = "Strelka2~VG[pangenome]",
    parse = TRUE,
    size = 6,
    colour = "black",
    hjust = -0.1,
    vjust = 2,
    family = "Helvetica"
  ) +
  # Headroom above the points so the title does not overlap them
  scale_y_continuous(breaks = seq(0, 20, by = 10), expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    x = "Mean Number of Reads",
    y = "Median Germline (REF - ALT) Read Count"
  ) +
  scale_fill_manual(
    values = ancestry_colours,
    name = NULL,
    labels = ancestry_labels
  ) +
  theme(
    text = element_text(family = "Helvetica"),
    legend.position = "inside",
    legend.position.inside = c(0.88, 0.80),
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

plot_depth

# Save plots and Source Data
save_plot(plot_depth, "supp10e_germline_read_depth_ref_alt", "Supp10", width = 6, height = 6)
write_source_data(read_depth, "Supplementary Figure 10e")
