# Supplementary Figure 5h-i
# Strelka2 VG pangenome - BWA GRCh38 F1 per donor against the African, European and East Asian admixture proportions
# (Carrot et al.), coloured by self-reported ancestry.
#
# Input:  Supplementary Table 1 (BLCA; Strelka_F1_Score and the admixture proportions)
# Output: <plots_dir>/Supp5/supp5h-i_*.{pdf,png}; Source Data
# Prints the donors per ancestry and the range of the F1 difference.

# Setup
library(ggplot2)
library(dplyr)
library(readxl)

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
st1 <- read_supp_table_1("BLCA")

strelka <- st1 %>% transmute(Donor_ID, Ancestry, Reference, F1_Score = Strelka_F1_Score)

# Admixture is the same on every alignment row
admixture <- st1 %>%
  filter(Reference == "T2T_pangenome") %>%
  transmute(Donor_ID, AFR = AA_Admixture_Proportion, EUR = EA_Admixture_Proportion, EAS = EAA_Admixture_Proportion)

# VG pangenome minus BWA GRCh38 F1, matched by donor
pan <- strelka %>%
  filter(Reference == "T2T_pangenome") %>%
  select(Donor_ID, Ancestry, F1_Score)

bwa <- strelka %>%
  filter(Reference == "GRCh38_bwa") %>%
  select(Donor_ID, F1_Score)

delta_strelka <- inner_join(pan, bwa, by = "Donor_ID", suffix = c("_pangenome", "_bwa")) %>%
  filter(Donor_ID %in% patients) %>%
  mutate(Delta_F1_Score = F1_Score_pangenome - F1_Score_bwa) %>%
  inner_join(admixture, by = "Donor_ID")

if (nrow(delta_strelka) != length(patients)) {
  stop("Delta F1 and admixture found for ", nrow(delta_strelka), " of the ", length(patients), " patients")
}

delta_strelka$Ancestry <- factor(delta_strelka$Ancestry, levels = ancestry_levels)

cat("Patients per ancestry:\n")
print(table(delta_strelka$Ancestry))
cat("Delta F1 range:", round(range(delta_strelka$Delta_F1_Score), 2), "\n")

# Useful functions
admixture_scatterplot <- function(x_var, x_label, show_legend) {

  ggplot(delta_strelka, aes(.data[[x_var]], Delta_F1_Score, fill = Ancestry)) +
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
      label = "Strelka2~VG[pangenome] - BWA[GRCh38]",
      parse = TRUE,
      size = 6,
      colour = "black",
      hjust = -0.1,
      vjust = 2,
      family = "Helvetica"
    ) +
    labs(
      x = x_label,
      y = expression(Delta ~ "F1 Score (%)")
    ) +
    scale_x_continuous(limits = c(-0.1, 1)) +
    scale_y_continuous(
      limits = c(-45, 45),
      breaks = seq(-40, 40, 20)
    ) +
    scale_fill_manual(
      values = ancestry_colours,
      name = NULL,
      labels = ancestry_labels
    ) +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = if (show_legend) "inside" else "none",
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
plots <- list(
  "supp5h-i_delta_f1_by_african_ancestry_scatterplot" = admixture_scatterplot("AFR", "Admixture % African Ancestry", TRUE),
  "supp5h-i_delta_f1_by_european_ancestry_scatterplot" = admixture_scatterplot("EUR", "Admixture % European Ancestry", FALSE),
  "supp5h-i_delta_f1_by_east_asian_ancestry_scatterplot" = admixture_scatterplot("EAS", "Admixture % East Asian Ancestry", FALSE)
)

for (p in plots) print(p)

# Save plots and Source Data
for (plot_name in names(plots)) {
  save_plot(plots[[plot_name]], plot_name, "Supp5", width = 6, height = 6)
}
write_source_data(delta_strelka, "Supplementary Figure 5h-i")
