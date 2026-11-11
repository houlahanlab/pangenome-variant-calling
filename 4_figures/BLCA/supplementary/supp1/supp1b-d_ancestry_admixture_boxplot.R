# Supplementary Figure 1b-d
# Genetic ancestry admixture (Carrot et al.) of the 45 BLCA donors by self-reported ancestry: African (1b),
# European (1c) and East Asian (1d) admixture proportions.
#
# Input:  Supplementary Table 1 (BLCA; AA / EA / EAA_Admixture_Proportion)
# Output: <plots_dir>/Supp1/supp1{b,c,d}_*.{pdf,png}; Source Data
# Prints the donors per ancestry.

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
# Admixture is the same on every alignment row, so one row per donor is kept
admixture <- read_supp_table_1("BLCA") %>%
  filter(Reference == "T2T_pangenome") %>%
  transmute(Donor_ID, Ancestry, AFR = AA_Admixture_Proportion, EUR = EA_Admixture_Proportion,
            EAS = EAA_Admixture_Proportion)

if (nrow(admixture) != length(read_patients())) {
  stop("Admixture found for ", nrow(admixture), " of the ", length(read_patients()), " patients")
}

admixture$Ancestry <- factor(admixture$Ancestry, levels = ancestry_levels)

cat("Patients per ancestry:\n")
print(table(admixture$Ancestry))

# Useful functions
# legend_x: horizontal position of the legend, placed above the lowest box (0.97 = top right, 0.03 = top left)
admixture_boxplot <- function(y_var, y_label, legend_x = 0.97) {

  set.seed(1)

  ggplot(admixture, aes(Ancestry, .data[[y_var]], fill = Ancestry)) +
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
      alpha = 0.9,
      show.legend = FALSE
    ) +
    labs(
      x = NULL,
      y = y_label
    ) +
    scale_x_discrete(labels = ancestry_labels) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = ancestry_colours, labels = ancestry_labels, name = "Ancestry") +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "inside",
      legend.position.inside = c(legend_x, 0.97),
      legend.justification = c(legend_x, 1),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 14, colour = "black"),
      axis.text.y = element_text(size = 15, colour = "black"),
      axis.title.y = element_text(size = 15),
      panel.border = element_blank(),
      axis.line.x.bottom = element_line(colour = "black"),
      axis.line.y.left = element_line(colour = "black")
    )
}

# Plots
plots <- list(
  "supp1b_blca_ancestry_score_boxplot_AA" = admixture_boxplot("AFR", "Admixture % African Ancestry"),
  "supp1c_blca_ancestry_score_boxplot_EA" = admixture_boxplot("EUR", "Admixture % European Ancestry"),
  "supp1d_blca_ancestry_score_boxplot_EAA" = admixture_boxplot("EAS", "Admixture % East Asian Ancestry", legend_x = 0.03)
)

for (p in plots) print(p)

# Save plots and Source Data
for (plot_name in names(plots)) {
  save_plot(plots[[plot_name]], plot_name, "Supp1", width = 5, height = 6)
}
write_source_data(admixture, "Supplementary Figure 1b-d")
