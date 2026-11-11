# Supplementary Figure 1e-g
# Ancestry proportions of the 45 BLCA donors from Yuan et al. against Carrot et al., with Spearman correlation:
# African (1e), European (1f) and East Asian (1g). Also saves the three panels side by side.
#
# Input:  Supplementary Table 1 (BLCA; Carrot et al. AA / EA / EAA_Admixture_Proportion, Yuan_* columns)
# Output: <plots_dir>/Supp1/supp1{e,f,g}_*.{pdf,png} and supp1e-g_*_combined.{pdf,png}; Source Data
# Prints the Spearman rho and p-value for each ancestry.

# Setup
library(ggplot2)
library(dplyr)
library(patchwork)
library(readxl)

source("4_figures/figure_helpers.R")

ancestry_colours <- c(
  "AA" = "#add8e6",
  "EA" = "#FEF5C3",
  "EAA" = "#fdc0cb"
)

# Read in files
# Ancestry proportions are the same on every alignment row, so one row per donor is kept
patients <- read_patients()

merged <- read_supp_table_1("BLCA") %>%
  filter(Reference == "T2T_pangenome") %>%
  transmute(
    Patient = Donor_ID,
    Carrot_AA = AA_Admixture_Proportion,
    Carrot_EA = EA_Admixture_Proportion,
    Carrot_EAA = EAA_Admixture_Proportion,
    Yuan_AA = Yuan_AA_Admixture_Proportion,
    Yuan_EA = Yuan_EA_Admixture_Proportion,
    Yuan_EAA = Yuan_EAA_Admixture_Proportion
  ) %>%
  filter(!is.na(Carrot_AA), !is.na(Yuan_AA))

cat("Patients with both Yuan and Carrot ancestry:", nrow(merged), "of", length(patients), "\n")
if (nrow(merged) != length(patients)) {
  stop("Missing from Yuan et al. or Carrot et al.: ", paste(setdiff(patients, merged$Patient), collapse = ", "))
}

# Useful functions
# p-value as plotmath, e.g. 1.23 x 10^-5
format_p <- function(p) {
  parts <- strsplit(formatC(p, format = "e", digits = 2), "e")[[1]]
  paste0(parts[1], " %*% 10^{", as.integer(parts[2]), "}")
}

# Yuan proportion on x, Carrot-Zhang proportion on y, for one ancestry
corr_plot <- function(ancestry) {
  x_col <- paste0("Yuan_", ancestry)
  y_col <- paste0("Carrot_", ancestry)

  corr <- cor.test(merged[[x_col]], merged[[y_col]], method = "spearman", exact = FALSE)

  cat("\n", ancestry, ": rho =", round(corr$estimate, 3), ", p =", signif(corr$p.value, 3), "\n")

  rho_label <- paste0(
    "atop(Spearman~rho == ", round(as.numeric(corr$estimate), 3),
    ", italic(p) == ", format_p(corr$p.value), ")"
  )

  ggplot(merged, aes(.data[[x_col]], .data[[y_col]])) +
    geom_point(
      size = 3,
      shape = 21,
      fill = ancestry_colours[[ancestry]],
      colour = "black"
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      colour = "black",
      linewidth = 0.8
    ) +
    annotate(
      "text",
      x = 0, y = 1,
      label = rho_label,
      parse = TRUE,
      size = 5,
      colour = "black",
      hjust = 0,
      vjust = 1.2,
      family = "Helvetica"
    ) +
    labs(
      x = paste0("Yuan et al. ", ancestry, " ancestry proportion"),
      y = paste0("Carrot-Zhang et al. ", ancestry, " ancestry proportion")
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "none",
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
  "supp1e_BLCA_AA_Yuan_vs_Carrot_correlation" = corr_plot("AA"),
  "supp1f_BLCA_EA_Yuan_vs_Carrot_correlation" = corr_plot("EA"),
  "supp1g_BLCA_EAA_Yuan_vs_Carrot_correlation" = corr_plot("EAA")
)

combined_plot <- wrap_plots(plots, ncol = 3)
combined_plot

# Save plots and Source Data
for (plot_name in names(plots)) {
  save_plot(plots[[plot_name]], plot_name, "Supp1", width = 5, height = 6)
}
save_plot(combined_plot, "supp1e-g_BLCA_Yuan_vs_Carrot_correlation_combined", "Supp1", width = 15, height = 6)

write_source_data(merged, "Supplementary Figure 1e-g")
