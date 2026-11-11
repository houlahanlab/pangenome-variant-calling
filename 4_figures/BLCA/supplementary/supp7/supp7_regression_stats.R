# Supplementary Figure 7 (regression models, printed only)
# Linear models of the Strelka2 VG pangenome - BWA GRCh38 F1, with African and East Asian donors grouped as
# Non-European and European as the reference group:
#   0. Delta F1 ~ ancestry group (45 donors, unadjusted)
#   1. Delta F1 ~ ancestry group + somatic mutation burden (45 donors)
#   2. Delta F1 ~ ancestry group + tumour purity (44 donors; TCGA-4Z-AA7N has no purity estimate)
# Burden and purity are modelled separately; together (beta = 5.03, P = 0.166) the model is underpowered.
#
# Input:  Supplementary Table 1 (BLCA; Strelka_F1_Score, Somatic_Mutation_Burden, Tumour_Purity)
# Output: printed only

# Setup
library(dplyr)

source("4_figures/figure_helpers.R")

# European first so the Non-European coefficient is compared with European
ancestry_group_levels <- c("European", "Non_European")

# Read in files
covariates <- strelka_delta_f1_covariates("BLCA")

with_group <- function(df, covariate) {
  df$Ancestry_Group <- factor(ifelse(df$Ancestry == "EA_Patients", "European", "Non_European"),
                              levels = ancestry_group_levels)
  df[complete.cases(df[, c("Delta_F1_Score", "Ancestry_Group", covariate)]), ]
}

mutation <- with_group(select(covariates, Donor_ID, Ancestry, Delta_F1_Score, Somatic), "Somatic")
purity <- with_group(select(covariates, Donor_ID, Ancestry, Delta_F1_Score, Purity), "Purity")

# Useful functions
fit_model <- function(formula, df, title) {
  cat("\n========================================\n")
  cat(title, "\n")
  cat("========================================\n")
  cat("Patients:", nrow(df), "\n")
  print(table(df$Ancestry, df$Ancestry_Group))

  fit <- lm(formula, data = df)
  print(summary(fit))
  fit
}

# Models
fit_unadjusted <- fit_model(Delta_F1_Score ~ Ancestry_Group, mutation, "Model 0: Delta F1 ~ Ancestry group")
fit_somatic <- fit_model(Delta_F1_Score ~ Ancestry_Group + Somatic, mutation, "Model 1: Delta F1 ~ Ancestry group + Somatic")
fit_purity <- fit_model(Delta_F1_Score ~ Ancestry_Group + Purity, purity, "Model 2: Delta F1 ~ Ancestry group + Purity")

# Ancestry coefficients side by side
coefficient_table <- function(fit, model) {
  coefs <- summary(fit)$coefficients
  ancestry_rows <- grep("^Ancestry_Group", rownames(coefs))

  data.frame(
    Model = model,
    Term = rownames(coefs)[ancestry_rows],
    Beta = round(coefs[ancestry_rows, "Estimate"], 2),
    P = round(coefs[ancestry_rows, "Pr(>|t|)"], 3)
  )
}

cat("\n========================================\n")
cat("Non-European coefficients (vs European)\n")
cat("========================================\n")
print(
  rbind(
    coefficient_table(fit_unadjusted, "Ancestry group"),
    coefficient_table(fit_somatic, "Ancestry group + Somatic"),
    coefficient_table(fit_purity, "Ancestry group + Purity")
  ),
  row.names = FALSE
)
