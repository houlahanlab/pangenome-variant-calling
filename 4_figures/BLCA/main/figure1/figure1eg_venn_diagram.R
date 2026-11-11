# Figure 1e (BWA GRCh38), Figure 1g (VG pangenome), Supplementary Figure 3 (VG GRCh38 Venn diagram)
# Overlap of MC3 with the Strelka2, Mutect2, SomaticSniper and DeepSomatic calls in the target panel, 45 BLCA donors.
# A call counts when it passes tumour depth >= 5 and normal depth > 8.
#
# Input:  <variant_tables_dir>/2026-10-02_<alignment>_BLCA_combined_variants.tsv (from pangenome-pipeline)
# Output: <plots_dir>/Figure1/figure1{e,g}_*_venn_diagram.{pdf,png}, <plots_dir>/Supp3/supp3_*_venn_diagram.{pdf,png},
#         *_percent_MC3_found.tsv (share of MC3 SNVs found by 0-4 callers, also printed); Source Data: count per region

# Setup
library(ggplot2)
library(ggVennDiagram)

source("4_figures/figure_helpers.R")

tool_set <- c("Strelka", "Mutect2", "SomaticSniper", "DeepSomatic")

# Read in files
variant_dir <- config$variant_tables_dir

variant_files <- c(
  "figure1e_GRCh38_bwa" = "2026-10-02_GRCh38bwa_BLCA_combined_variants.tsv",
  "figure1g_T2T_pangenome" = "2026-10-02_T2Thprc1.1_BLCA_combined_variants.tsv",
  "supp3_GRCh38_linear" = "2026-10-02_GRCh38linear_BLCA_combined_variants.tsv"
)

patients <- read_patients()

source_data_panels <- c(
  "figure1e_GRCh38_bwa" = "Figure 1e",
  "figure1g_T2T_pangenome" = "Figure 1g",
  "supp3_GRCh38_linear" = "Supplementary Figure 3 VG GRCh38 Venn"
)

# Useful functions
# Target panel variants of the 45 patients, with missing calls set to 0
read_variants <- function(file) {
  df <- read.delim(file.path(variant_dir, file))

  df <- df[df$Donor %in% patients, ]
  if (length(unique(df$Donor)) != length(patients)) {
    stop(file, " has ", length(unique(df$Donor)), " of the ", length(patients), " patients")
  }

  for (tool in c(tool_set, "MC3")) {
    df[is.na(df[, tool]), tool] <- 0
  }

  # TARGETPANEL is NA off target, so which() rather than a logical index
  df[which(df$TARGETPANEL == "Yes"), ]
}

# Fraction of MC3 variants found by 0, 1, 2, 3 or 4 callers
mc3_found <- function(df) {
  mc3 <- df[df$MC3 == 1, ]
  n_callers <- rowSums(mc3[, tool_set])

  data.frame(
    Callers = 0:4,
    Fraction = as.numeric(table(factor(n_callers, levels = 0:4))) / nrow(mc3)
  )
}

# Plots
for (plot_name in names(variant_files)) {
  df <- read_variants(variant_files[[plot_name]])
  out_dir <- file.path(config$plots_dir, ifelse(startsWith(plot_name, "supp3"), "Supp3", "Figure1"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  cat("\n========================================\n")
  cat(plot_name, ":", length(unique(df$Donor)), "patients,", nrow(df), "target panel variants\n")
  cat("========================================\n")
  found <- mc3_found(df)
  print(found)

  write.table(
    found,
    file.path(out_dir, paste0(plot_name, "_percent_MC3_found.tsv")),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )

  plot_data <- list(
    MC3 = which(df$MC3 == 1),
    Strelka = which(df$Strelka == 1),
    Mutect2 = which(df$Mutect2 == 1),
    SomSnip = which(df$SomaticSniper == 1),
    DeepSom = which(df$DeepSomatic == 1)
  )

  p <- ggVennDiagram(
    plot_data,
    label_alpha = 0,
    label = "count"
  ) +
    scale_fill_distiller(palette = "YlGnBu", direction = 1)

  print(p)

  for (ext in c("pdf", "png")) {
    ggsave(file.path(out_dir, paste0(plot_name, "_venn_diagram.", ext)),
           p, width = 6, height = 5, dpi = 300)
  }

  # Source Data: SNVs per combination of sets (1 = in the set), as counted in each Venn region
  membership <- data.frame(lapply(plot_data, function(rows) as.integer(seq_len(nrow(df)) %in% rows)))
  membership <- membership[rowSums(membership) > 0, ]
  regions <- aggregate(list(SNVs = rep(1, nrow(membership))), membership, sum)
  write_source_data(regions[do.call(order, -regions[names(plot_data)]), ], source_data_panels[[plot_name]])
}
