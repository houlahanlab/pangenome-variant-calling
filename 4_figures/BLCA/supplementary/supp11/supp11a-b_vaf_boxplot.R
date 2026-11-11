# Supplementary Figure 11a-b: VAF, Strelka2 (11a) and Mutect2 (11b)
# Median VAF per donor of the Strelka2 and Mutect2 variants found by both BWA GRCh38 and VG pangenome, missed by both, or
# missed by VG pangenome only, measured in each alignment, with paired Wilcoxon tests. The VG pangenome values are read
# from the VG GRCh38 pangenome BAMs (see measure_vaf.py).
#
# Input:  <data_dir>/False_Negatives_Analysis/VAF_results.tsv (measure_vaf.py)
# Output: <plots_dir>/Supp11/supp11a_Strelka_*, supp11b_Mutect2_* (*.pdf, *.png); Source Data

source("4_figures/figure_helpers.R")
source("4_figures/BLCA/fn_plot_helpers.R")

out_dir <- file.path(config$plots_dir, "Supp11")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Same dot jitter on every rerun
set.seed(1)

vaf <- read_vaf_results()

panel_letters <- c("Strelka" = "a", "Mutect2" = "b")

for (caller in c("Strelka", "Mutect2")) {

  vaf_long <- vaf %>%
    filter(TOOL == caller) %>%
    transmute(Variant_ID, DONOR_ID, Counted_In, Value = VAF)

  medians <- patient_medians(variant_groups(vaf, caller), vaf_long, caller, "VAF")

  p <- make_plot(
    medians,
    caller_label = ifelse(caller == "Strelka", "Strelka2", caller),
    y_label = "Median Variant Allele Frequency",
    y_max = 1.05,
    y_breaks = seq(0, 1, by = 0.2),
    output_file = paste0("supp11", panel_letters[[caller]], "_", caller, "_BWA_vs_VG_pangenome_found_and_missed_groups_three_BAMs_VAF")
  )
  write_source_data(p$data, paste0("Supplementary Figure 11", panel_letters[[caller]]))
}
