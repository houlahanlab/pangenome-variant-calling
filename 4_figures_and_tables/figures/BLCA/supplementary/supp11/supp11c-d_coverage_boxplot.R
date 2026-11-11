# Supplementary Figure 11c-d: coverage depth, Strelka2 (11c) and Mutect2 (11d)
# Median coverage depth (samtools mpileup depth at the variant) per donor of the Strelka2 and Mutect2 variants found by
# both BWA GRCh38 and VG pangenome, missed by both, or missed by VG pangenome only, measured in each alignment, with
# paired Wilcoxon tests. The VG pangenome values are read from the VG GRCh38 pangenome BAMs (see measure_vaf.py).
#
# Input:  <data_dir>/False_Negatives_Analysis/VAF_results.tsv (measure_vaf.py)
# Output: <plots_dir>/Supp11/supp11c_Strelka_*, supp11d_Mutect2_* (*.pdf, *.png); Source Data

source("4_figures/figure_helpers.R")
source("4_figures/BLCA/fn_plot_helpers.R")

out_dir <- file.path(config$plots_dir, "Supp11")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Same dot jitter on every rerun
set.seed(1)

vaf <- read_vaf_results()

panel_letters <- c("Strelka" = "c", "Mutect2" = "d")

for (caller in c("Strelka", "Mutect2")) {

  depth_long <- vaf %>%
    filter(TOOL == caller) %>%
    transmute(Variant_ID, DONOR_ID, Counted_In, Value = MPILEUP_DEPTH)

  medians <- patient_medians(variant_groups(vaf, caller), depth_long, caller, "coverage depth")

  p <- make_plot(
    medians,
    caller_label = ifelse(caller == "Strelka", "Strelka2", caller),
    y_label = "Median Coverage Depth",
    legend_position = c(0.55, 0.72),
    output_file = paste0("supp11", panel_letters[[caller]], "_", caller, "_BWA_vs_VG_pangenome_found_and_missed_groups_three_BAMs_coverage")
  )
  write_source_data(p$data, paste0("Supplementary Figure 11", panel_letters[[caller]]))
}
