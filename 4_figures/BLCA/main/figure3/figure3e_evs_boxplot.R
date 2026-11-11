# Figure 3e
# Strelka2 somatic EVS (SomaticEVS) of the variants found by both BWA GRCh38 and VG pangenome, missed by both, or missed by VG pangenome only;
# each point is one donor's median in each alignment, with paired Wilcoxon tests.
#
# Input:  <data_dir>/False_Negatives_Analysis/Strelka_EVS/strelka_evs_filter_results.tsv (evs_measure_strelka.py)
#         <data_dir>/False_Negatives_Analysis/VAF_results.tsv (measure_vaf.py, for the TP / FN status)
# Output: <plots_dir>/Figure3/figure3e_*.{pdf,png}; Source Data

source("4_figures/figure_helpers.R")
source("4_figures/BLCA/fn_plot_helpers.R")

out_dir <- file.path(config$plots_dir, "Figure3")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Same dot jitter on every rerun
set.seed(1)

evs_long <- read.delim(
  file.path(config$data_dir, "False_Negatives_Analysis", "Strelka_EVS", "strelka_evs_filter_results.tsv"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    Variant_ID = trimws(variant_id),
    DONOR_ID = trimws(donor),
    Counted_In = trimws(alignment),
    Value = suppressWarnings(as.numeric(somatic_evs))
  )

vaf <- read_vaf_results()

medians <- patient_medians(variant_groups(vaf, "Strelka"), evs_long, "Strelka", "EVS")

p <- make_plot(
  medians,
  caller_label = "Strelka2",
  y_label = "Median Somatic Empirical Variant Score",
  y_max = 27,
  y_breaks = seq(0, 20, by = 5),
  pvalues = paired_pvalues(medians),
  legend_position = c(0.55, 0.60),
  output_file = "figure3e_Strelka_BWA_vs_VG_pangenome_found_and_missed_groups_three_alignments_EVS"
)

write_source_data(p$data, "Figure 3e")
