# Figure 3f
# Mutect2 TLOD of the variants found by both BWA GRCh38 and VG pangenome, missed by both, or missed by VG pangenome only;
# each point is one donor's median in each alignment, with paired Wilcoxon tests.
#
# Input:  <data_dir>/False_Negatives_Analysis/Mutect2_TLOD/mutect2_tlod_filter_results.tsv (tlod_measure_mutect2.py)
#         <data_dir>/False_Negatives_Analysis/VAF_results.tsv (measure_vaf.py, for the TP / FN status)
# Output: <plots_dir>/Figure3/figure3f_*.{pdf,png}; Source Data

source("4_figures/figure_helpers.R")
source("4_figures/BLCA/fn_plot_helpers.R")

out_dir <- file.path(config$plots_dir, "Figure3")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Same dot jitter on every rerun
set.seed(1)

tlod_long <- read.delim(
  file.path(config$data_dir, "False_Negatives_Analysis", "Mutect2_TLOD", "mutect2_tlod_filter_results.tsv"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    Variant_ID = trimws(variant_id),
    DONOR_ID = trimws(donor),
    Counted_In = trimws(alignment),
    Value = suppressWarnings(as.numeric(tlod))
  )

vaf <- read_vaf_results()

medians <- patient_medians(variant_groups(vaf, "Mutect2"), tlod_long, "Mutect2", "TLOD")

p <- make_plot(
  medians,
  caller_label = "Mutect2",
  y_label = "Median Somatic Log-odds Score",
  y_max = 270,
  legend_position = c(0.82, 0.55),
  y_breaks = seq(0, 200, by = 50),
  pvalues = paired_pvalues(medians),
  output_file = "figure3f_Mutect2_BWA_vs_VG_pangenome_found_and_missed_groups_three_alignments_TLOD"
)

write_source_data(p$data, "Figure 3f")
