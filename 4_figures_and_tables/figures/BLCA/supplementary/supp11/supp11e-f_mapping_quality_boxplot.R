# Supplementary Figure 11e-f: mapping quality, Strelka2 (11e) and Mutect2 (11f)
# Median mapping quality of the ALT-supporting reads (mean_alt_mq) per donor of the Strelka2 and Mutect2 variants found
# by both BWA GRCh38 and VG pangenome, missed by both, or missed by VG pangenome only, with paired Wilcoxon tests.
# The VG pangenome values come from the VG GRCh38 pangenome table (GRCh38hprc1.1).
#
# Input:  <quality_scores_dir>/2026-09-30_*__variant_read_qualities.tsv (from pangenome-pipeline)
#         <data_dir>/False_Negatives_Analysis/VAF_results.tsv (measure_vaf.py, for the TP / FN status)
# Output: <plots_dir>/Supp11/supp11e_Strelka_*, supp11f_Mutect2_* (*.pdf, *.png); Source Data

source("4_figures/figure_helpers.R")
source("4_figures/BLCA/fn_plot_helpers.R")

out_dir <- file.path(config$plots_dir, "Supp11")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Same dot jitter on every rerun
set.seed(1)

quality_scores_dir <- config$quality_scores_dir

quality_file_map <- c(
  BWA_GRCh38_linear = "2026-09-30_GRCh38bwa__variant_read_qualities.tsv",
  VG_GRCh38_linear = "2026-09-30_GRCh38linear__variant_read_qualities.tsv",
  VG_T2T_pangenome = "2026-09-30_GRCh38hprc1.1__variant_read_qualities.tsv"
)

mq_long <- bind_rows(
  lapply(names(quality_file_map), function(label) {
    read.delim(
      file.path(quality_scores_dir, quality_file_map[[label]]),
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE
    ) %>%
      transmute(
        Variant_ID = paste(sub("^chr", "", chrom), pos, ref, alt, sep = ":"),
        DONOR_ID = trimws(as.character(Donor)),
        Counted_In = label,
        Value = suppressWarnings(as.numeric(mean_alt_mq))
      )
  })
)

cat("Patients in mapping-quality tables:", n_distinct(mq_long$DONOR_ID), "\n")

vaf <- read_vaf_results()

panel_letters <- c("Strelka" = "e", "Mutect2" = "f")

for (caller in c("Strelka", "Mutect2")) {

  medians <- patient_medians(variant_groups(vaf, caller), mq_long, caller, "mapping quality")

  p <- make_plot(
    medians,
    caller_label = ifelse(caller == "Strelka", "Strelka2", caller),
    y_label = "Median Mapping Quality",
    y_max = 70,
    y_breaks = seq(0, 60, by = 10),
    legend_position = c(0.81, 0.30),
    output_file = paste0("supp11", panel_letters[[caller]], "_", caller, "_BWA_vs_VG_pangenome_found_and_missed_groups_three_BAMs_MQ")
  )
  write_source_data(p$data, paste0("Supplementary Figure 11", panel_letters[[caller]]))
}
