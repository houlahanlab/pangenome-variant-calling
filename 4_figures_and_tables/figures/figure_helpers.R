# Shared setup for the figure scripts (sourced at the top of each one; scripts are run from the repository root):
# paths from config.yaml, Supplementary Table 1, saving plots and writing Source Data.
#
# Source Data: write_source_data() saves the values behind one panel as source_data/<panel>.csv in this repository;
# 4_figures/build_source_data.R can combine them into one Excel file (one sheet per panel).

config <- yaml::read_yaml("config.yaml")
source_data_dir <- "source_data"

# Alignment names used in Supplementary Table 1 -> names used in the stats tables and figure scripts
st1_alignments <- c(BWA_GRCh38 = "GRCh38_bwa", VG_GRCh38 = "GRCh38_linear", VG_Pangenome = "T2T_pangenome")

# Caller column prefix in Supplementary Table 1 -> caller name used in the figure scripts
st1_callers <- c(SomaticSniper = "SomaticSniper", Strelka = "Strelka", Mutect = "Mutect2", DeepSomatic = "DeepSomatic",
                 Consensus = "Consensus")

# Supplementary Table 1 for one cohort, one row per donor and alignment (Reference uses the stats table names)
read_supp_table_1 <- function(cohort = "BLCA") {
  st1 <- readxl::read_xlsx(
    file.path(config$data_dir, "Supplementary_Tables", "Supplementary_Table_1", "Supplementary_Table_1.xlsx"),
    sheet = cohort,
    na = "NA"
  )
  st1$Reference <- unname(st1_alignments[st1$Reference])
  as.data.frame(st1)
}

# Recall, precision and F1 from Supplementary Table 1 in long form: one row per donor, alignment and caller,
# callers stacked in the order given
supp_table_1_metrics <- function(cohort = "BLCA", callers = c("SomaticSniper", "Strelka", "Mutect2", "DeepSomatic")) {
  st1 <- read_supp_table_1(cohort)
  dplyr::bind_rows(lapply(callers, function(caller) {
    prefix <- names(st1_callers)[st1_callers == caller]
    data.frame(
      Donor_ID = st1$Donor_ID,
      Ancestry = st1$Ancestry,
      Variant_Caller = caller,
      Reference = st1$Reference,
      Recall = st1[[paste0(prefix, "_Recall")]],
      Precision = st1[[paste0(prefix, "_Precision")]],
      F1_Score = st1[[paste0(prefix, "_F1_Score")]],
      stringsAsFactors = FALSE
    )
  }))
}

# One row per donor: Strelka2 VG pangenome - BWA GRCh38 F1 (2 decimals, as the F1 values), somatic mutation burden
# and tumour purity (NA when not available)
strelka_delta_f1_covariates <- function(cohort = "BLCA") {
  st1 <- read_supp_table_1(cohort)
  f1 <- function(ref) st1$Strelka_F1_Score[st1$Reference == ref]
  donors <- st1[st1$Reference == "T2T_pangenome", ]
  data.frame(
    Donor_ID = donors$Donor_ID,
    Ancestry = donors$Ancestry,
    Delta_F1_Score = round(f1("T2T_pangenome") - f1("GRCh38_bwa"), 2),
    Somatic = donors$Somatic_Mutation_Burden,
    Purity = donors$Tumour_Purity,
    stringsAsFactors = FALSE
  )
}

# Patient list: donor, tumour barcode, normal barcode, ancestry
read_patients <- function(file = "metadata/BLCA_Patients.txt") {
  read.table(file, header = FALSE, stringsAsFactors = FALSE)$V1
}

# Saves one plot as PDF and PNG in <plots_dir>/<folder>
save_plot <- function(plot, name, folder = "", width = 11, height = 5) {
  out_dir <- file.path(config$plots_dir, folder)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (ext in c("pdf", "png")) {
    ggplot2::ggsave(file.path(out_dir, paste0(name, ".", ext)), plot, width = width, height = height, dpi = 300)
  }
}

# Saves the values behind one figure panel, e.g. write_source_data(df, "Figure 1b")
write_source_data <- function(data, panel) {
  dir.create(source_data_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(data, file.path(source_data_dir, paste0(gsub(" ", "_", panel), ".csv")), row.names = FALSE)
}
