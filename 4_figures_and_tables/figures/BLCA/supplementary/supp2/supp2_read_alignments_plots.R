# Supplementary Figure 2a-b and 2e-h
# VG pangenome - BWA GRCh38 difference in the proportion of unaligned reads and of properly paired reads per sample,
# by ancestry, for blood and tumour samples (45 BLCA donors). Proportion = reads / total reads.
#   2a-b: VG pangenome counts from vg (alignment to the graph, before projection to T2T-CHM13), "GAM vs. BAM"
#   2e-h: VG pangenome counts from GATK CollectAlignmentSummaryMetrics on the projected BAM, "BAM vs. BAM"
# BWA GRCh38 counts are from GATK in both. All plots share one y axis per metric, so the 0 line sits at the same height.
#
# Input:  <alignments_dir>/T2T_pangenome/2026-10-08_T2Thprc1.1_{vg,gatk}_alignment_statistics.tsv,
#         <alignments_dir>/GRCh38_bwa/2026-10-09_GRCh38bwa_gatk_alignment_statistics.tsv (from pangenome-pipeline)
# Output: <plots_dir>/Supp2/supp2_*_delta_{vg_gatk,gatk}_boxplot.{pdf,png}; Source Data
# Prints the samples per tissue and ancestry and the range of each difference.

# Setup
library(ggplot2)
library(dplyr)

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
patients <- read.table("metadata/BLCA_Patients.txt", header = FALSE, stringsAsFactors = FALSE)
patients <- data.frame(Donor_ID = patients$V1, Ancestry = patients$V4)

bwa <- read.delim(
  file.path(config$alignments_dir, "GRCh38_bwa", "2026-10-09_GRCh38bwa_gatk_alignment_statistics.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# VG pangenome counts: from vg (before projection) for 2a-b, from GATK (projected BAM) for 2e-h
comparisons <- list(
  "vg_gatk" = list(file = "2026-10-08_T2Thprc1.1_vg_alignment_statistics.tsv", title_tag = "GAM vs. BAM",
                   panels = "Supplementary Figure 2a-b"),
  "gatk" = list(file = "2026-10-08_T2Thprc1.1_gatk_alignment_statistics.tsv", title_tag = "BAM vs. BAM",
                panels = "Supplementary Figure 2e-h")
)

# Useful functions
# Proportions per sample, TCGA sample type 01 = tumour, 10/11 = blood
prepare_alignment_data <- function(data) {
  data %>%
    mutate(
      Donor_ID = sub("^(TCGA-[^-]+-[^-]+)-.*", "\\1", SAMPLE),
      Sample_Type = substr(sub("^TCGA-[^-]+-[^-]+-", "", SAMPLE), 1, 2),
      Tissue = case_when(
        Sample_Type == "01" ~ "Tumour",
        Sample_Type %in% c("10", "11") ~ "Blood"
      ),
      Proportion_Unaligned = UNALIGNED_READS / TOTAL_READS,
      Proportion_Paired = PROPER_PAIRS / TOTAL_READS
    ) %>%
    filter(Donor_ID %in% patients$Donor_ID, !is.na(Tissue)) %>%
    select(SAMPLE, Donor_ID, Tissue, Proportion_Unaligned, Proportion_Paired)
}

delta_boxplot <- function(data, y_var, title, y_title, limits) {
  set.seed(1)

  ggplot(data, aes(Ancestry, .data[[y_var]], fill = Ancestry)) +
    stat_boxplot(
      geom = "errorbar",
      width = 0.3,
      colour = "black",
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 1) +
    geom_boxplot(
      width = 0.6,
      outlier.shape = NA,
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
      colour = "black",
      alpha = 0.9,
      show.legend = FALSE
    ) +
    # Title centred over the middle box, leaving room on the left for arrows added in Adobe
    annotate(
      "text",
      x = 2, y = Inf,
      label = title,
      size = 5.5,
      colour = "black",
      hjust = 0.5,
      vjust = 1.5,
      family = "Helvetica"
    ) +
    labs(
      x = NULL,
      y = y_title,
      fill = "Ancestry"
    ) +
    scale_x_discrete(labels = ancestry_labels) +
    # Ticks every 0.02 on every plot
    scale_y_continuous(limits = limits, breaks = seq(limits[1], limits[2], by = 0.02)) +
    scale_fill_manual(values = ancestry_colours, labels = ancestry_labels) +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "inside",
      legend.position.inside = c(0.75, 0.16),
      legend.background = element_rect(fill = NA, colour = NA),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.title = element_text(size = 12.5),
      legend.text = element_text(size = 12.5),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 11, colour = "black", hjust = c(0.5, 0.5, 0.4)),
      axis.text.y = element_text(size = 11.5, colour = "black"),
      axis.title.y = element_text(size = 12.5, margin = margin(r = 8)),
      panel.border = element_blank(),
      axis.line.x.bottom = element_line(colour = "black"),
      axis.line.y.left = element_line(colour = "black"),
      plot.margin = margin(t = 5.5, r = 8, b = 5.5, l = 5.5)
    )
}

save_supp2_plot <- function(p, output_name, width, height) {
  out_dir <- file.path(config$plots_dir, "Supp2")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  print(p)
  for (ext in c("pdf", "png")) {
    ggsave(file.path(out_dir, paste0(output_name, ".", ext)),
           p, width = width, height = height, dpi = 300, bg = "white")
  }
}

metrics <- list(
  "unaligned_reads" = list(
    y_var = "Delta_Unaligned",
    y_title = expression(Delta ~ "Proportion of Unaligned Reads"),
    limits = c(-0.02, 0.06)
  ),
  "properly_paired_reads" = list(
    y_var = "Delta_Paired",
    y_title = expression(Delta ~ "Proportion of Properly Paired Reads"),
    limits = c(-0.06, 0.06)
  )
)

# Plots and Source Data
plots <- list()
for (file_tag in names(comparisons)) {
  title_tag <- comparisons[[file_tag]]$title_tag

  pangenome <- read.delim(
    file.path(config$alignments_dir, "T2T_pangenome", comparisons[[file_tag]]$file),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  # VG statistics use "Sample", GATK statistics use "SAMPLE"
  colnames(pangenome)[colnames(pangenome) == "Sample"] <- "SAMPLE"

  alignment_data <- inner_join(
    prepare_alignment_data(pangenome),
    prepare_alignment_data(bwa),
    by = c("SAMPLE", "Donor_ID", "Tissue"),
    suffix = c("_Pangenome", "_BWA")
  ) %>%
    left_join(patients, by = "Donor_ID") %>%
    mutate(
      Delta_Unaligned = Proportion_Unaligned_Pangenome - Proportion_Unaligned_BWA,
      Delta_Paired = Proportion_Paired_Pangenome - Proportion_Paired_BWA,
      Ancestry = factor(Ancestry, levels = ancestry_levels)
    ) %>%
    arrange(Ancestry, Donor_ID)

  cat("Samples per tissue and ancestry:\n")
  print(table(alignment_data$Tissue, alignment_data$Ancestry))

  for (tissue in c("Blood", "Tumour")) {
    n_tissue <- sum(alignment_data$Tissue == tissue)
    if (n_tissue != nrow(patients) || anyDuplicated(alignment_data$Donor_ID[alignment_data$Tissue == tissue]) > 0) {
      stop(tissue, ": ", n_tissue, " samples for the ", nrow(patients), " patients")
    }
  }

  cat("\nDelta ranges:\n")
  print(
    alignment_data %>%
      group_by(Tissue) %>%
      summarise(
        across(c(Delta_Unaligned, Delta_Paired), list(min = min, median = median, max = max)),
        .groups = "drop"
      )
  )


  for (metric in names(metrics)) {
    y_var <- metrics[[metric]]$y_var
    y_title <- metrics[[metric]]$y_title

    # Same axis for blood and tumour and for supp2ab and supp2e-h, so the 0 line sits at the same height
    limits <- metrics[[metric]]$limits
    if (min(alignment_data[[y_var]]) < limits[1] || max(alignment_data[[y_var]]) > limits[2]) {
      stop(metric, " deltas (", paste(signif(range(alignment_data[[y_var]]), 3), collapse = " to "),
           ") fall outside the y axis ", limits[1], " to ", limits[2])
    }

    for (tissue in c("Blood", "Tumour")) {
      tissue_df <- alignment_data %>% filter(Tissue == tissue)
      output_name <- paste0("supp2_", metric, "_", tolower(tissue), "_delta_", file_tag, "_boxplot")

      box <- delta_boxplot(tissue_df, y_var, paste0(tissue, " (", title_tag, ")"), y_title, limits)
      save_supp2_plot(box, output_name, 3.5, 4.7)
      plots[[output_name]] <- box
    }
  }

  write_source_data(
    select(alignment_data, SAMPLE, Donor_ID, Ancestry, Tissue, Proportion_Unaligned_Pangenome, Proportion_Unaligned_BWA,
           Delta_Unaligned, Proportion_Paired_Pangenome, Proportion_Paired_BWA, Delta_Paired),
    comparisons[[file_tag]]$panels
  )
}
