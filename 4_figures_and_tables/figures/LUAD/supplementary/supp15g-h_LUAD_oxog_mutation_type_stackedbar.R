# Supplementary Figure 15g-h (LUAD)
# Mutation types of the Strelka2 TP, FP and FN SNVs of TCGA-05-4396, a LUAD tumour flagged for OxoG oxidative damage,
# for BWA GRCh38 (15g) and VG pangenome (15h). Strand-complement types are side by side so the C>A / G>T OxoG pattern
# is easy to see.
#
# Input:  <data_dir>/LUAD/TCGA_Comparison/strelka_vs_tcga/<alignment>/EA_Patients/TCGA-05-4396/TCGA-05-4396-01A/
#             common_variants_full.txt (TP), strelka_unique_pass.txt (FP), tcga_unique_pass.txt (FN) (compare_to_mc3.py)
# Output: <plots_dir>/LUAD/Supp15/TCGA-05-4396_{BWA_GRCh38,VG_pangenome}_oxog_mutation_type_stackedbar.{pdf,png};
#         Source Data
# Prints the % of FPs that are C>A / G>T for each alignment.

# Setup
library(ggplot2)
library(dplyr)

source("4_figures/figure_helpers.R")

donor_id <- "TCGA-05-4396"

call_colours <- c(
  "False_Positive" = "#FCE4EC",
  "True_Positive" = "#F06292",
  "False_Negative" = "#F48FB1"
)

call_labels <- c(
  "False_Positive" = "False Positive",
  "True_Positive" = "True Positive",
  "False_Negative" = "False Negative"
)

# Strand-complement pairs side by side
mutation_type_levels <- c("C>A", "G>T", "C>G", "G>C", "C>T", "G>A", "T>A", "A>T", "T>C", "A>G", "T>G", "A>C")

# Useful functions
# TP, FP and FN SNVs of the donor for one alignment, counted per mutation type
count_mutation_types <- function(alignment) {
  comparison_dir <- file.path(config$data_dir, "LUAD", "TCGA_Comparison", "strelka_vs_tcga", alignment, "EA_Patients",
                              donor_id, paste0(donor_id, "-01A"))
  files <- c(
    True_Positive = "common_variants_full.txt",
    False_Positive = "strelka_unique_pass.txt",
    False_Negative = "tcga_unique_pass.txt"
  )
  variants <- bind_rows(lapply(names(files), function(status) {
    df <- read.table(file.path(comparison_dir, files[[status]]), sep = "\t",
                     col.names = c("Chrom", "Start", "Pos", "Ref", "Alt"), stringsAsFactors = FALSE)
    mutate(df, Call_Status = status)
  }))

  cat("\n========================================\n")
  cat("Processing:", alignment, "\n")
  cat("========================================\n")
  cat("True Positives: ", sum(variants$Call_Status == "True_Positive"), "\n")
  cat("False Positives:", sum(variants$Call_Status == "False_Positive"), "\n")
  cat("False Negatives:", sum(variants$Call_Status == "False_Negative"), "\n")

  plot_data <- variants %>%
    mutate(Mutation_Type = paste0(Ref, ">", Alt)) %>%
    filter(Mutation_Type %in% mutation_type_levels) %>%
    mutate(
      Mutation_Type = factor(Mutation_Type, levels = mutation_type_levels),
      Call_Status = factor(Call_Status, levels = names(call_colours))
    ) %>%
    count(Mutation_Type, Call_Status, name = "Count")

  cat("\nTotal SNVs plotted:", sum(plot_data$Count), "\n")
  cat("\nCounts by mutation type and call status:\n")
  print(arrange(plot_data, Mutation_Type, Call_Status))
  plot_data
}

make_oxog_plot <- function(plot_data, plot_title) {
  ggplot(plot_data, aes(Mutation_Type, Count, fill = Call_Status)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.4) +
    annotate(
      "text",
      x = 0.5, y = Inf,
      label = plot_title,
      parse = TRUE,
      size = 7,
      colour = "black",
      hjust = 0,
      vjust = 1.5,
      family = "Helvetica"
    ) +
    labs(x = NULL, y = "Number of Variants") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    scale_fill_manual(values = call_colours, labels = call_labels, name = "Call Type") +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "inside",
      legend.position.inside = c(0.85, 0.80),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_blank(),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 16, colour = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 18, colour = "black"),
      axis.title.y = element_text(size = 20),
      panel.border = element_blank(),
      axis.line.x.bottom = element_line(colour = "black"),
      axis.line.y.left = element_line(colour = "black")
    )
}

# Plots
alignments <- data.frame(
  Alignment = c("GRCh38_bwa", "T2T_pangenome"),
  Label = c("BWA GRCh38", "VG pangenome"),
  File_Label = c("BWA_GRCh38", "VG_pangenome"),
  Plot_Title = c("BWA[GRCh38]", "VG[pangenome]"),
  Panel = c("Supplementary Figure 15g", "Supplementary Figure 15h"),
  stringsAsFactors = FALSE
)

oxog_summary <- bind_rows(lapply(seq_len(nrow(alignments)), function(i) {
  a <- alignments[i, ]
  plot_data <- count_mutation_types(a$Alignment)
  plot_oxog <- make_oxog_plot(plot_data, a$Plot_Title)
  print(plot_oxog)

  # Save plots and Source Data
  save_plot(plot_oxog, paste0(donor_id, "_", a$File_Label, "_oxog_mutation_type_stackedbar"), file.path("LUAD", "Supp15"),
            width = 9, height = 6)
  write_source_data(plot_data, a$Panel)

  fp <- filter(plot_data, Call_Status == "False_Positive")
  oxog_fp <- sum(fp$Count[fp$Mutation_Type %in% c("C>A", "G>T")])
  data.frame(
    Alignment = a$Label,
    OxoG_FP = oxog_fp,
    Total_FP = sum(fp$Count),
    OxoG_FP_Percent = if (sum(fp$Count) > 0) 100 * oxog_fp / sum(fp$Count) else NA_real_
  )
}))

cat("\n========================================\n")
cat("C>A / G>T share of the false positives\n")
cat("========================================\n")
print(oxog_summary, row.names = FALSE)
