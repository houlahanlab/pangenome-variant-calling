# Supplementary Figure 9i
# F1 score of Strelka2 and Mutect2 calls from VG pangenome alignment before and after removing calls
# found in gnomAD (exact CHROM, POS, REF and ALT match). The MC3 truth set is not filtered.
#
# Input:  <data_dir>/Stats_Summary/<caller>/mega_stats_table.txt (before filtering)
#         <data_dir>/Gnomad_Filtering/<caller>/Stats_Summary/stats_table_pangenome_hard_gnomad.txt (after filtering)
# Output: <plots_dir>/Supp9/supp9i_gnomad_pre_post_f1_boxplot.{pdf,png}; Source Data
# Prints the calls removed per patient and a paired Wilcoxon signed-rank test of F1 before vs after filtering.

# Setup
library(ggplot2)
library(dplyr)

source("4_figures/figure_helpers.R")

caller_colours <- c(
  "Strelka" = "#f1a1ff",
  "Mutect2" = "#9c8fff"
)

caller_labels <- c(
  "Strelka" = "Strelka2",
  "Mutect2" = "Mutect2"
)

callers <- c("Strelka", "Mutect2")

# Read in files
patients <- read_patients()

read_stats <- function(file) {
  read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
}

pre <- bind_rows(lapply(callers, function(caller) {
  read_stats(file.path(config$data_dir, "Stats_Summary", caller, "mega_stats_table.txt")) %>%
    filter(Reference == "T2T_pangenome", Donor_ID %in% patients) %>%
    transmute(Donor_ID, Variant_Caller = caller, F1 = as.numeric(`F1_Score(%)`))
}))

post <- bind_rows(lapply(callers, function(caller) {
  read_stats(file.path(config$data_dir, "Gnomad_Filtering", caller, "Stats_Summary", "stats_table_pangenome_hard_gnomad.txt")) %>%
    filter(Donor_ID %in% patients) %>%
    transmute(Donor_ID, Variant_Caller = caller, F1 = as.numeric(`F1_Score(%)`),
              Removed_by_gnomAD, Pre_gnomAD_Calls)
}))

f1 <- bind_rows(
  pre %>% mutate(Filter = "Pre"),
  post %>% select(Donor_ID, Variant_Caller, F1) %>% mutate(Filter = "Post")
) %>%
  mutate(
    Filter = factor(Filter, levels = c("Pre", "Post")),
    Variant_Caller = factor(Variant_Caller, levels = callers)
  )

# Statistics
cat("\nCalls removed by gnomAD per patient (median, % of pre-gnomAD calls):\n")
print(as.data.frame(post %>% group_by(Variant_Caller) %>%
  summarise(Median_Removed = median(Removed_by_gnomAD),
            Median_Percent_Removed = round(median(100 * Removed_by_gnomAD / Pre_gnomAD_Calls), 1), .groups = "drop")))

cat("\nF1 before vs after gnomAD filtering:\n")
for (caller in callers) {
  wide <- inner_join(pre %>% filter(Variant_Caller == caller) %>% select(Donor_ID, Pre = F1),
                     post %>% filter(Variant_Caller == caller) %>% select(Donor_ID, Post = F1), by = "Donor_ID")
  test <- wilcox.test(wide$Post, wide$Pre, paired = TRUE, exact = FALSE)
  cat(sprintf("%-8s n = %d  median pre = %.2f  median post = %.2f  median change = %.2f  lower after: %d  p = %.3g\n",
              caller, nrow(wide), median(wide$Pre), median(wide$Post), median(wide$Post - wide$Pre),
              sum(wide$Post < wide$Pre), test$p.value))
}

# Plot
dodge <- position_dodge(width = 0.9)
jitterdodge <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9, seed = 1)

set.seed(1)

plot_f1 <- ggplot(f1, aes(Filter, F1, fill = Variant_Caller)) +
  annotate(
    "rect",
    xmin = 1.5, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    alpha = 0.2,
    fill = "grey"
  ) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.36,
    colour = "black",
    linewidth = 0.6,
    position = dodge
  ) +
  geom_boxplot(
    width = 0.72,
    outlier.shape = NA,
    alpha = 1,
    colour = "black",
    linewidth = 0.6,
    position = dodge
  ) +
  geom_point(
    aes(colour = Variant_Caller),
    size = 2,
    shape = 21,
    fill = NA,
    alpha = 0.9,
    position = jitterdodge
  ) +
  labs(
    x = NULL,
    y = "F1 Score (%)"
  ) +
  scale_x_discrete(labels = c("Pre" = "Pre-gnomAD", "Post" = "Post-gnomAD")) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  scale_fill_manual(values = caller_colours, guide = "none") +
  scale_colour_manual(
    values = c("Strelka" = "black", "Mutect2" = "black"),
    labels = caller_labels,
    name = "Variant Caller"
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        shape = 22,
        size = 6,
        fill = unname(caller_colours),
        colour = "black"
      )
    )
  ) +
  annotate(
    "text",
    x = 0.5, y = Inf,
    label = "VG[pangenome]",
    parse = TRUE,
    size = 6,
    colour = "black",
    hjust = 0,
    vjust = 2,
    family = "Helvetica"
  ) +
  geom_vline(
    xintercept = 1.5,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.6
  ) +
  theme(
    text = element_text(family = "Helvetica"),
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.15),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 15, colour = "black"),
    axis.title.y = element_text(size = 15),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black")
  )

plot_f1

# Save plots and Source Data
save_plot(plot_f1, "supp9i_gnomad_pre_post_f1_boxplot", "Supp9", width = 8, height = 6)
write_source_data(left_join(select(f1, Donor_ID, Variant_Caller, Filter, F1),
                            select(post, Donor_ID, Variant_Caller, Pre_gnomAD_Calls, Removed_by_gnomAD),
                            by = c("Donor_ID", "Variant_Caller")),
                  "Supplementary Figure 9i")
