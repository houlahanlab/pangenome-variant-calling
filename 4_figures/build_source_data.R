# Source Data: combines the per-panel CSVs written by the figure scripts (write_source_data() in figure_helpers.R)
# into one Excel file with one sheet per panel, ordered by figure, plus a Contents sheet with the full panel names
# (Excel sheet names are limited to 31 characters).
#
# Input:  source_data/*.csv (written by the figure scripts)
# Output: <data_dir>/Source_Data/Source_Data.xlsx
#
# Usage:  Rscript 4_figures/build_source_data.R

source("4_figures/figure_helpers.R")

files <- list.files(source_data_dir, pattern = "\\.csv$", full.names = TRUE)
panels <- gsub("_", " ", sub("\\.csv$", "", basename(files)))

# Main figures first, then supplementary, each in figure and panel order
is_supp <- grepl("^Supplementary", panels)
figure_number <- as.numeric(sub("^(Supplementary )?Figure ([0-9]+).*", "\\2", panels))
order_panels <- order(is_supp, figure_number, panels)
files <- files[order_panels]
panels <- panels[order_panels]

sheet_names <- sub("^Supplementary Figure", "SFig", sub("^Figure", "Fig", panels))
too_long <- nchar(sheet_names) > 31
sheet_names[too_long] <- paste0(substr(sheet_names[too_long], 1, 27), " ", seq_len(sum(too_long)))

sheets <- c(list(Contents = data.frame(Sheet = sheet_names, Panel = panels)),
            setNames(lapply(files, read.csv, check.names = FALSE), sheet_names))

out_dir <- file.path(config$data_dir, "Source_Data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, "Source_Data.xlsx")
openxlsx::write.xlsx(sheets, out_file, overwrite = TRUE)
cat(length(files), "panels ->", out_file, "\n")
