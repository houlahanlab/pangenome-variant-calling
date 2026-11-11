#!/usr/bin/env python3

# Supplementary Table 3: one row per Strelka2 SNV in a cancer driver gene (Bailey et al.; BLCA and PANCAN genes) that is
# shared with or exclusive to MC3, for VG pangenome and BWA GRCh38. Uses the first SnpEff annotation of each SNV.
#
# Input:  <data_dir>/Variant_Caller_Overlap/<alignment>/Shared_Variants_Annotations/<set>/
#             <status>_filtered_shared_by_1_caller_in_capture_annotated.vcf (3_analyses/impact_variants)
#         driver gene list from config.yaml
# Output: <data_dir>/Supplementary_Tables/Supplementary_Table_3/Supplementary_Table_3.xlsx
#         Prints the SNV, tumour and gene counts quoted in the text
#
# Usage:  python3 2_benchmarking/supplementary_tables/build_supplementary_table_3.py

import argparse
import csv
import re
import statistics
import subprocess
from collections import Counter
from pathlib import Path

import yaml

ALIGNMENT_LABELS = {"T2T_pangenome": "VG_Pangenome", "GRCh38_bwa": "BWA_GRCh38"}
HIGH_IMPACT = {"HIGH", "MODERATE"}

# Alignment, status, annotation folder, set name
SNV_SETS = [
    ("T2T_pangenome", "TP", "True_Positives", "Shared"),
    ("T2T_pangenome", "FN", "False_Negatives", "MC3-exclusive"),
    ("T2T_pangenome", "FP", "False_Positives", "VG_Pangenome-exclusive"),
    ("GRCh38_bwa", "FN", "False_Negatives", "MC3-exclusive"),
    ("GRCh38_bwa", "FP", "False_Positives", "BWA_GRCh38-exclusive"),
]

COLUMNS = ["Donor_ID", "Ancestry", "Reference", "SNV_Set", "Chromosome", "Position", "Gene", "Predicted_Impact"]


# BLCA and PANCAN driver genes
def read_driver_genes(path, rscript):
    text = subprocess.run(
        [rscript, "-e", f"genes <- readxl::read_xlsx('{path}'); "
         "write.table(genes[, c('Gene', 'Cancer')], stdout(), sep = '\\t', row.names = FALSE, quote = FALSE)"],
        capture_output=True, text=True, check=True).stdout
    blca, pancan = set(), set()
    for line in text.splitlines()[1:]:
        gene, cancer = line.split("\t")[:2]
        if cancer == "BLCA":
            blca.add(gene)
        elif cancer == "PANCAN":
            pancan.add(gene)
    return blca, pancan


# Strelka2 driver gene SNVs from an annotated VCF, using the first SnpEff annotation
def read_annotated_vcf(path, driver_genes):
    variants = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt, _, _, info = line.rstrip("\n").split("\t")[:8]
            sample = re.search(r"SAMPLE=([^;]+)", info)
            tool = re.search(r"TOOL=([^;]+)", info)
            ann = re.search(r"ANN=([^\t]+)", info)
            if len(ref) != 1 or len(alt) != 1 or not (sample and tool and ann) or tool.group(1) != "Strelka":
                continue

            ann = ann.group(1).split(",")[0].split("|")
            if len(ann) < 11 or ann[3] not in driver_genes:
                continue

            variants.append({"chrom": chrom, "pos": int(pos), "donor": sample.group(1), "gene": ann[3], "impact": ann[2]})
    return variants


def main():
    parser = argparse.ArgumentParser(description="Build Supplementary Table 3.")
    parser.add_argument("--patients", default="metadata/BLCA_Patients.txt")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    data_dir = Path(config["data_dir"])
    out_dir = data_dir / "Supplementary_Tables" / "Supplementary_Table_3"
    out_dir.mkdir(parents=True, exist_ok=True)

    ancestry = {}
    with open(args.patients) as f:
        for line in f:
            fields = line.split()
            if len(fields) >= 4:
                ancestry[fields[0]] = fields[3]

    blca, pancan = read_driver_genes(config["resources"]["driver_genes"], config["tools"]["rscript"])

    snvs = []
    for alignment, status, folder, set_name in SNV_SETS:
        vcf = (data_dir / "Variant_Caller_Overlap" / alignment / "Shared_Variants_Annotations" / folder /
               f"{status}_filtered_shared_by_1_caller_in_capture_annotated.vcf")
        for v in read_annotated_vcf(vcf, blca | pancan):
            v.update(alignment=alignment, set=set_name)
            snvs.append(v)

    rows = []
    for v in sorted(snvs, key=lambda v: (list(ALIGNMENT_LABELS).index(v["alignment"]), v["set"], v["gene"], v["donor"], v["pos"])):
        rows.append([v["donor"], ancestry[v["donor"]], ALIGNMENT_LABELS[v["alignment"]], v["set"], f"chr{v['chrom']}",
                     v["pos"], v["gene"], v["impact"]])

    tsv = out_dir / ".Supplementary_Table_3.tsv"
    xlsx = out_dir / "Supplementary_Table_3.xlsx"
    with open(tsv, "w") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        writer.writerows(rows)
    subprocess.run([config["tools"]["rscript"], "-e",
                    f"x <- read.delim('{tsv}', check.names = FALSE, stringsAsFactors = FALSE, quote = '', na.strings = ''); "
                    f"openxlsx::write.xlsx(list('Supplementary Table 3' = x), '{xlsx}', "
                    f"keepNA = TRUE, na.string = 'NA', overwrite = TRUE)"], check=True)
    tsv.unlink()

    # Counts quoted in the text; discordant = MC3-exclusive + alignment-exclusive
    groups = [(a, s, [s]) for a, _, _, s in SNV_SETS]
    groups.append(("T2T_pangenome", "Discordant", ["MC3-exclusive", "VG_Pangenome-exclusive"]))
    groups.append(("GRCh38_bwa", "Discordant", ["MC3-exclusive", "BWA_GRCh38-exclusive"]))

    print(f"{len(rows)} SNVs -> {xlsx}")
    print("alignment\tset\timpact\tSNVs\ttumours\tmedian_per_tumour\tbladder_genes\tpancan_only_genes")
    for alignment, set_name, members in groups:
        for impact_name, impacts in [("all", None), ("high/moderate", HIGH_IMPACT)]:
            group = [v for v in snvs if v["alignment"] == alignment and v["set"] in members
                     and (not impacts or v["impact"] in impacts)]
            per_tumour = Counter(v["donor"] for v in group)
            gene_set = {v["gene"] for v in group}
            print(f"{ALIGNMENT_LABELS[alignment]}\t{set_name}\t{impact_name}\t{len(group)}\t{len(per_tumour)}\t"
                  f"{statistics.median(per_tumour.get(d, 0) for d in ancestry)}\t{len(gene_set & blca)}\t{len(gene_set - blca)}")


if __name__ == "__main__":
    main()
