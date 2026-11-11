#!/usr/bin/env python3

# Figure 4, Supplementary Figure 12, Supplementary Table 3: SnpEff annotation of the Strelka2 SNV sets compared with MC3:
#   VG pangenome TP (shared with MC3), FN (MC3-exclusive) and FP (VG pangenome-exclusive)
#   BWA GRCh38 FN (MC3-exclusive) and FP (BWA GRCh38-exclusive)
# Each set is de-duplicated per donor, kept if it lies in the Agilent SureSelect V5 targets, written as a minimal VCF
# and annotated with SnpEff (GRCh38.86)
#
# Input:  <data_dir>/MegaTables/Megatable_all_alignments_all_callers.txt (build_variant_table.py)
# Output: <data_dir>/Variant_Caller_Overlap/<alignment>/Shared_Variants_Annotations/<True|False>_<Positives|Negatives>/
#             <status>_filtered_shared_by_1_caller_in_capture_annotated.vcf (and the files before annotation)
#
# Usage:  python3 3_analyses/impact_variants/annotate_variant_sets.py

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path

import yaml

SNPEFF_DB = "GRCh38.86"
CALLERS = ["Strelka"]
COLUMNS = ["CHROM", "START", "END", "REF", "ALT", "SAMPLE", "TOOL"]

# Alignment, status, folder
SETS = [
    ("GRCh38_bwa", "FP", "False_Positives"),
    ("GRCh38_bwa", "FN", "False_Negatives"),
    ("T2T_pangenome", "FP", "False_Positives"),
    ("T2T_pangenome", "FN", "False_Negatives"),
    ("T2T_pangenome", "TP", "True_Positives"),
]


# Variants for one status and alignment; the same variant called by several callers is kept once
def get_variants(variant_table, status, alignment):
    variants = defaultdict(set)
    for row in variant_table:
        if row["STATUS"] != status or row["ALIGNMENT"] != alignment or row["TOOL"] not in CALLERS:
            continue
        key = tuple(row[c].strip() for c in ["CHROM", "START", "END", "REF", "ALT", "DONOR_ID"])
        if "" not in key:
            variants[key].add(row["TOOL"])
    return [dict(zip(COLUMNS, key + ("+".join(sorted(variants[key])),)))
            for key in sorted(variants, key=lambda k: (k[0], int(k[1]), k[5], k[3], k[4]))]


def write_table(rows, path):
    with open(path, "w") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


# Keeps the variants inside the Agilent targets (bedtools intersect); the Agilent BED uses chr names
def capture_filter(rows, bedtools, target_bed, bed, captured_bed):
    with open(bed, "w") as f:
        for row in rows:
            f.write(f"chr{row['CHROM']}\t{row['START']}\t{row['END']}\n")
    with open(captured_bed, "w") as f:
        subprocess.run([bedtools, "intersect", "-a", bed, "-b", target_bed, "-u"], stdout=f, check=True)
    with open(captured_bed) as f:
        captured = {(c.removeprefix("chr"), s, e) for c, s, e in (line.split("\t")[:3] for line in f.read().splitlines())}
    return [row for row in rows if (row["CHROM"], row["START"], row["END"]) in captured]


# START is 0-based (BED), so the 1-based VCF POS is END
def write_vcf(rows, path):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in rows:
            f.write(f"{row['CHROM']}\t{row['END']}\t.\t{row['REF']}\t{row['ALT']}\t.\tPASS\t"
                    f"SAMPLE={row['SAMPLE']};TOOL={row['TOOL']}\n")


# Run from the VCF's folder so snpEff_summary.html and snpEff_genes.txt land next to it
def run_snpeff(snpeff_jar, vcf, annotated_vcf):
    with open(annotated_vcf, "w") as f:
        subprocess.run(["java", "-Xmx4g", "-jar", snpeff_jar, SNPEFF_DB, vcf.name], cwd=vcf.parent, stdout=f, check=True)


def main():
    parser = argparse.ArgumentParser(description="SnpEff annotation of the Strelka2 TP, FN and FP SNV sets.")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    data_dir = Path(config["data_dir"])
    with open(data_dir / "MegaTables" / "Megatable_all_alignments_all_callers.txt") as f:
        variant_table = list(csv.DictReader(f, delimiter="\t"))

    for alignment, status, folder in SETS:
        out_dir = data_dir / "Variant_Caller_Overlap" / alignment / "Shared_Variants_Annotations" / folder
        out_dir.mkdir(parents=True, exist_ok=True)

        rows = get_variants(variant_table, status, alignment)
        write_table(rows, out_dir / f"{status}_shared_by_1_caller_noduplicates.txt")
        kept = capture_filter(rows, config["tools"]["bedtools"], config["resources"]["target_bed_grch38"],
                              out_dir / f"{status}_shared_by_1_caller_noduplicates.bed",
                              out_dir / f"{status}_shared_by_1_caller_in_capture.bed")
        write_table(kept, out_dir / f"{status}_filtered_shared_by_1_caller_in_capture.txt")
        vcf = out_dir / f"{status}_filtered_shared_by_1_caller_in_capture.vcf"
        write_vcf(kept, vcf)
        run_snpeff(config["tools"]["snpeff"], vcf, out_dir / f"{status}_filtered_shared_by_1_caller_in_capture_annotated.vcf")
        print(f"{alignment} {status}: {len(rows)} variants, {len(kept)} in the capture targets")


if __name__ == "__main__":
    main()
