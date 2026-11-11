#!/usr/bin/env python3

# Variant table: combines the TP, FP and FN SNVs of every caller and alignment into one table (one row per SNV, donor, caller
# and alignment). Used by the false negative, false positive and variant impact analyses.
#
# Input:  <data_dir>/TCGA_Comparison/<caller>_vs_tcga/<alignment>/<ancestry>/<donor>/<tumour>/ from compare_to_mc3.py
# Output: <data_dir>/MegaTables/Megatable_all_alignments_all_callers.txt
#
# Usage:  python3 2_benchmarking/variant_table/build_variant_table.py

import argparse
import os

import yaml

ALIGNMENTS = ["GRCh38_bwa", "GRCh38_linear", "GRCh38_linear_plus", "T2T_pangenome", "GRCh38_hprc_v1.1"]

# Caller -> file prefix from compare_to_mc3.py
CALLERS = {"Mutect2": "mutect", "Strelka": "strelka", "SomaticSniper": "somaticsniper", "DeepSomatic": "deepsomatic"}

HEADER = ["CHROM", "START", "END", "REF", "ALT", "DONOR_ID", "TUMOUR_ID", "STATUS", "TOOL", "ALIGNMENT", "ANCESTRY"]


def short_id(sample_id):
    """TCGA-FD-A3B5-01A-11D-A20D-08 -> TCGA-FD-A3B5-01A"""
    return "-".join(sample_id.split("-")[:4])


# Writes each CHROM START END REF ALT row of one TP, FP or FN file with its labels
def write_variants(path, out, labels):
    if not os.path.exists(path):
        print(f"Warning: missing file, skipping: {path}")
        return
    with open(path) as f:
        for line in f:
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            if len(fields) < 5:
                print(f"Skipping malformed line in {path}: {line.strip()}")
                continue
            out.write("\t".join(fields[:5] + labels) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Combine TP, FP and FN SNVs of every caller and alignment.")
    parser.add_argument("--patients", default="metadata/BLCA_Patients.txt")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    out_dir = os.path.join(config["data_dir"], "MegaTables")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "Megatable_all_alignments_all_callers.txt")

    with open(args.patients) as f:
        patients = [line.split() for line in f if line.strip()]

    with open(out_path, "w") as out:
        out.write("\t".join(HEADER) + "\n")
        for entry in patients:
            if len(entry) < 4:
                print(f"Skipping incomplete line: {' '.join(entry)}")
                continue
            donor, tumour, ancestry = entry[0], entry[1], entry[3]
            for alignment in ALIGNMENTS:
                for caller, prefix in CALLERS.items():
                    comparison_dir = os.path.join(config["data_dir"], "TCGA_Comparison", f"{prefix}_vs_tcga",
                                                  alignment, ancestry, donor, short_id(tumour))
                    files = {"TP": "common_variants_full.txt", "FN": "tcga_unique_pass.txt", "FP": f"{prefix}_unique_pass.txt"}
                    for status, name in files.items():
                        write_variants(os.path.join(comparison_dir, name), out,
                                       [donor, tumour, status, caller, alignment, ancestry])

    print(f"\nMega variant table created:\n{out_path}")


if __name__ == "__main__":
    main()
