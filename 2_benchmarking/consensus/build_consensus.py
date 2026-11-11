#!/usr/bin/env python3

# Builds per-donor consensus SNV sets (called by >= 2 and >= 3 of Strelka2, Mutect2, SomaticSniper and DeepSomatic)
# for one alignment and compares them with MC3 by CHROM + POS.
#
# Input:  <data_dir>/TCGA_Comparison/<caller>_vs_tcga/<alignment>/<ancestry>/<donor>/<tumour>/ from compare_to_mc3.py
#         (calls in GRCh38 coordinates, and the MC3 truth set from the DeepSomatic folder)
# Output: <data_dir>/Consensus/MC3_Consensus/<alignment>/<ancestry>/<donor>/<tumour>/
#             consensus_ge2_positions.txt, consensus_ge3_positions.txt, consensus_sources.txt
#         <data_dir>/Stats_Summary/MC3_Consensus/stats_table_<alignment>.txt (one row per donor and consensus, appended)
#
# Usage:  python3 2_benchmarking/consensus/build_consensus.py --alignment T2T_pangenome

import argparse
import os
from collections import defaultdict

import yaml

# Caller label used in consensus_sources.txt -> file prefix from compare_to_mc3.py
CALLERS = {
    "Mutect2": "mutect",
    "Strelka2": "strelka",
    "SomaticSniper": "somaticsniper",
    "DeepSomatic": "deepsomatic",
}

ALIGNMENTS = ["GRCh38_bwa", "GRCh38_linear", "GRCh38_linear_plus", "T2T_pangenome"]

STATS_HEADER = [
    "Donor_ID", "Tumour_ID", "Ancestry", "Variant_Caller", "Reference",
    "True_Positive", "False_Positive", "False_Negative", "Total_Variant_Called", "Total_TCGA",
    "Recall(%)", "Precision(%)", "F1_Score(%)",
]


def short_id(sample_id):
    """TCGA-BT-A42C-01A-11D-A23M-08 -> TCGA-BT-A42C-01A"""
    return "-".join(sample_id.split("-")[:4])


def safe_div(n, d):
    return n / d if d else 0.0


def read_patients(path):
    """Patient list columns: donor, tumour barcode, normal barcode, ancestry."""
    patients = []
    with open(path) as f:
        for line in f:
            entry = line.split()
            if len(entry) < 4:
                if entry:
                    print(f"Skipping malformed line: {line.strip()}")
                continue
            patients.append((entry[0], entry[1], entry[3]))
    return patients


def load_positions(path):
    if not os.path.exists(path):
        print(f"WARNING: missing file: {path}")
        return set()
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())


# For T2T_pangenome the calls and the truth set have the non-callable MC3 sites removed
def calls_file(comparison_dir, prefix, alignment):
    if alignment == "T2T_pangenome":
        return os.path.join(comparison_dir, f"{prefix}_grch38_eval_positions.txt")
    return os.path.join(comparison_dir, f"{prefix}_grch38_positions.txt")


def truth_file(comparison_dir, alignment):
    if alignment == "T2T_pangenome":
        return os.path.join(comparison_dir, "tcga_patient_callable_positions.txt")
    return os.path.join(comparison_dir, "tcga_patient_target_positions.txt")


def process_patient(donor, tumour, ancestry, alignment, data_dir, stats_file):
    tumour_short = short_id(tumour)
    print(f"\nProcessing: {donor} | {tumour_short} | {ancestry} | {alignment}")

    def comparison_dir(prefix):
        return os.path.join(data_dir, "TCGA_Comparison", f"{prefix}_vs_tcga", alignment, ancestry, donor, tumour_short)

    truth = load_positions(truth_file(comparison_dir("deepsomatic"), alignment))
    if not truth:
        print(f"WARNING: truth set is empty for {donor}")

    variant_callers = defaultdict(set)
    for caller, prefix in CALLERS.items():
        variants = load_positions(calls_file(comparison_dir(prefix), prefix, alignment))
        print(f"  {caller}: {len(variants)}")
        for variant in variants:
            variant_callers[variant].add(caller)

    out_dir = os.path.join(data_dir, "Consensus", "MC3_Consensus", alignment, ancestry, donor, tumour_short)
    os.makedirs(out_dir, exist_ok=True)

    consensus = {}
    for n in [2, 3]:
        consensus[n] = {v for v, callers in variant_callers.items() if len(callers) >= n}
        with open(os.path.join(out_dir, f"consensus_ge{n}_positions.txt"), "w") as f:
            for variant in sorted(consensus[n]):
                f.write(variant + "\n")

    with open(os.path.join(out_dir, "consensus_sources.txt"), "w") as f:
        f.write("Variant\tNum_Callers\tCallers\n")
        for variant in sorted(variant_callers):
            callers = sorted(variant_callers[variant])
            f.write(f"{variant}\t{len(callers)}\t{'+'.join(callers)}\n")

    for n in [2, 3]:
        tp = len(consensus[n] & truth)
        fp = len(consensus[n] - truth)
        fn = len(truth - consensus[n])
        precision = safe_div(tp, tp + fp)
        recall = safe_div(tp, tp + fn)
        f1 = safe_div(2 * precision * recall, precision + recall)

        write_header = not os.path.exists(stats_file) or os.path.getsize(stats_file) == 0
        with open(stats_file, "a") as f:
            if write_header:
                f.write("\t".join(STATS_HEADER) + "\n")
            f.write("\t".join([
                donor, tumour, ancestry, f"Consensus_ge{n}", alignment,
                str(tp), str(fp), str(fn), str(tp + fp), str(tp + fn),
                f"{recall * 100:.2f}", f"{precision * 100:.2f}", f"{f1 * 100:.2f}",
            ]) + "\n")
        print(f"Consensus_ge{n}: TP={tp}, FP={fp}, FN={fn}, "
              f"Precision={precision:.3f}, Recall={recall:.3f}, F1={f1:.3f}")


def main():
    parser = argparse.ArgumentParser(description="Build consensus SNV sets and compare them with MC3.")
    parser.add_argument("--alignment", choices=ALIGNMENTS, required=True)
    parser.add_argument("--patients", default="metadata/BLCA_Patients.txt")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    stats_dir = os.path.join(config["data_dir"], "Stats_Summary", "MC3_Consensus")
    os.makedirs(stats_dir, exist_ok=True)
    stats_file = os.path.join(stats_dir, f"stats_table_{args.alignment}.txt")

    for donor, tumour, ancestry in read_patients(args.patients):
        process_patient(donor, tumour, ancestry, args.alignment, config["data_dir"], stats_file)
    print(f"\nAppended metrics to: {stats_file}")


if __name__ == "__main__":
    main()
