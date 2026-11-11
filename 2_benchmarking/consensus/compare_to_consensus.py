#!/usr/bin/env python3

# Compares Strelka2 and Mutect2 calls with a consensus truth set instead of MC3: the SNVs called by >= 3 of the
# four callers in one alignment (from build_consensus.py). Matching is by CHROM + POS in GRCh38 coordinates.
#   TP = calls in the consensus, FP = calls not in the consensus, FN = consensus SNVs not called
#
# Input:  <data_dir>/Consensus/MC3_Consensus/<truth>/<ancestry>/<donor>/<tumour>/consensus_ge3_positions.txt
#         <data_dir>/TCGA_Comparison/<caller>_vs_tcga/<alignment>/... from compare_to_mc3.py
# Output: <data_dir>/Stats_Summary/<BWA_GRCh38|VG_Pangenome>_Consensus/
#             stats_table_<alignment>_vs_consensus_ge3.txt (rewritten on each run)
#             Per_Patient/<alignment>/<ancestry>/<donor>/<tumour>/<caller>/ (shared, unique and missed SNVs)
#
# Usage:  python3 2_benchmarking/consensus/compare_to_consensus.py --truth T2T_pangenome --alignment GRCh38_bwa

import argparse
import os

import yaml

CALLERS = {"Mutect2": "mutect", "Strelka": "strelka"}

TRUTH_DIRS = {"GRCh38_bwa": "BWA_GRCh38_Consensus", "T2T_pangenome": "VG_Pangenome_Consensus"}

ALIGNMENTS = ["GRCh38_bwa", "GRCh38_linear", "GRCh38_linear_plus", "T2T_pangenome"]

STATS_HEADER = [
    "Donor_ID", "Tumour_ID", "Ancestry", "Variant_Caller", "Reference", "Ground_Truth",
    "True_Positive", "False_Positive", "False_Negative", "Total_Variant_Called", "Total_Consensus_ge3",
    "Recall(%)", "Precision(%)", "F1_Score(%)",
]


def short_id(sample_id):
    """TCGA-FD-A3B5-01A-11D-A20D-08 -> TCGA-FD-A3B5-01A"""
    return "-".join(sample_id.split("-")[:4])


def safe_div(n, d):
    return n / d if d else 0.0


def read_patients(path):
    """Patient list columns: donor, tumour barcode, normal barcode, ancestry."""
    patients = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            entry = line.split()
            if len(entry) < 4:
                print(f"Skipping malformed line: {line.strip()}")
                continue
            patients.append((entry[0], entry[1], entry[3]))
    return patients


def load_positions(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing position file: {path}")
    with open(path) as f:
        return {line.strip() for line in f if line.strip() and not line.startswith("#")}


def write_positions(variants, path):
    with open(path, "w") as f:
        for variant in sorted(variants):
            f.write(variant + "\n")


def process_patient(donor, tumour, ancestry, args, data_dir, stats_dir, stats_file):
    tumour_short = short_id(tumour)
    print(f"\nProcessing {donor} | {tumour_short} | {ancestry} | {args.alignment}")

    consensus = load_positions(os.path.join(data_dir, "Consensus", "MC3_Consensus", args.truth, ancestry,
                                            donor, tumour_short, "consensus_ge3_positions.txt"))
    if not consensus:
        raise ValueError(f"Consensus >=3 is empty for {donor}")

    patient_dir = os.path.join(stats_dir, "Per_Patient", args.alignment, ancestry, donor, tumour_short)

    for caller, prefix in CALLERS.items():
        # For T2T_pangenome the calls have the non-callable MC3 sites removed
        name = f"{prefix}_grch38_eval_positions.txt" if args.alignment == "T2T_pangenome" else f"{prefix}_grch38_positions.txt"
        calls = load_positions(os.path.join(data_dir, "TCGA_Comparison", f"{prefix}_vs_tcga", args.alignment,
                                            ancestry, donor, tumour_short, name))

        shared, unique, missed = calls & consensus, calls - consensus, consensus - calls
        tp, fp, fn = len(shared), len(unique), len(missed)
        precision = safe_div(tp, tp + fp)
        recall = safe_div(tp, tp + fn)
        f1 = safe_div(2 * precision * recall, precision + recall)

        caller_dir = os.path.join(patient_dir, caller)
        os.makedirs(caller_dir, exist_ok=True)
        write_positions(shared, os.path.join(caller_dir, f"{caller}_shared_with_consensus_ge3.txt"))
        write_positions(unique, os.path.join(caller_dir, f"{caller}_unique_vs_consensus_ge3.txt"))
        write_positions(missed, os.path.join(caller_dir, "consensus_ge3_not_called_by_caller.txt"))

        with open(stats_file, "a") as f:
            f.write("\t".join([
                donor, tumour, ancestry, caller, args.alignment, f"{args.truth}_Consensus_ge3",
                str(tp), str(fp), str(fn), str(tp + fp), str(tp + fn),
                f"{recall * 100:.2f}", f"{precision * 100:.2f}", f"{f1 * 100:.2f}",
            ]) + "\n")
        print(f"{caller}: TP={tp}, FP={fp}, FN={fn}, Precision={precision:.3f}, Recall={recall:.3f}, F1={f1:.3f}")


def main():
    parser = argparse.ArgumentParser(description="Compare Strelka2 and Mutect2 calls with a consensus truth set.")
    parser.add_argument("--truth", choices=TRUTH_DIRS, required=True, help="alignment the consensus was built from")
    parser.add_argument("--alignment", choices=ALIGNMENTS, required=True, help="alignment of the calls")
    parser.add_argument("--patients", default="metadata/BLCA_Patients.txt")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    stats_dir = os.path.join(config["data_dir"], "Stats_Summary", TRUTH_DIRS[args.truth])
    os.makedirs(stats_dir, exist_ok=True)
    stats_file = os.path.join(stats_dir, f"stats_table_{args.alignment}_vs_consensus_ge3.txt")
    with open(stats_file, "w") as f:
        f.write("\t".join(STATS_HEADER) + "\n")

    for donor, tumour, ancestry in read_patients(args.patients):
        try:
            process_patient(donor, tumour, ancestry, args, config["data_dir"], stats_dir, stats_file)
        except (FileNotFoundError, ValueError) as error:
            print(f"Skipping {donor}: {error}")

    print(f"\nStats table saved to: {stats_file}")


if __name__ == "__main__":
    main()
