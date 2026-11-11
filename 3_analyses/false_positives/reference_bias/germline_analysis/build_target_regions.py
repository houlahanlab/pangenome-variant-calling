#!/usr/bin/env python3

# Reference bias, step 1: HaplotypeCaller target regions around the Strelka2 VG pangenome somatic calls.
#
# For each donor, the Strelka2 VG pangenome TP + FP SNVs (GRCh38 coordinates) are lifted to T2T-CHM13, padded by
# 74 bp on each side (a germline and a somatic base 74 bp apart still fit on one read), sorted and merged.
# HaplotypeCaller is then run on the VG pangenome normal BAM in these regions (run_haplotypecaller_array.sh),
# and its heterozygous calls become the germline sites in paired_sites/build_paired_sites.py.
#
# Input:  <data_dir>/TCGA_Comparison/strelka_vs_tcga/T2T_pangenome/<ancestry>/<donor>/<tumour>/
#             common_variants_full.txt (TP), strelka_unique_pass.txt (FP) from compare_to_mc3.py
# Output: <data_dir>/Reference_Bias/germline_analysis/
#             regions/vg_pangenome/<donor>/target_regions_merged.bed (and the lifted sites)
#             manifest.tsv (one HaplotypeCaller task per donor), qc_partial.tsv
#         Prints the sbatch command for run_haplotypecaller_array.sh
#
# Usage:  python3 3_analyses/false_positives/reference_bias/germline_analysis/build_target_regions.py

import argparse
import os
import subprocess

import yaml

DISTANCE = 74


def short_id(sample_id):
    """TCGA-FD-A3N6-01A-11D-A21A-08 -> TCGA-FD-A3N6-01A"""
    return "-".join(sample_id.split("-")[:4])


def run_command(cmd):
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def read_patients(path):
    """Patient list columns: donor, tumour barcode, normal barcode, ancestry."""
    patients = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 4:
                print(f"Skipping malformed patient line: {line.strip()}")
                continue
            patients.append({"donor": fields[0], "tumour": fields[1], "normal": fields[2], "ancestry": fields[3]})
    return patients


# TP + FP SNVs as {chrom, pos, status}, one per position (TP wins); also the TP and FP line counts
def load_tp_fp_sites(comparison_dir):
    site_status, counts = {}, {"TP": 0, "FP": 0}
    for name, status in [("common_variants_full.txt", "TP"), ("strelka_unique_pass.txt", "FP")]:
        path = os.path.join(comparison_dir, name)
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing {status} file: {path}")
        with open(path) as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                counts[status] += 1
                key = (fields[0], int(fields[2]))
                if key not in site_status or status == "TP":
                    site_status[key] = status
    sites = [{"chrom": chrom, "pos": pos, "status": status} for (chrom, pos), status in site_status.items()]
    return sites, counts["TP"], counts["FP"]


# GRCh38 -> T2T. The GRCh38 side of this chain is also named CHM13#0#chr<N>, so the input BED uses that name
def liftover_to_t2t(sites, work_dir, liftover, chain):
    input_bed = os.path.join(work_dir, "somatic_sites_grch38_for_liftover.bed")
    output_bed = os.path.join(work_dir, "somatic_sites_t2t.bed")
    unmapped_bed = os.path.join(work_dir, "liftover_unmapped.bed")

    with open(input_bed, "w") as out:
        for i, site in enumerate(sites):
            out.write(f"CHM13#0#chr{site['chrom']}\t{site['pos'] - 1}\t{site['pos']}\tVARIANT_{i}\n")
    run_command([liftover, input_bed, chain, output_bed, unmapped_bed])

    lifted, seen = [], set()
    with open(output_bed) as f:
        for line in f:
            if not line.strip():
                continue
            chrom, _, end, variant_id = line.rstrip("\n").split("\t")[:4]
            if variant_id in seen:
                print(f"WARNING: liftOver produced >1 output for {variant_id}; keeping first.")
                continue
            seen.add(variant_id)
            site = sites[int(variant_id.split("_")[1])]
            lifted.append({"chrom": chrom, "pos": int(end), "status": site["status"],
                           "orig_grch38_chrom": site["chrom"], "orig_grch38_pos": site["pos"]})

    with open(unmapped_bed) as f:
        n_unmapped = sum(1 for line in f if line.strip() and not line.startswith("#"))
    return lifted, n_unmapped


# +/- 74 bp around each site, sorted and merged with bedtools
def build_target_bed(sites, work_dir, bedtools):
    padded = os.path.join(work_dir, "target_regions_padded.bed")
    sorted_bed = os.path.join(work_dir, "target_regions_sorted.bed")
    merged = os.path.join(work_dir, "target_regions_merged.bed")
    with open(padded, "w") as out:
        for site in sites:
            out.write(f"{site['chrom']}\t{max(0, site['pos'] - 1 - DISTANCE)}\t{site['pos'] + DISTANCE}\n")
    run_command(["bash", "-c", f"sort -k1,1 -k2,2n {padded} > {sorted_bed}"])
    run_command(["bash", "-c", f"{bedtools} merge -i {sorted_bed} > {merged}"])
    with open(merged) as f:
        n_merged = sum(1 for line in f if line.strip())
    return merged, n_merged


def write_table(rows, path, columns):
    with open(path, "w") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(str(row[c]) for c in columns) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Build HaplotypeCaller target regions for the reference bias analysis.")
    parser.add_argument("--patients", default="metadata/BLCA_Patients.txt")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    tools, resources = config["tools"], config["resources"]

    germline_dir = os.path.join(config["data_dir"], "Reference_Bias", "germline_analysis")
    manifest_file = os.path.join(germline_dir, "manifest.tsv")
    qc_file = os.path.join(germline_dir, "qc_partial.tsv")

    manifest, qc_rows = [], []
    for patient in read_patients(args.patients):
        donor, tumour_short, normal_short = patient["donor"], short_id(patient["tumour"]), short_id(patient["normal"])
        work_dir = os.path.join(germline_dir, "regions", "vg_pangenome", donor)
        os.makedirs(work_dir, exist_ok=True)

        comparison_dir = os.path.join(config["data_dir"], "TCGA_Comparison", "strelka_vs_tcga", "T2T_pangenome",
                                      patient["ancestry"], donor, tumour_short)
        try:
            sites, n_tp, n_fp = load_tp_fp_sites(comparison_dir)
        except FileNotFoundError as error:
            print(f"WARNING: {error}\nSkipping {donor}.")
            continue

        lifted, n_unmapped = liftover_to_t2t(sites, work_dir, tools["liftover"], resources["chain_grch38_to_t2t"])
        exact_sites_file = os.path.join(work_dir, "exact_somatic_sites.tsv")
        write_table(lifted, exact_sites_file, ["chrom", "pos", "status", "orig_grch38_chrom", "orig_grch38_pos"])
        target_bed, n_merged = build_target_bed(lifted, work_dir, tools["bedtools"])

        normal_bam = os.path.join(config["alignments_dir"], "T2T_pangenome", donor, normal_short,
                                  f"{normal_short}_recalibrated.bam")
        if not os.path.exists(normal_bam):
            print(f"WARNING: normal BAM not found for {donor}: {normal_bam}")

        manifest.append({
            "array_index": len(manifest), "donor": donor, "alignment": "vg_pangenome", "ancestry": patient["ancestry"],
            "normal_bam": normal_bam, "reference_fasta": resources["fasta_t2t"], "target_bed": target_bed,
            "out_prefix": os.path.join(germline_dir, "vg_pangenome", donor, f"{donor}_vg_pangenome"),
            "exact_sites_file": exact_sites_file,
        })
        qc_rows.append({
            "Donor": donor, "Ancestry": patient["ancestry"], "Alignment": "vg_pangenome",
            "Number_of_TP_sites": n_tp, "Number_of_FP_sites": n_fp, "Number_of_TP_FP_sites": len(sites),
            "Number_of_target_regions_before_merge": len(lifted), "Number_of_target_regions_after_merge": n_merged,
            "Number_of_somatic_sites_lifted_to_T2T": len(lifted), "Number_of_somatic_sites_unmapped_by_liftover": n_unmapped,
        })
        print(f"{donor}: TP={n_tp} FP={n_fp} sites={len(sites)} lifted={len(lifted)} regions={n_merged}")

    write_table(manifest, manifest_file, ["array_index", "donor", "alignment", "ancestry", "normal_bam",
                                          "reference_fasta", "target_bed", "out_prefix", "exact_sites_file"])
    write_table(qc_rows, qc_file, list(qc_rows[0]))

    log_dir = os.path.join(germline_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    print(f"\nManifest: {manifest_file} ({len(manifest)} tasks)")
    print("Next, run HaplotypeCaller:")
    print(f"  sbatch {config.get('slurm_options', '')} --array=0-{len(manifest) - 1} --output={log_dir}/%x_%A_%a.out --error={log_dir}/%x_%A_%a.err "
          f"--export=ALL,MANIFEST={manifest_file},GATK={tools['gatk']},BCFTOOLS={tools['bcftools']} "
          f"3_analyses/false_positives/reference_bias/germline_analysis/run_haplotypecaller_array.sh")


if __name__ == "__main__":
    main()
