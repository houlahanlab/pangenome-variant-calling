#!/usr/bin/env python3

# Converts the TCGA tumour and normal BAMs (GDC, BWA GRCh38) back to paired FASTQs for alignment with vg giraffe.
# Each BAM is name-sorted (samtools sort -n), written as R1 / R2 FASTQs (samtools fastq -n), and the sorted BAM is
# deleted. Samples whose FASTQs already exist are skipped.
#
# Input:  <tcga_bams_dir>/<donor>/<sample>/<barcode>.bam
# Output: <fastq_dir>/<barcode>_R1.fastq, <barcode>_R2.fastq
#
# Usage:  python3 1_pipeline/bam_to_fastq.py --patients metadata/BLCA_Patients.txt

import argparse
import os
import subprocess

import yaml


def run_command(cmd):
    print(f"Running: {cmd}")
    if subprocess.run(cmd, shell=True).returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}")


def short_id(sample_id):
    """TCGA-BT-A2LD-01A-12D-A20D-08 -> TCGA-BT-A2LD-01A"""
    return "-".join(sample_id.split("-")[:4])


def convert_sample(donor, barcode, config):
    samtools = config["tools"]["samtools"]
    fastq_dir = config["fastq_dir"]
    bam = os.path.join(config["tcga_bams_dir"], donor, short_id(barcode), f"{barcode}.bam")
    r1 = os.path.join(fastq_dir, f"{barcode}_R1.fastq")
    r2 = os.path.join(fastq_dir, f"{barcode}_R2.fastq")

    if not os.path.exists(bam):
        print(f"Missing BAM: {bam}")
        return
    if os.path.exists(r1) and os.path.exists(r2):
        print(f"FASTQs already exist for {barcode}")
        return

    sorted_bam = os.path.join(fastq_dir, f"{barcode}_sorted.bam")
    run_command(f"{samtools} sort -n {bam} -o {sorted_bam}")
    run_command(f"{samtools} fastq -1 {r1} -2 {r2} -n {sorted_bam}")
    os.remove(sorted_bam)


def main():
    parser = argparse.ArgumentParser(description="Convert the TCGA BAMs to paired FASTQs.")
    parser.add_argument("--patients", default="metadata/BLCA_Patients.txt",
                        help="patient list: donor, tumour barcode, normal barcode, ...")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    os.makedirs(config["fastq_dir"], exist_ok=True)

    with open(args.patients) as f:
        for line in f:
            fields = line.split()
            if len(fields) < 3:
                continue
            donor, tumour, normal = fields[:3]
            print(f"\nDonor: {donor}")
            for barcode in (tumour, normal):
                convert_sample(donor, barcode, config)


if __name__ == "__main__":
    main()
