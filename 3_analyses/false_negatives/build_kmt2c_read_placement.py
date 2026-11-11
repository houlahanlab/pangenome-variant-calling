#!/usr/bin/env python3

# Supplementary Figure 12 (KMT2C): where the VG pangenome places the reads that carry the ALT allele of four KMT2C
# FN SNVs (in MC3, not called by Strelka2 on either alignment).
# The ALT-carrying reads at each site are taken from the BWA GRCh38 tumour BAM, then found by name in the VG pangenome
# tumour BAM, recording their chromosome and MAPQ. Duplicate, secondary, supplementary, QC-fail and unmapped reads are
# skipped in BWA; secondary and supplementary records are skipped in the VG pangenome BAM.
#
# Input:  tumour BAMs in <alignments_dir>/GRCh38_bwa and <alignments_dir>/T2T_pangenome
# Output: <data_dir>/False_Negatives_Analysis/kmt2c_read_placement/kmt2c_alt_read_placement.txt
#
# Usage:  python3 3_analyses/false_negatives/build_kmt2c_read_placement.py   (needs pysam)

import argparse
import csv
import os
import subprocess
from collections import Counter

import pysam
import yaml

# GRCh38 position, T2T-CHM13 position (liftOver), ALT base, donor, tumour barcode
SITES = [
    (152265063, 153437953, "T", "TCGA-FD-A43P", "TCGA-FD-A43P-01A-31D-A23U-08"),
    (152265080, 153437970, "T", "TCGA-E5-A4TZ", "TCGA-E5-A4TZ-01A-11D-A31L-08"),
    (152265120, 153438010, "G", "TCGA-CF-A3MF", "TCGA-CF-A3MF-01A-12D-A21A-08"),
    (152265156, 153438046, "C", "TCGA-CF-A3MI", "TCGA-CF-A3MI-01A-11D-A20D-08"),
]


def keep(read):
    return not (read.is_duplicate or read.is_secondary or read.is_supplementary or read.is_qcfail or read.is_unmapped)


def base_at(read, pos):
    for q, p in read.get_aligned_pairs(matches_only=True):
        if p == pos - 1:
            return read.query_sequence[q]


# BWA reads carrying the ALT base: {(name, is_read1): MAPQ}
def bwa_alt_reads(bam, pos, alt):
    reads = {}
    for read in pysam.AlignmentFile(bam).fetch("chr7", pos - 1, pos):
        if keep(read) and base_at(read, pos) == alt:
            reads[(read.query_name, read.is_read1)] = read.mapping_quality
    return reads


# VG pangenome placement of the same reads, searched by read name over the whole BAM: {(name, is_read1): (chrom, MAPQ)}
def vg_placement(samtools, bam, reads, chm13_pos, out_dir):
    names = os.path.join(out_dir, "read_names.tmp")
    with open(names, "w") as f:
        f.write("\n".join({name for name, _ in reads}) + "\n")
    sam = subprocess.run([samtools, "view", "-@", "8", "-N", names, bam], capture_output=True, text=True, check=True).stdout
    os.remove(names)

    placed = {}
    for line in sam.splitlines():
        fields = line.split("\t")
        flag = int(fields[1])
        key = (fields[0], bool(flag & 64))
        if key not in reads or flag & 0x900:
            continue
        chrom = fields[2].replace("CHM13#0#", "")
        if flag & 4:
            chrom = "unmapped"
        elif chrom == "chr7" and abs(int(fields[3]) - chm13_pos) < 300:
            chrom = "chr7 (KMT2C)"
        placed[key] = (chrom, int(fields[4]))
    return placed


def mapq_class(mapq):
    return "MAPQ 0" if mapq == 0 else "MAPQ 1-19" if mapq < 20 else "MAPQ 20+"


def main():
    parser = argparse.ArgumentParser(description="VG pangenome placement of the ALT reads at four KMT2C FN SNVs.")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    alignments_dir = config["alignments_dir"]
    out_dir = os.path.join(config["data_dir"], "False_Negatives_Analysis", "kmt2c_read_placement")
    os.makedirs(out_dir, exist_ok=True)

    rows = []
    for pos, chm13_pos, alt, donor, bwa_sample in SITES:
        tumour = "-".join(bwa_sample.split("-")[:4])
        reads = bwa_alt_reads(os.path.join(alignments_dir, "GRCh38_bwa", donor, tumour, f"{bwa_sample}.bam"), pos, alt)
        placed = vg_placement(config["tools"]["samtools"],
                              os.path.join(alignments_dir, "T2T_pangenome", donor, tumour, f"{tumour}_recalibrated.bam"),
                              reads, chm13_pos, out_dir)

        bwa = Counter(("chr7 (KMT2C)", mapq_class(mapq)) for mapq in reads.values())
        vg = Counter((chrom, mapq_class(mapq)) for chrom, mapq in placed.values())
        vg["not found", "MAPQ 0"] += len(reads) - len(placed)

        for alignment, counts in [("GRCh38_bwa", bwa), ("T2T_pangenome", vg)]:
            for (chrom, mapq), n in counts.items():
                if n > 0:
                    rows.append({"site": f"chr7:{pos}", "donor": donor, "alignment": alignment,
                                 "chrom": chrom, "mapq": mapq, "reads": n})
        print(f"chr7:{pos} {donor}: {len(reads)} ALT reads, {len(placed)} found in the pangenome BAM")

    with open(os.path.join(out_dir, "kmt2c_alt_read_placement.txt"), "w") as f:
        writer = csv.DictWriter(f, fieldnames=["site", "donor", "alignment", "chrom", "mapq", "reads"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
