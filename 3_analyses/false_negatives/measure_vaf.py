#!/usr/bin/env python3

# Supplementary Figure 11: tumour REF / ALT read counts, VAF and depth at every Strelka2 and Mutect2 TP and FN SNV.
# The sites come from the TP and FN SNVs of BWA GRCh38, VG GRCh38 and VG pangenome (Defined_By). Every site is
# measured in three tumour BAMs in GRCh38 coordinates (Counted_In), with samtools mpileup at the site:
#   BWA_GRCh38_linear  GRCh38_bwa BAM
#   VG_GRCh38_linear   GRCh38_linear BAM (VG GRCh38)
#   VG_T2T_pangenome   T2T_pangenome BAM
# Reads are counted as REF, ALT or another base; VAF = ALT / (REF + ALT).
#
# Input:  <data_dir>/MegaTables/Megatable_all_alignments_all_callers.txt (build_variant_table.py)
#         tumour BAMs in <alignments_dir>
# Output: <data_dir>/False_Negatives_Analysis/VAF_results.tsv (one row per site and BAM)
#
# Usage:  python3 3_analyses/false_negatives/measure_vaf.py

import argparse
import csv
import os
import subprocess

import yaml

STATUSES = {"TP", "FN"}
TOOLS = {"Strelka", "Mutect2"}

# Alignment in the variant table -> Defined_By label
ALIGNMENTS = {"T2T_pangenome": "VG_T2T_pangenome", "GRCh38_bwa": "BWA_GRCh38_linear", "GRCh38_linear": "VG_GRCh38_linear"}

COLUMNS = ["Variant_ID", "CHROM", "START", "END", "REF", "ALT", "DONOR_ID", "TUMOUR_ID", "STATUS", "TOOL", "Defined_By",
           "ANCESTRY", "Counted_In", "BAM", "REF_COUNT", "ALT_COUNT", "OTHER_COUNT", "TOTAL_DEPTH", "MPILEUP_DEPTH",
           "REF_ALT_DEPTH", "VAF", "ALT_FRACTION_DEPTH", "REF_FRACTION_DEPTH", "OTHER_FRACTION_DEPTH"]


def short_id(sample_id):
    """TCGA-FD-A3N5-01A-11D-A21A-08 -> TCGA-FD-A3N5-01A"""
    return "-".join(sample_id.split("-")[:4])


# The three tumour BAMs each site is measured in (all GRCh38 coordinates)
def bam_paths(alignments_dir, donor, tumour):
    short = short_id(tumour)
    return {
        "VG_GRCh38_linear": os.path.join(alignments_dir, "GRCh38_linear", donor, short, f"{short}_recalibrated.bam"),
        "VG_T2T_pangenome": os.path.join(alignments_dir, "GRCh38_hprc_v1.1", donor, short, f"{short}_recalibrated.bam"),
        "BWA_GRCh38_linear": os.path.join(alignments_dir, "GRCh38_bwa", donor, short, f"{tumour}.bam"),
    }


# mpileup bases without read starts (^ + MAPQ), read ends ($) and indels; keeps . , A C G T N (either case)
def clean_pileup_bases(bases):
    cleaned, i = [], 0
    while i < len(bases):
        char = bases[i]
        if char == "^":
            i += 2
            continue
        if char == "$":
            i += 1
            continue
        if char in "+-":
            i += 1
            number = ""
            while i < len(bases) and bases[i].isdigit():
                number += bases[i]
                i += 1
            if number:
                i += int(number)
            continue
        if char in ".,ACGTNacgtn":
            cleaned.append(char)
        i += 1
    return cleaned


def count_ref_alt(pileup_bases, ref, alt, mpileup_ref):
    ref, alt, mpileup_ref = ref.upper(), alt.upper(), mpileup_ref.upper()
    ref_count = alt_count = other_count = 0
    for base in clean_pileup_bases(pileup_bases):
        # . and , are matches to the reference base in the pileup
        observed = mpileup_ref if base in ".," else base.upper()
        if observed == ref:
            ref_count += 1
        elif observed == alt:
            alt_count += 1
        elif observed in "ACGT":
            other_count += 1

    total = ref_count + alt_count + other_count
    ref_alt = ref_count + alt_count
    return {
        "REF_COUNT": ref_count, "ALT_COUNT": alt_count, "OTHER_COUNT": other_count,
        "TOTAL_DEPTH": total, "REF_ALT_DEPTH": ref_alt,
        "VAF": alt_count / ref_alt if ref_alt > 0 else 0,
        "ALT_FRACTION_DEPTH": alt_count / total if total > 0 else 0,
        "REF_FRACTION_DEPTH": ref_count / total if total > 0 else 0,
        "OTHER_FRACTION_DEPTH": other_count / total if total > 0 else 0,
    }


def measure_site(samtools, bam, chrom, pos, ref, alt):
    chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    result = subprocess.run(f"{samtools} mpileup -r {chrom}:{pos}-{pos} {bam}", shell=True,
                            capture_output=True, text=True)
    output = result.stdout.strip() if result.returncode == 0 else ""
    if result.returncode != 0:
        print(f"Command failed: samtools mpileup {chrom}:{pos} {bam}\n{result.stderr.strip()}")

    fields = output.splitlines()[0].split("\t") if output else []
    if len(fields) < 5:
        counts = count_ref_alt("", ref, alt, "N")
        counts["MPILEUP_DEPTH"] = 0
        return counts
    counts = count_ref_alt(fields[4], ref, alt, fields[2])
    try:
        counts["MPILEUP_DEPTH"] = int(fields[3])
    except ValueError:
        counts["MPILEUP_DEPTH"] = 0
    return counts


def main():
    parser = argparse.ArgumentParser(description="REF / ALT read counts and VAF at the TP and FN SNVs.")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    input_table = os.path.join(config["data_dir"], "MegaTables", "Megatable_all_alignments_all_callers.txt")
    out_dir = os.path.join(config["data_dir"], "False_Negatives_Analysis")
    os.makedirs(out_dir, exist_ok=True)
    output_file = os.path.join(out_dir, "VAF_results.tsv")

    n_rows = n_written = n_non_snv = n_missing_bams = 0
    with open(input_table) as infile, open(output_file, "w") as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        writer = csv.DictWriter(outfile, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        for row in reader:
            status, alignment = row["STATUS"].strip(), row["ALIGNMENT"].strip()
            tool = "Strelka" if row["TOOL"].strip() in {"Strelka", "Strelka2"} else row["TOOL"].strip()
            if status not in STATUSES or tool not in TOOLS or alignment not in ALIGNMENTS:
                continue
            n_rows += 1

            chrom, ref, alt = row["CHROM"].strip(), row["REF"].strip(), row["ALT"].strip()
            try:
                start, end = int(row["START"]), int(row["END"])
            except ValueError:
                print(f"Skipping row with invalid START or END: {row}")
                continue
            if len(ref) != 1 or len(alt) != 1:
                n_non_snv += 1
                continue

            donor, tumour = row["DONOR_ID"].strip(), row["TUMOUR_ID"].strip()
            for counted_in, bam in bam_paths(config["alignments_dir"], donor, tumour).items():
                if not os.path.exists(bam):
                    n_missing_bams += 1
                    print(f"Missing BAM file: {bam}")
                    continue
                counts = measure_site(config["tools"]["samtools"], bam, chrom, end, ref, alt)
                for key in ["VAF", "ALT_FRACTION_DEPTH", "REF_FRACTION_DEPTH", "OTHER_FRACTION_DEPTH"]:
                    counts[key] = f"{counts[key]:.6f}"
                writer.writerow({
                    "Variant_ID": f"{chrom}:{end}:{ref}:{alt}", "CHROM": chrom, "START": start, "END": end,
                    "REF": ref, "ALT": alt, "DONOR_ID": donor, "TUMOUR_ID": tumour, "STATUS": status, "TOOL": tool,
                    "Defined_By": ALIGNMENTS[alignment], "ANCESTRY": row["ANCESTRY"].strip(),
                    "Counted_In": counted_in, "BAM": bam, **counts,
                })
                n_written += 1

    print(f"TP / FN rows of Strelka2 and Mutect2: {n_rows:,}")
    print(f"Skipped non-SNVs: {n_non_snv:,}")
    print(f"Output rows written: {n_written:,}")
    print(f"Missing BAMs: {n_missing_bams:,}")
    print(f"Output: {output_file}")


if __name__ == "__main__":
    main()
