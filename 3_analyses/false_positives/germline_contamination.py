#!/usr/bin/env python3

# Supplementary Figures 8 and 9, Supplementary Table 1: which somatic SNVs are known germline variants in gnomAD.
# Every SNV in the variant table is compared with the gnomAD af-only biallelic SNVs (GRCh38) using bedtools closest:
#   GNOMAD_POSITION_OVERLAP  a gnomAD SNV at the same position
#   GNOMAD_EXACT_MATCH       the same position, REF and ALT
#   COMMON_GERMLINE(_EXACT)  ... with gnomAD AF >= 0.01. The paper calls this MAF >= 0.01: AF and MAF only differ
#                            when AF > 0.99, and no Strelka2 or Mutect2 call has AF above 0.989
# The gnomAD VCF is converted to a BED once (9.4 GB) and reused.
#
# Input:  <data_dir>/MegaTables/Megatable_all_alignments_all_callers.txt (build_variant_table.py), gnomAD from config.yaml
# Output: <data_dir>/False_Positive_Analysis/Germline_Contamination/
#             Megatable_all_alignments_all_callers_germline_annotated.tsv (the variant table plus the columns above)
#             germline_annotation_QC.tsv, unique_megatable_variants.bed, unique_megatable_variants_gnomad_closest.tsv,
#             af-only-gnomad.hg38.biallelic_snvs.bed
#
# Usage:  python3 3_analyses/false_positives/germline_contamination.py

import argparse
import csv
import os
import subprocess
from collections import defaultdict

import yaml

COMMON_AF_THRESHOLD = 0.01

REQUIRED_COLUMNS = ["CHROM", "START", "END", "REF", "ALT", "DONOR_ID", "TUMOUR_ID", "STATUS", "TOOL", "ALIGNMENT", "ANCESTRY"]

NEW_COLUMNS = ["DISTANCE_TO_CLOSEST_GNOMAD", "GNOMAD_POSITION_OVERLAP", "GNOMAD_POSITION_VARIANTS", "COMMON_GERMLINE",
               "GNOMAD_EXACT_MATCH", "GNOMAD_AF_EXACT", "COMMON_GERMLINE_EXACT"]


def run_command(command, output_file=None):
    print("Running:", " ".join(command))
    if output_file is None:
        subprocess.run(command, check=True)
    else:
        with open(output_file, "w") as handle:
            subprocess.run(command, stdout=handle, text=True, check=True)


def standardize_chromosome(chrom):
    chrom = str(chrom).strip()
    if chrom.startswith("chr"):
        return chrom
    if chrom in {"M", "MT"}:
        return "chrM"
    return f"chr{chrom}"


def chromosome_order(chrom):
    value = standardize_chromosome(chrom).replace("chr", "", 1)
    if value.isdigit() and 1 <= int(value) <= 22:
        return int(value)
    return {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(value, 999)


def is_valid_snv(row):
    try:
        start, end = int(row["START"]), int(row["END"])
    except (ValueError, TypeError):
        return False
    ref, alt = row["REF"].strip().upper(), row["ALT"].strip().upper()
    return (start >= 0 and end - start == 1 and len(ref) == 1 and len(alt) == 1 and "," not in alt
            and chromosome_order(row["CHROM"]) != 999)


def parse_af(value):
    if value in {"", ".", "NA"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def format_af(value):
    return "NA" if value is None else f"{value:.10g}"


def bool_text(value):
    return "TRUE" if value else "FALSE"


def read_variant_table(path):
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise RuntimeError("Variant table is missing columns: " + ", ".join(missing))
        rows = []
        for row in reader:
            row["_ELIGIBLE"] = is_valid_snv(row)
            row["_VARIANT_ID"] = (f"{standardize_chromosome(row['CHROM'])}:{int(row['END'])}:"
                                  f"{row['REF'].strip().upper()}:{row['ALT'].strip().upper()}") if row["_ELIGIBLE"] else None
            rows.append(row)
    return rows, reader.fieldnames


# One sorted BED row per unique SNV: chrom, start, end, ref, alt, variant id
def write_unique_variant_bed(rows, path):
    unique = {}
    for row in rows:
        if row["_ELIGIBLE"] and row["_VARIANT_ID"] not in unique:
            unique[row["_VARIANT_ID"]] = [standardize_chromosome(row["CHROM"]), int(row["START"]), int(row["END"]),
                                          row["REF"].strip().upper(), row["ALT"].strip().upper(), row["_VARIANT_ID"]]
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for record in sorted(unique.values(), key=lambda r: (chromosome_order(r[0]), r[1], r[2], r[3], r[4])):
            writer.writerow(record)
    return unique


# Per SNV: distance to the closest gnomAD SNV, position overlap, exact match and their gnomAD AF
def annotate(closest_file):
    matches = defaultdict(list)
    with open(closest_file) as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if len(fields) != 13:
                raise RuntimeError(f"Unexpected bedtools output: expected 13 columns, found {len(fields)}")
            try:
                distance = int(fields[12])
            except ValueError:
                distance = None
            overlap = distance == 0
            matches[fields[5]].append({
                "distance": distance, "position_overlap": overlap,
                "exact_match": overlap and fields[0] == fields[6] and fields[3].upper() == fields[9].upper()
                               and fields[4].upper() == fields[10].upper(),
                "gnomad_ref": fields[9].upper(), "gnomad_alt": fields[10].upper(), "gnomad_af": parse_af(fields[11]),
            })

    annotations = {}
    for variant_id, records in matches.items():
        distances = [r["distance"] for r in records if r["distance"] is not None]
        position = [r for r in records if r["position_overlap"]]
        exact = [r for r in records if r["exact_match"]]

        position_variants = []
        for r in position:
            text = f"{r['gnomad_ref']}>{r['gnomad_alt']}:{format_af(r['gnomad_af'])}"
            if text not in position_variants:
                position_variants.append(text)
        exact_afs = [r["gnomad_af"] for r in exact if r["gnomad_af"] is not None]
        af_exact = max(exact_afs) if exact_afs else None

        annotations[variant_id] = {
            "DISTANCE_TO_CLOSEST_GNOMAD": min(distances) if distances else "NA",
            "GNOMAD_POSITION_OVERLAP": bool_text(position),
            "GNOMAD_POSITION_VARIANTS": ";".join(position_variants) if position_variants else "NA",
            "COMMON_GERMLINE": bool_text(any(r["gnomad_af"] is not None and r["gnomad_af"] >= COMMON_AF_THRESHOLD
                                             for r in position)),
            "GNOMAD_EXACT_MATCH": bool_text(exact),
            "GNOMAD_AF_EXACT": format_af(af_exact),
            "COMMON_GERMLINE_EXACT": bool_text(af_exact is not None and af_exact >= COMMON_AF_THRESHOLD),
        }
    return annotations


def write_annotated_table(rows, fieldnames, annotations, path):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames + NEW_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            out = {c: row.get(c, "") for c in fieldnames}
            if row["_ELIGIBLE"] and row["_VARIANT_ID"] in annotations:
                out.update(annotations[row["_VARIANT_ID"]])
            else:
                out.update({c: "NA" for c in NEW_COLUMNS})
            writer.writerow(out)


def write_qc(rows, unique, annotations, path):
    eligible = [annotations.get(r["_VARIANT_ID"]) for r in rows if r["_ELIGIBLE"]]
    eligible_annotated = [a for a in eligible if a is not None]

    def count(column):
        return sum(1 for a in eligible_annotated if a[column] == "TRUE")

    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["METRIC", "VALUE"])
        writer.writerows([
            ["TOTAL_MEGATABLE_ROWS", len(rows)],
            ["ELIGIBLE_SNV_ROWS", len(eligible)],
            ["NON_ELIGIBLE_ROWS", len(rows) - len(eligible)],
            ["UNIQUE_VARIANTS_TESTED", len(unique)],
            ["POSITION_OVERLAP_ROWS", count("GNOMAD_POSITION_OVERLAP")],
            ["EXACT_MATCH_ROWS", count("GNOMAD_EXACT_MATCH")],
            ["COMMON_POSITION_OVERLAP_ROWS", count("COMMON_GERMLINE")],
            ["COMMON_EXACT_MATCH_ROWS", count("COMMON_GERMLINE_EXACT")],
        ])


def main():
    parser = argparse.ArgumentParser(description="Annotate the somatic SNVs with gnomAD.")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    tools = config["tools"]
    variant_table = os.path.join(config["data_dir"], "MegaTables", "Megatable_all_alignments_all_callers.txt")
    out_dir = os.path.join(config["data_dir"], "False_Positive_Analysis", "Germline_Contamination")
    os.makedirs(out_dir, exist_ok=True)
    gnomad_bed = os.path.join(out_dir, "af-only-gnomad.hg38.biallelic_snvs.bed")
    unique_bed = os.path.join(out_dir, "unique_megatable_variants.bed")
    closest_file = os.path.join(out_dir, "unique_megatable_variants_gnomad_closest.tsv")
    annotated_file = os.path.join(out_dir, "Megatable_all_alignments_all_callers_germline_annotated.tsv")
    qc_file = os.path.join(out_dir, "germline_annotation_QC.tsv")

    if os.path.isfile(gnomad_bed) and os.path.getsize(gnomad_bed) > 0:
        print(f"Using existing gnomAD BED: {gnomad_bed}")
    else:
        run_command([tools["bcftools"], "query", "-f", "%CHROM\\t%POS0\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\n",
                     config["resources"]["gnomad_grch38"]], output_file=gnomad_bed)

    rows, fieldnames = read_variant_table(variant_table)
    unique = write_unique_variant_bed(rows, unique_bed)
    run_command([tools["bedtools"], "closest", "-sorted", "-d", "-t", "all", "-a", unique_bed, "-b", gnomad_bed],
                output_file=closest_file)
    annotations = annotate(closest_file)
    write_annotated_table(rows, fieldnames, annotations, annotated_file)
    write_qc(rows, unique, annotations, qc_file)
    print(f"Annotated variant table: {annotated_file}\nQC: {qc_file}")


if __name__ == "__main__":
    main()
