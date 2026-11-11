#!/usr/bin/env python3

# Supplementary Table 1: one row per donor and alignment, with a BLCA and a LUAD tab (LUAD: Strelka2 and Mutect2 only).
#
# Input:  <data_dir>/[LUAD/]Stats_Summary/ (compare_to_mc3.py, build_consensus.py; ABSOLUTE purity and PanCancer Atlas
#             mutation burden in Strelka/)
#         <data_dir>/False_Positive_Analysis/Germline_Contamination/ (germline_contamination.py)
#         <data_dir>/Reference_Bias/paired_sites/ (combine_paired_counts.py)
#         metadata/admixture/, GATK alignment summary tables in <alignments_dir>
# Output: <data_dir>/Supplementary_Tables/Supplementary_Table_1/Supplementary_Table_1.xlsx
#         Prints the number of donors filled per column and alignment
#
# Usage:  python3 2_benchmarking/supplementary_tables/build_supplementary_table_1.py

import argparse
import csv
import glob
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path

import yaml

ALIGNMENT_LABELS = {"GRCh38_bwa": "BWA_GRCh38", "GRCh38_linear": "VG_GRCh38", "T2T_pangenome": "VG_Pangenome"}
STATS_FILES = {"GRCh38_bwa": "grch38_bwa", "GRCh38_linear": "grch38_linear", "T2T_pangenome": "pangenome"}

# Column prefix -> Stats_Summary caller folder (None = consensus)
CALLERS = {"Consensus": None, "Mutect": "Mutect2", "SomaticSniper": "SomaticSniper", "Strelka": "Strelka",
           "DeepSomatic": "DeepSomatic"}

# Alignment names used in the paired-site reference bias tables
PAIRED_ALIGNMENTS = {"GRCh38_bwa": "bwa_grch38", "GRCh38_linear": "vg_grch38", "T2T_pangenome": "vg_pangenome"}

# GATK alignment summary tables (each holds both cohorts)
QC_FILES = {
    "GRCh38_bwa": "2026-10-09_GRCh38bwa_gatk_alignment_statistics.tsv",
    "GRCh38_linear": "2026-09-16_GRCh38_linear_gatk_alignment_statistics.tsv",
    "T2T_pangenome": "2026-10-08_T2Thprc1.1_gatk_alignment_statistics.tsv",
}

# The VG GRCh38 alignment of this normal is labelled 10A; it is the same 10B sample
SAMPLE_OVERRIDE = {"TCGA-E5-A4TZ-10B": "TCGA-E5-A4TZ-10A"}

COHORTS = {
    "BLCA": {"patients": "metadata/BLCA_Patients.txt", "subdir": "",
             "admixture": "metadata/admixture/2026_10_06_Carrot_et_al_BLCA_ancestry_admixture_proportion.xlsx",
             "yuan": "metadata/admixture/2026_10_06_BLCA_Yuan_et_al_ancestry_proportion.xlsx"},
    "LUAD": {"patients": "metadata/LUAD_Patients.txt", "subdir": "LUAD",
             "admixture": "metadata/admixture/2026_06_02_Carrot_et_al_LUAD_ancestry_admixture.xlsx", "yuan": None},
}

LUAD_DROP = ("Yuan_", "Consensus_", "SomaticSniper_", "DeepSomatic_", "Mutect_Germline", "Strelka_Germline", "Median_REF")

COLUMNS = ["Donor_ID", "Reference", "Tumour_ID", "Blood_ID", "Ancestry"]
for prefix in CALLERS:
    COLUMNS += [f"{prefix}_Recall", f"{prefix}_Precision", f"{prefix}_F1_Score"]
for prefix in ["Mutect", "SomaticSniper", "Strelka", "DeepSomatic"]:
    COLUMNS += [f"{prefix}_Germline_Somatic_Overlap", f"{prefix}_Germline_Somatic_Overlap_MAF>=0.01"]
COLUMNS += ["Median_REF_Frequency_Germline", "Median_REF_Frequency_Somatics", "Tumour_Purity", "Somatic_Mutation_Burden",
            "AA_Admixture_Proportion", "EA_Admixture_Proportion", "EAA_Admixture_Proportion",
            "Yuan_AA_Admixture_Proportion", "Yuan_EA_Admixture_Proportion", "Yuan_EAA_Admixture_Proportion",
            "Tumour_Total_Reads", "Normal_Total_Reads", "Tumour_Properly_Paired", "Normal_Properly_Paired",
            "Tumour_Unaligned_Reads", "Normal_Unaligned_Reads"]


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_xlsx(path, rscript):
    text = subprocess.run([rscript, "-e", f"write.table(readxl::read_xlsx('{path}'), stdout(), sep = '\\t', "
                           "row.names = FALSE, quote = FALSE)"], capture_output=True, text=True, check=True).stdout
    return {r["Patient"]: r for r in csv.DictReader(text.splitlines(), delimiter="\t")}


def short_id(sample_id):
    """TCGA-FD-A3B5-01A-11D-A20D-08 -> TCGA-FD-A3B5-01A"""
    return "-".join(sample_id.split("-")[:4])


# Full TCGA barcode from the BWA BAM name (the LUAD patient list only has short IDs)
def full_barcode(alignments_dir, donor, sample_id):
    bams = glob.glob(str(alignments_dir / "GRCh38_bwa" / donor / short_id(sample_id) / f"{short_id(sample_id)}-*.bam"))
    return Path(bams[0]).stem if len(bams) == 1 else sample_id


# Median REF frequency per donor and alignment at the paired sites that pass the depth filter in all three alignments
# and have the same reference base in GRCh38 and T2T (as in Figure 3d)
def median_ref_frequency(path, somatic):
    sites = defaultdict(lambda: defaultdict(list))
    for r in read_tsv(path):
        if somatic:
            covered = all(int(r[f"{a}_tumour_depth"]) >= 5 and int(r[f"{a}_normal_depth"]) > 8 for a in PAIRED_ALIGNMENTS.values())
        else:
            covered = all(int(r[f"{a}_normal_depth"]) > 8 for a in PAIRED_ALIGNMENTS.values())
        if covered and r["grch38_ref"] == r["t2t_ref"]:
            for alignment, a in PAIRED_ALIGNMENTS.items():
                sites[r["donor"]][alignment].append(float(r[f"{a}_ref_freq"]))
    return {(d, a): round(statistics.median(v), 4) for d, by in sites.items() for a, v in by.items()}


def build_cohort(name, cohort, config):
    data_dir = Path(config["data_dir"])
    alignments_dir = Path(config["alignments_dir"])
    stats_dir = data_dir / cohort["subdir"] / "Stats_Summary"
    rscript = config["tools"]["rscript"]

    patients = []
    with open(cohort["patients"]) as f:
        for line in f:
            fields = line.split()
            if len(fields) >= 4:
                patients.append({"donor": fields[0], "tumour": fields[1], "normal": fields[2], "ancestry": fields[3]})
    donors = {p["donor"] for p in patients}
    values = defaultdict(dict)

    # Recall / precision / F1 vs MC3, per caller and alignment
    for prefix, caller in CALLERS.items():
        for alignment, stem in STATS_FILES.items():
            if caller:
                path = stats_dir / caller / f"stats_table_{stem}.txt"
            else:
                path = stats_dir / "MC3_Consensus" / f"stats_table_{alignment}.txt"
            if not path.exists():
                continue
            for r in read_tsv(path):
                if r["Donor_ID"] in donors and (caller or r["Variant_Caller"] == "Consensus_ge3"):
                    row = values[(r["Donor_ID"], alignment)]
                    row[f"{prefix}_Recall"], row[f"{prefix}_Precision"], row[f"{prefix}_F1_Score"] = (
                        r["Recall(%)"], r["Precision(%)"], r["F1_Score(%)"])

    # Germline-somatic overlap and REF frequency (BLCA only)
    if name == "BLCA":
        prefixes = {"Mutect2": "Mutect", "SomaticSniper": "SomaticSniper", "Strelka": "Strelka", "DeepSomatic": "DeepSomatic"}
        # Every caller ran on every alignment, so a donor with no somatic calls has 0 overlaps (not NA)
        overlap = {(d, a, prefix): [0, 0] for d in donors for a in ALIGNMENT_LABELS for prefix in prefixes.values()}
        for r in read_tsv(data_dir / "False_Positive_Analysis/Germline_Contamination/"
                          "Megatable_all_alignments_all_callers_germline_annotated.tsv"):
            if (r["STATUS"].strip().upper() in ("TP", "FP") and r["DONOR_ID"] in donors and r["TOOL"] in prefixes
                    and r["ALIGNMENT"] in ALIGNMENT_LABELS):
                counts = overlap[(r["DONOR_ID"], r["ALIGNMENT"], prefixes[r["TOOL"]])]
                counts[0] += r["GNOMAD_EXACT_MATCH"].strip().upper() == "TRUE"
                counts[1] += r["COMMON_GERMLINE_EXACT"].strip().upper() == "TRUE"
        for (donor, alignment, prefix), (n_all, n_common) in overlap.items():
            values[(donor, alignment)][f"{prefix}_Germline_Somatic_Overlap"] = n_all
            values[(donor, alignment)][f"{prefix}_Germline_Somatic_Overlap_MAF>=0.01"] = n_common

        paired_dir = data_dir / "Reference_Bias" / "paired_sites"
        for column, file, somatic in [("Median_REF_Frequency_Germline", "germline_paired_counts.tsv", False),
                                      ("Median_REF_Frequency_Somatics", "somatic_paired_counts.tsv", True)]:
            for (donor, alignment), median in median_ref_frequency(paired_dir / file, somatic).items():
                if donor in donors:
                    values[(donor, alignment)][column] = median

    # Donor-level values, the same on every alignment row
    purity = {r["Donor_ID"]: r["Purity"] for r in read_tsv(stats_dir / "Strelka" / "tumour_purity.txt")}
    burden = {r["Donor_ID"]: r["Somatic"] for r in read_tsv(stats_dir / "Strelka" / "mutation_burden.txt")}
    admixture = read_xlsx(cohort["admixture"], rscript)
    yuan = read_xlsx(cohort["yuan"], rscript) if cohort["yuan"] else {}

    # Reads per tumour / normal BAM
    for alignment, gatk_file in QC_FILES.items():
        gatk = {r["SAMPLE"][:16]: r for r in read_tsv(alignments_dir / alignment / gatk_file) if r["CATEGORY"] == "PAIR"}
        for p in patients:
            row = values[(p["donor"], alignment)]
            for kind, label in [("tumour", "Tumour"), ("normal", "Normal")]:
                sample = short_id(p[kind])
                g = gatk.get(sample) or gatk.get(SAMPLE_OVERRIDE.get(sample, sample))
                if g:
                    row[f"{label}_Total_Reads"] = g["TOTAL_READS"]
                    row[f"{label}_Properly_Paired"] = g["PROPER_PAIRS"]
                    row[f"{label}_Unaligned_Reads"] = g["UNALIGNED_READS"]

    # One row per donor and alignment; alignments with nothing for this cohort are left out
    rows = []
    for alignment in ALIGNMENT_LABELS:
        if not any(values.get((p["donor"], alignment)) for p in patients):
            continue
        for p in patients:
            row = values.get((p["donor"], alignment), {})
            a, y = admixture.get(p["donor"], {}), yuan.get(p["donor"], {})
            row.update({
                "Donor_ID": p["donor"], "Reference": ALIGNMENT_LABELS[alignment],
                "Tumour_ID": full_barcode(alignments_dir, p["donor"], p["tumour"]),
                "Blood_ID": full_barcode(alignments_dir, p["donor"], p["normal"]),
                "Ancestry": p["ancestry"], "Tumour_Purity": "" if purity.get(p["donor"]) in (None, "NA") else purity[p["donor"]],
                "Somatic_Mutation_Burden": burden.get(p["donor"], ""),
                "AA_Admixture_Proportion": a.get("AA", ""), "EA_Admixture_Proportion": a.get("EA", ""),
                "EAA_Admixture_Proportion": a.get("EAA", ""), "Yuan_AA_Admixture_Proportion": y.get("AA", ""),
                "Yuan_EA_Admixture_Proportion": y.get("EA", ""), "Yuan_EAA_Admixture_Proportion": y.get("EAA", ""),
            })
            rows.append([row.get(c, "") for c in COLUMNS])

    keep = [i for i, c in enumerate(COLUMNS) if not (name == "LUAD" and c.startswith(LUAD_DROP))]
    return [COLUMNS[i] for i in keep], [[r[i] for i in keep] for r in rows]


def main():
    parser = argparse.ArgumentParser(description="Build Supplementary Table 1.")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    out_dir = Path(config["data_dir"]) / "Supplementary_Tables" / "Supplementary_Table_1"
    out_dir.mkdir(parents=True, exist_ok=True)

    tsvs = []
    for name, cohort in COHORTS.items():
        header, rows = build_cohort(name, cohort, config)
        tsv = out_dir / f".{name}_tab.tsv"
        with open(tsv, "w") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerows(rows)
        tsvs.append((name, tsv))

        # Donors filled per column and alignment
        print(f"== {name}: {len(rows)} rows")
        references = sorted({r[1] for r in rows}, key=list(ALIGNMENT_LABELS.values()).index)
        for i, column in enumerate(header[5:], start=5):
            filled = [f"{ref} {sum(1 for r in rows if r[1] == ref and r[i] != '')}" for ref in references]
            print(f"{column}: " + ", ".join(filled))

    # Write the xlsx through R (missing values as NA)
    xlsx = out_dir / "Supplementary_Table_1.xlsx"
    subprocess.run([config["tools"]["rscript"], "-e",
                    "read <- function(x) read.delim(x, check.names = FALSE, stringsAsFactors = FALSE, na.strings = ''); "
                    "openxlsx::write.xlsx(list(" + ", ".join(f"'{n}' = read('{t}')" for n, t in tsvs) +
                    f"), '{xlsx}', keepNA = TRUE, na.string = 'NA', overwrite = TRUE)"], check=True)
    for _, tsv in tsvs:
        tsv.unlink()
    print(f"-> {xlsx}")


if __name__ == "__main__":
    main()
