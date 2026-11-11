#!/usr/bin/env python3

# Figure 3f: Mutect2 TLOD and FILTER of every Mutect2 TP and FN site (Mutect2 counterpart of evs_measure_strelka.py).
# Uses <tumour>_Mutect2_variants.vcf.gz, the FilterMutectCalls output that keeps every candidate with its FILTER tags.
# The TP / FN sites come from VAF_results.tsv (measure_vaf.py). VG pangenome VCFs are in T2T-CHM13 coordinates, so
# those sites are lifted GRCh38 -> T2T first; a site needs the same position, REF and ALT in the VCF.
#
# Input:  <data_dir>/False_Negatives_Analysis/VAF_results.tsv
#         <calls_dir>/Mutect2/<alignment>/<donor>/<tumour>/<tumour>_Mutect2_variants.vcf.gz
# Output: <data_dir>/False_Negatives_Analysis/Mutect2_TLOD/mutect2_tlod_filter_results.tsv: one row per site, alignment
#         and status with FILTER (or NO_RECORD / NO_VCF / LIFTOVER_FAILED), TLOD, GERMQ, NALOD, NLOD, tumour AF and DP
#
# Usage:  python3 3_analyses/false_negatives/tlod_measure_mutect2.py

import argparse
import os
import subprocess

import yaml

# Alignment name in VAF_results.tsv -> calls folder
ALIGNMENTS = {"VG_T2T_pangenome": "T2T_pangenome", "VG_GRCh38_linear": "GRCh38_linear", "BWA_GRCh38_linear": "GRCh38_bwa"}


def short_id(sample_id):
    return "-".join(sample_id.split("-")[:4])


def read_tumours(path):
    """donor -> full tumour barcode"""
    tumours = {}
    with open(path) as f:
        for line in f:
            fields = line.split()
            if len(fields) >= 4:
                tumours[fields[0]] = fields[1]
    return tumours


# One row per (variant, donor, status, alignment) for one caller; each site is in VAF_results.tsv once per BAM
def load_tp_fn_sites(vaf_file, tool):
    seen, sites = set(), []
    with open(vaf_file) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if cols[idx["TOOL"]] != tool or cols[idx["Defined_By"]] not in ALIGNMENTS:
                continue
            key = (cols[idx["Variant_ID"]], cols[idx["DONOR_ID"]], cols[idx["STATUS"]], cols[idx["Defined_By"]])
            if key in seen:
                continue
            seen.add(key)
            sites.append({"variant_id": key[0], "donor": key[1], "status": key[2], "alignment": key[3],
                          "chrom": cols[idx["CHROM"]], "pos": int(cols[idx["END"]]),
                          "ref": cols[idx["REF"]], "alt": cols[idx["ALT"]]})
    return sites


# Adds query_chrom / query_pos: lifted to T2T for VG pangenome (None if it does not lift), unchanged otherwise
def liftover_sites(sites, out_dir, suffix, liftover, chain):
    to_lift = [s for s in sites if s["alignment"] == "VG_T2T_pangenome"]
    in_bed = os.path.join(out_dir, f"grch38_positions_{suffix}.bed")
    out_bed = os.path.join(out_dir, f"chm13_positions_{suffix}.bed")
    unmapped_bed = os.path.join(out_dir, f"unmapped_positions_{suffix}.bed")
    with open(in_bed, "w") as f:
        for i, s in enumerate(to_lift):
            f.write(f"{s['chrom']}\t{s['pos'] - 1}\t{s['pos']}\tVAR_{i}\n")
    subprocess.run([liftover, in_bed, chain, out_bed, unmapped_bed], check=True)

    lifted = {}
    with open(out_bed) as f:
        for line in f:
            chrom, _, end, name = line.strip().split("\t")
            lifted[name] = (chrom, int(end))
    for i, s in enumerate(to_lift):
        s["query_chrom"], s["query_pos"] = lifted.get(f"VAR_{i}", (None, None))
    for s in sites:
        if s["alignment"] != "VG_T2T_pangenome":
            s["query_chrom"], s["query_pos"] = s["chrom"], s["pos"]


# Mutect2 sample columns are TCGA barcodes, full or short depending on the donor:
# exact full barcode, then short ID inside the name, then the first of two sample columns
def find_tumour_column(header_fields, tumour_full_id):
    samples = header_fields[9:]
    for i, sample in enumerate(samples):
        if sample == tumour_full_id:
            return 9 + i
    for i, sample in enumerate(samples):
        if short_id(tumour_full_id) in sample:
            return 9 + i
    if len(samples) == 2:
        return 9
    return None


def parse_info(info_str):
    info = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            info[k] = v
        else:
            info[field] = True
    return info


def query_mutect2(tabix, header_cache, vcf_path, tumour_full_id, site):
    prefix = "CHM13#0#chr" if site["alignment"] == "VG_T2T_pangenome" else "chr"
    region = f"{prefix}{site['query_chrom']}:{site['query_pos']}-{site['query_pos']}"

    if vcf_path not in header_cache:
        header = subprocess.run([tabix, "-H", vcf_path], capture_output=True, text=True, check=False).stdout
        header_cache[vcf_path] = next((find_tumour_column(line.split("\t"), tumour_full_id)
                                       for line in header.splitlines() if line.startswith("#CHROM")), None)
    tumour_col = header_cache[vcf_path]
    if tumour_col is None:
        return None

    result = subprocess.run([tabix, vcf_path, region], capture_output=True, text=True, check=False)
    if result.returncode != 0:
        return None
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if int(fields[1]) != site["query_pos"] or fields[3] != site["ref"] or site["alt"] not in fields[4].split(","):
            continue
        info = parse_info(fields[7])
        tumour_af, tumour_dp = "", ""
        if len(fields) > tumour_col:
            tags, values = fields[8].split(":"), fields[tumour_col].split(":")
            if len(tags) == len(values):
                tag_map = dict(zip(tags, values))
                tumour_af, tumour_dp = tag_map.get("AF", ""), tag_map.get("DP", "")
        return {"filter": fields[6] if fields[6] not in (".", "") else "PASS", "tlod": info.get("TLOD", ""),
                "germq": info.get("GERMQ", ""), "nalod": info.get("NALOD", ""), "nlod": info.get("NLOD", ""),
                "tumour_af": tumour_af, "tumour_dp": tumour_dp}
    return None


def main():
    parser = argparse.ArgumentParser(description="Mutect2 TLOD and FILTER at the TP and FN sites.")
    parser.add_argument("--patients", default="metadata/BLCA_Patients.txt")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    out_dir = os.path.join(config["data_dir"], "False_Negatives_Analysis", "Mutect2_TLOD")
    os.makedirs(out_dir, exist_ok=True)

    tumours = read_tumours(args.patients)
    sites = load_tp_fn_sites(os.path.join(config["data_dir"], "False_Negatives_Analysis", "VAF_results.tsv"), "Mutect2")
    print(f"Loaded {len(tumours)} samples and {len(sites)} Mutect2 sites")
    liftover_sites(sites, out_dir, "tlod", config["tools"]["liftover"], config["resources"]["chain_grch38_to_t2t_bare"])

    header_cache, n_record, n_no_record, n_skipped = {}, 0, 0, 0
    output_path = os.path.join(out_dir, "mutect2_tlod_filter_results.tsv")
    with open(output_path, "w") as out:
        out.write("variant_id\tdonor\talignment\tstatus\tfilter\ttlod\tgermq\tnalod\tnlod\ttumour_af\ttumour_dp\n")
        for s in sites:
            label = f"{s['variant_id']}\t{s['donor']}\t{s['alignment']}\t{s['status']}"
            if s["query_pos"] is None:
                n_skipped += 1
                out.write(f"{label}\tLIFTOVER_FAILED\t\t\t\t\t\n")
                continue
            tumour_full = tumours.get(s["donor"])
            if tumour_full is None:
                n_skipped += 1
                continue

            tumour = short_id(tumour_full)
            vcf_path = os.path.join(config["calls_dir"], "Mutect2", ALIGNMENTS[s["alignment"]], s["donor"], tumour,
                                    f"{tumour}_Mutect2_variants.vcf.gz")
            if not os.path.exists(vcf_path):
                n_skipped += 1
                out.write(f"{label}\tNO_VCF\t\t\t\t\t\n")
                continue

            rec = query_mutect2(config["tools"]["tabix"], header_cache, vcf_path, tumour_full, s)
            if rec is None:
                n_no_record += 1
                out.write(f"{label}\tNO_RECORD\t\t\t\t\t\n")
                continue
            n_record += 1
            out.write(f"{label}\t{rec['filter']}\t{rec['tlod']}\t{rec['germq']}\t{rec['nalod']}\t{rec['nlod']}\t"
                      f"{rec['tumour_af']}\t{rec['tumour_dp']}\n")

    print(f"Sites with a raw Mutect2 record: {n_record}")
    print(f"Sites with no raw Mutect2 record: {n_no_record}")
    print(f"Sites skipped (liftover failure / missing sample / missing VCF): {n_skipped}")
    print(f"Output: {output_path}")


if __name__ == "__main__":
    main()
