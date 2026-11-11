#!/usr/bin/env python3

# Figure 3e: Strelka2 SomaticEVS and FILTER of every Strelka2 TP and FN site, from Strelka2's raw output
# (results/variants/somatic.snvs.vcf.gz, all candidates, PASS or not).
# The TP / FN sites come from VAF_results.tsv (measure_vaf.py). VG pangenome VCFs are in T2T-CHM13 coordinates, so
# those sites are lifted GRCh38 -> T2T first; a site needs the same position, REF and ALT in the VCF.
#
# Input:  <data_dir>/False_Negatives_Analysis/VAF_results.tsv
#         <calls_dir>/Strelka/<alignment>/<donor>/<tumour>/results/variants/somatic.snvs.vcf.gz
# Output: <data_dir>/False_Negatives_Analysis/Strelka_EVS/strelka_evs_filter_results.tsv: one row per site, alignment
#         and status with FILTER (or NO_RECORD / LIFTOVER_FAILED), SomaticEVS, MQ, SNVSB, QSS_NT, ReadPosRankSum,
#         tumour DP and FDP
#
# Usage:  python3 3_analyses/false_negatives/evs_measure_strelka.py   (needs pysam)

import argparse
import os
import subprocess

import pysam
import yaml

# Alignment name in VAF_results.tsv -> calls folder
ALIGNMENTS = {"VG_T2T_pangenome": "T2T_pangenome", "VG_GRCh38_linear": "GRCh38_linear", "BWA_GRCh38_linear": "GRCh38_bwa"}


def read_tumours(path):
    """donor -> short tumour ID"""
    tumours = {}
    with open(path) as f:
        for line in f:
            fields = line.split()
            if len(fields) >= 4:
                tumours[fields[0]] = "-".join(fields[1].split("-")[:4])
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


def query_strelka(vcf_cache, vcf_path, site):
    if vcf_path not in vcf_cache:
        vcf_cache[vcf_path] = pysam.VariantFile(vcf_path) if os.path.exists(vcf_path) else None
    vcf = vcf_cache[vcf_path]
    if vcf is None:
        return None
    prefix = "CHM13#0#chr" if site["alignment"] == "VG_T2T_pangenome" else "chr"
    try:
        for rec in vcf.fetch(f"{prefix}{site['query_chrom']}", site["query_pos"] - 1, site["query_pos"]):
            if rec.pos == site["query_pos"] and rec.ref == site["ref"] and site["alt"] in rec.alts:
                return rec
    except ValueError:
        return None
    return None


def main():
    parser = argparse.ArgumentParser(description="Strelka2 SomaticEVS and FILTER at the TP and FN sites.")
    parser.add_argument("--patients", default="metadata/BLCA_Patients.txt")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    out_dir = os.path.join(config["data_dir"], "False_Negatives_Analysis", "Strelka_EVS")
    os.makedirs(out_dir, exist_ok=True)

    tumours = read_tumours(args.patients)
    sites = load_tp_fn_sites(os.path.join(config["data_dir"], "False_Negatives_Analysis", "VAF_results.tsv"), "Strelka")
    print(f"Loaded {len(tumours)} samples and {len(sites)} Strelka2 sites")
    liftover_sites(sites, out_dir, "evs", config["tools"]["liftover"], config["resources"]["chain_grch38_to_t2t_bare"])

    vcf_cache, n_record, n_no_record, n_skipped = {}, 0, 0, 0
    output_path = os.path.join(out_dir, "strelka_evs_filter_results.tsv")
    with open(output_path, "w") as out:
        out.write("variant_id\tdonor\talignment\tstatus\tfilter\tsomatic_evs\tmq\tsnvsb\tqss_nt\treadposranksum"
                  "\ttumor_dp\ttumor_fdp\n")
        for s in sites:
            label = f"{s['variant_id']}\t{s['donor']}\t{s['alignment']}\t{s['status']}"
            if s["query_pos"] is None:
                n_skipped += 1
                out.write(f"{label}\tLIFTOVER_FAILED\t\t\t\t\t\t\n")
                continue
            tumour = tumours.get(s["donor"])
            if tumour is None:
                n_skipped += 1
                continue

            vcf_path = os.path.join(config["calls_dir"], "Strelka", ALIGNMENTS[s["alignment"]], s["donor"], tumour,
                                    "results", "variants", "somatic.snvs.vcf.gz")
            rec = query_strelka(vcf_cache, vcf_path, s)
            if rec is None:
                n_no_record += 1
                out.write(f"{label}\tNO_RECORD\t\t\t\t\t\t\n")
                continue

            n_record += 1
            filt = ";".join(rec.filter.keys()) if rec.filter.keys() else "PASS"
            tumour_sample = rec.samples.get("TUMOR")
            tumour_dp = tumour_sample.get("DP", "") if tumour_sample is not None else ""
            tumour_fdp = tumour_sample.get("FDP", "") if tumour_sample is not None else ""
            info = rec.info
            out.write(f"{label}\t{filt}\t{info.get('SomaticEVS', '')}\t{info.get('MQ', '')}\t{info.get('SNVSB', '')}\t"
                      f"{info.get('QSS_NT', '')}\t{info.get('ReadPosRankSum', '')}\t{tumour_dp}\t{tumour_fdp}\n")

    print(f"Sites with a raw Strelka2 record: {n_record}")
    print(f"Sites with no raw Strelka2 record: {n_no_record}")
    print(f"Sites skipped (liftover failure / missing sample): {n_skipped}")
    print(f"Output: {output_path}")


if __name__ == "__main__":
    main()
