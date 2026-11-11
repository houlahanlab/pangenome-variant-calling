#!/usr/bin/env python3

# Reference bias, step 5: combines the per-BAM read counts into one row per site with the three alignments side by side.
# 
# REF is the reference base of the genome each alignment used (grch38_ref for BWA GRCh38 and VG GRCh38, t2t_ref for
# VG pangenome), so ref_freq = ref_count / (ref_count + alt_count) is the share of reads matching that alignment's
# reference (NA with no reads). Depth = all reads at the site, to mirror the depth filter used for calling
# (normal depth > 8, tumour depth >= 5).
#
# Input:  <data_dir>/Reference_Bias/paired_sites[_mutect2]/ manifest.tsv, somatic_sites.tsv, germline_sites.tsv, counts/
# Output: the same folder: somatic_paired_counts.tsv, germline_paired_counts.tsv,
#         counts_QC.tsv (sites passing the depth filter per donor and alignment)
#
# Usage:  python3 3_analyses/false_positives/reference_bias/paired_sites/combine_paired_counts.py --caller strelka

import argparse
import os

import yaml

ALIGNMENTS = ["bwa_grch38", "vg_grch38", "vg_pangenome"]
SAMPLE_TYPE = {"somatic": "tumour", "germline": "normal"}
MIN_TUMOUR_DEPTH = 5   # tumour depth >= 5
MIN_NORMAL_DEPTH = 9   # normal depth > 8


def read_table(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f if line.strip()]


def read_pileup_counts(path):
    """{(chrom, pos): {base: reads}} from the bcftools query output."""
    counts = {}
    for row in read_table(path):
        alleles = [row["REF"].upper()] + row["ALT"].upper().split(",")
        depths = row["AD"].split(",")
        base_counts = {}
        for allele, depth in zip(alleles, depths):
            if len(allele) == 1 and depth != ".":
                base_counts[allele] = base_counts.get(allele, 0) + int(depth)
        counts[(row["CHROM"], int(row["POS"]))] = base_counts
    return counts


def main():
    parser = argparse.ArgumentParser(description="Combine the paired-site read counts.")
    parser.add_argument("--caller", choices=["strelka", "mutect2"], default="strelka")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    root = os.path.join(config["data_dir"], "Reference_Bias", "paired_sites" if args.caller == "strelka" else "paired_sites_mutect2")

    # Read counts per (donor, alignment, sample type)
    pileups = {}
    for row in read_table(os.path.join(root, "manifest.tsv")):
        if not os.path.exists(row["out_file"]):
            raise FileNotFoundError(f"Missing counts for {row['donor']}/{row['alignment']}/{row['sample_type']}: "
                                    f"{row['out_file']} -- did every Slurm task finish?")
        pileups[(row["donor"], row["alignment"], row["sample_type"])] = read_pileup_counts(row["out_file"])

    qc = {}
    for kind in ["somatic", "germline"]:
        sites = read_table(os.path.join(root, f"{kind}_sites.tsv"))
        out_rows = []
        for site in sites:
            out = dict(site)
            for alignment in ALIGNMENTS:
                if alignment == "vg_pangenome":
                    key, ref_base = (site["t2t_chrom"], int(site["t2t_pos"])), site["t2t_ref"]
                else:
                    key, ref_base = (site["grch38_chrom"], int(site["grch38_pos"])), site["grch38_ref"]
                alt_base = site["allele_2"] if ref_base == site["allele_1"] else site["allele_1"]

                base_counts = pileups[(site["donor"], alignment, SAMPLE_TYPE[kind])].get(key, {})
                ref_count = base_counts.get(ref_base, 0)
                alt_count = base_counts.get(alt_base, 0)
                depth = ref_count + alt_count
                out[f"{alignment}_ref_count"] = ref_count
                out[f"{alignment}_alt_count"] = alt_count
                out[f"{alignment}_ref_freq"] = ref_count / depth if depth > 0 else "NA"

                if kind == "somatic":
                    tumour_depth = sum(base_counts.values())
                    normal_depth = sum(pileups[(site["donor"], alignment, "normal_somatic")].get(key, {}).values())
                    out[f"{alignment}_tumour_depth"] = tumour_depth
                    out[f"{alignment}_normal_depth"] = normal_depth
                    passes = tumour_depth >= MIN_TUMOUR_DEPTH and normal_depth >= MIN_NORMAL_DEPTH
                else:
                    normal_depth = sum(base_counts.values())
                    out[f"{alignment}_normal_depth"] = normal_depth
                    passes = normal_depth >= MIN_NORMAL_DEPTH

                qc_key = (site["donor"], site["ancestry"], kind, alignment)
                qc.setdefault(qc_key, [0, 0])
                qc[qc_key][0] += 1
                qc[qc_key][1] += passes
            out_rows.append(out)

        depth_columns = ("tumour_depth", "normal_depth") if kind == "somatic" else ("normal_depth",)
        columns = list(sites[0].keys()) + [f"{a}_{x}" for a in ALIGNMENTS
                                           for x in ("ref_count", "alt_count", "ref_freq") + depth_columns]
        out_path = os.path.join(root, f"{kind}_paired_counts.tsv")
        with open(out_path, "w") as f:
            f.write("\t".join(columns) + "\n")
            for r in out_rows:
                f.write("\t".join(str(r[c]) for c in columns) + "\n")
        print(f"Wrote {out_path} ({len(out_rows)} sites)")

    qc_path = os.path.join(root, "counts_QC.tsv")
    with open(qc_path, "w") as f:
        f.write("Donor\tAncestry\tSite_Type\tAlignment\tSites\tSites_Passing_Depth_Filter\n")
        for (donor, ancestry, kind, alignment), (n, n_covered) in sorted(qc.items()):
            f.write(f"{donor}\t{ancestry}\t{kind}\t{alignment}\t{n}\t{n_covered}\n")
    print(f"Wrote {qc_path}")


if __name__ == "__main__":
    main()
