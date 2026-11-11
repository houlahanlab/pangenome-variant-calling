#!/usr/bin/env python3

# Figure 4: tumour and normal REF / ALT read counts of the Strelka2 FN and FP SNVs in cancer driver genes
# (Bailey et al.; BLCA and PANCAN genes), for VG pangenome or BWA GRCh38. Uses the first SnpEff annotation of each SNV.
#   FP sites (Strelka2 call): tier 1 counts (AU/CU/GU/TU) from the donor's filtered Strelka2 VCF
#   FN sites (no call):       samtools mpileup on the tumour and normal BAMs (-aa -Q 0)
# Only the REF and ALT bases are counted, freq = count / (REF + ALT). VG pangenome is in T2T-CHM13 coordinates: FP
# sites take their T2T position from the comparison files, FN sites are lifted GRCh38 -> T2T
#
# Input:  <data_dir>/Variant_Caller_Overlap/<alignment>/Shared_Variants_Annotations/{False_Negatives,False_Positives}/
#             *_filtered_shared_by_1_caller_in_capture_annotated.vcf (impact_variants/fp_fn_variant_caller_overlap_annotations.py)
#         <filtered_calls_dir>/Strelka/<alignment>/.../*_Strelka_pass_variants_DP_target.vcf.gz (compare_to_mc3.py)
#         BAMs in <alignments_dir>; driver gene list from config.yaml
# Output: <data_dir>/False_Positive_Analysis/Allele_Freq/{tumour,normal}_{fp,fn}[_bwa_linear][_high_impact]_plot_data.txt
#         (all impacts, and HIGH + MODERATE only)
#
# Usage:  python3 3_analyses/impact_variants/allele_frequency/build_ref_alt_freq.py --alignment T2T_pangenome   (needs pysam)

import argparse
import os
import re
import subprocess

import pysam
import yaml

# Output name suffix per alignment
SUFFIX = {"T2T_pangenome": "", "GRCh38_bwa": "_bwa_linear"}


def read_samples(path):
    """donor -> ancestry, short tumour / normal IDs (folder names) and full IDs (BWA BAM names)"""
    samples = {}
    with open(path) as f:
        for line in f:
            fields = line.split()
            if len(fields) >= 4:
                donor, tumour, normal, ancestry = fields[:4]
                samples[donor] = {"ancestry": ancestry,
                                  "tumour": "-".join(tumour.split("-")[:4]), "tumour_full": tumour,
                                  "normal": "-".join(normal.split("-")[:4]), "normal_full": normal}
    return samples


# BLCA and PANCAN driver genes, read through R since the list is an xlsx
def read_driver_genes(path, rscript):
    text = subprocess.run(
        [rscript, "-e", f"genes <- readxl::read_xlsx('{path}'); "
         "write.table(genes[, c('Gene', 'Cancer')], stdout(), sep = '\\t', row.names = FALSE, quote = FALSE)"],
        capture_output=True, text=True, check=True).stdout
    return {line.split("\t")[0] for line in text.splitlines()[1:] if line.split("\t")[1] in ("BLCA", "PANCAN")}


# Strelka2 driver gene SNVs from an annotated VCF, using the first SnpEff annotation
def read_annotated_vcf(path, driver_genes, impacts=None):
    variants = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt, _, _, info = line.rstrip("\n").split("\t")[:8]
            sample = re.search(r"SAMPLE=([^;]+)", info)
            tool = re.search(r"TOOL=([^;]+)", info)
            ann = re.search(r"ANN=([^\t]+)", info)
            if len(ref) != 1 or len(alt) != 1 or not (sample and tool and ann) or tool.group(1) != "Strelka":
                continue
            ann = ann.group(1).split(",")[0].split("|")
            if len(ann) < 4 or ann[3] not in driver_genes or (impacts and ann[2] not in impacts):
                continue
            variants.append({"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt,
                             "donor": sample.group(1), "gene": ann[3], "impact": ann[2]})
    return variants


# REF and ALT rows for one sample
def freq_rows(v, ancestry, ref_count, alt_count):
    total = ref_count + alt_count
    ref_freq = round(ref_count / total, 4) if total else 0.0
    alt_freq = round(alt_count / total, 4) if total else 0.0
    gene_pos = f"{v['gene']}:{v['pos']}"
    return [(gene_pos, ancestry, v["impact"], "REF", ref_freq, ref_count),
            (gene_pos, ancestry, v["impact"], "ALT", alt_freq, alt_count)]


def write_rows(rows, path):
    with open(path, "w") as f:
        f.write("gene_pos\tancestry\timpact\tpart\tfreq\tcount\n")
        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")
    print(f"Wrote {len(rows)} rows -> {os.path.basename(path)}")


# GRCh38 -> T2T positions for all variants in one liftOver call (v["chm13"] = (chrom, pos) or None)
def lift_to_chm13(variants, out_dir, liftover, chain):
    in_bed = os.path.join(out_dir, "grch38_positions.bed")
    out_bed = os.path.join(out_dir, "chm13_positions.bed")
    with open(in_bed, "w") as f:
        for i, v in enumerate(variants):
            f.write(f"{v['chrom']}\t{v['pos'] - 1}\t{v['pos']}\tVAR_{i}\n")
    subprocess.run([liftover, in_bed, chain, out_bed, os.path.join(out_dir, "unmapped_positions.bed")], check=True)
    with open(out_bed) as f:
        lifted = {name: (chrom, int(end)) for chrom, _, end, name in (line.split() for line in f)}
    for i, v in enumerate(variants):
        v["chm13"] = lifted.get(f"VAR_{i}")
    return variants


# REF and ALT read counts at one position from the samtools mpileup read bases
def mpileup_counts(samtools, reference, bam, chrom, pos, ref, alt):
    out = subprocess.run([samtools, "mpileup", "-f", reference, "-r", f"{chrom}:{pos}-{pos}", "-aa", "-Q", "0", bam],
                         capture_output=True, text=True).stdout.split("\t")
    if len(out) < 5:
        return 0, 0
    ref_base, pileup = out[2].upper(), out[4]

    ref_count = alt_count = i = 0
    while i < len(pileup):
        c = pileup[i]
        if c == "^":            # read start, next character is its mapping quality
            i += 2
            continue
        if c in "+-":           # indel after this base, skip its length and sequence
            j = i + 1
            while pileup[j].isdigit():
                j += 1
            i = j + int(pileup[i + 1:j])
            continue
        base = ref_base if c in ".," else c.upper()
        if base == ref:
            ref_count += 1
        elif base == alt:
            alt_count += 1
        i += 1
    return ref_count, alt_count


# GRCh38 -> T2T positions of one donor's Strelka2 VG pangenome calls (strelka_patient.bed and strelka_grch38.bed are row matched)
def strelka_chm13_lookup(comparison_dir, sample):
    folder = os.path.join(comparison_dir, sample["ancestry"], sample["donor"], sample["tumour"])
    patient_bed, grch38_bed = os.path.join(folder, "strelka_patient.bed"), os.path.join(folder, "strelka_grch38.bed")
    if not (os.path.exists(patient_bed) and os.path.exists(grch38_bed)):
        return {}
    lookup = {}
    with open(patient_bed) as fp, open(grch38_bed) as fg:
        for line_p, line_g in zip(fp, fg):
            p, g = line_p.split("\t"), line_g.rstrip("\n").split("\t")
            lookup[(g[0], int(g[2]), g[3], g[4])] = (p[0].replace("CHM13#0#chr", ""), int(p[2]))
    return lookup


def fp_counts(variants, samples, args, config):
    t2t = args.alignment == "T2T_pangenome"
    strelka_dir = os.path.join(config["filtered_calls_dir"], "Strelka", args.alignment)
    comparison_dir = os.path.join(config["data_dir"], "TCGA_Comparison", "strelka_vs_tcga", args.alignment)
    tumour_rows, normal_rows, lookups, vcfs = [], [], {}, {}
    for v in variants:
        donor = v["donor"]
        if donor not in samples:
            continue

        if t2t:
            if donor not in lookups:
                lookups[donor] = strelka_chm13_lookup(comparison_dir, {**samples[donor], "donor": donor})
            mapped = lookups[donor].get((v["chrom"], v["pos"], v["ref"], v["alt"]))
            if mapped is None:
                print(f"WARNING: no CHM13 position for FP {v['gene']} {v['chrom']}:{v['pos']} ({donor})")
                continue
            chrom, pos = f"CHM13#0#chr{mapped[0]}", mapped[1]

        tumour = samples[donor]["tumour"]
        vcf = os.path.join(strelka_dir, donor, tumour, f"{tumour}_Strelka_pass_variants_DP_target.vcf.gz")
        if vcf not in vcfs:
            if t2t:
                vcfs[vcf] = pysam.VariantFile(vcf)
            else:
                vcfs[vcf] = pysam.VariantFile(vcf) if os.path.exists(vcf) else None
                if vcfs[vcf] is None:
                    print(f"WARNING: missing Strelka VCF for {donor}: {vcf}")
        if vcfs[vcf] is None:
            continue
        if not t2t:
            chrom, pos = f"chr{v['chrom']}", v["pos"]

        try:
            record = next((r for r in vcfs[vcf].fetch(chrom, pos - 1, pos)
                           if r.pos == pos and r.ref == v["ref"] and v["alt"] in r.alts), None)
        except ValueError:
            record = None
        if record is None:
            print(f"WARNING: no Strelka record for FP {v['gene']} {chrom}:{pos} ({donor})")
            continue

        for sample, rows in [("TUMOR", tumour_rows), ("NORMAL", normal_rows)]:
            counts = {base: record.samples[sample][f"{base}U"][0] for base in "ACGT"}
            rows += freq_rows(v, samples[donor]["ancestry"], counts[v["ref"]], counts[v["alt"]])
    return tumour_rows, normal_rows


def fn_counts(variants, samples, args, config, out_dir):
    t2t = args.alignment == "T2T_pangenome"
    tools, resources = config["tools"], config["resources"]
    reference = resources["fasta_t2t"] if t2t else resources["fasta_grch38"]
    alignment_dir = os.path.join(config["alignments_dir"], args.alignment)
    if t2t:
        variants = lift_to_chm13(variants, out_dir, tools["liftover"], resources["chain_grch38_to_t2t_bare"])

    tumour_rows, normal_rows = [], []
    for v in variants:
        donor = v["donor"]
        if donor not in samples:
            continue
        if t2t and v["chm13"] is None:
            print(f"WARNING: liftOver failed for FN {v['gene']} {v['chrom']}:{v['pos']} ({donor})")
            continue

        s = samples[donor]
        if t2t:
            bams = [os.path.join(alignment_dir, donor, s[k], f"{s[k]}_recalibrated.bam") for k in ["tumour", "normal"]]
            chrom, pos = f"CHM13#0#chr{v['chm13'][0]}", v["chm13"][1]
        else:
            bams = [os.path.join(alignment_dir, donor, s[k], f"{s[k + '_full']}.bam") for k in ["tumour", "normal"]]
            chrom, pos = f"chr{v['chrom']}", v["pos"]
        if not all(os.path.exists(bam) for bam in bams):
            print(f"WARNING: missing BAM(s) for FN {v['gene']} ({donor})")
            continue

        for bam, rows in zip(bams, [tumour_rows, normal_rows]):
            ref_count, alt_count = mpileup_counts(tools["samtools"], reference, bam, chrom, pos, v["ref"], v["alt"])
            rows += freq_rows(v, s["ancestry"], ref_count, alt_count)
    return tumour_rows, normal_rows


def main():
    parser = argparse.ArgumentParser(description="REF / ALT read counts of the Strelka2 FN and FP driver gene SNVs.")
    parser.add_argument("--alignment", choices=SUFFIX, required=True)
    parser.add_argument("--patients", default="metadata/BLCA_Patients.txt")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    out_dir = os.path.join(config["data_dir"], "False_Positive_Analysis", "Allele_Freq")
    os.makedirs(out_dir, exist_ok=True)
    annotation_dir = os.path.join(config["data_dir"], "Variant_Caller_Overlap", args.alignment, "Shared_Variants_Annotations")

    samples = read_samples(args.patients)
    driver_genes = read_driver_genes(config["resources"]["driver_genes"], config["tools"]["rscript"])
    print(f"{len(samples)} samples, {len(driver_genes)} driver genes")

    for impacts, impact_suffix in [(None, ""), ({"HIGH", "MODERATE"}, "_high_impact")]:
        suffix = SUFFIX[args.alignment] + impact_suffix
        for status, folder in [("fp", "False_Positives"), ("fn", "False_Negatives")]:
            vcf = os.path.join(annotation_dir, folder, f"{status.upper()}_filtered_shared_by_1_caller_in_capture_annotated.vcf")
            variants = read_annotated_vcf(vcf, driver_genes, impacts)
            print(f"{status.upper()}{suffix}: {len(variants)} driver gene variants")
            if status == "fp":
                tumour_rows, normal_rows = fp_counts(variants, samples, args, config)
            else:
                tumour_rows, normal_rows = fn_counts(variants, samples, args, config, out_dir)
            write_rows(tumour_rows, os.path.join(out_dir, f"tumour_{status}{suffix}_plot_data.txt"))
            write_rows(normal_rows, os.path.join(out_dir, f"normal_{status}{suffix}_plot_data.txt"))


if __name__ == "__main__":
    main()
