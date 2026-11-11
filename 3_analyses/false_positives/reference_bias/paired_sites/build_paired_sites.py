#!/usr/bin/env python3

# Reference bias, step 3: one list of somatic and germline sites per donor, counted later in every alignment.
#
# Every alignment is measured at the same loci:
#   Somatic sites:  the caller's VG pangenome TP + FP SNVs (GRCh38), lifted GRCh38 -> T2T
#   Germline sites: heterozygous SNVs from HaplotypeCaller on the VG pangenome normal BAM (T2T; GT 0/1, 0|1 or 1/2,
#                   REF + ALT >= 5 reads) within 1-74 bp of a somatic site on the same chromosome, lifted T2T -> GRCh38.
#                   Each germline SNV is kept once, with its nearest somatic site (ties: TP before FP, then lower position)
# A site is kept when it lifts over and both reference genomes carry one of its two alleles. The reference base of
# each genome is saved (grch38_ref, t2t_ref), since reference bias is the share of reads carrying the reference base
# of the genome each alignment used.
#
# Strelka2 uses the HaplotypeCaller run from germline_analysis/. Mutect2 needs its own run around its calls:
#   python3 .../build_paired_sites.py --caller mutect2 --hc-manifest   (then run_haplotypecaller_array.sh)
#   python3 .../build_paired_sites.py --caller mutect2
#
# Input:  <data_dir>/TCGA_Comparison/<caller>_vs_tcga/T2T_pangenome/.../ (compare_to_mc3.py)
#         Strelka2: <data_dir>/Reference_Bias/germline_analysis/vg_pangenome/<donor>/<donor>_vg_pangenome_germline_summary.tsv
#         Mutect2:  <data_dir>/Reference_Bias/paired_sites_mutect2/haplotypecaller/vg_pangenome/<donor>/...
# Output: <data_dir>/Reference_Bias/paired_sites[_mutect2]/
#             somatic_sites.tsv, germline_sites.tsv, sites_QC.tsv, sites/<donor>/ (bcftools -R region files)
#             manifest.tsv (9 rows per donor: 3 alignments x tumour at somatic sites, normal at somatic sites,
#             normal at germline sites); prints the sbatch command for run_paired_pileup_array.sh
#
# Usage:  python3 3_analyses/false_positives/reference_bias/paired_sites/build_paired_sites.py --caller strelka

import argparse
import os
import subprocess

import yaml

ALIGNMENTS = ["bwa_grch38", "vg_grch38", "vg_pangenome"]

# The VG GRCh38 alignment of this normal is labelled 10A; it is the same 10B sample
VG_GRCH38_SAMPLE_OVERRIDE = {"TCGA-E5-A4TZ-10B": "TCGA-E5-A4TZ-10A"}

DISTANCE = 74
MIN_DEPTH = 5
HET_GENOTYPES = {"0/1", "0|1", "1/2"}


def short_id(sample_id):
    """TCGA-FD-A3N6-01A-11D-A21A-08 -> TCGA-FD-A3N6-01A"""
    return "-".join(sample_id.split("-")[:4])


def read_patients(path):
    patients = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            donor, tumour, normal, ancestry = line.split()[:4]
            patients.append({"donor": donor, "tumour": tumour, "normal": normal, "ancestry": ancestry})
    return patients


class Fasta:
    """Single-base lookup from a .fai-indexed FASTA."""

    def __init__(self, path):
        self.handle = open(path, "rb")
        self.index = {}
        with open(path + ".fai") as f:
            for line in f:
                name, length, offset, line_bases, line_bytes = line.split("\t")[:5]
                self.index[name] = (int(length), int(offset), int(line_bases), int(line_bytes))

    def base(self, chrom, pos):
        if chrom not in self.index:
            return None
        length, offset, line_bases, line_bytes = self.index[chrom]
        if pos < 1 or pos > length:
            return None
        i = pos - 1
        self.handle.seek(offset + (i // line_bases) * line_bytes + i % line_bases)
        return self.handle.read(1).decode().upper()


def bam_path(alignments_dir, alignment, donor, sample_id):
    short = short_id(sample_id)
    if alignment == "bwa_grch38":
        return os.path.join(alignments_dir, "GRCh38_bwa", donor, short, f"{sample_id}.bam")
    if alignment == "vg_grch38":
        short = VG_GRCH38_SAMPLE_OVERRIDE.get(short, short)
        return os.path.join(alignments_dir, "GRCh38_linear", donor, short, f"{short}_recalibrated.bam")
    return os.path.join(alignments_dir, "T2T_pangenome", donor, short, f"{short}_recalibrated.bam")


def bam_index(bam):
    for candidate in (bam + ".bai", bam[:-4] + ".bai"):
        if os.path.exists(candidate):
            return candidate
    return None


# VG pangenome TP + FP SNVs in GRCh38 (chrom without chr); TP wins on a duplicated position
def load_pangenome_calls(calls_dir, fp_file, patient):
    call_dir = os.path.join(calls_dir, patient["ancestry"], patient["donor"], short_id(patient["tumour"]))
    sites = {}
    for filename, status in [("common_variants_full.txt", "TP"), (fp_file, "FP")]:
        with open(os.path.join(call_dir, filename)) as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                chrom, pos, ref, alt = fields[0], int(fields[2]), fields[3].upper(), fields[4].upper()
                if len(ref) != 1 or len(alt) != 1:
                    continue
                key = (chrom, pos)
                if key not in sites or status == "TP":
                    sites[key] = {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "status": status}
    return list(sites.values())


# records: dicts with id, chrom, pos. Returns {id: (chrom, pos)} for the ones that lift over
def run_liftover(records, liftover, chain, work_dir, prefix, in_chrom, out_chrom):
    in_bed = os.path.join(work_dir, f"{prefix}_in.bed")
    out_bed = os.path.join(work_dir, f"{prefix}_out.bed")
    unmapped_bed = os.path.join(work_dir, f"{prefix}_unmapped.bed")
    with open(in_bed, "w") as out:
        for r in records:
            out.write(f"{in_chrom(r['chrom'])}\t{r['pos'] - 1}\t{r['pos']}\t{r['id']}\n")
    subprocess.run([liftover, in_bed, chain, out_bed, unmapped_bed], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    lifted = {}
    with open(out_bed) as f:
        for line in f:
            chrom, _, end, site_id = line.rstrip("\n").split("\t")[:4]
            if site_id not in lifted:
                lifted[site_id] = (out_chrom(chrom), int(end))
    return lifted


# Heterozygous SNVs with REF + ALT >= 5 reads from the HaplotypeCaller summary, or None if there is no summary
def read_pangenome_germline(germline_dir, donor):
    path = os.path.join(germline_dir, donor, f"{donor}_vg_pangenome_germline_summary.tsv")
    if not os.path.exists(path):
        return None
    calls = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            ad = row["AD"].split(",")
            if row["GT"] not in HET_GENOTYPES or len(row["REF"]) != 1 or len(row["ALT"]) != 1 or len(ad) != 2:
                continue
            try:
                depth = int(ad[0]) + int(ad[1])
            except ValueError:
                continue
            if depth < MIN_DEPTH:
                continue
            calls.append({"chrom": row["CHROM"], "pos": int(row["POS"]), "ref": row["REF"].upper(), "alt": row["ALT"].upper()})
    return calls


# bcftools -R regions file (chrom, pos), sorted, one line per position
def write_regions(rows, path, chrom_key, pos_key):
    positions = sorted({(r[chrom_key], r[pos_key]) for r in rows})
    with open(path, "w") as out:
        for chrom, pos in positions:
            out.write(f"{chrom}\t{pos}\n")


# +/- 74 bp around each T2T somatic site, sorted and merged (HaplotypeCaller -L)
def write_hc_target_bed(sites, path):
    intervals = sorted((s["t2t_chrom"], max(0, s["t2t_pos"] - 1 - DISTANCE), s["t2t_pos"] + DISTANCE) for s in sites)
    merged = []
    for chrom, start, end in intervals:
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], end)
        else:
            merged.append([chrom, start, end])
    with open(path, "w") as out:
        for chrom, start, end in merged:
            out.write(f"{chrom}\t{start}\t{end}\n")


def write_table(rows, path, columns):
    with open(path, "w") as out:
        out.write("\t".join(columns) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in columns) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Build the paired somatic and germline sites for the reference bias analysis.")
    parser.add_argument("--caller", choices=["strelka", "mutect2"], default="strelka")
    parser.add_argument("--hc-manifest", action="store_true",
                        help="only write HaplotypeCaller target regions and hc_manifest.tsv, then stop")
    parser.add_argument("--patients", default="metadata/BLCA_Patients.txt")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    tools, resources = config["tools"], config["resources"]
    data_dir = config["data_dir"]
    reference_bias_dir = os.path.join(data_dir, "Reference_Bias")
    fasta = {"bwa_grch38": resources["fasta_grch38"], "vg_grch38": resources["fasta_grch38_primary"],
             "vg_pangenome": resources["fasta_t2t"]}

    # Per caller: VG pangenome TP / FP files, the HaplotypeCaller calls and the output folder
    if args.caller == "strelka":
        calls_dir = os.path.join(data_dir, "TCGA_Comparison", "strelka_vs_tcga", "T2T_pangenome")
        fp_file = "strelka_unique_pass.txt"
        germline_dir = os.path.join(reference_bias_dir, "germline_analysis", "vg_pangenome")
        output_root = os.path.join(reference_bias_dir, "paired_sites")
    else:
        calls_dir = os.path.join(data_dir, "TCGA_Comparison", "mutect_vs_tcga", "T2T_pangenome")
        fp_file = "mutect_unique_pass.txt"
        germline_dir = os.path.join(reference_bias_dir, "paired_sites_mutect2", "haplotypecaller", "vg_pangenome")
        output_root = os.path.join(reference_bias_dir, "paired_sites_mutect2")
    sites_dir = os.path.join(output_root, "sites")
    os.makedirs(sites_dir, exist_ok=True)
    print(f"Caller: {args.caller}  output: {output_root}")

    grch38_fasta = Fasta(fasta["bwa_grch38"])
    t2t_fasta = Fasta(fasta["vg_pangenome"])
    patients = read_patients(args.patients)
    print(f"Loaded {len(patients)} patients")

    all_somatic, all_germline, manifest, qc_rows, hc_manifest = [], [], [], [], []

    for patient in patients:
        donor = patient["donor"]
        work_dir = os.path.join(sites_dir, donor)
        os.makedirs(work_dir, exist_ok=True)

        # Somatic: VG pangenome TP + FP, GRCh38 -> T2T
        calls = load_pangenome_calls(calls_dir, fp_file, patient)
        for i, c in enumerate(calls):
            c["id"] = f"{donor}_S{i}"
        lifted = run_liftover(calls, tools["liftover"], resources["chain_grch38_to_t2t"], work_dir,
                              "somatic_grch38_to_t2t", in_chrom=lambda c: f"CHM13#0#chr{c}", out_chrom=lambda c: c)

        somatic, n_unlifted, n_ref_mismatch = [], 0, 0
        for c in calls:
            if c["id"] not in lifted:
                n_unlifted += 1
                continue
            t2t_chrom, t2t_pos = lifted[c["id"]]
            grch38_chrom = f"chr{c['chrom']}"
            grch38_ref = grch38_fasta.base(grch38_chrom, c["pos"])
            t2t_ref = t2t_fasta.base(t2t_chrom, t2t_pos)
            if grch38_ref not in (c["ref"], c["alt"]) or t2t_ref not in (c["ref"], c["alt"]):
                n_ref_mismatch += 1
                continue
            somatic.append({
                "site_id": c["id"], "donor": donor, "ancestry": patient["ancestry"], "status": c["status"],
                "grch38_chrom": grch38_chrom, "grch38_pos": c["pos"],
                "t2t_chrom": t2t_chrom, "t2t_pos": t2t_pos, "allele_1": c["ref"], "allele_2": c["alt"],
                "grch38_ref": grch38_ref, "t2t_ref": t2t_ref,
            })

        # --hc-manifest: HaplotypeCaller target regions only
        if args.hc_manifest:
            target_bed = os.path.join(work_dir, "hc_target_regions_t2t.bed")
            write_hc_target_bed(somatic, target_bed)
            normal_bam = bam_path(config["alignments_dir"], "vg_pangenome", donor, patient["normal"])
            if not os.path.exists(normal_bam):
                print(f"WARNING: missing pangenome normal BAM for {donor}: {normal_bam}")
            hc_manifest.append({
                "donor": donor, "alignment": "vg_pangenome", "ancestry": patient["ancestry"],
                "normal_bam": normal_bam, "reference_fasta": fasta["vg_pangenome"], "target_bed": target_bed,
                "out_prefix": os.path.join(germline_dir, donor, f"{donor}_vg_pangenome"), "exact_sites_file": target_bed,
            })
            print(f"{donor}: {len(somatic)} somatic sites -> HaplotypeCaller target regions")
            continue

        # Germline: heterozygous SNVs near those sites, T2T -> GRCh38
        germline_calls = read_pangenome_germline(germline_dir, donor)
        if germline_calls is None:
            print(f"WARNING: no pangenome germline summary for {donor}")
            germline_calls = []

        sites_by_chrom = {}
        for s in somatic:
            sites_by_chrom.setdefault(s["t2t_chrom"], []).append(s)

        near = []
        for g in germline_calls:
            nearby = [s for s in sites_by_chrom.get(g["chrom"], []) if 1 <= abs(g["pos"] - s["t2t_pos"]) <= DISTANCE]
            if not nearby:
                continue
            g["nearest"] = min(nearby, key=lambda s: (abs(g["pos"] - s["t2t_pos"]), s["status"] != "TP", s["t2t_pos"]))
            g["id"] = f"{donor}_G{len(near)}"
            near.append(g)

        lifted_g = run_liftover(near, tools["liftover"], resources["chain_t2t_to_grch38"], work_dir,
                                "germline_t2t_to_grch38", in_chrom=lambda c: c, out_chrom=lambda c: c.replace("CHM13#0#", ""))
        germline, n_g_unlifted, n_g_mismatch = [], 0, 0
        for g in near:
            if g["id"] not in lifted_g:
                n_g_unlifted += 1
                continue
            grch38_chrom, grch38_pos = lifted_g[g["id"]]
            t2t_ref = t2t_fasta.base(g["chrom"], g["pos"])
            grch38_ref = grch38_fasta.base(grch38_chrom, grch38_pos)
            if t2t_ref not in (g["ref"], g["alt"]) or grch38_ref not in (g["ref"], g["alt"]):
                n_g_mismatch += 1
                continue
            site = g["nearest"]
            germline.append({
                "site_id": g["id"], "donor": donor, "ancestry": patient["ancestry"],
                "somatic_site_id": site["site_id"], "somatic_status": site["status"],
                "distance": abs(g["pos"] - site["t2t_pos"]),
                "grch38_chrom": grch38_chrom, "grch38_pos": grch38_pos,
                "t2t_chrom": g["chrom"], "t2t_pos": g["pos"], "allele_1": g["ref"], "allele_2": g["alt"],
                "grch38_ref": grch38_ref, "t2t_ref": t2t_ref,
            })

        all_somatic += somatic
        all_germline += germline

        # Region files and manifest rows
        regions = {}
        for kind, rows in [("somatic", somatic), ("germline", germline)]:
            for coords, chrom_key, pos_key in [("grch38", "grch38_chrom", "grch38_pos"), ("t2t", "t2t_chrom", "t2t_pos")]:
                path = os.path.join(work_dir, f"{kind}_regions_{coords}.tsv")
                write_regions(rows, path, chrom_key, pos_key)
                regions[(kind, coords)] = path

        for alignment in ALIGNMENTS:
            coords = "t2t" if alignment == "vg_pangenome" else "grch38"
            for sample_type, sample_id, kind in [("tumour", patient["tumour"], "somatic"),
                                                 ("normal_somatic", patient["normal"], "somatic"),
                                                 ("normal", patient["normal"], "germline")]:
                bam = bam_path(config["alignments_dir"], alignment, donor, sample_id)
                index = bam_index(bam)
                if not os.path.exists(bam) or index is None:
                    print(f"WARNING: missing BAM or index for {donor}/{alignment}/{sample_type}: {bam}")
                manifest.append({
                    "donor": donor, "alignment": alignment, "sample_type": sample_type,
                    "bam": bam, "bam_index": index, "reference_fasta": fasta[alignment],
                    "regions": regions[(kind, coords)],
                    "out_file": os.path.join(output_root, "counts", alignment, donor, f"{donor}_{alignment}_{sample_type}_counts.tsv"),
                })

        qc_rows.append({
            "Donor": donor, "Ancestry": patient["ancestry"],
            "Pangenome_TP_FP_sites": len(calls), "Somatic_unlifted": n_unlifted,
            "Somatic_allele_mismatch": n_ref_mismatch, "Somatic_sites_kept": len(somatic),
            "Germline_hets_within_74bp": len(near), "Germline_unlifted": n_g_unlifted,
            "Germline_allele_mismatch": n_g_mismatch, "Germline_sites_kept": len(germline),
        })
        print(f"{donor}: somatic {len(somatic)}/{len(calls)} kept, germline {len(germline)}/{len(near)} kept")

    sbatch = f"sbatch {config.get('slurm_options', '')}"
    if args.hc_manifest:
        hc_file = os.path.join(output_root, "hc_manifest.tsv")
        for i, row in enumerate(hc_manifest):
            row["array_index"] = i
        write_table(hc_manifest, hc_file, ["array_index", "donor", "alignment", "ancestry", "normal_bam",
                                           "reference_fasta", "target_bed", "out_prefix", "exact_sites_file"])
        log_dir = os.path.join(output_root, "haplotypecaller", "logs")
        os.makedirs(log_dir, exist_ok=True)
        print(f"\nHaplotypeCaller manifest: {hc_file} ({len(hc_manifest)} tasks). Next:")
        print(f"  {sbatch} --array=0-{len(hc_manifest) - 1} --output={log_dir}/%x_%A_%a.out --error={log_dir}/%x_%A_%a.err "
              f"--export=ALL,MANIFEST={hc_file},GATK={tools['gatk']},BCFTOOLS={tools['bcftools']} "
              f"3_analyses/false_positives/reference_bias/germline_analysis/run_haplotypecaller_array.sh")
        return

    write_table(all_somatic, os.path.join(output_root, "somatic_sites.tsv"),
                ["site_id", "donor", "ancestry", "status", "grch38_chrom", "grch38_pos", "t2t_chrom", "t2t_pos",
                 "allele_1", "allele_2", "grch38_ref", "t2t_ref"])
    write_table(all_germline, os.path.join(output_root, "germline_sites.tsv"),
                ["site_id", "donor", "ancestry", "somatic_site_id", "somatic_status", "distance",
                 "grch38_chrom", "grch38_pos", "t2t_chrom", "t2t_pos", "allele_1", "allele_2", "grch38_ref", "t2t_ref"])
    for i, row in enumerate(manifest):
        row["array_index"] = i
    manifest_file = os.path.join(output_root, "manifest.tsv")
    write_table(manifest, manifest_file, ["array_index", "donor", "alignment", "sample_type", "bam", "bam_index",
                                          "reference_fasta", "regions", "out_file"])
    write_table(qc_rows, os.path.join(output_root, "sites_QC.tsv"), list(qc_rows[0].keys()))

    log_dir = os.path.join(output_root, "logs")
    os.makedirs(log_dir, exist_ok=True)
    print(f"\nSomatic sites: {len(all_somatic)}  Germline sites: {len(all_germline)}")
    print(f"Manifest: {manifest_file} ({len(manifest)} rows, {len(patients)} array tasks). Next:")
    print(f"  {sbatch} --array=0-{len(patients) - 1} --output={log_dir}/%x_%A_%a.out --error={log_dir}/%x_%A_%a.err "
          f"--export=ALL,MANIFEST={manifest_file},BCFTOOLS={tools['bcftools']} "
          f"3_analyses/false_positives/reference_bias/paired_sites/run_paired_pileup_array.sh")


if __name__ == "__main__":
    main()
