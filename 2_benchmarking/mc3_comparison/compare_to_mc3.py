#!/usr/bin/env python3

# Compares somatic SNV calls with the MC3 truth set, per donor, by CHROM and POS in GRCh38 coordinates.
#
# For each donor:
#   1. Keep PASS SNVs with normal depth > 8 and tumour depth >= 5 that lie in the Agilent targets
#   2. VG pangenome (T2T_pangenome) calls are lifted from T2T-CHM13 to GRCh38; GRCh38 alignments are used as is
#   3. MC3 truth set: PASS SNVs for the same tumour-normal pair with tumour depth >= 5 and normal depth > 8,
#      in the Agilent targets. For T2T_pangenome, MC3 sites that are not callable in the pangenome
#      alignments are removed from both the truth set and the calls
#   4. TP, FP and FN by CHROM + POS (REF/ALT are kept in the output files but not compared)
#
# With --gnomad (Supplementary Figure 9i), calls found in gnomAD (exact CHROM, POS, REF and ALT match) are
# removed after step 1. The MC3 truth set is not filtered.
#
# Input:  <calls_dir>/<caller>/<alignment>/<donor>/<tumour>/<tumour>_<caller>_pass_variants.vcf.gz
#         (SomaticSniper: _hc_filtered.vcf.gz)
#         MC3 MAF, target BEDs, liftover chain, non-callable sites and gnomAD from config.yaml
# Output: <filtered_calls_dir>/<caller>/<alignment>/<donor>/<tumour>/<tumour>_<caller>_pass_variants_DP_target.vcf.gz
#         <data_dir>/[LUAD/]TCGA_Comparison/<caller>_vs_tcga/<alignment>/<ancestry>/<donor>/<tumour>/
#             common_variants_full.txt (TP), <caller>_unique_pass.txt (FP), tcga_unique_pass.txt (FN)
#         <data_dir>/[LUAD/]Stats_Summary/<caller>/stats_table_<alignment>.txt (one row per donor, appended)
#         With --gnomad, everything is written under <data_dir>/[LUAD/]Gnomad_Filtering/<caller>/:
#             <alignment>/<ancestry>/<donor>/<tumour>/ (filtered VCFs, calls removed by gnomAD, TP/FP/FN files)
#             Stats_Summary/stats_table_<alignment>_hard_gnomad.txt
#
# Usage:  python3 2_benchmarking/mc3_comparison/compare_to_mc3.py --cohort BLCA --caller Strelka --alignment T2T_pangenome
#         python3 2_benchmarking/mc3_comparison/compare_to_mc3.py --caller Mutect2 --alignment T2T_pangenome --gnomad

import argparse
import os
import subprocess

import yaml

# Differences between callers:
#   vcf_suffix:   name of the PASS VCF from pangenome-pipeline
#   prefix:       name used for the comparison folder and output files
#   depth_format: per-sample fields given to bcftools query
#   depth_filter: awk test on those fields ($5, $6) for normal depth > 8 and tumour depth >= 5.
#                 Strelka2 and SomaticSniper list NORMAL then TUMOUR; Mutect2 lists the tumour first.
#                 DeepSomatic has one sample column with NDP (normal) then DP (tumour)
CALLERS = {
    "Strelka": {
        "vcf_suffix": "Strelka_pass_variants",
        "prefix": "strelka",
        "depth_format": "[%DP\\t]",
        "depth_filter": "$5>8 && $6>=5",
    },
    "Mutect2": {
        "vcf_suffix": "Mutect2_pass_variants",
        "prefix": "mutect",
        "depth_format": "[%DP\\t]",
        "depth_filter": "$6>8 && $5>=5",
    },
    "SomaticSniper": {
        "vcf_suffix": "SomaticSniper_hc_filtered",
        "prefix": "somaticsniper",
        "depth_format": "[%DP\\t]",
        "depth_filter": "$5>8 && $6>=5",
    },
    "DeepSomatic": {
        "vcf_suffix": "DeepSomatic_pass_variants",
        "prefix": "deepsomatic",
        "depth_format": "[%NDP\\t%DP\\t]",
        "depth_filter": "$5>8 && $6>=5",
    },
}

# Differences between cohorts:
#   subdir:          output folder inside data_dir
#   short_mc3_ids:   match MC3 barcodes on the first four fields (TCGA-XX-XXXX-01A) instead of the full barcode
#   default_patients
COHORTS = {
    "BLCA": {"subdir": "", "short_mc3_ids": False, "default_patients": "metadata/BLCA_Patients.txt"},
    "LUAD": {"subdir": "LUAD", "short_mc3_ids": True, "default_patients": "metadata/LUAD_All_Patients_Full_ID.txt"},
}

ALIGNMENTS = {
    "GRCh38_bwa": "stats_table_grch38_bwa.txt",
    "GRCh38_linear": "stats_table_grch38_linear.txt",
    "GRCh38_linear_plus": "stats_table_grch38_linear_plus.txt",
    "GRCh38_hprc_v1.1": "stats_table_grch38_hprc_v1.1.txt",
    "T2T_pangenome": "stats_table_pangenome.txt",
}

STATS_HEADER = [
    "Donor_ID", "Ancestry", "Variant_Caller", "Reference",
    "Pre_gnomAD_Calls", "Removed_by_gnomAD", "Post_gnomAD_Calls",
    "True_Positive", "False_Positive", "False_Negative", "Non_Callable_TCGA",
    "Total_Variant_Called", "Total_TCGA",
    "Recall(%)", "Precision(%)", "F1_Score(%)",
]


def run_command(cmd):
    print(f"Running:\n{cmd}\n")
    subprocess.run(cmd, shell=True, check=True)


def count_vcf_variants(vcf, bcftools):
    result = subprocess.run(f"{bcftools} view -H {vcf} | wc -l", shell=True, check=True,
                            capture_output=True, text=True)
    return int(result.stdout.strip())


def count_lines(path):
    if not os.path.exists(path):
        return 0
    with open(path) as f:
        return sum(1 for line in f if line.strip())


def safe_div(n, d):
    return n / d if d else 0.0


def short_id(sample_id):
    """TCGA-2F-A9KP-01A-11D-A38G-08 -> TCGA-2F-A9KP-01A"""
    return "-".join(sample_id.split("-")[:4])


def read_patients(path):
    """Patient list columns: donor, tumour barcode, normal barcode, ancestry."""
    patients = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            entry = line.split()
            if len(entry) < 4:
                print(f"Skipping bad line: {line.strip()}")
                continue
            patients.append(entry[:4])
    return patients


# Keeps PASS SNVs with normal depth > 8 and tumour depth >= 5 that lie in the target BED
def filter_vcf(donor, tumour, caller, alignment, config, out_dir):
    tools = config["tools"]
    settings = CALLERS[caller]
    tumour_short = short_id(tumour)

    input_vcf = os.path.join(config["calls_dir"], caller, alignment, donor, tumour_short,
                             f"{tumour_short}_{settings['vcf_suffix']}.vcf.gz")
    if not os.path.exists(input_vcf):
        print(f"Missing {caller} VCF: {input_vcf}")
        return None

    depth_txt = os.path.join(out_dir, f"{tumour_short}_depths.txt")
    pass_sites = os.path.join(out_dir, f"{tumour_short}_DP_sites.bed")
    filtered_vcf = os.path.join(out_dir, f"{tumour_short}_{caller}_pass_variants_DP_target.vcf.gz")

    if alignment == "T2T_pangenome":
        target_bed = config["resources"]["target_bed_t2t"]
    else:
        target_bed = config["resources"]["target_bed_grch38"]

    run_command(
        f"{tools['bcftools']} query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t{settings['depth_format']}\\n' "
        f"{input_vcf} > {depth_txt}"
    )
    run_command(
        f"""awk 'BEGIN{{OFS="\\t"}}
        length($3)==1 && length($4)==1 && {settings['depth_filter']} {{
            print $1, $2-1, $2
        }}' {depth_txt} > {pass_sites}"""
    )
    run_command(
        f"{tools['bedtools']} intersect -u -header -a {input_vcf} -b {pass_sites} | "
        f"{tools['bedtools']} intersect -u -header -a - -b {target_bed} | "
        f"{tools['bgzip']} -c > {filtered_vcf}"
    )
    run_command(f"{tools['tabix']} -f -p vcf {filtered_vcf}")

    return filtered_vcf


# Splits the calls into those absent from gnomAD (kept) and those present (removed); -c none needs REF and ALT to match
def remove_gnomad(filtered_vcf, tumour_short, caller, alignment, config, out_dir):
    tools = config["tools"]
    if alignment == "T2T_pangenome":
        gnomad_vcf = config["resources"]["gnomad_t2t"]
    else:
        gnomad_vcf = config["resources"]["gnomad_grch38"]

    kept_vcf = os.path.join(out_dir, f"{tumour_short}_{caller}_hard_gnomad_filtered.vcf.gz")
    removed_vcf = os.path.join(out_dir, f"{tumour_short}_{caller}_removed_by_gnomad.vcf.gz")

    run_command(f"{tools['bcftools']} isec -c none -C -w1 -Oz -o {kept_vcf} {filtered_vcf} {gnomad_vcf}")
    run_command(f"{tools['tabix']} -f -p vcf {kept_vcf}")
    run_command(f"{tools['bcftools']} isec -c none -n=2 -w1 -Oz -o {removed_vcf} {filtered_vcf} {gnomad_vcf}")
    run_command(f"{tools['tabix']} -f -p vcf {removed_vcf}")

    pre = count_vcf_variants(filtered_vcf, tools["bcftools"])
    kept = count_vcf_variants(kept_vcf, tools["bcftools"])
    removed = count_vcf_variants(removed_vcf, tools["bcftools"])
    if pre != kept + removed:
        print(f"WARNING: calls before gnomAD ({pre}) != kept ({kept}) + removed ({removed})")

    return kept_vcf, [str(pre), str(removed), str(kept)]


# MC3 SNVs for one tumour-normal pair, written as CHROM POS REF ALT FILTER tumour_GT tumour_DP normal_GT normal_DP
def extract_mc3(tumour, normal, short_ids, mc3_maf, output_file):
    if short_ids:
        tumour, normal = short_id(tumour), short_id(normal)
        shorten = """
            split(tumour_id, tumour_parts, "-")
            tumour_id = tumour_parts[1] "-" tumour_parts[2] "-" tumour_parts[3] "-" tumour_parts[4]
            split(normal_id, normal_parts, "-")
            normal_id = normal_parts[1] "-" normal_parts[2] "-" normal_parts[3] "-" normal_parts[4]
"""
    else:
        shorten = ""

    run_command(
        f"""gzip -dc {mc3_maf} | awk -F'\\t' -v OFS='\\t' \
        -v tumour='{tumour}' \
        -v normal='{normal}' '
        function clean(x) {{
            gsub(/^"|"$/, "", x)
            return x
        }}

        /^#/ {{
            next
        }}

        !header_found {{
            for (i = 1; i <= NF; i++) {{
                name = clean($i)

                if (name == "Chromosome")                  chrom_col = i
                if (name == "Start_Position")              pos_col = i
                if (name == "Reference_Allele")            ref_col = i
                if (name == "Tumor_Seq_Allele2")           alt_col = i
                if (name == "Variant_Type")                type_col = i
                if (name == "Tumor_Sample_Barcode")        tumour_col = i
                if (name == "Matched_Norm_Sample_Barcode") normal_col = i
                if (name == "t_depth")                     tumour_dp_col = i
                if (name == "n_depth")                     normal_dp_col = i
                if (name == "FILTER")                      filter_col = i
            }}

            header_found = 1
            next
        }}

        {{
            chrom       = clean($chrom_col)
            pos         = clean($pos_col)
            ref         = clean($ref_col)
            alt         = clean($alt_col)
            variant     = clean($type_col)
            tumour_id   = clean($tumour_col)
            normal_id   = clean($normal_col)
            tumour_dp   = clean($tumour_dp_col)
            normal_dp   = clean($normal_dp_col)
            filter      = clean($filter_col)
{shorten}
            if (tumour_id == tumour &&
                normal_id == normal &&
                variant == "SNP") {{
                print chrom, pos, ref, alt, filter, \
                    "0/1", tumour_dp, "0/0", normal_dp
            }}
        }}
        ' > {output_file}"""
    )


# Non-callable MC3 sites for one donor, as CHROM POS (the GRCh38 column looks like chr21:10485700T>C)
def make_noncallable_positions(donor, noncallable_file, output_file):
    run_command(
        f"""awk -F'\\t' -v donor='{donor}' 'BEGIN{{OFS="\\t"}}
        NR>1 {{
            gsub(/"/, "", $1);
            gsub(/"/, "", $4);

            if ($1 == donor) {{
                split($4, a, ":");
                chrom = a[1];
                sub(/^chr/, "", chrom);

                match(a[2], /^[0-9]+/);
                pos = substr(a[2], RSTART, RLENGTH);

                if (chrom != "" && pos != "") {{
                    print chrom, pos;
                }}
            }}
        }}' {noncallable_file} > {output_file}"""
    )


# Drops rows of input_file whose CHROM + POS (columns 1 and 3) are in the non-callable list
def remove_noncallable(noncallable_positions, input_file, output_file):
    run_command(
        f"""awk 'BEGIN{{OFS="\\t"}}
        FILENAME == badfile {{bad[$1 FS $2]=1; next}}
        !(($1 FS $3) in bad) {{print}}
        ' badfile={noncallable_positions} {noncallable_positions} {input_file} > {output_file}"""
    )


# Rows of query_file whose CHROM + POS are (keep_shared=True) or are not (False) in ref_file, once per position
def match_positions(ref_file, query_file, output_file, keep_shared):
    test = "(key in ref)" if keep_shared else "!(key in ref)"
    run_command(
        f"""awk 'BEGIN{{OFS="\t"}}
        NR==FNR {{
            ref[$1 SUBSEP $3]=1;
            next
        }}
        {{
            key=$1 SUBSEP $3;
            if ({test} && !(key in written)) {{
                print;
                written[key]=1;
            }}
        }}' {ref_file} {query_file} > {output_file}"""
    )


def write_stats_row(stats_file, values, gnomad):
    header = STATS_HEADER if gnomad else [c for c in STATS_HEADER if "gnomAD" not in c]
    write_header = not os.path.exists(stats_file) or os.path.getsize(stats_file) == 0
    with open(stats_file, "a") as f:
        if write_header:
            f.write("\t".join(header) + "\n")
        f.write("\t".join(values) + "\n")
    print(f"Appended metrics to: {stats_file}")


def process_patient(donor, tumour, normal, ancestry, args, config):
    print(f"\nStarting processing for Donor: {donor} [{args.alignment}]")
    tools = config["tools"]
    resources = config["resources"]
    cohort = COHORTS[args.cohort]
    prefix = CALLERS[args.caller]["prefix"]
    t2t = args.alignment == "T2T_pangenome"
    tumour_short = short_id(tumour)

    if args.gnomad:
        gnomad_dir = os.path.join(config["data_dir"], cohort["subdir"], "Gnomad_Filtering", args.caller)
        out_dir = os.path.join(gnomad_dir, args.alignment, ancestry, donor, tumour_short)
        filtered_dir = out_dir
    else:
        out_dir = os.path.join(config["data_dir"], cohort["subdir"], "TCGA_Comparison", f"{prefix}_vs_tcga",
                               args.alignment, ancestry, donor, tumour_short)
        filtered_dir = os.path.join(config["filtered_calls_dir"], args.caller, args.alignment, donor, tumour_short)
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(filtered_dir, exist_ok=True)

    calls_vcf = filter_vcf(donor, tumour, args.caller, args.alignment, config, filtered_dir)
    if calls_vcf is None:
        print(f"Skipping {donor} due to missing {args.caller} VCF.")
        return ""

    gnomad_counts = []
    if args.gnomad:
        calls_vcf, gnomad_counts = remove_gnomad(calls_vcf, tumour_short, args.caller, args.alignment, config, out_dir)

    calls_bed = f"{out_dir}/{prefix}_patient.bed"
    calls_liftover_bed = f"{out_dir}/{prefix}_grch38.bed"
    calls_positions = f"{out_dir}/{prefix}_grch38_positions.txt"
    calls_eval_positions = f"{out_dir}/{prefix}_grch38_eval_positions.txt"
    mc3_depths = f"{out_dir}/tcga_patient_depths.txt"
    mc3_unrestricted = f"{out_dir}/tcga_patient_positions_before_target.bed"
    mc3_positions = f"{out_dir}/tcga_patient_target_positions.txt"
    mc3_callable_positions = f"{out_dir}/tcga_patient_callable_positions.txt"
    noncallable_positions = f"{out_dir}/noncallable_tcga_positions.txt"
    tp_file = f"{out_dir}/common_variants_full.txt"
    fp_file = f"{out_dir}/{prefix}_unique_pass.txt"
    fn_file = f"{out_dir}/tcga_unique_pass.txt"

    # 1) Calls -> BED (SNVs only)
    run_command(
        f"{tools['bcftools']} query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {calls_vcf} | "
        f"awk 'BEGIN{{OFS=\"\\t\"}} "
        f"length($3)==1 && length($4)==1 "
        f"{{print $1, $2-1, $2, $3, $4}}' > {calls_bed}"
    )

    # 2) Calls in GRCh38 coordinates, without the chr prefix
    if t2t:
        run_command(
            f"{tools['liftover']} {calls_bed} {resources['chain_t2t_to_grch38']} "
            f"{calls_liftover_bed} {out_dir}/unmapped_variants.bed"
        )
        run_command(f"sed -i 's/CHM13#0#chr//' {calls_liftover_bed}")
        grch38_bed = calls_liftover_bed
    else:
        grch38_bed = calls_bed
    run_command(f"cut -f1-5 {grch38_bed} | sed -E 's/^chr([0-9XYM]+)/\\1/' > {calls_positions}")

    # 3) MC3 truth set: PASS heterozygous SNVs, tumour depth >= 5, normal depth > 8, in the GRCh38 targets
    extract_mc3(tumour, normal, cohort["short_mc3_ids"], resources["mc3_maf"], mc3_depths)
    run_command(
        f"""awk 'BEGIN{{OFS="\\t"}}
        length($3)==1 && length($4)==1 &&
        $5=="PASS" &&
        ($6=="0/1" || $6=="1/0") &&
        $7>=5 &&
        $9>8 {{
            chrom=$1;
            if (chrom !~ /^chr/) chrom="chr" chrom;
            print chrom, $2-1, $2, $3, $4
        }}' {mc3_depths} > {mc3_unrestricted}"""
    )
    run_command(
        f"{tools['bedtools']} intersect -u -a {mc3_unrestricted} -b {resources['target_bed_grch38']} | "
        f"sed -E 's/^chr([0-9XYM]+)/\\1/' > {mc3_positions}"
    )

    # 4) T2T_pangenome only: remove MC3 sites that are not callable in the pangenome alignments
    if t2t:
        make_noncallable_positions(donor, resources["noncallable_sites"], noncallable_positions)
        remove_noncallable(noncallable_positions, mc3_positions, mc3_callable_positions)
        remove_noncallable(noncallable_positions, calls_positions, calls_eval_positions)
        mc3_eval, calls_eval = mc3_callable_positions, calls_eval_positions
        non_callable = count_lines(noncallable_positions)
    else:
        mc3_eval, calls_eval = mc3_positions, calls_positions
        non_callable = 0

    # 5) TP, FP and FN by CHROM + POS
    match_positions(mc3_eval, calls_eval, tp_file, keep_shared=True)
    match_positions(mc3_eval, calls_eval, fp_file, keep_shared=False)
    match_positions(calls_eval, mc3_eval, fn_file, keep_shared=False)

    tp, fp, fn = count_lines(tp_file), count_lines(fp_file), count_lines(fn_file)
    recall = safe_div(tp, tp + fn)
    precision = safe_div(tp, tp + fp)
    f1 = safe_div(2 * precision * recall, precision + recall)

    if args.gnomad:
        stats_dir = os.path.join(gnomad_dir, "Stats_Summary")
        stats_name = ALIGNMENTS[args.alignment].replace(".txt", "_hard_gnomad.txt")
    else:
        stats_dir = os.path.join(config["data_dir"], cohort["subdir"], "Stats_Summary", args.caller)
        stats_name = ALIGNMENTS[args.alignment]
    os.makedirs(stats_dir, exist_ok=True)
    write_stats_row(os.path.join(stats_dir, stats_name), [
        donor, ancestry, args.caller, args.alignment, *gnomad_counts,
        str(tp), str(fp), str(fn), str(non_callable), str(tp + fp), str(tp + fn),
        f"{recall * 100:.2f}", f"{precision * 100:.2f}", f"{f1 * 100:.2f}",
    ], args.gnomad)

    return (
        f"{donor} [{ancestry}] {args.alignment} — TP={tp}, FP={fp}, FN={fn}, NonCallable={non_callable}, "
        f"Prec={precision:.3f}, Rec={recall:.3f}, F1={f1:.3f} -> {out_dir}"
    )


def main():
    parser = argparse.ArgumentParser(description="Compare somatic SNV calls with MC3.")
    parser.add_argument("--cohort", choices=COHORTS, default="BLCA")
    parser.add_argument("--caller", choices=CALLERS, required=True)
    parser.add_argument("--alignment", choices=ALIGNMENTS, required=True)
    parser.add_argument("--gnomad", action="store_true", help="remove calls found in gnomAD before comparing")
    parser.add_argument("--patients", help="patient list (default: the cohort's list in metadata/)")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    patients = read_patients(args.patients or COHORTS[args.cohort]["default_patients"])

    summaries = [process_patient(*patient, args, config) for patient in patients]

    print("\nFinal Summary:")
    for summary in summaries:
        if summary:
            print(summary)
    print("\nPipeline completed!")


if __name__ == "__main__":
    main()
