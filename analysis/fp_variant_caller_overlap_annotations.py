# Shared FP Overlap + Annotations Builder (no snpEff!)
# Part 1: Reads caller-specific "unique_pass" FP lists (vs TCGA) for Strelka, SomaticSniper, Mutect,
#         counts overlaps (2+ callers) per donor, and writes per-donor outputs + a summary.
# Part 2: From mega list, writes FP_shared_by_{1,2,3}_callers_noduplicates.txt to the
#         Variant_Caller_Overlap/<GRCh38_Linear|GRCh38_Pangenome>/Shared_Variants_Annotations/False_Positives/.
# Part 3 (NO snpEff, please just do it in bash following the notion guide): Post-process the three files:
#         - drop headers to BED
#         - intersect with TCGA capture kit
#         - filter back to TXT
#         - build VCFs to run Snpeff

import os
import subprocess
from collections import Counter, defaultdict

base = "/workspace/lab/houlahanlab/phamc20/data"
mega_base = "/workspace/lab/houlahanlab/data/data/mega_variant_summary"
capture_bed = "/workspace/lab/houlahanlab/resources/tcga_capture_kit/Filtering_TCGA_MC3_grch38.bed"

def run_command(cmd: str):
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    print("Command completed successfully!\n")

def map_ref_dirs(ref_choice_user: str):
    """
    Map user's choice ('GRCh38' or 'T2T-CHM13') to:
      - caller_ref: directory name used by strelka_vs_tcga2/somatic_vs_tcga2/mutect_vs_tcga inputs
      - overlap_ref: directory name for Variant_Caller_Overlap outputs
    """
    if ref_choice_user == "GRCh38":
        return "GRCh38", "GRCh38_Linear"
    else:
        return "T2T-CHM13", "GRCh38_Pangenome"

def load_variants(path: str):
    if not os.path.exists(path):
        print(f"File not found: {path}")
        return set()
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())

##### Part 1: Per-donor 2+ callers overlap #####

def compare_shared_false_positives(donor_id: str, ancestry: str, ref_choice_user: str):
    caller_ref, overlap_ref = map_ref_dirs(ref_choice_user)

    # Caller-vs-TCGA input paths
    strelka_file = os.path.join(
        base, f"strelka_vs_tcga2/{caller_ref}/{ancestry}/{donor_id}/strelka_unique_pass.txt"
    )
    somatic_file = os.path.join(
        base, f"somatic_vs_tcga2/{caller_ref}/{ancestry}/{donor_id}/somatic_unique_pass.txt"
    )
    mutect_file = os.path.join(
        base, f"mutect_vs_tcga/{caller_ref}/{ancestry}/{donor_id}/mutect_unique_pass.txt"
    )

    # Per-donor overlap outputs
    output_dir = os.path.join(
        base, f"Variant_Caller_Overlap/{overlap_ref}/Shared_Variants/False_Positives/{ancestry}/{donor_id}"
    )
    os.makedirs(output_dir, exist_ok=True)

    output_file_simple  = os.path.join(output_dir, "fp_shared_by_2_or_more_variant_callers.txt")
    output_file_sources = os.path.join(output_dir, "fp_shared_variants_with_sources.txt")

    # Load variants
    strelka = load_variants(strelka_file)
    somatic = load_variants(somatic_file)
    mutect  = load_variants(mutect_file)

    all_variants = list(strelka | somatic | mutect)
    variant_counter = Counter()
    variant_sources  = {}

    for var in all_variants:
        sources = []
        if var in strelka: sources.append("Strelka")
        if var in somatic: sources.append("SomaticSniper")
        if var in mutect:  sources.append("Mutect")
        variant_counter[var] = len(sources)
        variant_sources[var] = "+".join(sources)

    shared_by_2_or_more = [v for v, c in variant_counter.items() if c >= 2]

    strelka_mutect = strelka_somatic = somatic_mutect = all_three = 0

    with open(output_file_simple, "w") as f_simple, open(output_file_sources, "w") as f_sources:
        for var in sorted(shared_by_2_or_more):
            parts = var.split("\t")
            chrom = parts[0] if len(parts) > 0 else ""
            pos   = parts[1] if len(parts) > 1 else ""
            sources = variant_sources[var]

            f_simple.write(var + "\n")
            f_sources.write(f"{chrom}\t{pos}\t{sources}\n")

            if sources == "Strelka+Mutect": strelka_mutect += 1
            elif sources == "Strelka+SomaticSniper": strelka_somatic += 1
            elif sources == "SomaticSniper+Mutect": somatic_mutect += 1
            elif sources == "Strelka+SomaticSniper+Mutect": all_three += 1

    return {
        "summary_string": (f"{donor_id}: {len(shared_by_2_or_more)} FP shared by 2+ callers "
                           f"(Strelka+Mutect: {strelka_mutect}, "
                           f"Strelka+SomaticSniper: {strelka_somatic}, "
                           f"SomaticSniper+Mutect: {somatic_mutect}, "
                           f"All three: {all_three})"),
        "strelka_mutect": strelka_mutect,
        "strelka_somatic": strelka_somatic,
        "somatic_mutect": somatic_mutect,
        "all_three": all_three
    }

##### Part 2: Mega-list annotations (no snpEff here) #####

def build_fp_annotations_from_mega(ref_choice_user: str):
    _, overlap_ref = map_ref_dirs(ref_choice_user)

    input_file = os.path.join(mega_base, overlap_ref, "FP_variants.txt")

    output_dir = os.path.join(
        base, f"Variant_Caller_Overlap/{overlap_ref}/Shared_Variants_Annotations/False_Positives"
    )
    os.makedirs(output_dir, exist_ok=True)

    out_1 = os.path.join(output_dir, "FP_shared_by_1_caller_noduplicates.txt")
    out_2 = os.path.join(output_dir, "FP_shared_by_2_callers_noduplicates.txt")
    out_3 = os.path.join(output_dir, "FP_shared_by_3_callers_noduplicates.txt")

    variant_to_tools = defaultdict(set)
    variant_to_start = {}

    if not os.path.exists(input_file):
        print(f" Mega list not found: {input_file}")
        return (out_1, out_2, out_3)

    with open(input_file, "r") as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 8: continue
            status, tool = parts[6], parts[7]
            if status != "FP": continue
            chrom, start, end, ref, alt, sample = parts[:6]
            key = f"{chrom}\t{end}\t{ref}\t{alt}\t{sample}"
            variant_to_tools[key].add(tool)
            variant_to_start[key] = start

    with open(out_1, "w") as f1, open(out_2, "w") as f2, open(out_3, "w") as f3:
        header_line = "CHROM\tSTART\tEND\tREF\tALT\tSAMPLE\tTOOL\n"
        for h in (f1, f2, f3): h.write(header_line)

        for key, tools in variant_to_tools.items():
            chrom, end, ref, alt, sample = key.split("\t")
            start = variant_to_start.get(key, ".")
            tool_str = "+".join(sorted(tools))
            line = f"{chrom}\t{start}\t{end}\t{ref}\t{alt}\t{sample}\t{tool_str}\n"
            if len(tools) == 1: f1.write(line)
            elif len(tools) == 2: f2.write(line)
            elif len(tools) == 3: f3.write(line)

    print(" Annotation files written to:")
    print(f" - {out_1}\n - {out_2}\n - {out_3}")
    return (out_1, out_2, out_3)

###### Part 3: Post-processing (capture kit intersect + VCF build) ######

def postprocess_build_vcfs(ref_choice_user: str):
    """
    For each N in {1,2,3}:
      - create BED with CHROM START END (drop header)
      - intersect with TCGA capture kit
      - filter original TXT to captured positions
      - make minimal VCF
    Uses singular 'caller' for N=1 and plural 'callers' for N=2,3.
    """
    _, overlap_ref = map_ref_dirs(ref_choice_user)
    work_dir = os.path.join(
        base, f"Variant_Caller_Overlap2/{overlap_ref}/Shared_Variants_Annotations/False_Positives"
    )
    if not os.path.isdir(work_dir):
        print(f" Missing directory: {work_dir}")
        return
    print(f"Post-processing in: {work_dir}")
    cwd = os.getcwd()
    try:
        os.chdir(work_dir)
        for N in (1, 2, 3):
            suffix = "caller" if N == 1 else "callers"
            stem = f"FP_shared_by_{N}_{suffix}"

            txt          = f"{stem}_noduplicates.txt"
            bed          = f"{stem}_noduplicates.bed"
            in_capture   = f"{stem}_in_capture.bed"
            filtered_pre = f"FP_filtered_shared_by_{N}_{suffix}_in_capture"
            filtered_txt = f"{filtered_pre}.txt"
            vcf          = f"{filtered_pre}.vcf"

            if not os.path.exists(txt):
                print(f" Skipping (not found): {txt}")
                continue

            # 1) TXT -> BED (drop header)
            run_command(f"tail -n +2 {txt} | awk -F'\\t' '{{print $1, $2, $3}}' OFS='\\t' > {bed}")

            # 2) Intersect with capture kit
            run_command(f"bedtools intersect -a {bed} -b {capture_bed} -u > {in_capture}")

            # 3) Filter original TXT rows to captured positions
            run_command(f"bedtools intersect -a {txt} -b {in_capture} -wa > {filtered_txt}")

            # 4) Build VCF
            run_command(f"echo -e '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO' > {vcf}")
            run_command(
                "awk -F'\\t' 'NR>1{printf \"%s\\t%s\\t.\\t%s\\t%s\\t.\\tPASS\\tSAMPLE=%s;TOOL=%s\\n\",$1,$2,$4,$5,$6,$7}' "
                f"{filtered_txt} >> {vcf}"
            )

            print(f" Built VCF: {vcf}")
        print()
    finally:
        os.chdir(cwd)

def main():
    print("Shared FP Comparison + Mega-list Annotation Builder")

    while True:
        ref_choice = input("Which alignment would you like to use? (GRCh38 Linear/T2T-CHM13 Pangenome): ").strip()
        if ref_choice in ["GRCh38", "T2T-CHM13"]:
            break
        print("Invalid input. Please enter 'GRCh38' or 'T2T-CHM13'.")

    sample_list = input("Enter path to 5-column patient list: ").strip()
    if not os.path.exists(sample_list):
        print(f"File not found: {sample_list}")
        return

    totals = dict(strelka_mutect=0, strelka_somatic=0, somatic_mutect=0, all_three=0)
    summaries = []

    with open(sample_list) as f:
        for line in f:
            if not line.strip(): continue
            cols = line.strip().split()
            if len(cols) < 5: continue
            donor_id, ancestry = cols[0], cols[3]
            res = compare_shared_false_positives(donor_id, ancestry, ref_choice)
            summaries.append(res["summary_string"])
            for k in totals: totals[k] += res[k]

    print("\n Summary per sample:")
    for s in summaries: print(s)

    print("\n Grand Totals:")
    for k, v in totals.items(): print(f" {k}: {v}")

    print("\nBuilding annotation seed files")
    build_fp_annotations_from_mega(ref_choice)

    print("\nCapture-kit intersect + VCF build")
    postprocess_build_vcfs(ref_choice)

    print("\n Done.")

if __name__ == "__main__":
    main()
