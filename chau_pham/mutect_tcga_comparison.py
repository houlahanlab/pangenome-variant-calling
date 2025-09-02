
# Compares Mutect2 Tumour-Normal pair called SNVs to TCGA Gold Standard SNVs calls
# Compares just CHROM, START, and END positions
# Uses already lifted over (T2T-CHM13 -> GRCh38) vcf 
# filters mutect VCF and keeps only SNVs
# Thus all Stats Analysis and comparing is done in GRCh38 coordinates (Pangenome, Linear, and TCGA)
# Outputs the Stats Analysis (Recall, Precision, F1 Score) and raw data into text

import os
import subprocess

def run_command(cmd):
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    print("Command completed successfully!\n")

def count_lines(path):
    if not os.path.exists(path):
        return 0
    with open(path) as f:
        return sum(1 for ln in f if ln.strip())

def safe_div(n, d): return n / d if d else 0.0
def f1_score(p, r): return (2*p*r/(p+r)) if (p+r) else 0.0

def _ref_label_and_suffix(ref_choice: str):
    if ref_choice == "T2T-CHM13":
        return "Pangenome", "pangenome"
    else:  # "GRCh38"
        return "Linear", "linear"

# Mega stats tables
def write_mega_table_row(row, ref_choice):
    """Append a TSV row to Stats_Summary mega_table_<suffix>.tsv with % metrics."""
    label, suffix = _ref_label_and_suffix(ref_choice)
    stats_dir = "/workspace/lab/houlahanlab/phamc20/data/Stats_Summary/Mutect2"
    os.makedirs(stats_dir, exist_ok=True)
    mega_path = os.path.join(stats_dir, f"stats_table_{suffix}.txt")

    header = [
        "Donor_ID", "Ancestry", "Variant_Caller", "Reference",
        "True_Positive", "False_Positive", "False_Negative",
        "Total_Variant_Called", "Total_TCGA",
        "Recall(%)", "Precision(%)", "F1_Score(%)"
    ]
    write_header = (not os.path.exists(mega_path)) or (os.path.getsize(mega_path) == 0)
    with open(mega_path, "a") as f:
        if write_header:
            f.write("\t".join(header) + "\n")
        f.write("\t".join([
            row["Donor_ID"],
            row["Ancestry"],         
            row["Variant_Caller"],   
            label,                    # Reference label ("Pangenome"/"Linear")
            str(row["True_Positive"]),
            str(row["False_Positive"]),
            str(row["False_Negative"]),
            str(row["Total_Variant_Called"]),
            str(row["Total_TCGA"]),
            f'{row["Recall"]*100:.2f}',
            f'{row["Precision"]*100:.2f}',
            f'{row["F1_Score"]*100:.2f}',
        ]) + "\n")
    print(f"Appended metrics to: {mega_path}")

def filter_snv_vcf(input_vcf, output_vcf):
        with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
            for line in infile:
                if line.startswith('#'):
                    outfile.write(line)
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 5:
                    continue
                ref = fields[3]
                alt = fields[4]
                if len(ref) == 1 and len(alt) == 1:
                    outfile.write(line)


def filter_vcf(donor_id, normal_id, ancestry, ref_choice):
    print(f"\nStarting VCF filtering for donor: {donor_id}")

    base_path = f"/home/phamc20/houlahanlab/data/data/{ancestry}/{donor_id}/Variant-Calling"
    second_base_path = f"/home/phamc20/houlahanlab/data/data/{ancestry}/{donor_id}/Variant-Calling/Original-Bam-Files"

    if ref_choice == "T2T-CHM13":
        input_vcf = os.path.join(base_path, f"liftover_{normal_id}_Mutect_T2T-CHM13_passfilter_exonregion.vcf")
        output_vcf = os.path.join(base_path, f"liftover_{normal_id}_Mutect_T2T-CHM13_passfilter_exonregion_snvs.vcf")
    elif ref_choice == "GRCh38":
        input_vcf = os.path.join(second_base_path, f"{normal_id}_Mutect2_GRCh38_passfilter_exonregions.vcf")
        output_vcf = os.path.join(second_base_path, f"{normal_id}_Mutect2_GRCh38_passfilter_exonregions_snvs.vcf")
    else:
        raise ValueError(f"Unsupported reference choice: {ref_choice}")

    if not os.path.exists(input_vcf):
        print(f"Input VCF not found: {input_vcf}")
        return None

    filter_snv_vcf(input_vcf, output_vcf)
    return output_vcf


def process_patient(donor_id, sample_id, normal_id, ancestry, ref_choice, tcga_vcf, liftover_tool, chain_file):
    print(f"\n Starting processing for Donor: {donor_id}")

    # First filter the VCF
    filtered_vcf_path = filter_vcf(donor_id, normal_id, ancestry, ref_choice)
    if filtered_vcf_path is None:
        print(f" Skipping {donor_id} due to missing files.")
        return ""

    if ref_choice == "GRCh38":
        output_dir = f"/home/phamc20/houlahanlab/phamc20/data/TCGA_Comparison/mutect_vs_tcga/GRCh38/{ancestry}/{donor_id}"
    else:
        output_dir = f"/home/phamc20/houlahanlab/phamc20/data/TCGA_Comparison/mutect_vs_tcga/T2T-CHM13/{ancestry}/{donor_id}"

    os.makedirs(output_dir, exist_ok=True)

    # Output files
    mutect_bed = f"{output_dir}/mutect_patient.bed"
    mutect_liftover_bed = f"{output_dir}/mutect_grch38.bed"
    mutect_grch38_positions = f"{output_dir}/mutect_grch38_positions.txt"
    tcga_positions = f"{output_dir}/tcga_patient_positions.txt"
    common_variants_full = f"{output_dir}/common_variants_full.txt"

    run_command(f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {filtered_vcf_path} | "
                f"awk '{{print $1 \"\\t\" ($2-1) \"\\t\" $2 \"\\t\" $3 \"\\t\" $4}}' > {mutect_bed}")

    # 1) Mutect VCF -> BED (0-based)
    if ref_choice == "T2T-CHM13":
        run_command(f"sed 's/CHM13#0#chr//' {mutect_bed} > {mutect_grch38_positions}")

    # 2) Just PASS sites (not all heterozygous sites like Strelka)
    run_command(f"""bcftools view -s {sample_id} {tcga_vcf} | \
                  grep -v '\\./\\.:.:.' | \
                  awk '$7 == "PASS" && $10 ~ /0\\/1/ {{print $1 \"\\t\" ($2-1) \"\\t\" $2 \"\\t\" $4 \"\\t\" $5}}' > {tcga_positions}""") #forgot to use the PASS variants only

    # 4) Chomosome naming
    if ref_choice == "GRCh38":
        run_command(f"cat {mutect_bed} | sed -E 's/^chr([0-9XYM]+)/\\1/' > {mutect_grch38_positions}")

    # 6) TCGA intersections
    run_command(f"grep -F -f {mutect_grch38_positions} {tcga_positions} > {common_variants_full} || true")

    # 7) Metrices based on PASS
    tp = count_lines(common_variants_full)
    mutect_pass_cnt = count_lines(mutect_grch38_positions)
    tcga_pass_cnt = count_lines(tcga_positions)
    fp = max(mutect_pass_cnt - tp, 0)
    fn = max(tcga_pass_cnt - tp, 0)

    total_called = tp + fp
    total_tcga   = tp + fn
    recall   = safe_div(tp, total_tcga)
    precision= safe_div(tp, total_called)
    f1       = f1_score(precision, recall)

    # Write megatable columns
    write_mega_table_row({
        "Donor_ID": donor_id,
        "Ancestry": ancestry,      
        "Variant_Caller": "Mutect2",
        "Reference": ref_choice,
        "True_Positive": tp,
        "False_Positive": fp,
        "False_Negative": fn,
        "Total_Variant_Called": total_called,
        "Total_TCGA": total_tcga,
        "Recall": recall,
        "Precision": precision,
        "F1_Score": f1,
    }, ref_choice)

    # Return summary string
    return (f"{donor_id} [{ancestry}] {ref_choice} — TP={tp}, FP={fp}, FN={fn}, "
            f"Prec={precision:.3f}, Rec={recall:.3f}, F1={f1:.3f}  → {output_dir}")

# Main function for inputs
def main():
    print("Welcome to the Mutect2 Filtering + Statistic Analysis Pipeline!")

    # Reference choice
    while True:
        ref_choice = input("Which alignment would you like to use? (GRCh38/T2T-CHM13): ").strip()
        if ref_choice in ["GRCh38", "T2T-CHM13"]:
            break
        print(" Invalid choice. Please enter 'GRCh38' or 'T2T-CHM13'.")

    tcga_vcf = "/home/phamc20/houlahanlab/abdelf6/raw_data/TCGA_Mutations/liftover_mc3.v0.2.8.PUBLIC.vcf"
    liftover_tool = "/workspace/lab/houlahanlab/software/liftOver"
    chain_file = "/home/phamc20/houlahanlab/resources/chain_files/chm13v2-grch38_fixed.chain"

    while True:
        mode = input("Would you like to manually process patients (1) or enter a list of patients (2)? ").strip()
        if mode in ["1", "2"]:
            break
        print(" Invalid choice. Please enter '1' or '2'.")

    all_summaries = []

    # Mode choice
    if mode == "1":
        donor_id = input("Enter Donor ID (e.g., TCGA-BT-A20V): ").strip()
        sample_id = input("Enter Tumour Sample ID (e.g., TCGA-BT-A20V-01A-11D-A14W-08): ").strip()
        normal_id = input("Enter Normal Sample ID (e.g., TCGA-BT-A20V-11A-11D-A14W-08): ").strip()
        ancestry = input("Enter Patient Ancestry (e.g., AA_Patients): ").strip()
        summary = process_patient(donor_id, sample_id, normal_id, ancestry, ref_choice, tcga_vcf, liftover_tool, chain_file)
        all_summaries.append(summary)
    else:
        sample_file = input("Enter path to patient list file: ").strip()

        if not os.path.exists(sample_file):
            print(f" Error: File {sample_file} not found!")
            return

        with open(sample_file, "r") as f:
            patient_entries = [line.strip().split() for line in f if line.strip()]

        for entry in patient_entries:
            if len(entry) < 5:
                continue

            donor_id, sample_id, normal_id, ancestry, _ = entry
            summary = process_patient(donor_id, sample_id, normal_id, ancestry, ref_choice, tcga_vcf, liftover_tool, chain_file)
            all_summaries.append(summary)

    if all_summaries:
        print("\n Final Summary for All Patients: ")
        for summary in all_summaries:
            print(summary)

    print("\n Pipeline Completed!")

if __name__ == "__main__":
    main()