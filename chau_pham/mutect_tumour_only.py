# Mutect2 Tumour-Only Pipeline
# Runs Mutect2 on either the pangenome aligned to T2T-CHM13 or the linear reference 
# aligned to GRCh38
# Also filters Mutect2 results to just SNVs and PASS

import os
import subprocess

def run_command(cmd):
    """Runs a shell command and prints status."""
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    print(" Command completed successfully!\n")

def get_pass_count(vcf_file):
    result = subprocess.run(f"bcftools view {vcf_file} | grep 'PASS' | wc -l", shell=True, capture_output=True, text=True)
    return int(result.stdout.strip())

def process_sample(donor_id, tumor_sample_id, normal_sample_id, ancestry, date, ref_choice, reference_fasta, base_dir_template):
    # format base directory using ancestry and donor_id
    base_dir = base_dir_template.format(ancestry=ancestry, donor_id=donor_id)

    # set bam path depending on reference
    if ref_choice == "GRCh38":
        tumor_bam = os.path.join(base_dir, f"{tumor_sample_id}.bam")
        pon = "/workspace/lab/houlahanlab/phamc20/resources/PON/1000g_pon.hg38.vcf.gz"
        germline_resource = "/workspace/lab/houlahanlab/resources/germline_resource/af-only-gnomad.hg38.vcf.gz"
        agilent_bed = "/workspace/lab/houlahanlab/resources/agilent_sureselect_exon/Agilent_SureSelect_All_Exon_V5.bed.gz"
    else:
        tumor_bam = os.path.join(base_dir, f"{date}_{tumor_sample_id}_fragment_length_150_CHM13v1.1_RG_markduplicates_RG_recalibrated.bam")
        pon = "/workspace/lab/houlahanlab/phamc20/resources/PON/1000g_pon.t2t.renamed.sorted.vcf.gz"
        germline_resource = "/workspace/lab/houlahanlab/resources/germline_resource/af-only-gnomad.t2t.renamed.vcf.gz"
        agilent_bed = "/workspace/lab/houlahanlab/resources/agilent_sureselect_exon/Agilent_SureSelect_T2T_sorted.bed.gz"

    # make output directory
    output_dir = f"/home/phamc20/houlahanlab/phamc20/data/Mutect_Tumour_Only/{ref_choice}/{ancestry}/{donor_id}"
    os.makedirs(output_dir, exist_ok=True)

    # output files
    output_prefix = os.path.join(output_dir, f"{tumor_sample_id}")
    vcf_raw = f"{output_prefix}.vcf.gz"
    vcf_filtered = f"{output_prefix}_filtered.vcf.gz"
    vcf_agilent = f"{output_prefix}_filtered_agilent.vcf.gz"
    vcf_pass = f"{output_prefix}_filtered_agilent_PASS.vcf"
    vcf_dp8 = f"{output_prefix}_filtered_agilent_PASS_DP8.vcf.gz"

    # run Mutect2
    run_command(f"gatk Mutect2 -R {reference_fasta} -I {tumor_bam} --panel-of-normals {pon} --germline-resource {germline_resource} -O {vcf_raw}")

    # run FilterMutectCalls
    run_command(f"gatk FilterMutectCalls -R {reference_fasta} -V {vcf_raw} -O {vcf_filtered}")

    # filter with Agilent BED
    run_command(f"bcftools view -R {agilent_bed} {vcf_filtered} -Oz -o {vcf_agilent}")

    # extract only PASS variants
    run_command(f"bcftools view -f PASS {vcf_agilent} -Ov -o {vcf_pass}")

    # filter DP >= 8
    run_command(f"bcftools filter -i 'INFO/DP>=8' {vcf_pass} | bgzip -c > {vcf_dp8}")

    # count number of PASS variants
    pass_count = get_pass_count(vcf_dp8)
    return donor_id, ref_choice, pass_count

def main():
    print("Welcome to the Mutect Tumor-Only Variant Calling Pipeline!")

    summary_results = []

    # choose reference genome
    while True:
        ref_choice = input("Which reference genome would you like to use? (GRCh38/T2T-CHM13): ").strip()
        if ref_choice in ["GRCh38", "T2T-CHM13"]:
            break
        print("Invalid choice. Please enter 'GRCh38' or 'T2T-CHM13'.")

    # set reference fasta and base directory pattern
    if ref_choice == "GRCh38":
        reference_fasta = "/workspace/lab/houlahanlab/resources/GRCh38_Indexes/Reference-Fasta-Files/GRCh38.d1.vd1.fa"
        base_dir_template = "/workspace/lab/houlahanlab/abdelf6/raw_data/{ancestry}/{donor_id}/GRCh38_linear"
    else:
        reference_fasta = "/workspace/lab/houlahanlab/resources/T2T-CHM13_Indexes/version_1.1/Reference-Fasta-Files/CHM13v11Y_test_v1.1.fa"
        base_dir_template = "/workspace/lab/houlahanlab/data/data/{ancestry}/{donor_id}/T2T-CHM13_linear"

    # single or list mode
    while True:
        mode = input("Would you like to manually process patients (1) or enter a list of patients (2)? ").strip()
        if mode in ["1", "2"]:
            break
        print("Invalid choice. Please enter '1' or '2'.")

    if mode == "1":
        donor_id = input("Enter Donor ID: (e.g TCGA-BT-A20V): ").strip()
        tumor_sample_id = input("Enter Tumor Sample ID (e.g TCGA-BT-A20V-01A-11D-A14W-08): ").strip()
        normal_sample_id = input("Enter Normal Sample ID (e.g TCGA-BT-A20V-11A-11D-A14W-08): ").strip()
        ancestry = input("Enter Patient Ancestry (e.g AA_Patients): ").strip()
        date = input("Enter BAM File Creation Date (YYYY_MM_DD) Leave blank for GRCh38: ").strip()
        result = process_sample(donor_id, tumor_sample_id, normal_sample_id, ancestry, date, ref_choice, reference_fasta, base_dir_template)
        summary_results.append(result)
    else:
        sample_file = input("Enter path to sample list file (DonorID TumorID NormalID Ancestry Date): ").strip()
        if not os.path.exists(sample_file):
            print(f"Error: Sample list file {sample_file} not found!")
            return

        with open(sample_file, "r") as f:
            for line in f:
                if line.strip():
                    donor_id, tumor_sample_id, normal_sample_id, ancestry, date = line.strip().split()
                    result = process_sample(donor_id, tumor_sample_id, normal_sample_id, ancestry, date, ref_choice, reference_fasta, base_dir_template)
                    summary_results.append(result)

    print("\nMutect Tumor-Only Pipeline Completed!")
    print("\nSummary of PASS variants:")
    for donor_id, ref_used, count in summary_results:
        print(f"Donor: {donor_id} | Reference: {ref_used} | PASS Variants: {count}")

if __name__ == "__main__":
    main()


