# SomaticSniper Pipeline
# Runs SomaticSniper on either the pangenome aligned to T2T-CHM13 or the linear reference 
# aligned to GRCh38
# Filters the SNVs results using perl scripts
# Filters to just TCGA exons and DP of >8 for the normal, >=5 for the tumour

import os
import subprocess

def run_command(cmd):
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    print("Command completed successfully!\n")

def run_somaticsniper(donor_id, tumour_sample_id, normal_sample_id, ancestry, date, ref_choice, reference_fasta, base_dir_template):
    """Processes a single patient sample pair"""
    print(f"\nProcessing: Donor: {donor_id} | Tumour: {tumour_sample_id} | Normal: {normal_sample_id} | Ancestry: {ancestry} | Date: {date}")

    base_dir = base_dir_template.format(ancestry=ancestry, donor_id=donor_id)
    
    # GRCh38 = linear
    if ref_choice == "GRCh38":
        tumour_bam = os.path.join(base_dir, f"{tumour_sample_id}.bam")
        normal_bam = os.path.join(base_dir, f"{normal_sample_id}.bam")
    else:  # T2T-CHM13 = pangenome
        bam_suffix = "_fragment_length_150_CHM13v1.1_RG_markduplicates_RG_recalibrated.bam"
        tumour_bam = os.path.join(base_dir, f"{date}_{tumour_sample_id}{bam_suffix}")
        normal_bam = os.path.join(base_dir, f"{date}_{normal_sample_id}{bam_suffix}")

    # Ensure BAM files exist
    missing_files = [bam for bam in [tumour_bam, normal_bam] if not os.path.exists(bam)]
    if missing_files:
        print(" Error: The required BAM files are missing:")
        for file in missing_files:
            print(f"   - {file}")
        print("Skipping this sample pair!")
        return

    # Output directory
    sniper_output_dir = f"/workspace/lab/houlahanlab/phamc20/data/Variant_Callers_Results/SomaticSniper/{ancestry}/{donor_id}/{ref_choice}"
    os.makedirs(sniper_output_dir, exist_ok=True)

    print(f" Output Directory: {sniper_output_dir}")

    # Run SomaticSniper
    sniper_output_vcf = os.path.join(sniper_output_dir, "snp_output_file.vcf") 

    run_command(
    f"/workspace/lab/houlahanlab/software/somatic-sniper/build/bin/bam-somaticsniper "
    f"-Q 40 -G -L -f {reference_fasta} {tumour_bam} {normal_bam} -F vcf {sniper_output_vcf}"
    )

    print("Somatic Sniper Completed Successfully! All output files are stored in:")
    print(f"{sniper_output_vcf}")
    return sniper_output_vcf, sniper_output_dir

# Perl script filters
def run_sniper_filters(donor_id, sample_id, ancestry, date, ref_choice, reference_fasta, sniper_output_dir):
    """Apply NM, mpileup, snpfilter, prepare_for_readcount, bam-readcount, fpfilter, highconfidence"""
    scripts_path = "/workspace/lab/houlahanlab/software/somatic-sniper/src/scripts"

    # Locate input BAM
    if ref_choice == "GRCh38":
        bam_dir = f"/workspace/lab/houlahanlab/data/data/{ancestry}/{donor_id}/GRCh38_linear"
        input_bam = os.path.join(bam_dir, f"{sample_id}.bam")
        nm_bam    = os.path.join(bam_dir, f"{sample_id}.NM.bam")
    else:  # T2T-CHM13
        bam_dir = f"/workspace/lab/houlahanlab/data/data/{ancestry}/{donor_id}/T2T-CHM13_linear"
        input_bam = os.path.join(bam_dir, f"{date}_{sample_id}_fragment_length_150_CHM13v1.1_RG_markduplicates_RG_recalibrated.bam")
        nm_bam    = os.path.join(bam_dir, f"{date}_{sample_id}_fragment_length_150_CHM13v1.1_RG_markduplicates_RG_recalibrated_NM.bam")

    if not os.path.exists(input_bam):
        print(f" Error: BAM not found for NM tagging: {input_bam}")
        return None

    pileup_file = os.path.join(sniper_output_dir, "tumour_NM.pileup")
    somatic_vcf = os.path.join(sniper_output_dir, "snp_output_file.vcf")
    filtered_vcf = somatic_vcf + ".SNPfilter"
    readcounts_file = os.path.join(sniper_output_dir, "readcounts.rc")
    fp_filtered_vcf = somatic_vcf + ".SNPfilter.fp_pass"
    high_conf_vcf = somatic_vcf + ".SNPfilter.fp_pass.hc"

    # Step 1: NM tag
    run_command(f"samtools calmd -b {input_bam} {reference_fasta} > {nm_bam}")
    # Step 2: index
    run_command(f"samtools index {nm_bam}")
    # Step 3: mpileup (NM BAM)
    run_command(f"samtools mpileup -f {reference_fasta} {nm_bam} > {pileup_file}")
    # Step 4: SNP filter
    run_command(f"perl {scripts_path}/snpfilter.pl --snp-file {somatic_vcf} --indel-file {pileup_file}")
    # Step 5: prepare_for_readcount
    run_command(f"perl {scripts_path}/prepare_for_readcount.pl --snp-file {filtered_vcf}")
    # Step 6: bam-readcount
    readcount_positions = f"{filtered_vcf}.pos"
    run_command(f"/workspace/lab/houlahanlab/software/micromamba/bin/bam-readcount -b 15 -f {reference_fasta} -l {readcount_positions} {nm_bam} > {readcounts_file}")
    # Step 7: fpfilter
    run_command(f"perl {scripts_path}/fpfilter.pl --snp-file {filtered_vcf} --readcount-file {readcounts_file}")
    # Step 8: highconfidence
    run_command(f"perl {scripts_path}/highconfidence.pl --snp-file {fp_filtered_vcf}")

    print(f"\n Pipeline Completed Successfully for {donor_id} ({ref_choice}). High-confidence VCF: {high_conf_vcf}")
    return high_conf_vcf


# Exon filtering/read count filter
def run_exon_dp_filters(patient_id, ancestry, reference, high_conf_vcf, out_dir):
    """Intersect HC VCF with Agilent exons, then DP filters"""

    # Choose BED file depending on reference
    if reference == "T2T-CHM13":
        bed_file = "/workspace/lab/houlahanlab/resources/agilent_sureselect_exon/Agilent_SureSelect_T2T_sorted.bed.gz"
    else: # GRCh38
        bed_file = "/workspace/lab/houlahanlab/resources/agilent_sureselect_exon/Agilent_SureSelect_All_Exon_V5.bed.gz"

    # Output files
    filtered_vcf = os.path.join(out_dir, f"{patient_id}_in_capture_exons.vcf")
    normal_txt = os.path.join(out_dir, "normal_exons.txt")
    filtered_normal = os.path.join(out_dir, "filtered_normal_exons.txt")
    filtered_normal_tumour = os.path.join(out_dir, "filtered_normal_tumour_exons.txt")

    # Bedtools intersect VCF with Agilent
    run_command(f"bedtools intersect -header -a {high_conf_vcf} -b {bed_file} > {filtered_vcf}")
    print(f"Bedtools exon filtering done: {filtered_vcf}")

    # Step 1: bcftools query
    run_command(f"bcftools query -f '%CHROM\\t%POS\\t[%DP\\t]\\n' {filtered_vcf} > {normal_txt}")

    # Step 2: awk '$3 > 8'
    run_command(f"awk '$3 > 8 {{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4}}' {normal_txt} > {filtered_normal}")

    # Step 3: awk '$4 >= 5'
    run_command(f"awk '$4 >= 5 {{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4}}' {filtered_normal} > {filtered_normal_tumour}")


def main():
    print("Welcome to the Somatic Sniper Variant Calling Pipeline!")

    # Ask for alignment used
    while True:
        ref_choice = input("Which alignment would you like to use? (GRCh38/T2T-CHM13): ").strip()
        if ref_choice in ["GRCh38", "T2T-CHM13"]:
            break
        print(" Invalid choice. Please enter 'GRCh38' or 'T2T-CHM13'.")

    # Reference genome paths + base dirs
    if ref_choice == "GRCh38":
        reference_fasta = "/workspace/lab/houlahanlab/resources/GRCh38_Indexes/Reference-Fasta-Files/GRCh38.d1.vd1.fa"
        base_dir_template = "/workspace/lab/houlahanlab/data/data/{ancestry}/{donor_id}/GRCh38_linear"
    else:
        reference_fasta = "/workspace/lab/houlahanlab/resources/T2T-CHM13_Indexes/version_1.1/Reference-Fasta-Files/CHM13v11Y_test_v1.1.fa"
        base_dir_template = "/workspace/lab/houlahanlab/data/data/{ancestry}/{donor_id}/T2T-CHM13_linear"

    # Simple processor for one patient
    def process_one(donor_id, tumour_sample_id, normal_sample_id, ancestry, date):
        vcf_path, out_dir = run_somaticsniper(
            donor_id, tumour_sample_id, normal_sample_id,
            ancestry, date, ref_choice, reference_fasta, base_dir_template
        )
        if not vcf_path:
            print(" Skipping due to missing BAMs.")
            return
        hc_vcf = run_sniper_filters(
            donor_id, tumour_sample_id, ancestry, date,
            ref_choice, reference_fasta, out_dir
        )
        if not hc_vcf:
            print(" Skipping due to filtering error.")
            return
        run_exon_dp_filters(donor_id, ancestry, ref_choice, hc_vcf, out_dir)

    # Manual vs list
    while True:
        mode = input("Would you like to manually process patients (1) or enter a list of patients (2)? ").strip()
        if mode in ["1", "2"]:
            break
        print(" Invalid choice. Please enter '1' or '2'.")

    if mode == "1":
        donor_id = input("Enter Donor ID: (e.g TCGA-BT-A20V): ").strip()
        tumour_sample_id = input("Enter Tumour Sample ID (e.g TCGA-BT-A20V-01A-11D-A14W-08): ").strip()
        normal_sample_id = input("Enter Normal Sample ID (e.g TCGA-BT-A20V-11A-11D-A14W-08): ").strip()
        ancestry = input("Enter Patient Ancestry (e.g AA_Patients): ").strip()
        date = input("Enter BAM File Creation Date (YYYY_MM_DD) Leave blank for GRCh38: ").strip()
        process_one(donor_id, tumour_sample_id, normal_sample_id, ancestry, date)
    else:
        sample_file = input("Enter path to sample list file (DonorID TumourID NormalID Ancestry Date): ").strip()
        if not os.path.exists(sample_file):
            print(f" Sample list file {sample_file} not found!")
            return
        with open(sample_file, "r") as f:
            for line in f:
                if not line.strip():
                    continue
                parts = line.strip().split()
                if len(parts) < 5:
                    continue
                donor_id, tumour_sample_id, normal_sample_id, ancestry, date = parts[:5]
                process_one(donor_id, tumour_sample_id, normal_sample_id, ancestry, date)

    print("\n SomaticSniper Completed Successfully!")

if __name__ == "__main__":
    main()

