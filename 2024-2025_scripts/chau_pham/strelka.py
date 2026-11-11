# Strelka Pipeline
# Runs Strelka on either the pangenome aligned to T2T-CHM13 or the linear reference 
# aligned to GRCh38

import os
import subprocess

def run_command(cmd):
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    print("Command completed successfully!\n")

def process_sample(donor_id, tumour_sample_id, normal_sample_id, ancestry, date, ref_choice, reference_fasta, base_dir_template):
    """Processes a single patient sample pair"""
    print(f"\nProcessing: Donor: {donor_id} | Tumour: {tumour_sample_id} | Normal: {normal_sample_id} | Ancestry: {ancestry} | Date: {date}")

    base_dir = base_dir_template.format(ancestry=ancestry, donor_id=donor_id)
    
    # GRCh38 = linear
    if ref_choice == "GRCh38":
        tumour_bam = os.path.join(base_dir, f"{tumour_sample_id}.bam")
        normal_bam = os.path.join(base_dir, f"{normal_sample_id}.bam")
        call_regions = "/workspace/lab/houlahanlab/resources/agilent_sureselect_exon/Agilent_SureSelect_All_Exon_V5.bed.gz"
    else:  # T2T-CHM13 = pangenome
        bam_suffix = "_fragment_length_150_CHM13v1.1_RG_markduplicates_RG_recalibrated.bam"
        tumour_bam = os.path.join(base_dir, f"{date}_{tumour_sample_id}{bam_suffix}")
        normal_bam = os.path.join(base_dir, f"{date}_{normal_sample_id}{bam_suffix}")
        call_regions = "/workspace/lab/houlahanlab/resources/agilent_sureselect_exon/Agilent_SureSelect_T2T_sorted.bed.gz"

    # Ensure BAM files exist
    missing_files = [bam for bam in [tumour_bam, normal_bam] if not os.path.exists(bam)]
    if missing_files:
        print(" Error: The required BAM files are missing:")
        for file in missing_files:
            print(f"   - {file}")
        print("Skipping this sample pair!")
        return
    
    tumour_bai = tumour_bam +".bai"
    normal_bai = normal_bam +".bai"

    # Output directory
    strelka_output_dir = f"/workspace/lab/houlahanlab/phamc20/data/Variant_Callers_Results/Strelka/{ancestry}/{donor_id}/{ref_choice}"
    os.makedirs(strelka_output_dir, exist_ok=True)

    missing_files = []
    for file_path in [tumour_bam, normal_bam, tumour_bai, normal_bai, reference_fasta, call_regions]:
        if not os.path.exists(file_path):
            missing_files.append(file_path)

    if missing_files:
        print("Error: The following required files are missing:")
        for file in missing_files:
            print(f"   - {file}")
        print("Ensure BAM and BAI files are correctly indexed before running Strelka.")
        return 

    print(f"Detected Tumour BAM: {tumour_bam}")
    print(f"Detected Normal BAM: {normal_bam}")
    print(f"Output Directory: {strelka_output_dir}")

    # Run Strelka
    # This can also just be ran in the linux terminal
    run_command(
    f"python2 /workspace/lab/houlahanlab/software/strelka/bin/configureStrelkaSomaticWorkflow.py "
    f"--tumorBam {tumour_bam} "
    f"--normalBam {normal_bam} "
    f"--referenceFasta {reference_fasta} "
    f"--callRegions {call_regions} "
    f"--exome "
    f"--runDir {strelka_output_dir}"
    )

    run_command(f"python2 {strelka_output_dir}/runWorkflow.py -m local -j 8")

    print("Strelka Completed Successfully! All output files are stored in:")
    print(f"{strelka_output_dir}")

#User inputs 
def main():
    print("Welcome to the Strelka Variant Calling Pipeline!")

    # Ask for reference genome choice
    while True:
        ref_choice = input("Which alignment would you like to use? (GRCh38/T2T-CHM13): ").strip()
        if ref_choice in ["GRCh38", "T2T-CHM13"]:
            break
        print(" Invalid choice. Please enter 'GRCh38' or 'T2T-CHM13'.")

    # Reference genome paths
    if ref_choice == "GRCh38":
        reference_fasta = "/workspace/lab/houlahanlab/resources/GRCh38_Indexes/Reference-Fasta-Files/GRCh38.d1.vd1.fa"
        base_dir_template = "/workspace/lab/houlahanlab/data/data/{ancestry}/{donor_id}/GRCh38_linear"
        #base_dir_template = "/workspace/lab/houlahanlab/abdelf6/raw_data/{ancestry}/{donor_id}/GRCh38_linear" old base path used for first 30 patients..

    else:  # T2T-CHM13
        reference_fasta = "/workspace/lab/houlahanlab/resources/T2T-CHM13_Indexes/version_1.1/Reference-Fasta-Files/CHM13v11Y_test_v1.1.fa"
        base_dir_template = "/workspace/lab/houlahanlab/data/data/{ancestry}/{donor_id}/T2T-CHM13_linear"

    # Manually input donors or use a list
    # So far i've been making the list within the ~/houlahanlab/phamc20/scripts directory using nano
    while True:
        mode = input("Would you like to manually process patients (1) or enter a list of patients (2)? ").strip()
        if mode in ["1", "2"]:
            break
        print(" Invalid choice. Please enter '1' or '2'.")

    if mode == "1":
        # Manually enters patient details
        donor_id = input("Enter Donor ID: (e.g TCGA-BT-A20V): ").strip()
        tumour_sample_id = input("Enter Tumour Sample ID (e.g TCGA-BT-A20V-01A-11D-A14W-08): ").strip()
        normal_sample_id = input("Enter Normal Sample ID (e.g TCGA-BT-A20V-11A-11D-A14W-08): ").strip()
        ancestry = input("Enter Patient Ancestry (e.g AA_Patients): ").strip()
        date = input("Enter BAM File Creation Date (YYYY_MM_DD) Leave blank for GRCh38: ").strip()

        process_sample(donor_id, tumour_sample_id, normal_sample_id, ancestry, date, ref_choice, reference_fasta, base_dir_template)

    else:
        # Providing a sample list file
        sample_file = input("Enter path to sample list file (DonorID, TumourID, NormalID, Ancestry, Date): ").strip()
        
        if not os.path.exists(sample_file):
            print(f"Sample list file {sample_file} not found!")
            return
        
        with open(sample_file, "r") as f:
            sample_entries = [line.strip().split() for line in f.readlines()]

        for donor_id, tumour_sample_id, normal_sample_id, ancestry, date in sample_entries:
            process_sample(donor_id, tumour_sample_id, normal_sample_id, ancestry, date, ref_choice, reference_fasta, base_dir_template)

    print("\nStrelka Completed Successfully!")

if __name__ == "__main__":
    main()