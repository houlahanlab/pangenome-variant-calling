# GOOD ONE
# Need to have extra filtering to stop over calling variants 

import os
import subprocess

def run_command(command):
    """Runs a shell command and prints status."""
    print(f"Running: {' '.join(command)}")
    subprocess.run(command, check=True)
    print("✅ Command completed successfully!\n")

def process_sample(donor_id, tumor_sample_id, normal_sample_id, ancestry, date, ref_choice, reference_fasta, base_dir_template):
    """Processes a single patient sample pair"""
    print(f"\n🧬 Processing: Donor: {donor_id} | Tumor: {tumor_sample_id} | Normal: {normal_sample_id} | Ancestry: {ancestry} | Date: {date}")

    base_dir = base_dir_template.format(ancestry=ancestry, donor_id=donor_id)
    
    # GRCh38 = linear
    if ref_choice == "GRCh38":
        tumor_bam = os.path.join(base_dir, f"{tumor_sample_id}.bam")
        normal_bam = os.path.join(base_dir, f"{normal_sample_id}.bam")
    else:  # T2T-CHM13 = pangenome
        bam_suffix = "_fragment_length_150_CHM13v1.1_RG_markduplicates_RG_recalibrated.bam"
        tumor_bam = os.path.join(base_dir, f"{date}_{tumor_sample_id}{bam_suffix}")
        normal_bam = os.path.join(base_dir, f"{date}_{normal_sample_id}{bam_suffix}")

    # Ensure BAM files exist
    missing_files = [bam for bam in [tumor_bam, normal_bam] if not os.path.exists(bam)]
    if missing_files:
        print("❌ Error: The following required BAM files are missing:")
        for file in missing_files:
            print(f"   - {file}")
        print("⚠️ Skipping this sample pair!")
        return

    # Output directory
    sniper_output_dir = f"/workspace/lab/houlahanlab/phamc20/data/somatic_sniper/{ancestry}/{donor_id}/{ref_choice}"
    os.makedirs(sniper_output_dir, exist_ok=True)

    print(f"📂 Output Directory: {sniper_output_dir}")

    # Run SomaticSniper
    sniper_output_vcf = os.path.join(sniper_output_dir, "snp_output_file.vcf")

    somatic_sniper_command = [
        "/home/phamc20/houlahanlab/software/somatic-sniper/build/bin/bam-somaticsniper",
        "-Q", "40",  # Minimum base quality
        "-G",  
        "-L",  # Log scale for quality scores
        "-f", reference_fasta, tumor_bam, normal_bam,
        "-F", "vcf", sniper_output_vcf # Output format
    ]

    run_command(somatic_sniper_command)


def main():
    print("✨💻💖 Welcome to the Somatic Sniper Variant Calling Pipeline! 🔬💖✨")

    # Ask for reference genome choice
    while True:
        ref_choice = input("Which reference genome would you like to use? (GRCh38/T2T-CHM13): ").strip()
        if ref_choice in ["GRCh38", "T2T-CHM13"]:
            break
        print("❌ Invalid choice. Please enter 'GRCh38' or 'T2T-CHM13'.")

    # Reference genome paths
    if ref_choice == "GRCh38":
        reference_fasta = "/workspace/lab/houlahanlab/resources/GRCh38_Indexes/Reference-Fasta-Files/GRCh38.d1.vd1.fa"
        base_dir_template = "/home/phamc20/houlahanlab/abdelf6/raw_data/{ancestry}/{donor_id}/GRCh38_linear"
    else:  # T2T-CHM13
        reference_fasta = "/workspace/lab/houlahanlab/resources/T2T-CHM13_Indexes/version_1.1/Reference-Fasta-Files/CHM13v11Y_test_v1.1.fa"
        base_dir_template = "/home/phamc20/houlahanlab/data/data/{ancestry}/{donor_id}/T2T-CHM13_linear"

    # Manually input donors or use a list
    # So far i've been making the list within the ~/houlahanlab/phamc20/scripts directory using nano
    while True:
        mode = input("Would you like to manually process patients (1) or enter a list of patients (2)? ").strip()
        if mode in ["1", "2"]:
            break
        print("❌ Invalid choice. Please enter '1' or '2'.")

    if mode == "1":
        # Manually enters patient details
        donor_id = input("Enter Donor ID: (e.g TCGA-BT-A20V): ").strip()
        tumor_sample_id = input("Enter Tumor Sample ID (e.g TCGA-BT-A20V-01A-11D-A14W-08): ").strip()
        normal_sample_id = input("Enter Normal Sample ID (e.g TCGA-BT-A20V-11A-11D-A14W-08): ").strip()
        ancestry = input("Enter Patient Ancestry (e.g AA_Patients): ").strip()
        date = input("Enter BAM File Creation Date (YYYY_MM_DD) Leave blank for GRCh38: ").strip()

        process_sample(donor_id, tumor_sample_id, normal_sample_id, ancestry, date, ref_choice, reference_fasta, base_dir_template)

    else:
        # Providing a sample list file
        sample_file = input("Enter path to sample list file (DonorID, TumorID, NormalID, Ancestry, Date): ").strip()
        
        if not os.path.exists(sample_file):
            print(f"❌ Error: Sample list file {sample_file} not found!")
            return
        
        with open(sample_file, "r") as f:
            sample_entries = [line.strip().split() for line in f.readlines()]

        for donor_id, tumor_sample_id, normal_sample_id, ancestry, date in sample_entries:
            process_sample(donor_id, tumor_sample_id, normal_sample_id, ancestry, date, ref_choice, reference_fasta, base_dir_template)

    print("\n🎉✅ SomaticSniper Completed Successfully!")

if __name__ == "__main__":
    main()