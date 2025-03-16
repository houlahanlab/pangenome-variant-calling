import os
import subprocess

def run_command(command):
    """Runs a shell command and prints status."""
    print(f"Running: {' '.join(command)}")
    subprocess.run(command, check=True)
    print("✅ Command completed successfully!\n")

def main():
    print("✨💻💖 Welcome to the Strelka Somatic Variant Calling Pipeline!🔬💖✨")

    # Get user inputs
    donor_id = input("Enter Donor ID: ")  # e.g., TCGA-BT-A20V
    tumor_sample_id = input("Enter Tumor Sample ID: ")  # e.g., TCGA-BT-A20V-01A-11D-A14W-08
    normal_sample_id = input("Enter Normal Sample ID: ")  # e.g., TCGA-BT-A20V-11A-11D-A14W-08
    ancestry = input("Enter Patient Ancestry: ")  # e.g., AA_Patients
    date = input("Enter BAM File Creation Date (YYYY_MM_DD) (Leave blank if using GRCh38): ").strip()  

    # Linear or Pangenome reference
    # GRCh38 in this case is the linear reference and T2T-CHM13 is the pangenome which later gets lifted over to GRCh38
    while True:
        ref_choice = input("Which reference genome would you like to use? (GRCh38/T2T-CHM13): ").strip()
        if ref_choice in ["GRCh38", "T2T-CHM13"]:
            break
        print("❌ Invalid choice. Please enter 'GRCh38' or 'T2T-CHM13'.")

    # GRCh38
    if ref_choice == "GRCh38":
        reference_fasta = "/workspace/lab/houlahanlab/resources/GRCh38_Indexes/Reference-Fasta-Files/GRCh38.d1.vd1.fa"
        call_regions = "/workspace/lab/houlahanlab/resources/agilent_sureselect_exon/Agilent_SureSelect_All_Exon_V5.bed.gz"
        base_dir = f"/home/phamc20/houlahanlab/abdelf6/raw_data/{ancestry}/{donor_id}/GRCh38_linear"

        # GRCh38 BAM files do not include a date in their filename
        tumor_bam = os.path.join(base_dir,f"{tumor_sample_id}.bam")
        normal_bam = os.path.join(base_dir,f"{normal_sample_id}.bam")

    else:  # T2T-CHM13
        reference_fasta = "/workspace/lab/houlahanlab/resources/T2T-CHM13_Indexes/version_1.1/Reference-Fasta-Files/CHM13v11Y_test_v1.1.fa"
        call_regions = "/home/phamc20/houlahanlab/resources/agilent_sureselect_exon/Agilent_SureSelect_T2T_sorted.bed.gz"
        base_dir = f"/home/phamc20/houlahanlab/data/data/{ancestry}/{donor_id}/T2T-CHM13_linear"

        # T2T-CHM13 BAM files include the date in their filename
        bam_suffix = "_fragment_length_150_CHM13v1.1_RG_markduplicates_RG_recalibrated.bam"
        tumor_bam = os.path.join(base_dir, f"{date}_{tumor_sample_id}{bam_suffix}")
        normal_bam = os.path.join(base_dir, f"{date}_{normal_sample_id}{bam_suffix}")

    # Make sure .bai files are present
    # For GRCh38, some of the bam files are not index (have .bai) so I manually indexed them?
    # Maybe code to index them if missing?
    tumor_bai = tumor_bam +".bai"
    normal_bai = normal_bam +".bai"

    # Output directory 
    # Maybe rename instead of it saying GRCh38 it'll say linear?
    strelka_output_dir = f"/workspace/lab/houlahanlab/phamc20/data/strelka/{ancestry}/{donor_id}/{ref_choice}_{donor_id}"

    # Verify if required files exist before running Strelka
    missing_files = []
    for file_path in [tumor_bam, normal_bam, tumor_bai, normal_bai, reference_fasta, call_regions]:
        if not os.path.exists(file_path):
            missing_files.append(file_path)

    if missing_files:
        print("❌ Error: The following required files are missing:")
        for file in missing_files:
            print(f"   - {file}")
        print("⚠️ Ensure BAM and BAI files are correctly indexed before running Strelka.")
        return 

    print(f"Detected Tumor BAM: {tumor_bam}")
    print(f"Detected Normal BAM: {normal_bam}")
    print(f"Output Directory: {strelka_output_dir}")

    # Run Strelka
    # This can also just be ran in the linux terminal
    configure_command = [
        "python2", "/workspace/lab/houlahanlab/software/strelka/bin/configureStrelkaSomaticWorkflow.py",
        "--tumorBam", tumor_bam,
        "--normalBam", normal_bam,
        "--referenceFasta", reference_fasta,
        "--callRegions", call_regions,
        "--exome",
        "--runDir", strelka_output_dir
    ]
    
    run_command(configure_command)

    run_workflow_command = [
        "python2", f"{strelka_output_dir}/runWorkflow.py",
        "-m", "local",
        "-j", "8"
    ]

    run_command(run_workflow_command)

    print("🎉✅ Strelka Completed Successfully! All output files are stored in:")
    print(f"{strelka_output_dir}")

if __name__ == "__main__":
    main()

