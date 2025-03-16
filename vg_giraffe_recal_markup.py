import os
import subprocess

def run_command(command, output_file=None):
    """Runs a shell command and logs output if provided."""
    print("Executing:", " ".join(command))
    
    if output_file:
        with open(output_file, "wb") as output:
            subprocess.run(command, check=True, stdout=output)
        print(f"Command completed successfully. Output saved to: {output_file}!")
    else:
        subprocess.run(command, check=True)

def run_vg_giraffe():
    print("✨💻💖 Running VG Giraffe Alignment! 💖🔬✨")
    
    # For user inputs
    donor_id = input("Enter Donor ID: ")       
    sample_id = input("Enter Sample ID: ")       
    ancestry = input("Enter Patient Ancestry: ") 
    date = input("Please enter today's date (YYYY_MM_DD): ")
    
    # Base directory for input FASTQ files
    input_base_dir = "/home/phamc20/houlahanlab/abdelf6/raw_data"
    
    # Base directory for output files
    output_base_dir = os.path.expanduser("~/houlahanlab/data/data")
    
    # Input directory for FASTQ files
    base_directory = os.path.join(input_base_dir, ancestry, donor_id, "GRCh38_linear")
    if not os.path.exists(base_directory):
        print(f"Error: Input directory {base_directory} does not exist.")
        return

    # Locate the FASTQ files
    fastq_files = [
        os.path.join(base_directory, f"{sample_id}_R1.fastq"),
        os.path.join(base_directory, f"{sample_id}_R2.fastq")
    ]
    for file in fastq_files:
        if not os.path.exists(file):
            print(f"Error: FASTQ file {file} not found.")
            return
    print(f"FASTQ files found: {fastq_files}!")
    
    # Reference genome choice
    references = {
        "GRCh38" : "/home/phamc20/houlahanlab/resources/GRCh38_Indexes/VG-Index-Files/hprc-v1.0-mc-grch38.gbz",
        "T2T-CHM13": "/home/phamc20/houlahanlab/resources/T2T-CHM13_Indexes/version_1.1/VG-Index-Files/hprc-v1.1-mc-chm13.d9.gbz"
    }

    dist_index = {
        "GRCh38" : "/home/phamc20/houlahanlab/resources/GRCh38_Indexes/VG-Index-Files/hprc-v1.0-mc-grch38.dist",
        "T2T-CHM13" : "/home/phamc20/houlahanlab/resources/T2T-CHM13_Indexes/version_1.1/VG-Index-Files/hprc-v1.1-mc-chm13.d9.dist"}
    
    minimizer_index = {
        "GRCh38" : "/home/phamc20/houlahanlab/resources/GRCh38_Indexes/VG-Index-Files/hprc-v1.0-mc-grch38.min", 
        "T2T-CHM13" : "/home/phamc20/houlahanlab/resources/T2T-CHM13_Indexes/version_1.1/VG-Index-Files/hprc-v1.1-mc-chm13.d9.min"}
    
    # Message 
    ref_choice = input("Please choose the reference genome (GRCh38 or T2T-CHM13): ")
    reference_genome = references.get(ref_choice)
    if not reference_genome:
        print("❌ Invalid reference genome choice.")
        return

    # Construct the output directory for pangenome alignment results using an absolute path
    # Figure out how to make it for relative path
    # For example: ~/houlahanlab/phamc20/data/AA_Patients/TCGA-BT-A20V/GRCh38_pangenome
    if ref_choice == "GRCh38":
        alignment_directory = os.path.join(output_base_dir, ancestry, donor_id, "GRCh38_pangenome")
    elif ref_choice == "T2T-CHM13":
        alignment_directory = os.path.join(output_base_dir, ancestry, donor_id, "T2T-CHM13_pangenome")
    else:
        print("❌ Invalid reference genome choice.")
        return
    os.makedirs(alignment_directory, exist_ok=True)

    # Output GAM file name
    output_gam = os.path.join(alignment_directory, f"{date}_{sample_id}_{ref_choice}_pangenome_fragment_length-150_aligned.gam")
    
    # VG Giraffe alignment
    vg_executable = "/workspace/lab/houlahanlab/software/modules/apps/vg/1.61.0/bin/vg"
    align_command = [
        vg_executable, "giraffe",
        "-p",         # tracks progress
        "-t", "16",   # adjust thread count as needed
        "-Z", reference_genome,
        "-d", dist_index[ref_choice],
        "-m", minimizer_index[ref_choice],
        "-f", fastq_files[0],
        "-f", fastq_files[1],
        "-D", "150"
    ]
    align_command = [arg for arg in align_command if arg]  # Remove any empty arguments
    print("Running alignment using vg giraffe...")
    run_command(align_command, output_gam)

    # Output directory for linear (BAM) outputs
    # For example: ~/houlahanlab/phamc20/data/AA_Patients/TCGA-BT-A20V/T2T-CHM13_linear
    if ref_choice == "GRCh38":
        bam_directory = os.path.join(output_base_dir, ancestry, donor_id, "GRCh38_linear")
    elif ref_choice == "T2T-CHM13":
        bam_directory = os.path.join(output_base_dir, ancestry, donor_id, "T2T-CHM13_linear")
    else:
        print("❌ Invalid reference genome choice.")
        return
    os.makedirs(bam_directory, exist_ok=True)

    # Run VG Surject (convert the GAM file to a BAM file)
    surjected_bam = os.path.join(bam_directory, f"{date}_{sample_id}_{ref_choice}_pangenome_fragment_length-150_aligned_surjected.bam")
    surject_command = [
        vg_executable, "surject",
        "-x", reference_genome,
        "-b", output_gam
    ]
    # For T2T-CHM13, add the -F parameter with its file path.
    if ref_choice == "T2T-CHM13":
        surject_command.insert(2, "-F")
        surject_command.insert(3, "/home/phamc20/houlahanlab/resources/T2T-CHM13_Indexes/version_1.1/VG-Index-Files/chm13_paths.txt")
    
    print("Running vg surject...")
    run_command(surject_command, surjected_bam)

    # Run Samtools Stats on the surjected BAM file and output stats
    surjected_bam_stats = os.path.join(bam_directory, f"{date}_{sample_id}_{ref_choice}_pangenome_fragment_length-150_aligned_surjected_stats.txt")
    surject_stats_command = ["samtools", "stats", surjected_bam]
    print("💻 Getting surjection statistics...")
    run_command(surject_stats_command, surjected_bam_stats)

    print("✅All processes completed successfully! Have a nice day!")

    return sample_id, date, ancestry, donor_id, ref_choice

def run_recal_markdup(sample_id, date, ancestry, donor_id, ref_choice):
    print("✨💻💖 Running Recalibration and MarkDuplicates Processing! 💖🔬✨")
    
    # Define input and output directories
    base_dir = f"/home/phamc20/houlahanlab/data/data/{ancestry}/{donor_id}/T2T-CHM13_linear" # figure out relative path...
    
    # Detect input BAM file based on sample_id
    input_bam = None
    for file in os.listdir(base_dir):
        if sample_id in file and file.endswith("_T2T-CHM13_pangenome_fragment_length-150_aligned_surjected.bam"):
            input_bam = os.path.join(base_dir, file)
            break

    if not input_bam:
        print(f"❌ Error: Could not find a BAM file for sample ID {sample_id} in {base_dir}.")
        return

    print(f"✅ Found input BAM file: {input_bam}")

    # Ooutput filenames
    rg_bam = f"{base_dir}/{sample_id}_CHM13_RG.bam"
    sorted_bam = f"{base_dir}/{sample_id}_CHM13_RG_sorted.bam"
    markdup_bam = f"{base_dir}/{sample_id}_CHM13_RG_sorted_markduplicates.bam"
    metrics_file = f"{base_dir}/{sample_id}_CHM13_sorted_markduplicates_metrics.txt"
    recal_table = f"{base_dir}/{sample_id}_recalibrator.table"
    recalibrated_bam = f"{base_dir}/{date}_{sample_id}_fragment_length_150_CHM13v1.1_RG_markduplicates_RG_recalibrated.bam"

    # Reference files
    reference_fasta = "/home/phamc20/houlahanlab/resources/T2T-CHM13_Indexes/version_1.1/Reference-Fasta-Files/CHM13v11Y_test_v1.1.fa"
    known_sites = "/home/phamc20/houlahanlab/resources/GRCh38_Indexes/Known-Indel-Files/Homo_sapiens_assembly38.known_indels_T2T-CHM13_edited.vcf"

    # Run commands
    print("Running gatk 'AddOrReplaceReadGroups'")
    run_command([
        "gatk", "AddOrReplaceReadGroups",
        "-I", input_bam,
        "-O", rg_bam,
        "--RGLB", "lib1",
        "--RGPU", "unit1",
        "--RGSM", sample_id,
        "--RGPL", "ILLUMINA"
    ])

    print("Sorting BAM files")
    run_command(["samtools", "sort", rg_bam, "-o", sorted_bam])

    print("Indexing BAM files")
    run_command(["samtools", "index", sorted_bam])

    print("Running Mark Duplicates")
    run_command([
        "gatk", "MarkDuplicates",
        "-I", sorted_bam,
        "-O", markdup_bam,
        "-M", metrics_file
    ])

    print("Running Base Recalibration")
    run_command([
        "gatk", "BaseRecalibrator",
        "-I", markdup_bam,
        "-R", reference_fasta,
        "--known-sites", known_sites,
        "-O", recal_table
    ])

    print("Applying Base Recalibration to BAM files")
    run_command([
        "gatk", "ApplyBQSR",
        "-I", markdup_bam,
        "-R", reference_fasta,
        "--bqsr-recal-file", recal_table,
        "-O", recalibrated_bam
    ])

    print("Indexing Recalibrated BAM files")
    run_command(["samtools", "index", recalibrated_bam])

    print(f"🎉✅ Vg Giraffe and Recalibration/Markduplicates completed successfully! All output files are stored in:")
    print(f"📁 {recalibrated_bam}")

def main():
    sample_id, date, ancestry, donor_id, ref_choice = run_vg_giraffe()
    if sample_id:
        choice = input("Do you want to run Recalibration and MarkDuplicates? (yes/no): ")
        if choice.lower() == 'yes':
            run_recal_markdup(sample_id, date, ancestry, donor_id, ref_choice)

if __name__ == "__main__":
    main()