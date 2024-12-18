import os
import subprocess
 
def run_command(command, output_file):
    print("Executing:", " ".join(command))
    with open(output_file, "wb") as output:
        result = subprocess.run(command, check=True, stdout=output)
    print(f"Command completed successfully. Output file saved to: {output_file}!")
    return result

def main():
    print("Welcome to the VG Giraffe Alignment Tool!")
    # Get user inputs
    donor_id = input("Enter Donor ID: ")
    sample_id = input("Enter Sample ID: ")
    ancestry = input("Enter Patient Ancestry: ")
    date = input("Please enter today's date (YYYY_MM_DD): ")
    # Construct the directory path
    base_directory = f"../{ancestry}/{donor_id}/GRCh38_linear"
    if not os.path.exists(base_directory):
        print(f"Error: Directory {base_directory} does not exist.")
        return
   # Extract FASTQ files
    fastq_files = [
        os.path.join(base_directory, file)
        for file in [f"{sample_id}_R1.fastq", f"{sample_id}_R2.fastq"]
    ]
    for file in fastq_files:
        if not os.path.exists(file):
            print("No FASTQ files found in the specified directory.")
            return
    print(f"FASTQ files found: {fastq_files}!")

    # Get reference genome choice
    references = {
        "GRCh38" : "../GRCh38_Indexes/VG-Index-Files/hprc-v1.0-mc-grch38.gbz",
        "T2T-CHM13": "../T2T-CHM13_Indexes/version_1.1/VG-Index-Files/hprc-v1.1-mc-chm13.d9.gbz "
    }

    dist_index = {
        "GRCh38" : "../GRCh38_Indexes/VG-Index-Files/hprc-v1.0-mc-grch38.dist",
        "T2T-CHM13" : "../T2T-CHM13_Indexes/version_1.1/VG-Index-Files/hprc-v1.1-mc-chm13.d9.dist"}
    
    minimizer_index = {
        "GRCh38" : "../GRCh38_Indexes/VG-Index-Files/hprc-v1.0-mc-grch38.min", 
        "T2T-CHM13" : "../T2T-CHM13_Indexes/version_1.1/VG-Index-Files/hprc-v1.1-mc-chm13.d9.min"}
    
    ref_choice = input("Choose the reference genome (GRCh38 or T2T-CHM13): ")
    reference_genome = references.get(ref_choice)
    if not reference_genome:
        print("Invalid reference genome choice.")
        return
    # Align FASTQ files using vg giraffe
    if ref_choice == "GRCh38":
        alignment_directory = f"../{ancestry}/{donor_id}/GRCh38_pangenome"
    
    if ref_choice == "T2T-CHM13":
        alignment_directory = f"../{ancestry}/{donor_id}/T2T-CHM13_pangenome"

    output_gam = os.path.join(alignment_directory, f"{date}_{sample_id}_{ref_choice}_pangenome_fragment_length-150_aligned.gam")
    align_command = [
        "vg", "giraffe",
        "-p", #tracks progress
        "-t", "16",  # Adjust thread count as needed]
        "-Z", reference_genome,
        "-d", dist_index[ref_choice],
        "-m", minimizer_index[ref_choice],
        "-f", fastq_files[0],
        "-f", fastq_files[1], 
        "-D",  "150",
        # "-o", "gam",
    ]
    align_command = [arg for arg in align_command if arg]  # Remove empty arguments
    print("Running alignment using vg giraffe...")
    run_command(align_command, output_gam)

     # Get alignment statistics
    output_vgstats_gam = os.path.join(alignment_directory, f"{date}_{sample_id}_{ref_choice}_pangenome_fragment_length-150_aligned_stats.txt")
    stats_command = ["vg", "stats", "-a", output_gam]
    print("Getting alignment statistics...")
    run_command(stats_command, output_vgstats_gam)

    # Surject the file
    if ref_choice == "GRCh38":
        bam_directory = f"../{ancestry}/{donor_id}/GRCh38_linear"
    
    if ref_choice == "T2T-CM13":
        bam_directory = f"../{ancestry}/{donor_id}/T2T-CHM13_linear"

    surjected_bam = os.path.join(bam_directory, f"{date}_{sample_id}_{ref_choice}_pangenome_fragment_length-150_aligned_surjected.bam")
    surject_command = [
        "vg", "surject",
        "-x", reference_genome,
        "-b", output_gam
    ]
    print("Running vg surject...")
    run_command(surject_command, surjected_bam)

    # Get surjection statistics
    surjected_bam_stats = os.path.join(bam_directory, f"{date}_{sample_id}_{ref_choice}_pangenome_fragment_length-150_aligned_surjected_stats.txt")
    surject_stats_command = ["samtools", "stats", surjected_bam]
    print("Getting surjection statistics...")
    run_command(surject_stats_command, surjected_bam_stats)

    print("All processes completed successfully! Have a nice day!")  

if __name__ == "__main__": 
    main()
