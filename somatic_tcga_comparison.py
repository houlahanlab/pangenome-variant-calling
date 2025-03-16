import subprocess
import os

def run_command(cmd):
    """Run a shell command and print output."""
    subprocess.run(cmd, shell=True, check=True)

def process_patient(donor_id, sample_id, ancestry, ref_choice, tcga_vcf, liftover_tool, chain_file):
    """Processes a single patient sample and returns the summary string."""
    
    print(f"\n Processing: Donor: {donor_id} | Tumor: {sample_id} | Ancestry: {ancestry}")

    #GRCh38
    if ref_choice == "GRCh38":
        somatic_vcf = f"/home/phamc20/houlahanlab/phamc20/data/somatic_sniper/{ancestry}/{donor_id}/{ref_choice}/snp_output_file.vcf"
        output_dir = f"/home/phamc20/houlahanlab/phamc20/data/somatic_vs_tcga/GRCh38/{ancestry}/{donor_id}"
    else:  # T2T-CHM13
        somatic_vcf = f"/home/phamc20/houlahanlab/phamc20/data/somatic_sniper/{ancestry}/{donor_id}/{ref_choice}/snp_output_file.vcf"
        output_dir = f"/home/phamc20/houlahanlab/phamc20/data/somatic_vs_tcga/T2T-CHM13/{ancestry}/{donor_id}"

    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists

    # Output file paths
    somatic_bed = f"{output_dir}/somatic_patient.bed"
    somatic_liftover_bed = f"{output_dir}/somatic_grch38.bed"
    somatic_grch38_positions = f"{output_dir}/somatic_grch38_positions.txt"
    tcga_positions = f"{output_dir}/tcga_patient_positions.txt"
    common_variants_full = f"{output_dir}/common_variants_full.txt"

    # Step 1: Convert Strelka VCF to BED
    run_command(f"bcftools query -f '%CHROM\t%POS\n' {somatic_vcf} | awk '{{print $1\"\t\"$2-1\"\t\"$2}}' > {somatic_bed}")

    # Step 2: Perform liftover (Only for T2T-CHM13)
    if ref_choice == "T2T-CHM13":
        run_command(f"{liftover_tool} {somatic_bed} {chain_file} {somatic_liftover_bed} {output_dir}/unmapped_variants.bed")
        run_command(f"sed -i 's/CHM13#0#chr//' {somatic_liftover_bed}")

    # Step 3: Extract TCGA variants
    run_command(f"""bcftools view -s {sample_id} {tcga_vcf} | \
                  grep -v '\./\.:.:.' | \
                  awk '{{if ($10 ~ /0\/1/) print $1\"\t\"$2}}' > {tcga_positions}""")

    # Step 4: Extract somatic positions
    if ref_choice == "GRCh38":
        run_command(f"cat {somatic_bed} | sed -E 's/^chr([0-9XYM]+)/\\1/' | cut -f1,3 > {somatic_grch38_positions}")
    else:  # T2T-CHM13 (using lifted over file)
        run_command(f"cut -f1,3 {somatic_liftover_bed} > {somatic_grch38_positions}")

    # Step 5: Compare somatic and TCGA variant lists
    cmd = f"grep -F -f {somatic_grch38_positions} {tcga_positions} > {common_variants_full}"

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode == 1:  # grep exits with code 1 when no matches are found
       print(f"Warning: No matching variants found between {somatic_grch38_positions} and {tcga_positions}.")
       open(common_variants_full, 'w').close()  # Create an empty file to avoid errors later.
    elif result.returncode != 0:
       print(f"Error: grep failed with code {result.returncode}")
       print(f"STDERR: {result.stderr}")

    # Step 6: Get variant counts after processing is done
    somatic_count = int(subprocess.getoutput(f"wc -l < {somatic_grch38_positions}").strip() or 0)
    tcga_count = int(subprocess.getoutput(f"wc -l < {tcga_positions}").strip() or 0)
    common_count = int(subprocess.getoutput(f"wc -l < {common_variants_full}").strip() or 0)

    # Generate summary string for this patient
    summary = f"""
===== Full Variant Summary for {donor_id} =====
    {somatic_count} somatic_grch38_positions.txt
    {tcga_count} tcga_patient_positions.txt
    {common_count} common_variants_full.txt
"""

    return summary


def main():
    print("✨💻💖 Welcome to the Somatic vs TCGA Variant Comparison Tool! 🔬💖✨")

    # Ask for reference genome choice
    while True:
        ref_choice = input("Which reference genome would you like to use? (GRCh38/T2T-CHM13): ").strip()
        if ref_choice in ["GRCh38", "T2T-CHM13"]:
            break
        print("❌ Invalid choice. Please enter 'GRCh38' or 'T2T-CHM13'.")

    # TCGA VCF file path
    tcga_vcf = "/home/phamc20/houlahanlab/abdelf6/raw_data/TCGA_Mutations/liftover_mc3.v0.2.8.PUBLIC.vcf"

    # Liftover paths (only used for T2T-CHM13)
    liftover_tool = "/workspace/lab/houlahanlab/software/liftOver"
    chain_file = "/home/phamc20/houlahanlab/resources/chain_files/chm13v2-grch38_fixed.chain"

    # Ask user whether to process a single patient or a list
    while True:
        mode = input("Would you like to manually process patients (1) or enter a list of patients (2)? ").strip()
        if mode in ["1", "2"]:
            break
        print("❌ Invalid choice. Please enter '1' or '2'.")

    all_summaries = []  # List to store all patient summaries

    if mode == "1":
        # Manually enters patient details
        donor_id = input("Enter Donor ID (e.g., TCGA-BT-A20V): ").strip()
        sample_id = input("Enter Tumor Sample ID (e.g., TCGA-BT-A20V-01A-11D-A14W-08): ").strip()
        ancestry = input("Enter Patient Ancestry (e.g AA_Patients): ").strip()

        summary = process_patient(donor_id, sample_id, ancestry, ref_choice, tcga_vcf, liftover_tool, chain_file)
        all_summaries.append(summary)

    else:
        # Use existing SomaticSniper patient list
        sample_file = input("Enter path to SomaticSniper patient list file: ").strip()

        if not os.path.exists(sample_file):
            print(f"❌ Error: Patient list file {sample_file} not found!")
            return

        with open(sample_file, "r") as f:
            patient_entries = [line.strip().split() for line in f.readlines()]

        # Extract only DonorID, TumorID, and Ancestry (Columns 1, 2, and 4)
        for entry in patient_entries:
            if not entry or len(entry) < 5:
                print("\n🚨 Detected empty entry. Ending batch processing and printing results...")
                break

            donor_id, sample_id, _, ancestry, _ = entry 
            summary = process_patient(donor_id, sample_id, ancestry, ref_choice, tcga_vcf, liftover_tool, chain_file)
            all_summaries.append(summary)

    # Print all summaries at the end
    if all_summaries:
        print("\n📢 **Final Summary for All Patients:**")
        for summary in all_summaries:
            print(summary)
    else:
        print("\n No valid patients processed. Please check the input file!")

    print("\n🎉✅ Analysis Completed Successfully!")


if __name__ == "__main__":
    main()



