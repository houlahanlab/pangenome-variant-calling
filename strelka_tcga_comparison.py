import subprocess
import os

def run_command(cmd):
    """Run a shell command and print output."""
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def main():
    print("✨💻💖 Welcome to the Strelka vs TCGA Variant Comparison Tool! 🔬💖✨")

    # Get user inputs
    donor_id = input("Enter Donor ID (e.g., TCGA-BT-A20V): ").strip()
    sample_id = input("Enter Tumor Sample ID (e.g., TCGA-BT-A20V-01A-11D-A14W-08): ").strip()
    ancestry = input("Enter Patient Ancestry): ").strip()  

    # Ask user for reference genome choice
    # GRCh38 = linear, T2T = pangenome
    while True:
        ref_choice = input("Which reference genome would you like to use? (GRCh38/T2T-CHM13): ").strip()
        if ref_choice in ["GRCh38", "T2T-CHM13"]:
            break
        print("❌ Invalid choice. Please enter 'GRCh38' or 'T2T-CHM13'.")

    # Strelka VCF paths based on reference genome choice
    # GRCh38 
    if ref_choice == "GRCh38":
        strelka_vcf = f"/home/phamc20/houlahanlab/phamc20/data/strelka/{ancestry}/{donor_id}/GRCh38_{donor_id}/results/variants/somatic.snvs.vcf.gz" ## change this to "/home/phamc20/houlahanlab/phamc20/data/strelka/GRCh38/{donor_id}/results/variants/somatic.snvs.vcf.gz"
        output_dir = f"/home/phamc20/houlahanlab/phamc20/data/strelka_vs_tcga/GRCh38/{ancestry}/{donor_id}"
    else:  # T2T-CHM13
        strelka_vcf = f"/home/phamc20/houlahanlab/phamc20/data/strelka/{ancestry}/{donor_id}/T2T-CHM13_{donor_id}/results/variants/somatic.snvs.vcf.gz" ## change this to "/home/phamc20/houlahanlab/phamc20/data/strelka/T2T-CHM13/{donor_id}/results/variants/somatic.snvs.vcf.gz"
        output_dir = f"/home/phamc20/houlahanlab/phamc20/data/strelka_vs_tcga/T2T-CHM13/{ancestry}/{donor_id}"

    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists

    # TCGA VCF file
    tcga_vcf = "/home/phamc20/houlahanlab/abdelf6/raw_data/TCGA_Mutations/liftover_mc3.v0.2.8.PUBLIC.vcf"

    # Liftover paths (only used for T2T-CHM13 to lift from T2T-CHM13 to GRCh38)
    liftover_tool = "/workspace/lab/houlahanlab/software/liftOver"
    chain_file = "/home/phamc20/houlahanlab/resources/chain_files/chm13v2-grch38_fixed.chain"
    unmapped_bed = f"{output_dir}/unmapped_variants.bed"

    # Output file paths
    strelka_bed = f"{output_dir}/strelka_patient.bed"
    strelka_liftover_bed = f"{output_dir}/strelka_grch38.bed"
    strelka_grch38_positions = f"{output_dir}/strelka_grch38_positions.txt"
    tcga_positions = f"{output_dir}/tcga_patient_positions.txt"
    strelka_pass_positions = f"{output_dir}/strelka_pass_positions.txt"
    tcga_pass_positions = f"{output_dir}/tcga_pass_positions.txt"
    common_variants_full = f"{output_dir}/common_variants_full.txt"
    common_variants_pass = f"{output_dir}/common_variants_pass.txt"

    # Step 1: Convert Strelka VCF to BED
    def convert_strelka_vcf_to_bed():
        cmd = f"bcftools query -f '%CHROM\t%POS\n' {strelka_vcf} | awk '{{print $1\"\t\"$2-1\"\t\"$2}}' > {strelka_bed}"
        run_command(cmd)

    # Step 2: Perform liftover (Only for T2T-CHM13)
    # Trying to comapare everything in GRCh38
    def liftover_strelka():
        if ref_choice == "T2T-CHM13":
            cmd = f"{liftover_tool} {strelka_bed} {chain_file} {strelka_liftover_bed} {unmapped_bed}"
            run_command(cmd)
            cmd_cleanup = f"sed -i 's/CHM13#0#chr//' {strelka_liftover_bed}"
            run_command(cmd_cleanup)

    # Step 3: Extract TCGA variants
    def extract_tcga_positions():
        cmd = f"""bcftools view -s {sample_id} {tcga_vcf} | \
                  grep -v '\./\.:.:.' | \
                  awk '{{if ($10 ~ /0\/1/) print $1\"\t\"$2}}' > {tcga_positions}"""
        run_command(cmd)

    # Step 4: Extract strelka positions
    def extract_strelka_positions():
        if ref_choice == "GRCh38":
            # Remove "chr" prefix for GRCh38
            cmd = f"cat {strelka_bed} | sed -E 's/^chr([0-9XYM]+)/\\1/' | cut -f1,3 > {strelka_grch38_positions}"
            run_command(cmd)
        else:  # T2T-CHM13 (using lifted over file)
            cmd = f"cut -f1,3 {strelka_liftover_bed} > {strelka_grch38_positions}"
        run_command(cmd)

    # Step 5: Extract PASS variants
    # Only really looking at the PASS variants but i included the full list incase...
    def filter_pass_variants():
        strelka_pass_bed = f"{output_dir}/strelka_pass.bed"
        cmd = f"bcftools view -f PASS {strelka_vcf} | bcftools query -f '%CHROM\t%POS\n' | awk '{{print $1\"\t\"$2-1\"\t\"$2}}' > {strelka_pass_bed}"
        run_command(cmd)

        if ref_choice == "T2T-CHM13":
            strelka_liftover_pass_bed = f"{output_dir}/strelka_pass_grch38.bed"
            unmapped_pass_bed = f"{output_dir}/unmapped_pass_variants.bed"
            cmd = f"{liftover_tool} {strelka_pass_bed} {chain_file} {strelka_liftover_pass_bed} {unmapped_pass_bed}"
            run_command(cmd)

            strelka_pass_positions = f"{output_dir}/strelka_pass_positions.txt"
            cmd = f"awk '{{print $1\"\t\"$3}}' {strelka_liftover_pass_bed} | sed 's/CHM13#0#chr//' > {strelka_pass_positions}"
        else:
            strelka_pass_positions = f"{output_dir}/strelka_pass_positions.txt"
            cmd = f"cat {strelka_pass_bed} | sed -E 's/^chr([0-9XYM]+)/\\1/' | cut -f1,3 > {strelka_pass_positions}"
        run_command(cmd)

        cmd = f"""bcftools view -s {sample_id} {tcga_vcf} | \
                  awk '$7 == "PASS" && $10 ~ /0\/1/ {{print $1\"\t\"$2}}' > {tcga_pass_positions}"""
        run_command(cmd)

    # Step 6: Compare Strelka and TCGA variant lists
    def compare_variants():
        cmd = f"grep -F -f {strelka_grch38_positions} {tcga_positions} > {common_variants_full}"
        run_command(cmd)

        cmd = f"grep -F -f {strelka_pass_positions} {tcga_pass_positions} > {common_variants_pass}"
        run_command(cmd)

    # Step 7: Output summary statistics
    def output_summary():
        print("\n===== full variant summary =====")
        os.system(f"wc -l {strelka_grch38_positions}")
        os.system(f"wc -l {tcga_positions}")
        os.system(f"wc -l {common_variants_full}")

        print("\n===== pass variant summary =====")
        os.system(f"wc -l {strelka_pass_positions}")
        os.system(f"wc -l {tcga_pass_positions}")
        os.system(f"wc -l {common_variants_pass}")
    # Calulate TPR (recall) by dividing the common_variants_pass by tcga_pass_positions
    # Precision = common_variants_pass divided by the trelka_pass_positions

    # Run it all
    convert_strelka_vcf_to_bed()
    liftover_strelka()
    extract_tcga_positions()
    extract_strelka_positions()
    filter_pass_variants()
    compare_variants()
    output_summary()

    print("\n🎉✅ Analysis complete! All output files are stored in:")
    print(f"📁 {output_dir}")

if __name__ == "__main__":
    main()
