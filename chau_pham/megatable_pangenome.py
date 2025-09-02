# Mega pangenome variant file 
# What this does:
#   1) Read a 5-column donor list (donor_id, tumor_id, normal_id, ancestry, date).
#   2) For each donor and each caller (Mutect2, Strelka, SomaticSniper), collect:
#        - TP file  (sites called by both caller and TCGA)
#        - FP file  (sites unique to the caller)
#        - FN file  (sites unique to TCGA)
#   3) Write a single "Megatable.txt" with CHROM/START/END/REF/ALT + sample/status/tool/ancestry.
# (has duplicates of those that may have been called by each of the variant callers)

import os

def run_command(cmd):
    print(f"Running: {cmd}")
    os.system(cmd)
    print("Command completed successfully!\n")

# Writes #REF/#ALT/ANNOTATIONS as empty placeholders, more to be added later
def write_variants_with_status(input_file, output_handle, status, tool, sample_id, ancestry):
    if not os.path.exists(input_file):
        print(f"Warning: {input_file} not found, skipping...")
        return

    with open(input_file, "r") as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) < 5:
                continue  # Skip incomplete lines
            chrom, start, end, ref, alt = fields[:5]
            output_handle.write(f"{chrom}\t{start}\t{end}\t{ref}\t{alt}\t{sample_id}\t{status}\t{tool}\t{ancestry}\t\t\n")

def process_variants(donor_id, sample_id, ancestry, caller_name, output_handle):
    
    if caller_name == "Mutect2":
        base_dir = f"/workspace/lab/houlahanlab/phamc20/data/TCGA_Comparison/mutect_vs_tcga/T2T-CHM13/{ancestry}/{donor_id}" # <-- changed reference 
        tp_file = os.path.join(base_dir, "common_variants_full.txt")
        fp_file = os.path.join(base_dir, "mutect_unique_pass.txt")
        fn_file = os.path.join(base_dir, "tcga_unique_pass.txt")

    elif caller_name == "Strelka":
        base_dir = f"/workspace/lab/houlahanlab/phamc20/data/TCGA_Comparison/strelka_vs_tcga2/T2T-CHM13/{ancestry}/{donor_id}" # <-- changed reference 
        tp_file = os.path.join(base_dir, "common_variants_pass.txt")  # <- different for Strelka
        fp_file = os.path.join(base_dir, "strelka_unique_pass.txt")
        fn_file = os.path.join(base_dir, "tcga_unique_pass.txt")

    elif caller_name == "SomaticSniper":
        base_dir = f"/workspace/lab/houlahanlab/phamc20/data/TCGA_Comparison/somatic_vs_tcga2/T2T-CHM13/{ancestry}/{donor_id}" # <-- changed reference 
        tp_file = os.path.join(base_dir, "common_variants_full.txt") 
        fp_file = os.path.join(base_dir, "somatic_unique_pass.txt")
        fn_file = os.path.join(base_dir, "tcga_unique_pass.txt")

    else:
        print(f"Unknown caller: {caller_name}")
        return

    write_variants_with_status(tp_file, output_handle, "TP", caller_name, sample_id, ancestry)
    write_variants_with_status(fn_file, output_handle, "FN", caller_name, sample_id, ancestry)
    write_variants_with_status(fp_file, output_handle, "FP", caller_name, sample_id, ancestry)

def main():
    print("Creating MEGA variant file!")

    sample_file = input("Enter path to your 5-column donor list file: ").strip()
    if not os.path.exists(sample_file):
        print(f"Error: File {sample_file} not found!")
        return

    mega_output_path = "/workspace/lab/houlahanlab/data/data/mega_variant_summary/GRCh38_Pangenome/Megatable.txt"
    os.makedirs(os.path.dirname(mega_output_path), exist_ok=True)

    with open(mega_output_path, "w") as mega_out:
        # Write header
        mega_out.write("CHROM\tSTART\tEND\tREF\tALT\tSAMPLE\tSTATUS\tTOOL\tANCESTRY\t#REF\t#ALT\tANNOTATIONS\n")

        with open(sample_file, "r") as f:
            patient_entries = [line.strip().split() for line in f if line.strip()]

        for entry in patient_entries:
            if len(entry) < 5:
                continue

            donor_id = entry[0]
            sample_id = entry[1]  # <-- 2nd column is tumour_id
            normal_id = entry[2]
            ancestry = entry[3]

            process_variants(donor_id, sample_id, ancestry, "Mutect2", mega_out)
            process_variants(donor_id, sample_id, ancestry, "Strelka", mega_out)
            process_variants(donor_id, sample_id, ancestry, "SomaticSniper", mega_out)

    print(f"\n Mega variant file created at: {mega_output_path}")

if __name__ == "__main__":
    main()