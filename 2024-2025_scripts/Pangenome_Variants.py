# This code finds variants only found in the Pangenome Aligned Variant Calling (compares both the Pangenome Aligned Variant Calling Files and the Linear Aligned Variant Calling Files)

import pandas as pd

def read_vcf(file_path):
    data = []
    with open(file_path, 'r') as f:
       for line in f:
           if line.startswith("chr is"):
               parts = line.strip().split(", ")
               chrom = parts[0].split(" ")[-1].replace("chr", "")  # Remove "chr" prefix if present
               pos = int(parts[1].split(" ")[-1])  # Ensure POS is an integer
               ref = parts[2].split(" ")[-1]
               alt = parts[3].split(" ")[-1]
               data.append([chrom, pos, ref, alt])
    return pd.DataFrame(data, columns=['CHROM', 'POS', 'REF', 'ALT'])

def calculate_unmatched_ratio(original_bam_file, GAMtoBAM_file, output_file):
    original_bam = read_vcf(original_bam_file)
    GAMtoBAM = read_vcf(GAMtoBAM_file)

#merge files to find the matching and unmatched variants in both 
    merged = pd.merge(
        GAMtoBAM, 
        original_bam,
        on = ['CHROM', 'POS', 'REF', 'ALT'],
        how = 'left',
        indicator=True
    )

# filters variants in the left file only --> only in the GAMtoBAM file
    unmatched_variants = merged[merged['_merge'] == 'left_only']

    total_GAMtoBAM_variants = count_variants_in_file(GAMtoBAM_file)
    unmatched_variants_count = len(unmatched_variants)
    ratio = unmatched_variants_count/total_GAMtoBAM_variants if total_GAMtoBAM_variants > 0 else 0

    with open(output_file, 'w') as f:
        f.write(f"Total Variants in pangenome alignment file: {total_GAMtoBAM_variants}\n")
        f.write(f"Variants only in pangenome alignment file: {unmatched_variants_count}\n")
        f.write(f"Ratio (Unmatched/Total Pangenome Variants): {ratio:.4f}\n")
        f.write("\nUnmatched Variants:\n")
        
        #Writes the unmatched variants into a txt file

        for index, row in unmatched_variants.iterrows():
            f.write(f"CHROM: {row['CHROM']}, POS: {row['POS']}, Reference: {row['REF']}, Alternate: {row['ALT']}\n")

    print(f"Comparison Complete! The results are in written to: {output_file}!")
    
#counts the variants in the GAMtoBAM file
def count_variants_in_file(file_path):
    count = 0
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("chr is"):
                count +=1
    return count 

original_bam_file = 'OriginalBAMfiles/vcf_filtered_somatic_originalBAMfilesfiltered_passfilter.vcf'
GAMtoBAM_file = 'GAMtoBAM150v1.1/vcf_filtered_liftover_somatic_GAMtoBAM150v1.1_filtered_passfilter.vcf'
output_file = 'TCGA-FD-A3B4-Pangenome-Specific-Variants.txt'

calculate_unmatched_ratio(original_bam_file, GAMtoBAM_file, output_file)

