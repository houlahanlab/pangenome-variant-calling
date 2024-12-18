# This code will go through the VCF file and only extract the chrosomone number, chromosome position, reference allele, and alternate allele for each variant 

import vcf
vcf_input_file = 'liftover_somatic150DLv1.1_markedduplicates_recalibrated_filtered_passfilter.vcf'
vcf_output_file = 'vcf_filtered_liftover_somatic150DLv1.1_markedduplicates_recalibrated_filtered_passfilter.vcf'

records = vcf.Reader(open(vcf_input_file, 'r'))

# records is an iterable, from which you can get attributes such as REF, ALT, POS etc.

with open(vcf_output_file, 'w') as outfile:

    outfile.write(f"{vcf_output_file} output file after being filtered and read\n")
    for row in records:
        chr = row.CHROM
        pos = row.POS
        id = row.ID
        ref = row.REF

        # If there is more than one alternate allele groups, to add them with commas 
        if row.ALT:
            alt = ','.join([str(allele) for allele in row.ALT])
        else:
            alt = "."
            
        outfile.write(f"chr is {chr}, pos is {pos}, reference is {ref}, alternate allele is {alt}\n")

print(f"Output written to {vcf_output_file}")

    
