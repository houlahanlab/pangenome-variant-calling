# This script will go through the filtered somatic variant output VCF and will filter the file, outputting a new VCF file that only contains variants with a FILTER status of PASS

import os

def filter_vcf(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
            else:
                columns = line.strip().split('\t')

                # The Filter column is column number 6, this checks this column to see if the value there is "PASS"
                if len(columns) > 6 and columns[6] == "PASS":
                    outfile.write(line)

input_vcf = "somatic_150v1.1_recalibratedxmarkduplicates_filtered.vcf"
output_vcf = "somatic_150v1.1_recalibratedxmarkduplicates_filtered_passfilter.vcf"

if os.path.exists(input_vcf):
    filter_vcf(input_vcf, output_vcf)
    print(f"Filtered VCF saved to {output_vcf}!")
else:
    print(f"Input file {input_vcf} does not exist!")