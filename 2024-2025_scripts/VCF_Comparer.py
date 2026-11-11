# This script compares the somatic variant calling between files and can be used to see matching variants between files (you can compare up to 4 files at the moment)

import pandas as pd
from collections import defaultdict

def read_vcf(file_path):
   # Reads the specific output format from the filtered VCF file and converts it to a DataFrame
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

def chrom_key(chrom):
   # Goes through the chromosomes to make sure they are either numbers or X/Y
   try:
       return int(chrom)  # Convert numeric chromosomes to integers
   except ValueError:
       if chrom == 'X':
           return float('inf')  # X comes after numbers
       elif chrom == 'Y':
           return float('inf') + 1  # Y comes after X
       else:
           return float('inf') + 2  # Other non-numeric chromosomes
       
def sort_variants(variants):
   # Sorts variants by chromosome and position
   return sorted(variants, key=lambda x: (chrom_key(str(x[0])), int(x[1])))

def compare_variants(vcf_files, output_file):
   # Compares variants across multiple VCF files and organizes by categories and subcategories
   vcf_data = [read_vcf(file) for file in vcf_files]
   all_variants = pd.concat(
       [df.assign(File=file) for df, file in zip(vcf_data, vcf_files)],
       ignore_index=True
   )

   grouped = all_variants.groupby(['CHROM', 'POS'])
   # Main categories with subcategories
   results = {
       "All Files": {"Same Everything": [], "Same Ref": [], "Different Ref/Alt": []},
       "3 Files": {"Same Everything": [], "Same Ref": [], "Different Ref/Alt": []},
       "2 Files": {"Same Everything": [], "Same Ref": [], "Different Ref/Alt": []},
   }

   # Chromosome-specific count for each category
   chrom_counts = defaultdict(lambda: defaultdict(int))
   all_files_set = set(vcf_files)
   for (chrom, pos), group in grouped:
       unique_refs = group['REF'].unique()
       unique_alts = group['ALT'].unique()
       files_present = set(group['File'].unique())
       file_count = len(files_present)
       # Determine the main category
       if file_count == len(vcf_files):
           category = "All Files"
       elif file_count == 3:
           category = "3 Files"
       elif file_count == 2:
           category = "2 Files"
       else:
           continue
       
       # Counts how many variants per chromosome 
       chrom_counts[category][chrom] += 1

       # Determines the subcategory of the variant 
       if len(unique_refs) == 1 and len(unique_alts) == 1:
           if category == "3 Files":
               files_not_present = list(all_files_set - files_present)
               results[category]["Same Everything"].append((chrom, pos, unique_refs[0], unique_alts[0], files_not_present))
           elif category == "2 Files":
               results[category]["Same Everything"].append((chrom, pos, unique_refs[0], unique_alts[0], list(files_present)))
           else:
               results[category]["Same Everything"].append((chrom, pos, unique_refs[0], unique_alts[0]))
       elif len(unique_refs) == 1:
           detailed_info = group.groupby('File').apply(
               lambda x: x[['REF', 'ALT']].to_dict('records')
           ).to_dict()
           if category == "3 Files":
               files_not_present = list(all_files_set - files_present)
               results[category]["Same Ref"].append((chrom, pos, unique_refs[0], detailed_info, files_not_present))
           elif category == "2 Files":
               results[category]["Same Ref"].append((chrom, pos, unique_refs[0], detailed_info, list(files_present)))
       else:
           detailed_info = group.groupby('File').apply(
               lambda x: x[['REF', 'ALT']].to_dict('records')
           ).to_dict()
           if category == "3 Files":
               files_not_present = list(all_files_set - files_present)
               results[category]["Different Ref/Alt"].append((chrom, pos, detailed_info, files_not_present))
           elif category == "2 Files":
               results[category]["Different Ref/Alt"].append((chrom, pos, detailed_info, list(files_present)))
   # Sort results within each subcategory
   for category in results.keys():
       for subcategory in results[category].keys():
           results[category][subcategory] = sorted(
               results[category][subcategory],
               key=lambda x: (chrom_key(str(x[0])), int(x[1]))
           )

   # Write results to the output file
   with open(output_file, 'w') as f:
       f.write("Comparison Results:\n\n")
       # Overall counts
       for category, subcategories in results.items():
           total_variants = sum(len(variants) for variants in subcategories.values())
           f.write(f"{category} ({total_variants} total variants):\n")
           for subcategory, variants in subcategories.items():
               f.write(f"  {subcategory} ({len(variants)} variants):\n")
               for variant in variants:
                   if subcategory == "Same Everything":
                       if category == "3 Files":
                           chrom, pos, ref, alt, files_not_present = variant
                           f.write(f"    CHROM: {chrom}, POS: {pos}, Reference: {ref}, Alternate: {alt}, Files Not Present: {', '.join(files_not_present)}\n")
                       elif category == "2 Files":
                           chrom, pos, ref, alt, files_present = variant
                           f.write(f"    CHROM: {chrom}, POS: {pos}, Reference: {ref}, Alternate: {alt}, Files Present: {', '.join(files_present)}\n")
                       else:
                           chrom, pos, ref, alt = variant
                           f.write(f"    CHROM: {chrom}, POS: {pos}, Reference: {ref}, Alternate: {alt}\n")
                   elif subcategory == "Same Ref":
                       chrom, pos, ref, detailed_info, files_info = variant
                       if category == "3 Files":
                           f.write(f"    CHROM: {chrom}, POS: {pos}, Reference: {ref}, Files Not Present: {', '.join(files_info)}\n")
                       elif category == "2 Files":
                           f.write(f"    CHROM: {chrom}, POS: {pos}, Reference: {ref}, Files Present: {', '.join(files_info)}\n")
                       for file, records in detailed_info.items():
                           f.write(f"      File: {file}\n")
                           for record in records:
                               f.write(f"        Alternate: {record['ALT']}\n")
                   elif subcategory == "Different Ref/Alt":
                       chrom, pos, detailed_info, files_info = variant
                       if category == "3 Files":
                           f.write(f"    CHROM: {chrom}, POS: {pos}, Files Not Present: {', '.join(files_info)}\n")
                       elif category == "2 Files":
                           f.write(f"    CHROM: {chrom}, POS: {pos}, Files Present: {', '.join(files_info)}\n")
                       for file, records in detailed_info.items():
                           f.write(f"      File: {file}\n")
                           for record in records:
                               f.write(f"        Reference: {record['REF']}, Alternate: {record['ALT']}\n")
               if not variants:
                   f.write("    None\n")
               f.write("\n")

       # Write chromosome-specific counts in sorted order
       f.write("Chromosome Counts per Category:\n")
       for category, counts in chrom_counts.items():
           f.write(f"{category}:\n")
           sorted_counts = sorted(counts.items(), key=lambda x: chrom_key(x[0]))
           for chrom, count in sorted_counts:
               f.write(f"  Chromosome {chrom}: {count} variants\n")
               
   # Print overall and chromosome-specific counts to the user
   print("Comparison Summary:")
   for category, subcategories in results.items():
       total_variants = sum(len(variants) for variants in subcategories.values())
       print(f"{category}: {total_variants} variants")
       for subcategory, variants in subcategories.items():
           print(f"  {subcategory}: {len(variants)} variants")
   print("\nChromosome Counts per Category:")
   for category, counts in chrom_counts.items():
       print(f"{category}:")
       sorted_counts = sorted(counts.items(), key=lambda x: chrom_key(x[0]))
       for chrom, count in sorted_counts:
           print(f"  Chromosome {chrom}: {count} variants")
           
# Input and output
vcf_files = [
   'GAMtoBAM150v1.1/Recalibration-and-MarkDuplicates/vcf_filtered_liftover_somatic150DLv1.1_markedduplicates_recalibrated_filtered_passfilter.vcf',
   'OriginalBAMfiles/vcf_filtered_somatic_originalBAMfiles_filtered_passfilter.vcf'
]

output_file = 'VCF_Comparison_Results_with_MarkDuplicates.txt'
compare_variants(vcf_files, output_file)
print(f"\nComparison results have been written to {output_file}!")

