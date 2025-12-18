# This file will compare the donor specific txt files extracted from the VCF TCGA Mutation File to the Variant Callings for that specific donor conducted by us 

import pandas as pd

def compare_patient_variants(liftover_txt, patient_vcf, output_file, patient_id_tumor):
   
    liftover_txt_df = pd.read_csv(liftover_txt, sep = "\t", header=None)
    liftover_txt_df = liftover_txt_df.iloc[:,:2]
    liftover_txt_df.columns = ["CHROM", "POS"]
    print(liftover_txt_df.dtypes)
    print(liftover_txt_df.head())

   # Read the patient's VCF
    with open(patient_vcf, "r") as f:
       header_line_patient = [line for line in f if line.startswith("#CHROM")][0]
    patient_specific_pd = pd.read_csv(patient_vcf, sep="\t", comment="#", header=None)
    header_patient = header_line_patient.strip().lstrip("#").split("\t")
    patient_specific_pd.columns = header_patient
    print(patient_specific_pd.head())
    print(patient_specific_pd.dtypes)

    patient_specific_pd = patient_specific_pd[["CHROM", "POS"]]

   # Merge on `#CHROM` and `POS` to compare variants for the same patient
    merged = pd.merge(
       liftover_txt_df, 
       patient_specific_pd, 
       on=["CHROM", "POS"], 
       how="inner"
   )

   # Count the number of matching variants
    match_count = len(merged)

   # Write the header for the output VCF
    header = [
       f"## Comparison results for patient ' {patient_id_tumor}'",
       f"## Number of matching variants: {match_count}",
       "#CHROM\tPOS" + patient_id_tumor
   ]

   # Write the header and matching variants to the output VCF
    with open(output_file, "w") as out_vcf:
       out_vcf.write("\n".join(header) + "\n")
       merged.to_csv(out_vcf, sep="\t", index=False, header=False)


    print(f"Number of matching variants for patient '{patient_id_tumor}': {match_count}")

compare_patient_variants(
   liftover_txt="TCGA_Mutations/TCGA-BT-A20V-01A-11D-A14W-08.txt",
   patient_vcf="AA_Patients/TCGA-BT-A20V/Variant-Calling/OriginalBAMfiles/somatic_OriginalBAMfiles_filtered_passfilter_edited.vcf",
   output_file="TCGA-BT-A20V-Linear-Variant-Comparison.vcf",
   patient_id_tumor = "TCGA-BT-A20V-01A-11D-A14W-08.txt"
   )