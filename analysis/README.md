# pangenome-variant-calling
Scripts to benchmark somatic mutation calling using the human pangenome.

Scripts Uploaded: 
- Scripts to create plots:
    - Variant Calling Plots: Variant-Calling-Plot-TP.py and Variant-Calling-Plot-FN.py
    - Plots for Alignment Statistics: GraphSubPlots.py
      
- Scripts to Filter VCF Files and the MAF TCGA Mutation File:
    - MAF-Filter.py
    - VCF-Variant-Filter=PASS-Filtering.py
    - VCF_Reader.py
      
- Scripts to Compare the Variants in the MAF TCGA Mutation File and Donor Variant Calls and to Compare the Variant Calls between the Pangenome Alignment and Linear Alignment:
    - MAF-VCF-Comparer.py
    - Pangenome_Variants.py
    - VCF_Comparer.py

- Scripts to Run the VG Giraffe Command (as well as VG Stats, VG Surject, and Samtools Stats):
    - VG-Giraffe.py
