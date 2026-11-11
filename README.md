# pangenome-variant-calling

This repository has the code used to compare somatic SNV calling with pangenome alignment (vg giraffe) and linear alignment (BWA) in TCGA bladder cancer (BLCA, 45 donors) and lung adenocarcinoma (LUAD, 29 donors).

## Folders

| Folder | What it has |
|---|---|
| `1_pipeline/` | BAM to FASTQ conversion and the alignment / variant calling config |
| `2_benchmarking/` | Comparison of the calls with MC3, consensus calls |
| `3_analyses/` | False positive, false negative, reference bias, germline contamination and variant impact analyses |
| `4_figures_and_tables/` | One R script per figure (BLCA: Figures 1-4 and Supplementary Figures 1-12; LUAD: Figure 5 and Supplementary Figures 13-15) and Supplementary Tables 1-3 |
| `source_data/` | The numbers behind each figure panel, one CSV per panel |
| `metadata/` | Patient lists and ancestry admixture |

## Setup

1. Copy `config.example.yaml` to `config.yaml` and fill in the paths for your system.
2. Run every script from the top folder of this repository.

Every script reads its paths from `config.yaml`. The top of each script says what it reads and what it makes.

## Alignment names

The code uses these names for the four alignments:

| In the code | In the paper |
|---|---|
| `GRCh38_bwa` | BWA GRCh38 |
| `GRCh38_linear` | VG GRCh38 |
| `GRCh38_hprc_v1.1` | VG GRCh38 pangenome |
| `T2T_pangenome` | VG pangenome |

## Steps

### 1. Alignment and variant calling

Alignment and somatic variant calling (Strelka2, Mutect2, SomaticSniper, DeepSomatic) were run with [pangenome-pipeline](https://github.com/houlahanlab/pangenome-pipeline).
The VG pangenome alignment (HPRC v1.1 + T2T-CHM13) was run with an earlier version of the same pipeline, using the settings in `1_pipeline/configs/T2T_pangenome_config.yaml`.
BWA GRCh38 alignments are the TCGA BAMs from the GDC.

The TCGA BAMs were turned back into FASTQs for vg giraffe with:

```bash
python3 1_pipeline/bam_to_fastq.py --patients metadata/BLCA_Patients.txt
python3 1_pipeline/bam_to_fastq.py --patients metadata/LUAD_All_Patients_Full_ID.txt
```

### 2. Compare the calls with MC3

```bash
# Each caller and alignment (for LUAD add --cohort LUAD; LUAD has Strelka and Mutect2 only)
for caller in Strelka Mutect2 SomaticSniper DeepSomatic; do
  for alignment in GRCh38_bwa GRCh38_linear T2T_pangenome; do
    python3 2_benchmarking/mc3_comparison/compare_to_mc3.py --caller $caller --alignment $alignment
  done
done

# After removing calls found in gnomAD (Supplementary Figure 9i)
python3 2_benchmarking/mc3_comparison/compare_to_mc3.py --caller Strelka --alignment T2T_pangenome --gnomad
python3 2_benchmarking/mc3_comparison/compare_to_mc3.py --caller Mutect2 --alignment T2T_pangenome --gnomad
```

### 3. Consensus calls

SNVs called by at least 2 or at least 3 of the four callers, then Strelka2 and Mutect2 compared with the 3-caller consensus.

```bash
python3 2_benchmarking/consensus/build_consensus.py --alignment T2T_pangenome
python3 2_benchmarking/consensus/build_consensus.py --alignment GRCh38_bwa
for alignment in GRCh38_bwa GRCh38_linear T2T_pangenome; do
  python3 2_benchmarking/consensus/compare_to_consensus.py --truth T2T_pangenome --alignment $alignment
  python3 2_benchmarking/consensus/compare_to_consensus.py --truth GRCh38_bwa --alignment $alignment
done
```

Note: steps 2 and 3 add rows to their results tables every time they run, so delete the old tables before running them again.

### 4. Table of all TP, FP and FN SNVs

Used by the analyses in step 5.

```bash
python3 2_benchmarking/variant_table/build_variant_table.py
```

### 5. Analyses

Scripts marked (pysam) need pysam installed.
Some steps print an `sbatch` command. Run it, and wait for all the jobs to finish before the next step.

False positives that are known germline variants (Supplementary Figures 8-9):

```bash
python3 3_analyses/false_positives/germline_contamination.py
```

Reference bias (Figure 3d, Supplementary Figure 10). Strelka2 first, then Mutect2:

```bash
python3 3_analyses/false_positives/reference_bias/germline_analysis/build_target_regions.py   # then run the printed sbatch
python3 3_analyses/false_positives/reference_bias/paired_sites/build_paired_sites.py          # then run the printed sbatch
python3 3_analyses/false_positives/reference_bias/paired_sites/combine_paired_counts.py

python3 3_analyses/false_positives/reference_bias/paired_sites/build_paired_sites.py --caller mutect2 --hc-manifest   # then run the printed sbatch
python3 3_analyses/false_positives/reference_bias/paired_sites/build_paired_sites.py --caller mutect2                 # then run the printed sbatch
python3 3_analyses/false_positives/reference_bias/paired_sites/combine_paired_counts.py --caller mutect2
```

False negatives (Figure 3e-f, Supplementary Figures 11-12):

```bash
python3 3_analyses/false_negatives/measure_vaf.py
python3 3_analyses/false_negatives/evs_measure_strelka.py                # (pysam)
python3 3_analyses/false_negatives/tlod_measure_mutect2.py
python3 3_analyses/false_negatives/false_negative_filtered_fraction.py   # prints numbers quoted in the text
python3 3_analyses/false_negatives/build_kmt2c_read_placement.py         # (pysam)
```

Variant impact (Figure 4, Supplementary Figure 12, Supplementary Table 3):

```bash
python3 3_analyses/impact_variants/annotate_variant_sets.py
python3 3_analyses/impact_variants/allele_frequency/build_ref_alt_freq.py --alignment T2T_pangenome   # (pysam)
python3 3_analyses/impact_variants/allele_frequency/build_ref_alt_freq.py --alignment GRCh38_bwa      # (pysam)
```

### 6. Supplementary Tables

```bash
python3 2_benchmarking/supplementary_tables/build_supplementary_table_1.py
python3 2_benchmarking/supplementary_tables/build_supplementary_table_2.py
python3 2_benchmarking/supplementary_tables/build_supplementary_table_3.py
```

### 7. Figures and Source Data

Each script in `4_figures/` makes the figure panels in its name, and saves the numbers behind each panel to `source_data/` (for example `source_data/Figure_3b.csv`).

```bash
Rscript 4_figures/BLCA/main/figure1/figure1b-d_variant_callers_alignments.R   # and the other figure scripts
```

To combine all the CSVs into one Excel file (one sheet per panel), run `Rscript 4_figures/build_source_data.R`.
