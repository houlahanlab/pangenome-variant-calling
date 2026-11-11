#!/bin/bash
# Reference bias, step 2: germline calls on the VG pangenome normal BAM, in the target regions from step 1.
# One array task per donor (one row of manifest.tsv): HaplotypeCaller (-ERC GVCF), GenotypeGVCFs, then all calls
# written as CHROM POS REF ALT FILTER GT AD (no filtering here; build_paired_sites.py keeps the heterozygous SNVs).
#
# Input:  MANIFEST (from build_target_regions.py, or build_paired_sites.py --caller mutect2 --hc-manifest)
# Output: <out_prefix>.g.vcf.gz, <out_prefix>.vcf.gz, <out_prefix>_germline_summary.tsv
#
# Usage:  the sbatch command printed by build_target_regions.py, which sets the array size, logs, MANIFEST, GATK and BCFTOOLS

#SBATCH --job-name=germline_ref_bias
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

set -euo pipefail

GATK="${GATK:-gatk}"
BCFTOOLS="${BCFTOOLS:-bcftools}"

# Data row for this task = header (row 1) + SLURM_ARRAY_TASK_ID + 1
ROW=$(awk -F'\t' -v idx="$SLURM_ARRAY_TASK_ID" 'NR==idx+2' "$MANIFEST")
if [ -z "$ROW" ]; then
    echo "ERROR: no manifest row found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi
IFS=$'\t' read -r ARRAY_INDEX DONOR ALIGNMENT ANCESTRY NORMAL_BAM REFERENCE TARGET_BED OUT_PREFIX EXACT_SITES_FILE <<< "$ROW"

echo "Donor: $DONOR  Alignment: $ALIGNMENT"
echo "Normal BAM: $NORMAL_BAM"
echo "Target BED: $TARGET_BED"
for f in "$NORMAL_BAM" "$TARGET_BED"; do
    [ -f "$f" ] || { echo "ERROR: not found: $f"; exit 1; }
done

mkdir -p "$(dirname "$OUT_PREFIX")"
GVCF="${OUT_PREFIX}.g.vcf.gz"
VCF="${OUT_PREFIX}.vcf.gz"
SUMMARY="${OUT_PREFIX}_germline_summary.tsv"

"$GATK" --java-options "-Xmx3g" HaplotypeCaller \
    -R "$REFERENCE" \
    -I "$NORMAL_BAM" \
    -L "$TARGET_BED" \
    -O "$GVCF" \
    -ERC GVCF

"$GATK" GenotypeGVCFs \
    -R "$REFERENCE" \
    -V "$GVCF" \
    -O "$VCF"

{
    echo -e "CHROM\tPOS\tREF\tALT\tFILTER\tGT\tAD"
    "$BCFTOOLS" query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT\t%AD]\n' "$VCF"
} > "$SUMMARY"

echo "Done: $SUMMARY"
