#!/bin/bash
# Reference bias, step 4: read counts for each allele at the same sites in every alignment.
# One array task per donor = its 9 rows of manifest.tsv (3 alignments x tumour at the somatic sites, normal at the
# somatic sites, normal at the germline sites). A pileup is used instead of each caller's VCF because a VCF has no
# record at sites that alignment did not call.
# Same settings for every BAM: mapping and base quality >= 20, no BAQ, duplicate / secondary / QC-fail reads skipped
# (bcftools defaults), depth cap 100000, SNVs only; FORMAT/AD gives the reads per base.
# The VG BAM indexes are named <sample>_recalibrated.bai, which htslib does not find on its own, so each BAM and its
# index are linked as <name>.bam / <name>.bam.bai first.
#
# Input:  MANIFEST (from build_paired_sites.py)
# Output: <out_file> per manifest row (CHROM POS REF ALT AD)
#
# Usage:  the sbatch command printed by build_paired_sites.py, which sets the array size, logs, MANIFEST and BCFTOOLS

#SBATCH --job-name=paired_pileup
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00

set -euo pipefail

BCFTOOLS="${BCFTOOLS:-bcftools}"
ROWS_PER_TASK=9
FIRST_ROW=$((SLURM_ARRAY_TASK_ID * ROWS_PER_TASK))

for ROW_INDEX in $(seq "$FIRST_ROW" $((FIRST_ROW + ROWS_PER_TASK - 1))); do
    ROW=$(awk -F'\t' -v idx="$ROW_INDEX" 'NR==idx+2' "$MANIFEST")
    if [ -z "$ROW" ]; then
        echo "ERROR: no manifest row for row index $ROW_INDEX"
        exit 1
    fi
    IFS=$'\t' read -r ARRAY_INDEX DONOR ALIGNMENT SAMPLE_TYPE BAM BAM_INDEX REFERENCE REGIONS OUT_FILE <<< "$ROW"
    echo "Donor: $DONOR  Alignment: $ALIGNMENT  Sample type: $SAMPLE_TYPE"

    OUT_DIR="$(dirname "$OUT_FILE")"
    mkdir -p "$OUT_DIR"
    {
        echo -e "CHROM\tPOS\tREF\tALT\tAD"
        # No sites for this donor (e.g. no germline variants nearby): header only
        if [ -s "$REGIONS" ]; then
            LINK="$OUT_DIR/${DONOR}_${ALIGNMENT}_${SAMPLE_TYPE}"
            ln -sf "$BAM" "$LINK.bam"
            ln -sf "$BAM_INDEX" "$LINK.bam.bai"
            "$BCFTOOLS" mpileup \
                -f "$REFERENCE" \
                -R "$REGIONS" \
                -a FORMAT/AD \
                -q 20 -Q 20 -B \
                -d 100000 \
                -I \
                "$LINK.bam" \
            | "$BCFTOOLS" query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n'
            rm -f "$LINK.bam" "$LINK.bam.bai"
        fi
    } > "$OUT_FILE"
    echo "Sites counted: $(($(wc -l < "$OUT_FILE") - 1)) of $(wc -l < "$REGIONS")"
done

echo "Task complete"
