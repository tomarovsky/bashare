#!/bin/bash
# Usage:
# $TOOLS/bashare/mask_from_coverage.sh PER_BASE_BED_FILE COVERAGE MALES HEMI_REGIONS_BED
# $TOOLS/bashare/mask_from_coverage.sh S49.mzib.hic.purged.mkdup.mapq10.per-base.bed.gz $(cat S49/*.mapq10_whole_genome_stats.csv | sed -n 2p | awk '{print $2}') males.txt hemi_regions.bed

set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 PER_BASE_BED_FILE COVERAGE MALES 'HEMI_REGION_COORDS'"
    exit 1
fi

PER_BASE_BED_FILE=$1 # mosdepth.per-base.bed.gz
COVERAGE=$2 # $(cat *.mapq20_whole_genome_stats.csv | sed -n 2p | awk '{print $2}')
MALES=$3 # sample per line
HEMI_REGIONS_BED=$4 # BED with hemizygous regions

SAMPLE=$(basename "$PER_BASE_BED_FILE" | cut -d. -f1)

source $(conda info --base)/etc/profile.d/conda.sh

if grep -qw "$SAMPLE" "$MALES"; then
    echo "[INFO] $SAMPLE == MALE"

    HEMI_TMP="${SAMPLE}.hemi.tmp.per-base.bed.gz"
    DIP_TMP="${SAMPLE}.diploid.tmp.per-base.bed.gz"

    conda activate varcall
    bedtools intersect -a <(zcat "$PER_BASE_BED_FILE") -b "$HEMI_REGIONS_BED" | gzip > "$HEMI_TMP"
    bedtools intersect -v -a <(zcat "$PER_BASE_BED_FILE") -b "$HEMI_REGIONS_BED" | gzip > "$DIP_TMP"
    conda deactivate

    conda activate py38
    "$TOOLS/MAVR/scripts/alignment/coverage/generate_mask_from_coverage_bed.py" \
        -c ${SAMPLE}.diploid.tmp.per-base.bed.gz -m "$COVERAGE" \
        --max_coverage_threshold 2.5 --min_coverage_threshold 0.33 \
        -o ${SAMPLE}.diploid.per-base.mask.bed

    "$TOOLS/MAVR/scripts/alignment/coverage/generate_mask_from_coverage_bed.py" \
        -c ${SAMPLE}.hemi.tmp.per-base.bed.gz -m "$(awk -v cov=$COVERAGE 'BEGIN{print cov/2}')" \
        --max_coverage_threshold 2.5 --min_coverage_threshold 0.33 \
        -o ${SAMPLE}.hemi.per-base.mask.bed

    cat ${SAMPLE}.diploid.per-base.mask.bed ${SAMPLE}.hemi.per-base.mask.bed | sort -k1,1 -k2,2n > "${PER_BASE_BED_FILE%.*.*}.max250.min33.bed"
    # rm -f ${SAMPLE}.diploid.tmp.per-base.bed.gz ${SAMPLE}.hemi.tmp.per-base.bed.gz ${SAMPLE}.diploid.per-base.mask.bed ${SAMPLE}.hemi.per-base.mask.bed

else
    echo "[INFO] $SAMPLE == FEMALE"

    "$TOOLS/MAVR/scripts/alignment/coverage/generate_mask_from_coverage_bed.py" \
        -c "$PER_BASE_BED_FILE" -m "$COVERAGE" --max_coverage_threshold 2.5 --min_coverage_threshold 0.33 \
        -o "${PER_BASE_BED_FILE%.*.*}.max250.min33.bed"
fi

