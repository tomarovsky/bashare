#!/bin/bash
# Usage:
# $TOOLS/bashare/mask_from_coverage.sh PER_BASE_BED_FILE WHOLE_GENOME_STATS MALES HEMI_REGION_BED [MIN_THRESHOLD] [MAX_THRESHOLD]
# $TOOLS/bashare/mask_from_coverage.sh S49.mzib.hic.purged.mkdup.mapq20.per-base.bed.gz whole_genome_stats.csv males.txt hemi_regions.bed 0.33 2.5

set -euo pipefail

if [[ $# -lt 4 || $# -gt 6 ]]; then
    echo "Usage: $0 PER_BASE_BED_FILE WHOLE_GENOME_STATS MALES HEMI_REGION_BED [MIN_THRESHOLD] [MAX_THRESHOLD]"
    exit 1
fi

PER_BASE_BED_FILE=$1 # mosdepth.per-base.bed.gz
COVERAGE=$(cat $2 | sed -n 2p | awk '{print $2}') # whole_genome_stats.csv
MALES=$3 # sample per line
HEMI_REGION_BED=$4 # BED with hemizygous regions
MIN_THRESHOLD=${5:-0.33}
MAX_THRESHOLD=${6:-2.5}

SAMPLE=$(basename "$PER_BASE_BED_FILE" | cut -d. -f1)

source $(conda info --base)/etc/profile.d/conda.sh
conda activate py38

if grep -qw "$SAMPLE" "$MALES"; then
    echo "[INFO] $SAMPLE == MALE"

    bedtools intersect -a <(zcat "$PER_BASE_BED_FILE") -b "$HEMI_REGION_BED" | gzip > "${SAMPLE}.hemi.tmp.per-base.bed.gz" &
    bedtools intersect -v -a <(zcat "$PER_BASE_BED_FILE") -b "$HEMI_REGION_BED" | gzip > "${SAMPLE}.diploid.tmp.per-base.bed.gz" &

    wait

    "$TOOLS/MAVR/scripts/alignment/coverage/generate_mask_from_coverage_bed.py" \
        -c ${SAMPLE}.diploid.tmp.per-base.bed.gz \
        -m "$COVERAGE" \
        --max_coverage_threshold "$MAX_THRESHOLD" \
        --min_coverage_threshold "$MIN_THRESHOLD" \
        -o ${SAMPLE}.diploid.per-base.mask.bed &
    "$TOOLS/MAVR/scripts/alignment/coverage/generate_mask_from_coverage_bed.py" \
        -c ${SAMPLE}.hemi.tmp.per-base.bed.gz \
        -m "$(awk -v cov=$COVERAGE 'BEGIN{print cov/2}')" \
        --max_coverage_threshold "$MAX_THRESHOLD" \
        --min_coverage_threshold "$MIN_THRESHOLD" \
        -o ${SAMPLE}.hemi.per-base.mask.bed &

    wait

    cat ${SAMPLE}.diploid.per-base.mask.bed ${SAMPLE}.hemi.per-base.mask.bed | sort -S 20G -k1,1 -k2,2n > "${PER_BASE_BED_FILE%.*.*}.max${MAX_THRESHOLD}.min${MIN_THRESHOLD}.bed"
    rm ${SAMPLE}.diploid.tmp.per-base.bed.gz ${SAMPLE}.hemi.tmp.per-base.bed.gz ${SAMPLE}.diploid.per-base.mask.bed ${SAMPLE}.hemi.per-base.mask.bed

else
    echo "[INFO] $SAMPLE == FEMALE"

    "$TOOLS/MAVR/scripts/alignment/coverage/generate_mask_from_coverage_bed.py" \
        -c "$PER_BASE_BED_FILE" \
        -m "$COVERAGE" \
        --max_coverage_threshold "$MAX_THRESHOLD" \
        --min_coverage_threshold "$MIN_THRESHOLD" \
        -o "${PER_BASE_BED_FILE%.*.*}.max${MAX_THRESHOLD}.min${MIN_THRESHOLD}.bed"
fi

echo "[INFO] ${SAMPLE}: Done"
