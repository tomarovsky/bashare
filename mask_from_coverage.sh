#!/bin/bash
# Usage:
# $TOOLS/bashare/mask_from_coverage.sh PER_BASE_BED_FILE COVERAGE MALES "HEMI_REGION_COORDS"

set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 PER_BASE_BED_FILE COVERAGE MALES 'HEMI_REGION_COORDS'"
    exit 1
fi

PER_BASE_BED_FILE=$1 # mosdepth.per-base.bed.gz
COVERAGE=$2 # $(cat *.mapq20_whole_genome_stats.csv | sed -n 2p | awk '{print $2}')
MALES=$3 # txt file
HEMI_REGION_COORDS=$4 # HiC_scaffold_19:6680001-124421298

SAMPLE=$(basename "$PER_BASE_BED_FILE" | cut -d. -f1)

source $(conda info --base)/etc/profile.d/conda.sh
conda activate py38

if grep -qw "$SAMPLE" "$MALES"; then
    echo "[INFO] $SAMPLE == MALE"

    # scaffold:start-end
    CHR=$(echo "$HEMI_REGION_COORDS" | cut -d: -f1)
    START=$(echo "$HEMI_REGION_COORDS" | cut -d: -f2 | cut -d- -f1)
    END=$(echo "$HEMI_REGION_COORDS" | cut -d: -f2 | cut -d- -f2)

    # split per-base.bed.gz
    zcat "$PER_BASE_BED_FILE" | awk -v chr="$CHR" -v start="$START" -v end="$END" \
        -v hemi="${SAMPLE}.hemi.tmp.per-base.bed.gz" -v dip="${SAMPLE}.diploid.tmp.per-base.bed.gz" '
        ($1==chr && $2>=start && $3<=end) {print | "gzip -c > " hemi}
        !($1==chr && $2>=start && $3<=end) {print | "gzip -c > " dip}'

    "$TOOLS/MAVR/scripts/alignment/coverage/generate_mask_from_coverage_bed.py" \
        -c ${SAMPLE}.diploid.tmp.per-base.bed.gz -m "$COVERAGE" --max_coverage_threshold 2.5 --min_coverage_threshold 0.33 \
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

