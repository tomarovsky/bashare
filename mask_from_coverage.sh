#!/bin/bash
# Usage:
# $TOOLS/bashare/mask_from_coverage.sh PER_BASE_BED_FILE COVERAGE MALES "HEMI_REGION_COORDS"

set -euo pipefail

PER_BASE_BED_FILE=$1 # mosdepth
COVERAGE=$2 # $(cat *.mapq20_whole_genome_stats.csv | sed -n 2p | awk '{print $2}')
MALES=$3 # txt file
HEMI_REGION_COORDS=$4 # HiC_scaffold_19:6680001-124421298

SAMPLE=$(basename "$PER_BASE_BED_FILE" | cut -d. -f1)

source $(conda info --base)/etc/profile.d/conda.sh

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 PER_BASE_BED_FILE COVERAGE MALES 'HEMI_REGION_COORDS'"
    exit 1
fi

conda activate py38

# male of female
if grep -qw "$SAMPLE" "$MALES"; then
    echo "[INFO] $SAMPLE == FEMALE"

    # scaffold:start-end
    CHR=$(echo "$HEMI_REGION_COORDS" | cut -d: -f1)
    START=$(echo "$HEMI_REGION_COORDS" | cut -d: -f2 | cut -d- -f1)
    END=$(echo "$HEMI_REGION_COORDS" | cut -d: -f2 | cut -d- -f2)

    # split per-base.bed.gz
    zcat "$PER_BASE_BED_FILE" | awk -v chr="$CHR" -v start="$START" -v end="$END" \
        '($1==chr && $2>=start && $3<=end) {print > "${SAMPLE}.hemi.tmp.bed"} 
         !($1==chr && $2>=start && $3<=end) {print > "${SAMPLE}.diploid.tmp.bed"}'

    "$TOOLS/MAVR/scripts/alignment/coverage/generate_mask_from_coverage_bed.py" \
        -c ${SAMPLE}.diploid.tmp.bed -m "$COVERAGE" --max_coverage_threshold 2.5 --min_coverage_threshold 0.33 \
        -o ${SAMPLE}.diploid.mask.bed

    "$TOOLS/MAVR/scripts/alignment/coverage/generate_mask_from_coverage_bed.py" \
        -c ${SAMPLE}.hemi.tmp.bed -m "$(awk -v c=$COVERAGE 'BEGIN{print c/2}')" \
        --max_coverage_threshold 2.5 --min_coverage_threshold 0.33 \
        -o ${SAMPLE}.hemi.mask.bed

    cat ${SAMPLE}.diploid.mask.bed ${SAMPLE}.hemi.mask.bed | sort -k1,1 -k2,2n > "${SAMPLE}.mask.max250.min33.bed"
    # rm -f ${SAMPLE}.diploid.tmp.bed ${SAMPLE}.hemi.tmp.bed ${SAMPLE}.diploid.mask.bed ${SAMPLE}.hemi.mask.bed

else
    echo "[INFO] $SAMPLE == MALE"

    "$TOOLS/MAVR/scripts/alignment/coverage/generate_mask_from_coverage_bed.py" \
        -c "$PER_BASE_BED_FILE" -m "$COVERAGE" --max_coverage_threshold 2.5 --min_coverage_threshold 0.33 \
        -o "${SAMPLE}.mask.max250.min33.bed"
fi

