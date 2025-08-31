#!/bin/bash
# Usage:
# $TOOLS/bashare/coverage_stats.sh MOSDEPTH_BED
# parallel -v --progress -j 16 "$TOOLS/bashare/coverage_stats.sh {} > {}.coverage_stats.log 2>&1" ::: *.per-base.bed.gz
# find . -name *.per-base.bed.gz | parallel -v --progress -j 64 "$TOOLS/bashare/coverage_stats.sh {} > {}.coverage_stats.log 2>&1"

set -e -u

MOSDEPTH_BED=$1

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 MOSDEPTH_BED"
    exit 1
fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate python3.8

$TOOLS/Biocrutch/scripts/Coverage/coverage_statistics.py -i $MOSDEPTH_BED -g -o ${MOSDEPTH_BED%.*.*.*}
$TOOLS/Biocrutch/scripts/Coverage/coverage_statistics.py -i $MOSDEPTH_BED -n -f 1000000 -o ${MOSDEPTH_BED%.*.*.*}

