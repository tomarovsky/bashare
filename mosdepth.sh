#!/bin/bash
# Usage:
# $TOOLS/bashare/mosdepth.sh BAM_FILE OUTPUT THREADS
# parallel -v --progress -j 8 "$TOOLS/bashare/mosdepth.sh {} {.} 4 > {.}.mosdepth.log 2>&1" ::: *.bam
# find . -name *.bam | parallel -v --progress -j 16 "$TOOLS/bashare/mosdepth.sh {} {.} 4 > {.}.mosdepth.log 2>&1"

set -euo pipefail

source "$TOOLS/bashare/lib/log_functions.sh"

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 BAM_FILE OUTPUT THREADS"
    exit 1
fi

BAM_FILE=$1
OUTPUT=$2
THREADS=$3

if [[ ! -f "${BAM_FILE}.bai" ]]; then
    log_warning "BAM index not found. ${BAM_FILE} indexing..."
    samtools index "$BAM_FILE"
fi

mosdepth --mapq 20 --threads $THREADS ${OUTPUT}.mosdepth.mapq20 $BAM_FILE
