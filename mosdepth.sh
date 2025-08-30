#!/bin/bash
# Usage:
# ./mosdepth.sh BAM_FILE OUTPUT THREADS
# parallel -j 8 "./mosdepth.sh {} {.} 4 > {.}.mosdepth.log 2>&1" ::: *.bam
# parallel --results logs -j 8 ./mosdepth.sh {} {.} 4 ::: *.bam

BAM_FILE=$1
OUTPUT=$2
THREADS=$3

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 BAM_FILE OUTPUT THREADS"
    exit 1
fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate mosdepth

mosdepth --mapq 20 --threads $THREADS ${OUTPUT}.mosdepth.mapq20 $BAM_FILE

