#!/bin/bash
# Usage:
# $TOOLS/bashare/picard_downsampling.sh INPUT_BAM OUTPUT_BAM PROBABILITY
# S=F1_MoCAN; $TOOLS/bashare/picard_downsampling.sh  ${S}/${S}.SE-PE.RG.bam  ${S}/${S}.SE-PE.RG.10x.bam  0.4 > ${S}/${S}.SE-PE.RG.10x.down.log 2>&1 &
# parallel -v --progress -j 8 "$TOOLS/bashare/picard_downsampling.sh {} {.}.22x.bam 0.1 > {.}.picard_downsampling.log 2>&1" ::: *.bam
# find . -name "*.bam" | parallel -v --progress -j 16 "$TOOLS/bashare/picard_downsampling.sh {} {.}.22x.bam 0.2 > {.}.picard_downsampling.log 2>&1"

INPUT_BAM=$1
OUTPUT_BAM=$2
P=$3

source $(conda info --base)/etc/profile.d/conda.sh
conda activate gatk

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 BAM_FILE OUTPUT PROBABILITY"
    exit 1
fi

picard DownsampleSam -Xmx20g I="${INPUT_BAM}" O="${OUTPUT_BAM}" P="${P}" 

