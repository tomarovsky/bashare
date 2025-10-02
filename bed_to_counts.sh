#!/bin/bash
# Usage:
# $TOOLS/bashare/bed_to_window_counts.sh MASK_BED WINDOW_BED

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 MASK_BED WINDOW_BED"
    exit 1
fi

MASK_BED=$1
WINDOW_BED=$2 # w1mb.s100kb.features.bed

source $(conda info --base)/etc/profile.d/conda.sh
conda activate varcall

# step 1: intersect
bedtools intersect -a $WINDOW_BED -b $MASK_BED -wao > ${MASK_BED%.*}.w1mb.s100kb.intersect_wao.bed

# step 2: groupby
conda activate python3.8
$TOOLS/Biocrutch/scripts/Convert/intersect_wao_groupby.py ${MASK_BED%.*}.w1mb.s100kb.intersect_wao.bed ${MASK_BED%.*}.w1mb.s100kb.counts.bed
rm ${MASK_BED%.*}.w1mb.s100kb.intersect_wao.bed 

