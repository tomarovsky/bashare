#!/bin/bash
# Usage:
# $TOOLS/bashare/draw_bed_counts.sh COUNTS_BED ASSEMBLY_PATH [DENSITY_THRESHOLDS DENSITY_MULTIPLIER]
# COUNTS_BED: *.w1mb.s100kb.counts.bed
# HiC_scaffold_1    0       1000000 213412
# HiC_scaffold_1    100000  1100000 210689
# HiC_scaffold_1    200000  1200000 204796

if [[ $# -ne 2 && $# -ne 4 ]]; then
    echo "Usage: $0 COUNTS_BED ASSEMBLY_PATH [DENSITY_THRESHOLDS DENSITY_MULTIPLIER]"
    exit 1
fi

COUNTS_BED=$1
ASSEMBLY_PATH=$2
DENSITY_THRESHOLDS=${3:-"0.0,0.1,0.5,0.75,1.0,1.25,1.5,2.0,2.5"}  # Comma-separated list of thresholds(SNPs/kb) for SNP densities to use for window coloring
DENSITY_MULTIPLIER=${4:-1000}  # Default: 1000, i.e. densities will be calculated per kbp

source $(conda info --base)/etc/profile.d/conda.sh
conda activate py38

$TOOLS/MACE/scripts/draw_variant_window_densities.py \
    -i ${COUNTS_BED} \
    -o ${COUNTS_BED%.*} \
    -t bedgraph \
    --scaffold_ordered_list ${ASSEMBLY_PATH}/*.orderlist \
    --scaffold_white_list ${ASSEMBLY_PATH}/*.whitelist \
    --scaffold_length_file ${ASSEMBLY_PATH}/*.len \
    --scaffold_syn_file ${ASSEMBLY_PATH}/*.syn \
    --centromere_bed ${ASSEMBLY_PATH}/*.centromere.bed \
    --hide_track_label \
    --rounded \
    --density_thresholds="${DENSITY_THRESHOLDS}" \
    --density_multiplier ${DENSITY_MULTIPLIER} \
    --window_size 1000000 \
    --window_step 100000
