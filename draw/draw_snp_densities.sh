#!/bin/bash
# Usage:
# $TOOLS/bashare/draw_bed_counts.sh VCF ASSEMBLY_PATH [DENSITY_THRESHOLDS \"LABEL\"]

if [[ $# -ne 2 && $# -ne 4 ]]; then
    echo "Usage: $0 COUNTS_BED ASSEMBLY_PATH [DENSITY_THRESHOLDS DENSITY_MULTIPLIER]"
    exit 1
fi

VCF=$1
ASSEMBLY_PATH=$2
DENSITY_THRESHOLDS=${3:-"0,0.1,0.5,1,2,3,4,5,6,7"}  # Comma-separated list of thresholds(SNPs/kb) for SNP densities to use for window coloring
LABEL=${4:-" "}

source $(conda info --base)/etc/profile.d/conda.sh
conda activate py38

$TOOLS/MACE/scripts/draw_variant_window_densities.py \
    -i ${VCF} \
    -o ${VCF%.*.*}.w1mb.s100kb \
    -l "${LABEL}" \
    --density_thresholds "${DENSITY_THRESHOLDS}" \
    --window_size 1000000 \
    --window_step 100000 \
    --scaffold_ordered_list ${ASSEMBLY_PATH}/*.orderlist \
    --scaffold_white_list ${ASSEMBLY_PATH}/*.whitelist \
    --scaffold_length_file ${ASSEMBLY_PATH}/*.len \
    --scaffold_syn_file ${ASSEMBLY_PATH}/*.syn \
    --centromere_bed $ASSEMBLY/*.centromere.bed \
    --hide_track_label \
    --rounded

