#!/bin/bash

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 VCF ASSEMBLY_BASEPATH [DENSITY_THRESHOLDS] [LABEL] [EXTRA_OPTIONS...]"
    echo "Example: $0 data.vcf path/to/assembly_basepath -l 'MyLabel' --density_thresholds "0,0.1,0.5,1,2,3,4,5,6,7" --centromere_bed $ASSEMBLY/*.centromere.bed"
    exit 1
fi

VCF=$1
ASSEMBLY_BASEPATH=$2
shift 2 # Remove VCF and ASSEMBLY_BASEPATH from the list of arguments ($@)

# Add the remaining arguments to EXTRA_OPTIONS
EXTRA_OPTIONS="$@"

$TOOLS/MACE/scripts/draw_variant_window_densities.py \
    -i ${VCF} \
    -o ${VCF%.*.*}.w1mb.s100kb \
    --window_size 1000000 \
    --window_step 100000 \
    --scaffold_ordered_list ${ASSEMBLY_BASEPATH}.orderlist \
    --scaffold_white_list ${ASSEMBLY_BASEPATH}.whitelist \
    --scaffold_length_file ${ASSEMBLY_BASEPATH}.len \
    --scaffold_syn_file ${ASSEMBLY_BASEPATH}.syn \
    --hide_track_label \
    --rounded \
    ${EXTRA_OPTIONS}
