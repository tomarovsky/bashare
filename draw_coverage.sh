#!/bin/bash
# Usage:
# $TOOLS/bashare/draw_coverage.sh WINDOWS_STATS_CSV ASSEMBLY_PREFIX_PATH
# find . -name *_windows_stats.csv | parallel -v --progress -j 64 "$TOOLS/bashare/draw_coverage.sh {} ../../assembly/genome"

WINDOWS_STATS_CSV=$1
ASSEMBLY_PREFIX_PATH=$2

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 WINDOWS_STATS_CSV ASSEMBLY_PREFIX_PATH"
    exit 1
fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate py38

$TOOLS/latest_mahajrod/scripts/draw_coverage.py \
    -i ${WINDOWS_STATS_CSV} \
    -o ${WINDOWS_STATS_CSV%.*}.track \
    -l ${WINDOWS_STATS_CSV%.*} \
    -m $(cat ${WINDOWS_STATS_CSV%_*_*_*}_whole_genome_stats.csv | sed -n 2p | awk '{print $2}') \
    --scaffold_length_file ${ASSEMBLY_PREFIX_PATH}.len \
    --scaffold_white_list ${ASSEMBLY_PREFIX_PATH}.whitelist \
    --scaffold_syn_file ${ASSEMBLY_PREFIX_PATH}.syn \
    --scaffold_column_name '#scaffold' \
    --coverage_column_name 'median' \
    --window_column_name 'frame' \
    --subplots_adjust_left 0.12 \
    --figure_width 12 \
    --hide_track_label \
    --rounded;


