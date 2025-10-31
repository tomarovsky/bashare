#!/bin/bash
# Usage:
# $TOOLS/bashare/bed_to_window_counts.sh MASK_BED WINDOW_BED

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 MASK_BED WINDOW_BED"
    exit 1
fi

MASK_BED=$1
WINDOW_BED=$2 # e.g., w1mb.s100kb.features.bed

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate varcall

# Step 1: bedtools intersect
INTERSECT_FILE="${MASK_BED%.*}.w1mb.s100kb.intersect_wao.bed"
bedtools intersect -a "$WINDOW_BED" -b "$MASK_BED" -wao > "$INTERSECT_FILE"

# Step 2: Grouping and counting
conda deactivate && conda activate py38
OUTPUT_FILE="${MASK_BED%.*}.w1mb.s100kb.counts.bed"

python3 - "$INTERSECT_FILE" "$OUTPUT_FILE" <<EOF
import pandas as pd
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(input_file, sep='\t', header=None)
df.columns = ['scaffold1', 'start', 'end', 'scaffold2', 'start2', 'end2', 'value']
result = df.groupby(['scaffold1', 'start', 'end'])['value'].sum().reset_index()
result.to_csv(output_file, sep='\t', index=False, header=False)

EOF

rm "$INTERSECT_FILE"
