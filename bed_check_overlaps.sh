#!/bin/bash

BED_FILE="$1"

if [[ -z "$BED_FILE" ]]; then
    echo "Usage: $0 input.bed"
    exit 1
fi

# Number of lines in the original BED file
ORIG_LINES=$(wc -l < "$BED_FILE")

# Number of lines after merge
MERGED_LINES=$(bedtools merge -i <(sort -k1,1 -k2,2n "$BED_FILE") | wc -l)

if [[ "$ORIG_LINES" -eq "$MERGED_LINES" ]]; then
    echo "No overlaps detected"
    exit 0
else
    echo "Overlapping regions found!"
    exit 1
fi
