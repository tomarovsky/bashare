#!/bin/bash
# Total length of regions in a BED file

BED_FILE="$1"

if [[ -z "$BED_FILE" ]]; then
    echo "Usage: $0 input.bed"
    exit 1
fi

bedtools merge -i <(sort -k1,1 -k2,2n "$BED_FILE") | awk '{sum += $3 - $2} END {print sum}'
