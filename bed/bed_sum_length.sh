#!/bin/bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Total length of regions in a BED file"
    echo "Usage: $0 input.bed"
    exit 1
fi

BED_FILE="$1"

bedtools merge -i <(sort -k1,1 -k2,2n "$BED_FILE") \
    | awk '{sum += $3 - $2} END {print sum}'
