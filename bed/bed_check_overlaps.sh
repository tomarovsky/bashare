#!/bin/bash
set -euo pipefail

source "$TOOLS/bashare/lib/log_functions.sh"


if [ "$#" -ne 1 ]; then
    echo "Usage: $0 BED_FILE"
    exit 1
fi

BED_FILE="$1"

# Number of lines in the original BED file
ORIG_LINES=$(wc -l < "$BED_FILE")

# Number of lines after merge
MERGED_LINES=$(bedtools merge -i <(sort -k1,1 -k2,2n "$BED_FILE") | wc -l)

if [[ "$ORIG_LINES" -eq "$MERGED_LINES" ]]; then
    log_info "No overlaps detected"
else
    log_warning "Overlapping regions found!"
fi
