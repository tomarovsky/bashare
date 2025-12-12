#!/bin/bash
set -euo pipefail

source "$TOOLS/bashare/lib/log_functions.sh"

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <merged_output.bed> <bed1> [bed2 ...]"
    exit 1
fi

OUT="$1"
shift
BEDS=("$@")

# Check input files
for bed in "${BEDS[@]}"; do
    if [ ! -s "$bed" ]; then
        log_error "Input BED not found or empty: $bed"
        exit 1
    fi
done

cat "${BEDS[@]}" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "$OUT"

log_info "Done! Merged BED written to: $OUT"
