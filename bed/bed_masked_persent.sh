#!/bin/bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Calculate genome coverage by BED"
    echo "Usage: $0 input.bed genome.len [exclude_scaffold1 exclude_scaffold2 ...]"
    echo "Output: genome_length\tcovered_length\tmasked_percentage"
    exit 1
fi

BED_FILE="$1"
LEN_FILE="$2"
shift 2
EXCLUDE_LIST=("$@")

filter_scaffolds() {
    if [ ${#EXCLUDE_LIST[@]} -eq 0 ]; then
        cat
    else
        grep -v -F -w -f <(printf "%s\n" "${EXCLUDE_LIST[@]}")
    fi
}

GENOME_LEN=$(filter_scaffolds < "$LEN_FILE" | awk '{sum += $2} END {print sum+0}')

REGIONS_LEN=$(filter_scaffolds < "$BED_FILE" | \
              sort -k1,1 -k2,2n | \
              bedtools merge -i - | \
              awk '{sum += $3 - $2} END {print sum}')

MASKED_PERCENT=$(awk -v total="$GENOME_LEN" -v covered="$REGIONS_LEN" \
    'BEGIN { printf "%.2f\n", (covered/total)*100 }')

echo -e "$GENOME_LEN\t$REGIONS_LEN\t$MASKED_PERCENT"
