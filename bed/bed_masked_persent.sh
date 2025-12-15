#!/bin/bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Calculate masked percent of genome by BED regions"
    echo "Usage: $0 mask.bed genome.len [comma separated scaffolds for exclusion]"
    echo "Example scaffolds for exclusion: scaffold_19,chrX"
    exit 1
fi

BED="$1"
LEN="$2"
EXCLUDE="${3:-}"

# Filter by scaffolds (if a list is provided)
if [[ -n "$EXCLUDE" ]]; then
    EXCL_AWK='
        BEGIN {
            n = split(excl, a, ",")
            for (i = 1; i <= n; i++) excl_map[a[i]] = 1
        }
        !($1 in excl_map)
    '
else
    EXCL_AWK='{print}'
fi

# Total length of masked regions
TOTAL_MASKED=$(
    awk -v excl="$EXCLUDE" "$EXCL_AWK" "$BED" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    | awk '{sum += $3 - $2} END {print sum}'
)

# Total length of genome
GENOME_LENGTH=$(
    awk -v excl="$EXCLUDE" "$EXCL_AWK" "$LEN" \
    | awk '{sum += $2} END {print sum}'
)

# Percentage
MASKED_PERCENT=$(awk -v bed="$TOTAL_MASKED" -v genome="$GENOME_LENGTH" \
    'BEGIN { printf "%.2f", (bed / genome) * 100 }')

echo "${GENOME_LENGTH}\t${TOTAL_MASKED}\t${MASKED_PERCENT}"
