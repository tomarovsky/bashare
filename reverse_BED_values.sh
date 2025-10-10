#!/bin/bash

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 BED_FILE REVERSE_LIST"
    echo "Example: $0 mmar.localADMIXTURE.0.2.snp400.cov75kb.mask.bed HiC_scaffold_2,HiC_scaffold_4"
    echo "NB! '_rc' was added to the scaffold name to indicate that the scaffold values were reversed."
    exit 1
fi

BED_FILE="$1"
REVERSE_LIST="$2"

# comma separated list to array
IFS=',' read -r -a REVERSE_SCAFFOLDS <<< "$REVERSE_LIST"

# create temporary directory
TMP_DIR=$(mktemp -d)

# for each scaffold in BED
for SCAFFOLD in $(awk '{print $1}' "$BED_FILE" | sort -u); do
    awk -v s="$SCAFFOLD" '$1 == s' "$BED_FILE" > "$TMP_DIR/$SCAFFOLD.bed"

    # check if scaffold needs to be reversed
    if printf '%s\n' "${REVERSE_SCAFFOLDS[@]}" | grep -qx "$SCAFFOLD"; then
        awk '{print $0}' "$TMP_DIR/$SCAFFOLD.bed" | tac | \
        awk -v s="${SCAFFOLD}_rc" 'BEGIN{OFS="\t"}{$1=s; print $0}'
    else
        cat "$TMP_DIR/$SCAFFOLD.bed"
    fi
done

# remove temporary directory
rm -r "$TMP_DIR"
