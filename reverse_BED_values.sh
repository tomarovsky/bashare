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

# Comma separated list to array
IFS=',' read -r -a REVERSE_SCAFFOLDS <<< "$REVERSE_LIST"

# Create temporary directory
TMP_DIR=$(mktemp -d)

# for each scaffold in BED
for SCAFFOLD in $(awk '{print $1}' "$BED_FILE" | sort -u); do
    awk -v s="$SCAFFOLD" '$1 == s' "$BED_FILE" > "$TMP_DIR/$SCAFFOLD.bed"

    if printf '%s\n' "${REVERSE_SCAFFOLDS[@]}" | grep -qx "$SCAFFOLD"; then
        # Get 4 column values
        mapfile -t values < <(awk '{print $4}' "$TMP_DIR/$SCAFFOLD.bed" | tac)
        # Process the reversed values
        awk -v s="${SCAFFOLD}_rc" -v vals="$(printf '%s ' "${values[@]}")" '
            BEGIN{
                split(vals, v, " ");
                i = 1;
                OFS="\t";
            }
            {
                $1 = s;
                $4 = v[i++];
                print $0;
            }' "$TMP_DIR/$SCAFFOLD.bed"
    else
        cat "$TMP_DIR/$SCAFFOLD.bed"
    fi
done

# Remove temporary directory
rm -r "$TMP_DIR"
