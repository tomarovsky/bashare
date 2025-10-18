#!/bin/bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <FEATURES_BED_FILES> GENOME_BASE SCAFFOLD"
    echo "Example: $0 *.features.bed /path/to/assembly/mzib HiC_scaffold_19"
    exit 1
fi

# Parse arguments safely
BED_FILES=("${@:1:$#-2}")     # all except last two
GENOME_BASE="${@: -2:1}"      # argument before last
SCAFFOLD="${@: -1}"           # last argument

# Input files
GENOME_PREFIX=$(basename "$GENOME_BASE")
LEN_FILE="${GENOME_BASE}.len"
SYN_FILE="${GENOME_BASE}.syn"
ORDEREDLIST_FILE="${GENOME_BASE}.orderedlist"

# Check existence
for f in "$LEN_FILE" "$SYN_FILE" "$ORDEREDLIST_FILE"; do
    [[ -f "$f" ]] || { echo "Error: $f not found"; exit 1; }
done

# Outfiles
OUT_PREFIX="all_samples.${SCAFFOLD}"
OUT_BED="${OUT_PREFIX}.features.bed"
OUT_LEN="${GENOME_PREFIX}.${SCAFFOLD}.len"
OUT_WHITELIST="${GENOME_PREFIX}.${SCAFFOLD}.whitelist"
OUT_ORDERED="${GENOME_PREFIX}.${SCAFFOLD}.orderedlist"
OUT_SYN="${GENOME_PREFIX}.${SCAFFOLD}.syn"

# Clear previous files if exist
> "$OUT_BED"
> "$OUT_LEN"
> "$OUT_WHITELIST"
> "$OUT_ORDERED"
> "$OUT_SYN"

# Get chromosome name from genome.syn
CHR_NAME=$(awk -v s="$SCAFFOLD" '$1==s {print $2}' "$SYN_FILE" | head -n1)
if [[ -z "$CHR_NAME" ]]; then
    echo "Error: scaffold $SCAFFOLD not found in $SYN_FILE"
    exit 1
fi

# Get length of scaffold from genome.len
SCAFFOLD_LEN=$(awk -v s="$SCAFFOLD" '$1==s {print $2}' "$LEN_FILE" | head -n1)
if [[ -z "$SCAFFOLD_LEN" ]]; then
    echo "Error: scaffold $SCAFFOLD not found in $LEN_FILE"
    exit 1
fi

# Prepare scaffold-specific files
for file in "${BED_FILES[@]}"; do
    filename=$(basename "$file")
    sample="${filename%%.*}"

    awk -v s="$SCAFFOLD" -v id="$sample" '$1==s {print id"."$0}' "$file" >> "$OUT_BED"

    echo -e "${sample}.${SCAFFOLD}\t${SCAFFOLD_LEN}" >> "$OUT_LEN"
    echo "${sample}.${SCAFFOLD}" >> "$OUT_WHITELIST"
    echo "${sample}.${CHR_NAME}" >> "$OUT_ORDERED"
    echo -e "${sample}.${SCAFFOLD}\t${sample}.${CHR_NAME}" >> "$OUT_SYN"
done

echo "[INFO] Done!"
echo "Created files:"
echo "  $OUT_BED"
echo "  $OUT_LEN"
echo "  $OUT_WHITELIST"
echo "  $OUT_ORDERED"
echo "  $OUT_SYN"
