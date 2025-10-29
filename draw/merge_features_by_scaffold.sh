#!/bin/bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <FEATURES_BED_FILES> GENOME_BASE SCAFFOLD"
    echo "Example: $0 *.features.bed /path/to/assembly/mzib HiC_scaffold_19"
    echo "Example (with mask): $0 *.features.bed /path/to/assembly/mzib HiC_scaffold_19 mask.bed"
    exit 1
fi

# Parse arguments
if [[ $# -ge 4 ]]; then
    MASK="${@: -1}"                # last argument — MASK
    SCAFFOLD="${@: -2:1}"          # argument before last — SCAFFOLD
    GENOME_BASE="${@: -3:1}"       # third from last — GENOME_BASE
    BED_FILES=("${@:1:$#-3}")      # all before them — BED files
else
    MASK=""                        # no mask
    SCAFFOLD="${@: -1}"            # last — SCAFFOLD
    GENOME_BASE="${@: -2:1}"       # second last — GENOME_BASE
    BED_FILES=("${@:1:$#-2}")      # all before them — BED files
fi

# Input files
GENOME_PREFIX=$(basename "$GENOME_BASE")
LEN_FILE="${GENOME_BASE}.len"
SYN_FILE="${GENOME_BASE}.syn"
ORDEREDLIST_FILE="${GENOME_BASE}.orderedlist"
CENTROMERE_BED="${GENOME_BASE}.centromere.bed"

# Check existence
for f in "$LEN_FILE" "$SYN_FILE" "$ORDEREDLIST_FILE" "$CENTROMERE_BED"; do
    [[ -f "$f" ]] || { echo "Error: $f not found"; exit 1; }
done

if [[ -n "$MASK" && ! -f "$MASK" ]]; then
    echo "Error: mask.bed '$MASK' not found"
    exit 1
fi

# Outfiles
OUT_PREFIX="all_samples.${SCAFFOLD}"
OUT_BED="${OUT_PREFIX}.features.bed"
OUT_LEN="${GENOME_PREFIX}.${SCAFFOLD}.len"
OUT_WHITELIST="${GENOME_PREFIX}.${SCAFFOLD}.whitelist"
OUT_ORDERED="${GENOME_PREFIX}.${SCAFFOLD}.orderedlist"
OUT_SYN="${GENOME_PREFIX}.${SCAFFOLD}.syn"
OUT_CENTROMERE_BED="${GENOME_PREFIX}.${SCAFFOLD}.centromere.bed"

# Clear previous files if exist
> "$OUT_BED"
> "$OUT_LEN"
> "$OUT_WHITELIST"
> "$OUT_ORDERED"
> "$OUT_SYN"
> "$OUT_CENTROMERE_BED"

# Get chromosome name from genome.syn
CHR_NAME=$(awk -v s="$SCAFFOLD" '$1==s {print $2}' "$SYN_FILE" | head -n1)
if [[ -z "$CHR_NAME" ]]; then
    echo "Error: $SCAFFOLD not found in $SYN_FILE"
    exit 1
fi

# Get length of scaffold from genome.len
SCAFFOLD_LEN=$(awk -v s="$SCAFFOLD" '$1==s {print $2}' "$LEN_FILE")
if [[ -z "$SCAFFOLD_LEN" ]]; then
    echo "Error: $SCAFFOLD not found in $LEN_FILE"
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
    awk -v s="$CHR_NAME" -v id="$sample" '$1==s {print id"."$0}' "$CENTROMERE_BED" >> "$OUT_CENTROMERE_BED"
done

echo "Created files:"
echo "  $OUT_BED"
echo "  $OUT_LEN"
echo "  $OUT_WHITELIST"
echo "  $OUT_ORDERED"
echo "  $OUT_SYN"

# Prepare mask file if provided
if [[ -n "$MASK" ]]; then
    OUT_MASK="${GENOME_PREFIX}.${SCAFFOLD}.mask.bed"
    > "$OUT_MASK"
    for file in "${BED_FILES[@]}"; do
        filename=$(basename "$file")
        sample="${filename%%.*}"
        awk -v s="$SCAFFOLD" -v id="$sample" '$1==s {print id"."$0}' "$MASK" >> "$OUT_MASK"
    done
fi

echo "  $OUT_MASK"
