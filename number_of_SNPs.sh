#!/bin/bash
set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 VCF_FILE MASK_BED THRESHOLD EXCLUDE_SCAFFOLD1[,SCAFFOLD2,...]"
    exit 1
fi

VCF_FILE=$1
MASK_BED=$2
THRESHOLD=$3
EXCLUDE_SC=$4

# grep regex
EXCLUDE_REGEX=$(echo "$EXCLUDE_SC" | sed 's/,/|/g')

# Number of SNPs
zcat "$VCF_FILE" | \
grep -v '^#' | \
grep -Pv "^($EXCLUDE_REGEX)\t" | \
awk '{print $1"\t"$2-1"\t"$2}' | \
bedtools intersect -v -a - -b <(awk '$4>'"$THRESHOLD"'' "$MASK_BED") | \
wc -l
