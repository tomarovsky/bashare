#!/bin/bash
set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 VCF_FILE EXCLUDE_SCAFFOLDS[,SCAFFOLD2,...]"
    exit 1
fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate varcall

VCF_FILE=$1
EXCLUDE_SC=$2

# grep regex
EXCLUDE_REGEX=$(echo "$EXCLUDE_SC" | sed 's/,/|/g')

# Number of SNPs
zcat "$VCF_FILE" | grep -v '^#' | grep -Pvc "^($EXCLUDE_REGEX)\t"
