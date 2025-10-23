#!/bin/bash
set -euo pipefail

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 VCF_FILE VARIANT_COUNTS EXCLUDE_SCAF EXCLUDE_SYN"
    exit 1
fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate py38

VCF_FILE=$1
VARIANT_COUNTS=$2
EXCLUDE_SCAF=$3
EXCLUDE_SYN=$4

# Number of SNPs
NUM_SNPS=$(zcat "$VCF_FILE" | grep -v '^#' | grep -Pvc "^$EXCLUDE_SCAF\t")

# Mean and median based VARIANT COUNTS
read SAMPLE MEAN MEDIAN < <(
grep -Pv "^$EXCLUDE_SYN\t" "$VARIANT_COUNTS" | python3 - <<'PYTHON'
import sys
import pandas as pd

WINDOW_SIZE = 1000000
MULTIPLICATOR = 1000

df = pd.read_csv(sys.stdin, sep='\t', header=0)
sample = df.columns[2]

df[sample] = df[sample] / WINDOW_SIZE * MULTIPLICATOR

mean = df[sample].mean()
median = df[sample].median()

print(sample, mean, median)
PYTHON
)

echo -e "${SAMPLE}\t${NUM_SNPS}\t${MEAN}\t${MEDIAN}"
