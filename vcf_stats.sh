#!/bin/bash
set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 VCF_FILE TRACK_FILE EXCLUDE_SCAFFOLD EXCLUDE_TRACK_CHR" >&2
    exit 1
fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate py38

VCF_FILE=$1
TRACK_FILE=$2
EXCLUDE_SCAFFOLD=$3
EXCLUDE_TRACK_CHR=$4

NUMBER_OF_SNPS=$(zcat "${VCF_FILE}" | grep -v '^#' | grep -v "^${EXCLUDE_SCAFFOLD}\t" | wc -l)

MEAN_MEDIAN=$(python3 - <<END
import pandas as pd
import sys

track_file = "$TRACK_FILE"
exclude_chr = "$EXCLUDE_TRACK_CHR"

df = pd.read_csv(track_file, sep="\t", header=0)
sample_id = df.columns[4]

# masking
df_filtered = df[(df['scaffold'] != exclude_chr) & (df['color'] != 'gray')].copy()

mean = round(df_filtered[sample_id].mean(), 2)
median = round(df_filtered[sample_id].median(), 2)

print(f"{mean}\t{median}")
END
)

echo -e "${NUMBER_OF_SNPS}\t${MEAN_MEDIAN}"
