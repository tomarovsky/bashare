#!/bin/bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <stat.file1> [<stat.file2> ...]"
    echo "Example: $0 $(ls admixtools/F3stat.*.txt | tr '\n' ' ')"
    exit 1
fi

# Header
echo -e "Source_1,Source_2,Target,F3,stderr,Z-score,SNPs"
echo -e "Target,Source_1,Source_2,Outgroup,D/F4,stderr,Z-score,ABBA,BABA,SNPs"

for file in "$@"; do
    grep "result:" "$file" | while read -r line; do
        line_clean=$(echo "$line" | sed 's/result:[[:space:]]*//')
        echo -e "$(echo "$line_clean" | awk '{$1=$1; print}' | tr -s ' ' ',')"
    done
done
