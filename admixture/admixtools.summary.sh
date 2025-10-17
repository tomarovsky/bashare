#!/bin/bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <stat.file1> [<stat.file2> ...]"
    exit 1
fi

# Header
echo -e "Source_1\tSource_2\tTarget\tOutgroup\tStat\tZ-score\tABBA\tBABA\tSNPs"

for file in "$@"; do
    grep "^result:" "$file" | while read -r line; do
        line_clean=$(echo "$line" | sed 's/^result:[[:space:]]*//')
        echo -e "$(echo "$line_clean" | awk '{$1=$1; print}' | tr -s ' ' '\t')"
    done
done
