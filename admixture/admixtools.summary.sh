#!/bin/bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <stat.file1> [<stat.file2> ...]"
    echo "Example: $0 $(ls admixtools/F3stat.*.txt | tr '\n' ' ')"
    exit 1
fi

# Header
echo -e "F3:Source_1\tSource_2\tTarget\tF3\tstderr\tZ-score\tSNPs"
echo -e "D/F4:Target\tSource_1\tSource_2\tOutgroup\tD/F4\tZ-score\tABBA\tBABA\tSNPs"

for file in "$@"; do
    grep "result:" "$file" | while read -r line; do
        line_clean=$(echo "$line" | sed 's/result:[[:space:]]*//')
        echo -e "$(echo "$line_clean" | awk '{$1=$1; print}' | tr -s ' ' '\t')"
    done
done
