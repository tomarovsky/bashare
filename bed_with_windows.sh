#!/bin/bash

if [ "$#" -lt 3 ] || [ "$#" -gt 5 ]; then
    echo "Usage: $0 LEN_FILE WHITELIST_FILE OUTPUT_FILE WINDOW_SIZE STEP_SIZE"
    exit 1
fi

LEN_FILE="$1"
WHITELIST_FILE="$2"
OUTPUT_FILE="$3"
WINDOW_SIZE="${4:-1000000}"
STEP_SIZE="${5:-100000}"

# Clear old output file
> "$OUTPUT_FILE"

while read -r scaffold_name length_from_len_file; do
    length_line=$(grep -Pw "^${scaffold_name}\t" "$LEN_FILE")

    if [ -n "$length_line" ]; then
        scaffold_length=$(echo "$length_line" | awk '{print $2}')

        # Generate windows
        start_coord=0
        while [ "$start_coord" -lt "$scaffold_length" ]; do
            stop_coord=$((start_coord + WINDOW_SIZE))

            if [ "$stop_coord" -gt "$scaffold_length" ]; then
                break
            fi

            printf "%s\t%d\t%d\n" "$scaffold_name" "$start_coord" "$stop_coord"

            start_coord=$((start_coord + STEP_SIZE))
        done
    fi

done < "$WHITELIST_FILE"
