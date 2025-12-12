#!/bin/bash

if [ "$#" -lt 3 ] || [ "$#" -gt 5 ]; then
    echo "Usage: $0 LEN_FILE WHITELIST_FILE OUTPUT_FILE [WINDOW_SIZE] [STEP_SIZE]" >&2
    echo "  [WINDOW_SIZE] default: 1000000" >&2
    echo "  [STEP_SIZE] default: 100000" >&2
    exit 1
fi

LEN_FILE="$1"
WHITELIST_FILE="$2"
OUTPUT_FILE="$3"
WINDOW_SIZE=${4:-1000000}
STEP_SIZE=${5:-100000}

{
    while read -r scaffold_name; do
        length_line=$(grep -w "^${scaffold_name}" "$LEN_FILE")

        if [ -n "$length_line" ]; then
            scaffold_length=$(echo "$length_line" | awk '{print $2}')

            if [ "$scaffold_length" -lt "$WINDOW_SIZE" ]; then
                echo "warning: scaffold $scaffold_name (${scaffold_length}) less than (${WINDOW_SIZE}). skipped." >&2
                continue
            fi

            # generate windows
            start_coord=0
            while [ "$start_coord" -lt "$scaffold_length" ]; do
                stop_coord=$((start_coord + WINDOW_SIZE))

                if [ "$stop_coord" -gt "$scaffold_length" ]; then
                    break
                fi

                echo -e "${scaffold_name}\t${start_coord}\t${stop_coord}"

                start_coord=$((start_coord + STEP_SIZE))
            done

        else
            echo "warning: '$scaffold_name' not found in '$LEN_FILE'" >&2
        fi

    done < "$WHITELIST_FILE"

} > "$OUTPUT_FILE"
