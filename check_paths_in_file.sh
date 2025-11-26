#!/bin/bash

set -euo pipefail

if [[ "$#" -ne 1 ]]; then
    echo "Usage: $0 <file_with_paths>"
    exit 1
fi

PATHS_FILE="$1"

echo "----------------"
while IFS= read -r path
do
    if [[ -e "$path" ]]; then
        echo "Exist: $path"
    else
        echo "NOT exist: $path"
    fi
done < "$PATHS_FILE"
echo "----------------"
