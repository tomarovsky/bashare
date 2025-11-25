#!/bin/bash

set -euo pipefail

PATHS_FILE="$1"

if [[ -z "$PATHS_FILE" ]]; then
    echo "Usage: $0 <file_with_paths>"
    exit 1
fi

while IFS= read -r path
do
    if [[ -e "$path" ]]; then
        echo "✅ Exist: $path"
    else
        echo "❌ Not exist: $path"
    fi
done < "$PATHS_FILE"
