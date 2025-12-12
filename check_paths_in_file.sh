#!/bin/bash
set -euo pipefail

source "$TOOLS/bashare/lib/log_functions.sh"

if [[ "$#" -ne 1 ]]; then
    echo "Usage: $0 <file_with_paths>"
    exit 1
fi

PATHS_FILE="$1"

log_info "----------------"
while IFS= read -r path
do
    if [[ -e "$path" ]]; then
        log_info "Exist: $path"
    else
        log_warning "NOT exist: $path"
    fi
done < "$PATHS_FILE"
log_info "----------------"
