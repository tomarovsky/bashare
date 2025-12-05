#!/bin/bash
# Logging utilities

# --- Colors ---
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[0;33m'
readonly RED='\033[0;31m'
readonly NC='\033[0m'

# --- Universal logging function ---
# Usage: log_message <LEVEL> <MESSAGE> [LOG_FILE]
# Examples: log_message INFO "Genotyping in progress..." "pipeline.log"
#           log_message ERROR "PLINK failure"
log_message() {
    local level="$1"
    local message="$2"
    local log_file="$3"
    local color=""
    local output_stream="1" # 1 = stdout (by default), 2 = stderr (WARNING/ERROR)

    # 1. Determine color and output stream (stdout/stderr)
    case "$level" in
        INFO)
            color="${GREEN}"
            output_stream="1"
            ;;
        WARNING)
            color="${YELLOW}"
            output_stream="2"
            ;;
        ERROR)
            color="${RED}"
            output_stream="2"
            ;;
        *)
            color="${RED}"
            level="FATAL" # If unknown level, consider it a fatal error
            output_stream="2"
            ;;
    esac

    local timestamp=$(date +%d-%m-%Y\ %H:%M:%S)
    local formatted_message="[${level}] | ${timestamp} | ${message}"

    # 2. Output to terminal
    if [[ "$output_stream" == "1" ]]; then
        echo -e "${color}${formatted_message}${NC}"
    else
        echo -e "${color}${formatted_message}${NC}" >&2
    fi

    # 3. Write to log file (plain text)
    if [[ -n "$log_file" ]]; then
        echo "${formatted_message}" >> "$log_file"
    fi
}

# --- Helper functions ---
log_info() {
    log_message INFO "$1" "${2:-}"
}

log_warning() {
    log_message WARNING "$1" "${2:-}"
}

log_error() {
    log_message ERROR "$1" "${2:-}"
}
