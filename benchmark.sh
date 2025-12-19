#!/bin/bash

# Check arguments
if [ -z "$1" ]; then
    echo "Usage: $0 <command or script> [arguments...]"
    echo "Example: "
    echo "    $0 bash cmd.sh"
    echo "    $0 command arg1 arg2"
    echo "    $0 -c bash 'command1; command2 arg1 arg2;'"
    exit 1
fi

# Temporary file for metrics (remove it on exit)
TMP_METRICS=$(mktemp)
trap "rm -f $TMP_METRICS" EXIT

# --- Running script ---
/usr/bin/time -f "%e %M %P %S %U" -o "$TMP_METRICS" "$@"

read real_sec max_mem_kb cpu_percent sys_sec user_sec < "$TMP_METRICS"

# Clear memory usage (KB -> MB)
mem_mb=$(echo "scale=2; $max_mem_kb / 1024" | bc)

# Calculate approximate number of threads (cores)
cpu_percent=${cpu_percent//%/}
threads=$(echo "scale=2; ${cpu_percent} / 100" | bc)

# --- Output ---
echo "----------------------------------------"
echo -e "\n=== Benchmark Report ==="

# Format output using printf for readability
printf "%-20s : %s sec\n" "Wall Time" "$real_sec"
printf "%-20s : %s sec (User) + %s sec (Sys)\n" "CPU Time" "$user_sec" "$sys_sec"
printf "%-20s : %s MB\n" "Max RAM (RSS)" "$mem_mb"
printf "%-20s : %s%% (~%s cores)\n" "Avg CPU Load" "$cpu_percent" "$threads"
echo ""
