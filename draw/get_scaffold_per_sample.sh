#!/bin/bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <*.local_admixture.features.bed ...> <config_dir> <scaffold>"
    echo "Example: $0 *.local_admixture.features.bed /path/to/dir/ HiC_scaffold_1"
    exit 1
fi

files=("$@")
config_dir="${files[-2]}"
scaffold="${files[-1]}"
unset 'files[-1]'
unset 'files[-1]'

# Input files
len_file="${config_dir}/*.len"
syn_file="${config_dir}/*.syn"
ord_file="${config_dir}/*.orderedlist"
white_file="${config_dir}/*.whitelist"

# Check files
for f in "$len_file" "$syn_file" "$ord_file" "$white_file"; do
    [[ -f "$f" ]] || { echo "Error: config file not found: $f"; exit 1; }
done

# Temporary and output files
out_prefix="all_samples"
out_features="${out_prefix}.local_admixture.features.bed"
out_len="${config_dir}/${scaffold}.len"
out_syn="${config_dir}/${scaffold}.syn"
out_ord="${config_dir}/${scaffold}.orderedlist"
out_white="${config_dir}/${scaffold}.whitelist"

> "$out_features"
> "$out_len"
> "$out_syn"
> "$out_ord"
> "$out_white"

# Find the name for scaffold in syntax (for MZIB.chrX)
target_chr=$(awk -v scaf="$scaffold" '$1==scaf {print $2}' "$syn_file")
if [[ -z "$target_chr" ]]; then
    echo "Error: scaffold $scaffold not found in genome.syn"
    exit 1
fi

# Process each file
for file in "${files[@]}"; do
    [[ -f "$file" ]] || continue
    sample=$(basename "$file" .local_admixture.features.bed)

    echo "Processing $sample ..."

    # get only one scaffold and add it to the beginning
    awk -v scaf="$scaffold" -v id="$sample" '$1==scaf {print id"_"$0}' "$file" >> "$out_features"

    # genome.len
    len_value=$(awk -v scaf="$scaffold" '$1==scaf {print $2}' "$len_file")
    echo -e "${sample}_${scaffold}\t${len_value}" >> "$out_len"

    # genome.whitelist
    echo "${sample}_${scaffold}" >> "$out_white"

    # genome.orderedlist
    echo "${sample}_${target_chr}" >> "$out_ord"

    # genome.syn
    echo -e "${sample}_${scaffold}\t${sample}_${target_chr}" >> "$out_syn"
done
