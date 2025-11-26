#!/bin/bash

set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <fasta_file> <orthodb> <threads>"
    exit 1
fi

FASTA=$1
SPECIES_NAME=${FASTA%.*}
ORTHO_DB=$2
THREADS=$3

# Create directory for the species results
mkdir -p "${SPECIES_NAME}"
cd "${SPECIES_NAME}"

# Run BUSCO
busco \
    --metaeuk \
    -m genome \
    -i "../${FASTA}" \
    -c "${THREADS}" \
    -l "${ORTHO_DB}" \
    -o "${SPECIES_NAME}" \
    --offline

# Normalize BUSCO output structure
mv "${SPECIES_NAME}"/* ./
rm -r "${SPECIES_NAME}"
mv run_*/* ./
rm -r run_*

# Rename BUSCO summary outputs
mv full_table.tsv "full_table_${SPECIES_NAME}.tsv"
mv missing_busco_list.tsv "missing_busco_list_${SPECIES_NAME}.tsv"
mv short_summary.txt "short_summary_${SPECIES_NAME}.txt"
