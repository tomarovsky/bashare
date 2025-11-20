#!/bin/bash

set -euo pipefail

WINDOW_SIZE="${3:-10000}"
MOSDEPTH_BEDGZ="${1:-}"
SCAFFOLD_NAME="${2:-}"

if [[ -z "$MOSDEPTH_BEDGZ" || -z "$SCAFFOLD_NAME" ]]; then
    echo "Usage: $0 mosdepth.per-base.bed.gz scaffold_name [window_size (default: 10000)]"
    exit 1
fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate py38

DIR=$(dirname "$MOSDEPTH_BEDGZ")
BASENAME=$(basename "$MOSDEPTH_BEDGZ")
PREFIX="${BASENAME%.per-base.bed.gz}"
PAR_DIR="${DIR}/PAR"
WG_STATS="${DIR}/${PREFIX}_whole_genome_stats.csv"


# 1. Whole Genome Stats
if [[ ! -f "${WG_STATS}" ]]; then
    echo "[1/4] Computing whole-genome stats..."
    "$TOOLS/Biocrutch/scripts/Coverage/coverage_statistics.py" -i "$MOSDEPTH_BEDGZ" -g -o "${DIR}/${PREFIX}"
else
    echo "[1/4] Whole-genome stats found: ${WG_STATS}, skipping."
fi

mkdir -p "$PAR_DIR"

# 2. Extract Scaffold
echo "[2/4] Extracting scaffold ${SCAFFOLD_NAME} to ${CHR_BED}..."
CHR_BED="${PAR_DIR}/${PREFIX}.${SCAFFOLD_NAME}.bed.gz"
zgrep -w "^${SCAFFOLD_NAME}$" "${MOSDEPTH_BEDGZ}" | gzip > "${CHR_BED}"

# 3. Window Stats
echo "[3/4] Computing window stats..."
CHR_PREFIX="${PAR_DIR}/${PREFIX}.${SCAFFOLD_NAME}"
STATS_CSV="${CHR_PREFIX}_${WINDOW_SIZE}_windows_stats.csv"
"$TOOLS/Biocrutch/scripts/Coverage/coverage_statistics.py" -i "${CHR_BED}" -n -f "$WINDOW_SIZE" -o "${CHR_PREFIX}"

# Remove header
tail -n +2 "$STATS_CSV" > "${STATS_CSV}.tmp" && mv "${STATS_CSV}.tmp" "$STATS_CSV"

# 4. PAR Detection
echo "[4/4] Running PAR identification..."
MEAN_COV=$(awk 'NR==2 {print $2}' "${WG_STATS}")
"$TOOLS/Biocrutch/scripts/PAR/pseudoautosomal_region.py" \
    -f "$WINDOW_SIZE" \
    -i "$STATS_CSV" \
    -s "$SCAFFOLD_NAME" \
    -m "$MEAN_COV" \
    -o "$CHR_PREFIX" \
    | tee "${CHR_PREFIX}_pseudo.log"

echo "Done! Results are in: $PAR_DIR"
