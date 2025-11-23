#!/bin/bash

set -euo pipefail

usage() {
    echo "Usage: $0 PER_BASE_BED_FILE WHOLE_GENOME_STATS MALES HEMI_REGION_BED [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --min_relative FLOAT    Minimum coverage threshold (relative to mean)"
    echo "  --max_relative FLOAT    Maximum coverage threshold (relative to mean)"
    echo "  --min_absolute FLOAT    Minimum coverage threshold (absolute value)"
    echo "  --max_absolute FLOAT    Maximum coverage threshold (absolute value)"
    echo ""
    echo "Examples:"
    echo "  $0 file.bed stats.csv males.txt hemi.bed --min_relative 0.33 --max_relative 2.5"
    echo "  $0 file.bed stats.csv males.txt hemi.bed --min_absolute 5 --max_relative 2.5"
    echo "  $0 file.bed stats.csv males.txt hemi.bed --min_absolute 3 --max_absolute 20"
    exit 1
}

# Required arguments
PER_BASE_BED_FILE=""
WHOLE_GENOME_STATS=""
MALES=""
HEMI_REGION_BED=""
MIN_RELATIVE=""
MAX_RELATIVE=""
MIN_ABSOLUTE=""
MAX_ABSOLUTE=""

if [ $# -lt 4 ]; then
    usage
fi

PER_BASE_BED_FILE=$1
WHOLE_GENOME_STATS=$2
MALES=$3
HEMI_REGION_BED=$4

shift 4

# Parsing optional arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --min_relative)
            MIN_RELATIVE="$2"
            shift 2
            ;;
        --max_relative)
            MAX_RELATIVE="$2"
            shift 2
            ;;
        --min_absolute)
            MIN_ABSOLUTE="$2"
            shift 2
            ;;
        --max_absolute)
            MAX_ABSOLUTE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check required arguments
if [[ -z "$PER_BASE_BED_FILE" || -z "$WHOLE_GENOME_STATS" || -z "$MALES" || -z "$HEMI_REGION_BED" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check conflicting arguments
if [[ -n "$MIN_RELATIVE" && -n "$MIN_ABSOLUTE" ]]; then
    echo "ERROR: Cannot use both --min_relative and --min_absolute simultaneously"
    exit 1
fi

if [[ -n "$MAX_RELATIVE" && -n "$MAX_ABSOLUTE" ]]; then
    echo "ERROR: Cannot use both --max_relative and --max_absolute simultaneously"
    exit 1
fi

# Check that at least one threshold is specified
if [[ -z "$MIN_RELATIVE" && -z "$MIN_ABSOLUTE" && -z "$MAX_RELATIVE" && -z "$MAX_ABSOLUTE" ]]; then
    echo "ERROR: No threshold specified. Use --min_relative/--min_absolute and/or --max_relative/--max_absolute"
    exit 1
fi

COVERAGE=$(cat "$WHOLE_GENOME_STATS" | sed -n 2p | awk '{print $2}')
SAMPLE=$(basename "$PER_BASE_BED_FILE" | cut -d. -f1)

source $(conda info --base)/etc/profile.d/conda.sh
conda activate py38

# Function to build command with correct parameters
build_command() {
    local coverage=$1
    local input_file=$2
    local output_file=$3

    local cmd="$TOOLS/MAVR/scripts/alignment/coverage/generate_mask_from_coverage_bed.py"
    cmd+=" -c $input_file"
    cmd+=" -m $coverage"

    # Add threshold parameters
    if [[ -n "$MIN_RELATIVE" ]]; then
        cmd+=" --min_coverage_threshold $MIN_RELATIVE"
    elif [[ -n "$MIN_ABSOLUTE" ]]; then
        cmd+=" --min_coverage_absolute $MIN_ABSOLUTE"
    fi

    if [[ -n "$MAX_RELATIVE" ]]; then
        cmd+=" --max_coverage_threshold $MAX_RELATIVE"
    elif [[ -n "$MAX_ABSOLUTE" ]]; then
        cmd+=" --max_coverage_absolute $MAX_ABSOLUTE"
    fi

    cmd+=" -o $output_file"

    echo "$cmd"
}

# Функция для построения команды для гемизиготных регионов
build_hemi_command() {
    local input_file=$1
    local output_file=$2

    local cmd="$TOOLS/MAVR/scripts/alignment/coverage/generate_mask_from_coverage_bed.py"
    cmd+=" -c $input_file"
    cmd+=" -m $(awk -v cov=$COVERAGE 'BEGIN{print cov/2}')"  # покрытие для гемизиготных

    # Для гемизиготных регионов абсолютные пороги тоже делим на 2
    if [[ -n "$MIN_RELATIVE" ]]; then
        cmd+=" --min_coverage_threshold $MIN_RELATIVE"
    elif [[ -n "$MIN_ABSOLUTE" ]]; then
        hemi_min_absolute=$(awk -v min=$MIN_ABSOLUTE 'BEGIN{print min/2}')
        cmd+=" --min_coverage_absolute $hemi_min_absolute"
    fi

    if [[ -n "$MAX_RELATIVE" ]]; then
        cmd+=" --max_coverage_threshold $MAX_RELATIVE"
    elif [[ -n "$MAX_ABSOLUTE" ]]; then
        hemi_max_absolute=$(awk -v max=$MAX_ABSOLUTE 'BEGIN{print max/2}')
        cmd+=" --max_coverage_absolute $hemi_max_absolute"
    fi

    cmd+=" -o $output_file"

    echo "$cmd"
}

# Function to generate output filename
get_output_filename() {
    local base_input="$1"
    local min_label=""
    local max_label=""

    if [[ -n "$MIN_RELATIVE" ]]; then
        min_label="min${MIN_RELATIVE}"
    elif [[ -n "$MIN_ABSOLUTE" ]]; then
        min_label="minAbs${MIN_ABSOLUTE}"
    fi

    if [[ -n "$MAX_RELATIVE" ]]; then
        max_label="max${MAX_RELATIVE}"
    elif [[ -n "$MAX_ABSOLUTE" ]]; then
        max_label="maxAbs${MAX_ABSOLUTE}"
    fi

    # If both thresholds are specified
    if [[ -n "$min_label" && -n "$max_label" ]]; then
        echo "${base_input%.*.*}.${max_label}.${min_label}.bed"
    elif [[ -n "$min_label" ]]; then
        echo "${base_input%.*.*}.${min_label}.bed"
    elif [[ -n "$max_label" ]]; then
        echo "${base_input%.*.*}.${max_label}.bed"
    fi
}

if grep -qw "$SAMPLE" "$MALES"; then
    echo "[INFO] $SAMPLE == MALE"

    bedtools intersect -a <(zcat "$PER_BASE_BED_FILE") -b "$HEMI_REGION_BED" | gzip > "${SAMPLE}.hemi.tmp.per-base.bed.gz" &
    bedtools intersect -v -a <(zcat "$PER_BASE_BED_FILE") -b "$HEMI_REGION_BED" | gzip > "${SAMPLE}.diploid.tmp.per-base.bed.gz" &

    wait

    # Диплоидные регионы
    diploid_cmd=$(build_command "$COVERAGE" "${SAMPLE}.diploid.tmp.per-base.bed.gz" "${SAMPLE}.diploid.per-base.mask.bed")
    echo "[INFO] Running: $diploid_cmd"
    eval "$diploid_cmd" &

    # Гемизиготные регионы (покрытие / 2)
    hemi_cmd=$(build_hemi_command "${SAMPLE}.hemi.tmp.per-base.bed.gz" "${SAMPLE}.hemi.per-base.mask.bed")
    echo "[INFO] Running: $hemi_cmd"
    eval "$hemi_cmd" &

    wait

    OUTPUT_FILE=$(get_output_filename "$PER_BASE_BED_FILE")
    cat ${SAMPLE}.diploid.per-base.mask.bed ${SAMPLE}.hemi.per-base.mask.bed | sort -S 20G -k1,1 -k2,2n > "$OUTPUT_FILE"
    rm ${SAMPLE}.diploid.tmp.per-base.bed.gz ${SAMPLE}.hemi.tmp.per-base.bed.gz ${SAMPLE}.diploid.per-base.mask.bed ${SAMPLE}.hemi.per-base.mask.bed

else
    echo "[INFO] $SAMPLE == FEMALE"

    OUTPUT_FILE=$(get_output_filename "$PER_BASE_BED_FILE")
    female_cmd=$(build_command "$COVERAGE" "$PER_BASE_BED_FILE" "$OUTPUT_FILE")
    echo "[INFO] Running: $female_cmd"
    eval "$female_cmd"
fi

echo "[INFO] ${SAMPLE}: Done"
echo "[INFO] Output file: $OUTPUT_FILE"
