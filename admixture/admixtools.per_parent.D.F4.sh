#!/bin/bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
    echo "Usage: $0 <stat:D|F4> <plink_prefix> <group1:name:samples> <group2:name:samples> <outgroup:name:samples> output_dir"
    echo
    echo "Example:"
    echo "  $0 D mzib.mfoi.allsamples.filt.mask.auto.snp.plink \\"
    echo "     sables:10xmzib,CHN,S26,T104,T118,T148,T150,T194,T26,T50,T72,T8,T90 \\"
    echo "     pines:10xmmar,S44,S46,S49,T149,T24,T76,T77,T82 \\"
    echo "     outgroup:10xmfoi"
    echo "     D/"
    exit 1
fi

script_path="$TOOLS/bashare/admixture/admixtools.D.F4.sh"

STAT_TYPE=$1
PREFIX=$2
GROUP1=$3
GROUP2=$4
OUTGROUP=$5
OUTPUT_DIR=$6

# Check stat type
if [[ "$STAT_TYPE" != "D" && "$STAT_TYPE" != "F4" ]]; then
    echo "Error: <stat> must be either 'D' or 'F4'"
    exit 1
fi

# Parse input groups
parse_group() {
    local arg="$1"
    local name=$(echo "$arg" | cut -d':' -f1)
    local samples=$(echo "$arg" | cut -d':' -f2)
    echo "$name" "$samples"
}

read GROUP1_NAME GROUP1_SAMPLES <<< "$(parse_group "$GROUP1")"
read GROUP2_NAME GROUP2_SAMPLES <<< "$(parse_group "$GROUP2")"
read OUTGROUP_NAME OUTGROUP_SAMPLES <<< "$(parse_group "$OUTGROUP")"

IFS=',' read -ra GROUP1_ARR <<< "$GROUP1_SAMPLES"
IFS=',' read -ra GROUP2_ARR <<< "$GROUP2_SAMPLES"

# per sample in group 1
for sample in "${GROUP1_ARR[@]}"; do
    reduced_group1=()
    for s in "${GROUP1_ARR[@]}"; do
        [[ "$s" != "$sample" ]] && reduced_group1+=("$s")
    done

    group1_str="${GROUP1_NAME}:$(IFS=,; echo "${reduced_group1[*]}")"
    group2_str="${GROUP2_NAME}:$(IFS=,; echo "${GROUP2_ARR[*]}")"

    "$script_path" \
        "$STAT_TYPE" \
        "$PREFIX" \
        "${sample}:${sample}" \
        "$group1_str" \
        "$group2_str" \
        "${OUTGROUP_NAME}:${OUTGROUP_SAMPLES}" \
        "${outdir}/${STAT_TYPE}stat.${sample}_${GROUP1_NAME}_${GROUP2_NAME}.txt"
done

# per sample in group 2
for sample in "${GROUP2_ARR[@]}"; do
    reduced_group2=()
    for s in "${GROUP2_ARR[@]}"; do
        [[ "$s" != "$sample" ]] && reduced_group2+=("$s")
    done

    group1_str="${GROUP1_NAME}:$(IFS=,; echo "${GROUP1_ARR[*]}")"
    group2_str="${GROUP2_NAME}:$(IFS=,; echo "${reduced_group2[*]}")"

    "$script_path" \
        "$STAT_TYPE" \
        "$PREFIX" \
        "${sample}:${sample}" \
        "$group1_str" \
        "$group2_str" \
        "${OUTGROUP_NAME}:${OUTGROUP_SAMPLES}" \
        "${outdir}/${STAT_TYPE}stat.${sample}_${GROUP1_NAME}_${GROUP2_NAME}.txt"
done
