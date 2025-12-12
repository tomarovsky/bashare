#!/bin/bash
set -euo pipefail

source "$TOOLS/bashare/lib/log_functions.sh"


if [[ $# -ne 2 ]]; then
    echo "Usage: $0 FEATURES_ROH SYN_FILE"
    echo 'Example: parallel -v --progress -j 8 "$TOOLS/bashare/mask/mask_by_roh.sh {} genome.syn" ::: *.w100kb.s10kb.features.roh'
    exit 1
fi

FEATURES_ROH=$1
SYN_FILE=$2

awk -F "\t" 'NR==FNR {map[$2]=$1; next} {if ($1 in map) $1=map[$1]; print}' OFS="\t" ${SYN_FILE} ${FEATURES_ROH} > ${FEATURES_ROH}.mask.bed

log_info "Done!"
