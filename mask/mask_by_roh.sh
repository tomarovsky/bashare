#!/bin/bash
# Usage:
# $TOOLS/bashare/mask_by_roh.sh FEATURES_ROH SYN_FILE
# parallel -v --progress -j 8 "$TOOLS/bashare/mask_by_roh.sh {} genome.syn" ::: *.w100kb.s10kb.features.roh

FEATURES_ROH=$1
SYN_FILE=$2

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 FEATURES_ROH SYN_FILE"
    exit 1
fi

awk -F "\t" 'NR==FNR {map[$2]=$1; next} {if ($1 in map) $1=map[$1]; print}' OFS="\t" ${SYN_FILE} ${FEATURES_ROH} > ${FEATURES_ROH}.mask.bed
