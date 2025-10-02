#!/bin/bash
# Usage:
# $TOOLS/bashare/admixture.sh PLINK_BED K_LIST THREADS

PLINK_BED=$1 # absolute path
K_LIST=$2
THREADS=$3

source $(conda info --base)/etc/profile.d/conda.sh
conda activate admixture

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 PLINK_BED K_LIST THREADS"
    exit 1
fi

for SEED in {41..43}; do
    DIR="seed${SEED}"
    mkdir -p $DIR
    cd $DIR

    for K in $K_LIST; do
        admixture --seed=$SEED -j$THREADS --cv $PLINK_BED $K | tee admixture.K_${K}.log
    done

    cd ..
done

