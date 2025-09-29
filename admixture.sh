#!/bin/bash
# Usage:
# $TOOLS/bashare/admixture.sh PLINK_BED THREADS

PLINK_BED=$1 # absolute path
THREADS=$2

source $(conda info --base)/etc/profile.d/conda.sh
conda activate admixture

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 PLINK_BED THREADS"
    exit 1
fi

for SEED in {41..43}; do
    DIR="seed${SEED}"
    mkdir -p $DIR
    cd $DIR

    for K in {1..6}; do
        admixture --seed=$SEED -j$THREADS --cv $PLINK_BED $K | tee admixture.K_${K}.log
    done

    cd ..
done

