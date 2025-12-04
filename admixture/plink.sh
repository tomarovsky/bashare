#!/bin/bash

set -euo pipefail

VCF=$1
OUTPREFIX=$2
LD=$3
THREADS=$4

GREEN='\033[0;32m'
NC='\033[0m'

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 VCF OUTPREFIX LD THREADS"
    exit 1
fi

# Step 1: geno and maf
echo "${GREEN}[INFO] | $(date) | Step 1: geno and maf${NC}"
plink --vcf $VCF --out ${OUTPREFIX}.prefiltered \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --geno 0.05 \
    --maf 0.03 \
    --snps-only \
    --make-bed \
    --threads $THREADS |& tee -a ${OUTPREFIX}.1.log

# Step 2: Scaffold renaming
echo "${GREEN}[INFO] | $(date) | Step 2: Renaming scaffolds to integers${NC}"

cp ${OUTPREFIX}.prefiltered.bim ${OUTPREFIX}.prefiltered.bim.raw

awk -v synfile="${OUTPREFIX}.prefiltered.bim.syn" '
BEGIN {
    OFS="\t";
    count=0;
}
{
    # $1
    if (!($1 in map)) {
        count++;
        map[$1] = count;    # Запоминаем mapping
        print $1, count > synfile
    }

    $1 = map[$1];
    print $0;
}' ${OUTPREFIX}.prefiltered.bim.raw > ${OUTPREFIX}.prefiltered.bim

# Step 3: LD pruning
echo "${GREEN}[INFO] | $(date) | Step 3: LD pruning${NC}"

plink --bfile ${OUTPREFIX}.prefiltered --out ${OUTPREFIX} \
    --indep-pairwise 50 10 $LD \
    --allow-extra-chr \
    --threads $THREADS |& tee -a ${OUTPREFIX}.2.log

# Step 4: PCA
echo "${GREEN}[INFO] | $(date) | Step 4: PCA${NC}"

plink --bfile ${OUTPREFIX}.prefiltered --out ${OUTPREFIX} \
    --extract ${OUTPREFIX}.prune.in \
    --allow-extra-chr \
    --pca \
    --make-bed \
    --threads $THREADS |& tee -a ${OUTPREFIX}.3.log

echo "${GREEN}[INFO] | $(date) | Done! ${NC}"
