#!/bin/bash
# Usage:
# $TOOLS/bashare/plink.sh VCF OUTPREFIX LD THREADS

VCF=$1
OUTPREFIX=$2
LD=$3
THREADS=$4

source $(conda info --base)/etc/profile.d/conda.sh
conda activate plink1.9

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 VCF OUTPREFIX LD THREADS"
    exit 1
fi

# Step 1: geno and maf
plink --vcf $VCF --out ${OUTPREFIX}.prefiltered \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --geno 0 \
    --maf 0.03 \
    --snps-only \
    --make-bed \
    --threads $THREADS

# Step 2: LD pruning
plink --bfile ${OUTPREFIX}.prefiltered --out ${OUTPREFIX} \
    --indep-pairwise 50 10 $LD \
    --threads $THREADS

# Step 3: PCA
plink --bfile ${OUTPREFIX}.prefiltered --out ${OUTPREFIX} \
    --extract ${OUTPREFIX}.prune.in \
    --pca \
    --make-bed \
    --threads $THREADS

# to fix 'HiC_scaffold_' naming
cat ${OUTPREFIX}.bim | awk -F '_' '{print $3"_"$4"_"$5}' > ${OUTPREFIX}.bim.tmp
mv ${OUTPREFIX}.bim ${OUTPREFIX}.bim.raw
mv ${OUTPREFIX}.bim.tmp ${OUTPREFIX}.bim

