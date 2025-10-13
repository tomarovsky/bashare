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
    --geno 0.05 \
    --maf 0.03 \
    --snps-only \
    --make-bed \
    --threads $THREADS |& tee -a ${OUTPREFIX}.1.log

# to fix 'HiC_scaffold_' naming
cat ${OUTPREFIX}.prefiltered.bim | awk -F '_' '{print $3"_"$4"_"$5}' > ${OUTPREFIX}.prefiltered.bim.tmp
mv ${OUTPREFIX}.prefiltered.bim ${OUTPREFIX}.prefiltered.bim.raw
mv ${OUTPREFIX}.prefiltered.bim.tmp ${OUTPREFIX}.prefiltered.bim

# Step 2: LD pruning
plink --bfile ${OUTPREFIX}.prefiltered --out ${OUTPREFIX} \
    --indep-pairwise 50 10 $LD \
    --threads $THREADS |& tee -a ${OUTPREFIX}.2.log

# Step 3: PCA
plink --bfile ${OUTPREFIX}.prefiltered --out ${OUTPREFIX} \
    --extract ${OUTPREFIX}.prune.in \
    --pca \
    --make-bed \
    --threads $THREADS |& tee -a ${OUTPREFIX}.3.log
