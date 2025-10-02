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

plink --vcf $VCF --out ${OUTPREFIX} \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --indep-pairwise 50 10 $LD \
    --threads $THREADS

plink --vcf $VCF --out ${OUTPREFIX} \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --extract ${OUTPREFIX}.prune.in \
    --geno 0 \
    --maf 0.03 \
    --snps-only \
    --pca \
    --make-bed \
    --threads $THREADS

# to fix 'HiC_scaffold_' naming
cat ${OUTPREFIX}.bim | awk -F '_' '{print $3"_"$4"_"$5}' > ${OUTPREFIX}.bim.tmp
mv ${OUTPREFIX}.bim ${OUTPREFIX}.bim.raw
mv ${OUTPREFIX}.bim.tmp ${OUTPREFIX}.bim
