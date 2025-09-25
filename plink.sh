#!/bin/bash
# Usage:
# $TOOLS/bashare/plink.sh VCF OUTPREFIX THREADS

VCF=$1
OUTPREFIX=$2
THREADS=$3

source $(conda info --base)/etc/profile.d/conda.sh
conda activate plink1.9

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 VCF OUTPREFIX THREADS"
    exit 1
fi

plink --vcf $VCF --out ${OUTPREFIX} \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --indep-pairwise 50 10 0.7 \
    --threads $THREADS

plink --vcf $VCF --out ${OUTPREFIX} \
    --allow-extra-chr
    --set-missing-var-ids @:# \
    --extract ${VCF%.*.*}.plink.prune.in \
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
