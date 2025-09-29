#!/bin/bash
# Usage:
# $TOOLS/bashare/plink_to_vcf.sh PLINK_PREFIX OUTPREFIX

PLINK_PREFIX=$1
OUTPREFIX=$2

source $(conda info --base)/etc/profile.d/conda.sh
conda activate plink1.9

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 BFILE OUTPREFIX"
    exit 1
fi

# plink -> vcf.gz
plink --bfile $PLINK_PREFIX --recode vcf bgz --out $OUTPREFIX

# index
conda activate varcall
bcftools index ${OUTPREFIX}.vcf.gz

