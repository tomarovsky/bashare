#!/bin/bash
# Usage:
# $TOOLS/bashare/plink_to_vcf.sh PLINK_PREFIX

PLINK_PREFIX=$1

source $(conda info --base)/etc/profile.d/conda.sh
conda activate plink1.9

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 PLINK_PREFIX"
    exit 1
fi

# plink -> vcf.gz
plink --bfile $PLINK_PREFIX --recode vcf bgz --out $PLINK_PREFIX

# index
conda activate varcall
bcftools index $PLINK_PREFIX.vcf.gz

