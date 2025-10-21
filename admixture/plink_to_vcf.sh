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
plink --bfile $PLINK_PREFIX --recode vcf bgz --out $PLINK_PREFIX |& tee -a ${PLINK_PREFIX}.to_vcf.log

# fix VCF header (sample_sample -> sample)
mv $PLINK_PREFIX.vcf.gz $PLINK_PREFIX.vcf.gz.tmp
zcat $PLINK_PREFIX.vcf.gz.tmp | awk 'BEGIN{OFS="\t"} /^#CHROM/ {for(i=1;i<=NF;i++){if(i>9){split($i,a,"_"); $i=a[1]}}}1' > $PLINK_PREFIX.vcf.gz
rm $PLINK_PREFIX.vcf.gz.tmp

# index
conda activate varcall
bcftools index $PLINK_PREFIX.vcf.gz
