#!/bin/bash
# Usage:
# $TOOLS/bashare/bcftools_filtration.sh VCF MASK ASSEMBLY THREADS

set -ex

export VCF=$1 # vcf.gz
export PREFIX=${VCF%.*.*}
export MASK=$2 # common mask from common_mask.sh
export ASSEMBLY=$3 # absolute path to assembly files
export THREADS=$4

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 VCF MASK ASSEMBLY THREADS"
    exit 1
fi

#---- script ----
mkdir -p bcftools_filtration/ROH/
cd bcftools_filtration/

echo "$(date) | VCF filtration"
# QUAL < 20.0: Quality, Phred-scaled.
# FORMAT/SP > 60.0: Strand Bias P-value.
# FORMAT/DP < 5.0: Depth, Check depth in specific sample less than 5 reads.
# FORMAT/GQ < 20.0: Genotype Quality.
bcftools filter --threads 30 -S . -O z -o $PREFIX.filt.vcf.gz --exclude 'QUAL < 20.0 || (FORMAT/SP > 60.0 | FORMAT/DP < 5.0 | FORMAT/GQ < 20.0)' ../${VCF}

echo "$(date) | VCF masking"
bedtools intersect -header -v -a ${PREFIX}.filt.vcf.gz -b ${MASK} | bgzip -c > ${PREFIX}.filt.mask.vcf.gz

echo "$(date) | Sample separation"
bcftools query -l $PREFIX.filt.mask.vcf.gz | parallel -j 6 '
    SAMPLE={}
    echo "Sample: ${SAMPLE}"
    bcftools view --threads 10 --min-ac 1 --with-header -s ${SAMPLE} -O z -o ${SAMPLE}.$PREFIX.filt.mask.vcf.gz $PREFIX.filt.mask.vcf.gz
    bcftools filter --threads 10 -i "TYPE=\"snp\"" -O z -o ${SAMPLE}.$PREFIX.filt.mask.snp.vcf.gz ${SAMPLE}.$PREFIX.filt.mask.vcf.gz
    bcftools filter --threads 10 -i "FMT/GT = \"het\"" -O z -o ${SAMPLE}.$PREFIX.filt.mask.snp.hetero.vcf.gz ${SAMPLE}.$PREFIX.filt.mask.snp.vcf.gz
    bcftools filter --threads 10 -i "FMT/GT = \"hom\"" -O z -o ${SAMPLE}.$PREFIX.filt.mask.snp.homo.vcf.gz ${SAMPLE}.$PREFIX.filt.mask.snp.vcf.gz
    # bcftools filter --threads 10 -i "TYPE=\"indel\"" -O z -o ${SAMPLE}.$PREFIX.filt.mask.indel.vcf.gz ${SAMPLE}.$PREFIX.filt.mask.vcf.gz
    # bcftools filter --threads 10 -i "FMT/GT = \"het\"" -O z -o ${SAMPLE}.$PREFIX.filt.mask.indel.hetero.vcf.gz ${SAMPLE}.$PREFIX.filt.mask.indel.vcf.gz
    # bcftools filter --threads 10 -i "FMT/GT = \"hom\"" -O z -o ${SAMPLE}.$PREFIX.filt.mask.indel.homo.vcf.gz ${SAMPLE}.$PREFIX.filt.mask.indel.vcf.gz
'
