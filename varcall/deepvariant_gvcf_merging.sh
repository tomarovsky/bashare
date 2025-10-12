#!/bin/bash
# Usage:
# $TOOLS/bashare/deepvariant_gvcf_merging.sh gVCF_GZ_FILES OUTPUT_PREFIX

GVCF_GZ_FILES=$1 # space separated
OUTPUT_PREFIX=$2

source $(conda info --base)/etc/profile.d/conda.sh
conda activate varcall

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 gVCF_GZ_FILES OUTPUT_PREFIX"
    exit 1
fi

glnexus_cli --config DeepVariantWES ${GVCF_GZ_FILES} | bcftools view - -O z -o ${OUTPUT_PREFIX}.vcf.gz
