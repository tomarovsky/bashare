#!/bin/bash
# Usage: 64 cpu
# $TOOLS/bashare/gatk_varcall.sh ASSEMBLY BAM_FILES OUTPUT_PREFIX
# $TOOLS/bashare/gatk_varcall.sh genomic.fasta S1.bam,S2.bam,S3.bam gulo_gulo.gatk |& tee -a gatk.log

source $(conda info --base)/etc/profile.d/conda.sh
conda activate python3.8

# GATK
export PATH=$(conda info --base)/envs/gatk/bin/:${TOOLS}/gatk-4.6.2.0/:${PATH}

ASSEMBLY=$1
BAM_FILES=$2 # comma separated
OUTPUT_PREFIX=$3

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 ASSEMBLY BAM_FILES OUTPUT_PREFIX"
    exit 1
fi

if [[ ! -f "${ASSEMBLY%.*}.dict" ]]; then
    picard CreateSequenceDictionary -R "$ASSEMBLY"
fi

$TOOLS/MAVR/scripts/snpcall/parallel_vcf_call_gatk4.py \
    -r ${ASSEMBLY} \
    -o $(pwd)/ \
    -b ${BAM_FILES} \
    -p ${OUTPUT_PREFIX} \
    -t 64 \
    -m 250g


