#!/bin/bash
# Usage: 64 cpu
# $TOOLS/bashare/bcftools_varcall.sh ASSEMBLY "BAM_FILES" OUTPUT_PREFIX [PLOIDY_FILE SAMPLES_FILE]
# $TOOLS/bashare/bcftools_varcall.sh /path/to/genomic.fasta "/path/to/S1.bam /path/to/S2.bam /path/to/S3.bam" prefix

source $(conda info --base)/etc/profile.d/conda.sh
conda activate py38

ASSEMBLY=$1
BAM_FILES=$2 # space separated
OUTPUT_PREFIX=$3
PLOIDY_FILE=$4
SAMPLES_FILE=$5

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 ASSEMBLY BAM_FILES OUTPUT_PREFIX [PLOIDY_FILE SAMPLES_FILE]"
    exit 1
fi

if [[ ! -f "${ASSEMBLY}.fai" ]]; then
    samtools faidx "$ASSEMBLY"
fi

mkdir -p split/ split/mpileup/ split/bcf/;

BCFTOOLS_CALL_OPTS="-m -O u -v -a GQ,GP"
if [[ -n "$PLOIDY_FILE" && -n "$SAMPLES_FILE" ]]; then
    BCFTOOLS_CALL_OPTS="--ploidy-file $PLOIDY_FILE --samples-file $SAMPLES_FILE --group-samples - $BCFTOOLS_CALL_OPTS"
fi

${TOOLS}/MAVR/scripts/sequence/prepare_region_list.py -r ${ASSEMBLY} -s -m 1500000 -n 100 -g samtools -x 1000 2>/dev/null | \
    parallel -j 20 "bcftools mpileup -d 250 -q 30 -Q 30 --adjust-MQ 50 \
    -a AD,INFO/AD,ADF,INFO/ADF,ADR,INFO/ADR,DP,SP,SCR,INFO/SCR -O u -f ${ASSEMBLY} -r {} ${BAM_FILES} | \
    tee split/mpileup/tmp.{#}.mpileup.bcf | \
    bcftools call $BCFTOOLS_CALL_OPTS > split/bcf/tmp.{#}.bcf"


bcftools concat -O u --threads 32 `ls split//bcf//tmp.*.bcf | sort -V` | bcftools view -O z -o ${OUTPUT_PREFIX}.vcf.gz - &
# bcftools concat -O u --threads 32 `ls split//mpileup//tmp.*.mpileup.bcf | sort -V` | bcftools view -O z -o ${OUTPUT_PREFIX}.mpileup.vcf.gz - &

wait

