#!/bin/bash
# Usage: 64 cpu
# $TOOLS/bashare/deepvariant_varcall.sh ASSEMBLY "BAM_FILES"
# $TOOLS/bashare/deepvariant_varcall.sh genomic.fasta "S1.bam S2.bam S3.bam" |& tee -a deepvariant.log

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 ASSEMBLY BAM_FILES"
    exit 1
fi

ASSEMBLY=$1
BAM_FILES=$2 # space separated
SIF="$TOOLS/deepvariant-1.9.0.sif"

source $(conda info --base)/etc/profile.d/conda.sh
conda activate sing

mkdir -p "$PWD/tmp"
for var in SINGULARITYENV_TMPDIR TMPDIR TMP TEMP; do
    export $var="$PWD/tmp"
done

if [[ "$ASSEMBLY" == *.gz ]]; then
    unpigz -p 8 "$ASSEMBLY"
    ASSEMBLY="${ASSEMBLY%.gz}"
fi

if [[ ! -f "${ASSEMBLY}.fai" ]]; then
    conda run -n varcall samtools faidx "$ASSEMBLY"
fi

# get bing directories
BIND_DIRS=()
BIND_DIRS+=($PWD)
BIND_DIRS+=("$(dirname "$(readlink -f ${ASSEMBLY})")")
for BAM in $BAM_FILES; do
    BIND_DIRS+=("$(dirname "$(readlink -f "$BAM")")")
done
# uniq only
BIND_DIRS=$(printf "%s\n" "${BIND_DIRS[@]}" | sort -u | paste -sd, -)


for BAM in $BAM_FILES; do 

    SAMPLE=$(conda run -n varcall samtools view -H $BAM | grep '@RG' | sed 's/.*SM:\([^ \t]*\).*/\1/' | uniq)

    if [[ ! -f "${BAM}.bai" ]]; then
        conda run -n varcall samtools index "$BAM"
    fi

    singularity run --bind "$BIND_DIRS" --pwd "$(pwd)" $SIF \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=${ASSEMBLY} \
    --reads=${BAM} \
    --output_vcf=${SAMPLE}.deepvariant.vcf.gz \
    --output_gvcf=${SAMPLE}.deepvariant.g.vcf.gz \
    --vcf_stats_report=true \
    --num_shards=64
    # --regions "chr20:10,000,000-10,010,000" \
    # --intermediate_results_dir "${OUTPUT_DIR}/intermediate_results_dir" \
done

