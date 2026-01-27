#!/bin/bash
set -euo pipefail

source "$TOOLS/bashare/lib/log_functions.sh"

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 ASSEMBLY READ1 READ2 SAMPLE"
    echo "$TOOLS/bashare/alignment/bwa.sh ref.fa Sample_1.fq.gz Sample_2.fq.gz Sample"
    echo "32 CPU"
    exit 1
fi

ASSEMBLY=$1
READ1=$2
READ2=$3
SAMPLE=$4

BAM="${SAMPLE}/${SAMPLE}.bam"

mkdir -p "$SAMPLE"

# Check input files
for f in "$ASSEMBLY" "$READ1" "$READ2"; do
    if [[ ! -f "$f" ]]; then
        log_error "File not found: $f"
        exit 1
    fi
done

# Check BWA index
if [[ ! -f "${ASSEMBLY}.bwt" ]]; then
    log_warning "BWA index not found, building index..."
    bwa index "$ASSEMBLY"
fi

log_info "${SAMPLE} | Running bwa mem + samtools pipeline"

bwa mem \
    -R "@RG\tID:${SAMPLE}\tPU:x\tSM:${SAMPLE}\tPL:Illumina\tLB:x" \
    -t 20 \
    "$ASSEMBLY" \
    <(zcat "$READ1") \
    <(zcat "$READ2") \
| samtools fixmate -@ 2 -m - - \
| samtools sort -@ 8 -m 12G \
| samtools markdup -@ 2 - "$BAM"

log_info "Indexing BAM"
samtools index "$BAM"

log_info "${SAMPLE} | Done! BAM file: ${BAM}"
