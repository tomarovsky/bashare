#!/bin/bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <sample.fq.gz> <mask.bed>" >&2
    exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate py38

CONSENSUS_FASTQ="$1"
MASK_BED="$2"
OUTPUT_FASTQ_GZ="${CONSENSUS_FASTQ%.*.*}.masked.fq.gz"

# temporary files
BASE_NAME="${CONSENSUS_FASTQ%.*.*}"
FASTA_TMP="${BASE_NAME}.tmp.fasta.gz"
QUAL_TMP="${BASE_NAME}.tmp.qual.gz"
FASTA_MASKED_TMP="${BASE_NAME}.tmp.masked.fasta.gz"
QUAL_MASKED_TMP="${BASE_NAME}.tmp.masked.qual.gz"

# cleanup() {
#     echo "$(date) | Cleaning up temporary files..." >&2
#     rm -f "$FASTA_TMP" "$QUAL_TMP" "$FASTA_MASKED_TMP" "$QUAL_MASKED_TMP"
# }
# trap cleanup EXIT

echo "$(date) | Split $CONSENSUS_FASTQ into FASTA and QUAL..." >&2
$TOOLS/Biocrutch/scripts/psmc/split_consensus_fq.py "$CONSENSUS_FASTQ" "$FASTA_TMP" "$QUAL_TMP"

echo "$(date) | Mask sequences (replace with 'n') using bedtools..." >&2
conda deactivate && conda activate varcall
bedtools maskfasta -fi "$FASTA_TMP" -bed "$MASK_BED" -fo "$FASTA_MASKED_TMP" -mc 'n'
conda deactivate && conda activate py38

echo "$(date) | Masking quality (to '!') based on masked sequence file..." >&2
$TOOLS/Biocrutch/scripts/psmc/mask_consensus_qual.py "$FASTA_MASKED_TMP" "$QUAL_TMP" "$QUAL_MASKED_TMP"

echo "$(date) | Merge masked FASTA and QUAL files into ${OUTPUT_FASTQ_GZ}..." >&2
$TOOLS/Biocrutch/scripts/psmc/merge_into_consensus_fq.py "$FASTA_MASKED_TMP" "$QUAL_MASKED_TMP" "$OUTPUT_FASTQ_GZ"

echo "$(date) | DONE | ${CONSENSUS_FASTQ} -> ${OUTPUT_FASTQ_GZ}" >&2
