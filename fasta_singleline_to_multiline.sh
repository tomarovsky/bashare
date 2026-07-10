#!/bin/bash
set -euo pipefail

source "$TOOLS/bashare/lib/log_functions.sh"

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 FASTA [LINE_WIDTH]"
    echo "  Convert a single-line FASTA/FASTA.gz into a multi-line (wrapped) FASTA/FASTA.gz."
    echo "  FASTA        input file, plain or gzip-compressed (.gz)"
    echo "  LINE_WIDTH   sequence line width, default 60"
    exit 1
fi

FASTA="$1"
LINE_WIDTH="${2:-60}"

if [[ ! -f "$FASTA" ]]; then
    log_error "Input file not found: $FASTA"
    exit 1
fi

if [[ "$FASTA" == *.gz ]]; then
    STEM="${FASTA%.gz}"                            # genome.fasta.gz -> genome.fasta
    OUTPUT="${STEM%.*}.multiline.${STEM##*.}.gz"   # -> genome.multiline.fasta.gz
    log_info "Wrapping (width=${LINE_WIDTH}): $FASTA -> $OUTPUT"
    seqtk seq -l "$LINE_WIDTH" "$FASTA" | pigz -p 8 -c > "$OUTPUT"
else
    OUTPUT="${FASTA%.*}.multiline.${FASTA##*.}"    # -> genome.multiline.fasta
    log_info "Wrapping (width=${LINE_WIDTH}): $FASTA -> $OUTPUT"
    seqtk seq -l "$LINE_WIDTH" "$FASTA" > "$OUTPUT"
fi

log_info "Done!"
