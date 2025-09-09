#!/bin/bash
# Usage:
# $TOOLS/bashare/assembly_stats.sh FASTA NUMBER_OF_CHROMOSOME

FASTA=$1
NUMBER_OF_CHROMOSOME=$2

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 FASTA NUMBER_OF_CHROMOSOME"
    exit 1
fi

if [[ "$FASTA" == *.gz ]]; then
    unpigz -p 8 "$FASTA"
    FASTA="${FASTA%.gz}"
fi

if [[ ! -f "${FASTA}.fai" ]]; then
    samtools faidx "$FASTA"
fi

FASTA_PREFIX=${FASTA%.*}

cat ${FASTA}.fai | awk '{print $1"\t"$2}' | sort -nr -k2 > ${FASTA_PREFIX}.len
cat ${FASTA_PREFIX}.len | awk '{print $1}' | head -n ${NUMBER_OF_CHROMOSOME} > ${FASTA_PREFIX}.whitelist
cat ${FASTA_PREFIX}.whitelist | awk '{print $1"\t"$1}' | head -n ${NUMBER_OF_CHROMOSOME} > ${FASTA_PREFIX}.syn
cat ${FASTA_PREFIX}.syn | awk '{print $2}' | head -n ${NUMBER_OF_CHROMOSOME} > ${FASTA_PREFIX}.renamelist

