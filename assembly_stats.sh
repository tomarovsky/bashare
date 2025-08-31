#!/bin/bash
# Usage:
# $TOOLS/bashare/assembly_stats.sh FASTA NUMBER_OF_CHROMOSOME

FASTA=$1
NUMBER_OF_CHROMOSOME=$2

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 FASTA NUMBER_OF_CHROMOSOME"
    exit 1
fi

if [[ "$FASTA" == *.gz ]]; then
    gunzip "$FASTA"
    FASTA="${FASTA%.gz}"
fi

if [[ ! -f "${FASTA}.fai" ]]; then
    samtools index "$FASTA"
fi

FASTA_PREFIX=${FASTA%.*}

cat ${FASTA}.fai | awk '{print $1"\t"$2}' | sort -nr -k2 > ${FASTA_PREFIX}.len
cat ${FASTA_PREFIX}.lengths | awk '{print $1}' > ${FASTA_PREFIX}.whitelist
cat ${FASTA_PREFIX}.whitelist | awk '{print $1"\t"$1}' > ${FASTA_PREFIX}.syn
cat ${FASTA_PREFIX}.syn | awk '{print $2}' > ${FASTA_PREFIX}.renamelist

