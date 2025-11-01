#!/bin/bash
# Usage:
# $TOOLS/bashare/wtr.sh FASTA

source $(conda info --base)/etc/profile.d/conda.sh
SIF="$TOOLS/dfam-tetools-latest.sif"

for var in SINGULARITYENV_TMPDIR TMPDIR TMP TEMP; do
    export $var=/mnt/tank/scratch/atomarovsky/tmp/
done

# script
FASTA=$1

if [[ "$FASTA" == *.gz ]]; then
    unpigz -p 8 "$FASTA"
    FASTA="${FASTA%.gz}"
fi

FASTA_PREFIX=$(basename $FASTA .fasta)
DIR_PATH=$(dirname $FASTA)
OUTPUT_DIR_PATH=${DIR_PATH}/${FASTA_PREFIX}

mkdir -p ${OUTPUT_DIR_PATH}/{WindowMasker-1.0.0,TRF-4.09.1,RepeatMasker-4.1.6}
ln -s $(realpath $FASTA) ${OUTPUT_DIR_PATH}/

FASTA=${OUTPUT_DIR_PATH}/$(basename $FASTA)


# WindowMasker
$TOOLS/windowmasker-1.0.0/windowmasker \
    -in $FASTA \
    -infmt fasta \
    -mk_counts \
    -parse_seqids \
    -out ${OUTPUT_DIR_PATH}/WindowMasker-1.0.0/${FASTA_PREFIX}.counts \
    >> ${OUTPUT_DIR_PATH}/WindowMasker-1.0.0/WindowMasker-1.0.0.log 2>&1

$TOOLS/windowmasker-1.0.0/windowmasker \
    -in $FASTA \
    -infmt fasta \
    -ustat ${OUTPUT_DIR_PATH}/WindowMasker-1.0.0/${PREFIX}.counts \
    -outfmt interval \
    -parse_seqids \
    -out ${OUTPUT_DIR_PATH}/WindowMasker-1.0.0/${FASTA_PREFIX}.interval \
    >> ${OUTPUT_DIR_PATH}/WindowMasker-1.0.0/WindowMasker-1.0.0.log 2>&1


# TRF
conda activate sing
singularity run --bind $(pwd) --pwd $(pwd) ${SIF} trf ${FASTA} 2 7 7 80 10 50 500 -l 10 -d -h >> ${OUTPUT_DIR_PATH}/TRF-4.09.1/TRF-4.09.1.log 2>&1 &
mv ${FASTA}.2.7.7.80.10.50.500.dat ${OUTPUT_DIR_PATH}/TRF-4.09.1/${FASTA_PREFIX}.2.7.7.80.10.50.500.l10.dat


# RepeatMasker
singularity run --bind $(pwd) --pwd $(pwd) ${SIF} RepeatMasker -pa 18 -a -species carnivora $FASTA > ${OUTPUT_DIR_PATH}/RepeatMasker-4.1.6/RepeatMasker-4.1.6.log 2>&1
mv ${FASTA}.out ${OUTPUT_DIR_PATH}/RepeatMasker-4.1.6/${FASTA_PREFIX}.out
mv ${FASTA}.tbl ${OUTPUT_DIR_PATH}/RepeatMasker-4.1.6/${FASTA_PREFIX}.tbl
mv ${FASTA}.align ${OUTPUT_DIR_PATH}/RepeatMasker-4.1.6/${FASTA_PREFIX}.align
mv ${FASTA}.masked ${OUTPUT_DIR_PATH}/RepeatMasker-4.1.6/${FASTA_PREFIX}.masked.fasta

    singularity run --bind $(pwd) --pwd $(pwd) /mnt/tank/scratch/atomarovsky/three_martens/WTR/dfam-tetools-latest.sif calcDivergenceFromAlign.pl -s ${sample_prefix}.divsum -a ${sample_prefix}.divsum.align ${sample_prefix}.align
    singularity run --bind $(pwd) --pwd $(pwd) /mnt/tank/scratch/atomarovsky/three_martens/WTR/dfam-tetools-latest.sif createRepeatLandscape.pl -div ${sample_prefix}.divsum -g $(sed -n '2p' ${sample_prefix}.quast | cut -d ',' -f 3) > ${sample_prefix}.html


conda deactivate && conda activate py38
$TOOLS/Biocrutch/scripts/RepeatMasking/WindowMasker.py \
    -i ${OUTPUT_DIR_PATH}/WindowMasker-1.0.0/${PREFIX}.interval \
    -o ${OUTPUT_DIR_PATH}/WindowMasker-1.0.0/${PREFIX}.interval
$TOOLS/Biocrutch/scripts/RepeatMasking/TRF.py \
    -i ${OUTPUT_DIR_PATH}/TRF-4.09.1/${FASTA}.2.7.7.80.10.50.500.l10.dat \
    -o ${OUTPUT_DIR_PATH}/TRF-4.09.1/${FASTA}.2.7.7.80.10.50.500.l10.dat
$TOOLS/Biocrutch/scripts/RepeatMasking/RepeatMasker.py \
    -i ${OUTPUT_DIR_PATH}/RepeatMasker-4.1.6/${FASTA}.out \
    -o ${OUTPUT_DIR_PATH}/RepeatMasker-4.1.6/${FASTA}.out
conda deactivate

    # gzip ${sample_prefix}.counts
    # gzip ${sample_prefix}.interval

    # gzip ${sample_prefix}.2.7.7.80.10.50.500.l10.dat
    # $TOOLS/BerryTart/python_scr/GFF_to_TAB.py ${sample_prefix}.2.7.7.80.10.50.500.l10.dat.gff.gz ${sample_prefix}.2.7.7.80.10.50.500.l10.dat.gff.tab

    # gzip ${sample_prefix}.align
    # gzip ${sample_prefix}.divsum.align
    # gzip ${sample_prefix}.masked.fasta
    # gzip ${sample_prefix}.out

# Concat GFF
zcat ${OUTPUT_DIR_PATH}/TRF-4.09.1/${FASTA}.2.7.7.80.10.50.500.l10.dat.gff.gz \
    ${OUTPUT_DIR_PATH}/RepeatMasker-4.1.6/${FASTA}.out.gff.gz \
    ${OUTPUT_DIR_PATH}/WindowMasker-1.0.0/${PREFIX}.interval.gff.gz |\
    sort -k1,1 -k4,4n -k5,5n |\
    gzip > ${OUTPUT_DIR_PATH}.trf.repeatmasker.windowmasker.gff.gz

# FASTA masking
conda activate varcall
bedtools maskfasta -soft -bed ${PREFIX}.trf.repeatmasker.windowmasker.gff.gz -fi ${FASTA} -fo ${PREFIX}.masked.fasta && pigz -p 8 ${PREFIX}.masked.fasta
conda deactivate
