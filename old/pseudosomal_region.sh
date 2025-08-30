#!/bin/bash

WINDOW_SIZE='10000'

print_usage() {
	echo "Script for PAR calculation" 
	echo "'-s' sex scsffold name"
	echo "'-i' mosdepth.per-base.bed.gz file"
	echo "'-w' window size (default 10000)"

}

while getopts 'i:s:w:g:' flag; do
	case "${flag}" in
		s) SCAFFOLD_NAME="${OPTARG}" ;;
		i) MOSDEPTH_BEDGZ="${OPTARG}" ;;
		w) WINDOW_SIZE="${OPTARG}" ;;
		*) print_usage
			exit 1 ;;
	esac
done

# script
PWD=$(pwd)
echo "whole genome stats..."
python3 $TOOLS/Biocrutch/scripts/Coverage/coverage_statistics.py -i ${MOSDEPTH_BEDGZ} --tool-name mosdepth -g -o ${MOSDEPTH_BEDGZ%.*.*.*}

mkdir ${PWD}/${MOSDEPTH_BEDGZ%.*.*.*}_PAR/; 
cd ${PWD}/${MOSDEPTH_BEDGZ%.*.*.*}_PAR/;

zcat ../${MOSDEPTH_BEDGZ} | awk '{ if ($1 == "'${SCAFFOLD_NAME}'") print $0}' | gzip --stdout > ${MOSDEPTH_BEDGZ%.*.*}.chrX.bed.gz

echo "windows stats..."
python3 $TOOLS/Biocrutch/scripts/Coverage/coverage_statistics.py -i ${MOSDEPTH_BEDGZ%.*.*}.chrX.bed.gz --tool-name mosdepth -n -f ${WINDOW_SIZE} -o ${MOSDEPTH_BEDGZ%.*.*}.chrX

sed -i '1,1d' ${MOSDEPTH_BEDGZ%.*.*}.chrX_${WINDOW_SIZE}_windows_stats.csv; # remove header 
python3 $TOOLS/Biocrutch/scripts/PAR/pseudoautosomal_region.py -f ${WINDOW_SIZE} -i ${MOSDEPTH_BEDGZ%.*.*}.chrX_${WINDOW_SIZE}_windows_stats.csv -s ${SCAFFOLD_NAME} -m $(cat ../${MOSDEPTH_BEDGZ%.*.*.*}_whole_genome_stats.csv | sed -n 2p | awk '{print $2}') -o ${MOSDEPTH_BEDGZ%.*.*}.chrX | tee ${MOSDEPTH_BEDGZ%.*.*}.chrX_pseudo.log

cd ../


