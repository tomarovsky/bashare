#!/bin/bash

MOSDEPTH_BED_GZ=''
WHOLE_GENOME_STATS=''

print_usage() {
	echo "Creating mask file for Mosdepth coverage."
	echo "'-i' - per-base.bed.gz from Mosdepth"
	echo "'-w' - whole_genome_stats.csv file"
}

while getopts 'i:w:' flag; do
	case "${flag}" in
		i) MOSDEPTH_BED_GZ="${OPTARG}" ;;
		w) WHOLE_GENOME_STATS="${OPTARG}" ;;
		*) print_usage
			exit 1 ;;
	esac
done

# script
python3 $TOOLS/MAVR/scripts/alignment/coverage/generate_mask_from_coverage_bed.py -c ${MOSDEPTH_BED_GZ} -m $(cat ${WHOLE_GENOME_STATS} | sed -n 2p | awk '{print $2}') -x 2.5 -n 0.33 -o ${MOSDEPTH_BED_GZ%.*.*.*}.max250.min33.bed
