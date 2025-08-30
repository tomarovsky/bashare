#!/bin/bash

perbase=''

print_usage() {
	echo "Creating mask file for mosdepth coverage."
}

while getopts 'i:' flag; do
	case "${flag}" in
		i) perbase="${OPTARG}" ;;
		*) print_usage
			exit 1 ;;
	esac
done

# script
$TOOLS/Biocrutch/scripts/Coverage/coverage_masking.py -i ${perbase} -w $(cat *_whole_genome_stats.csv | sed -n 2p | awk '{print $2}') -o ${perbase%.*.*.*}
