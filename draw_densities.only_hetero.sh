#!/bin/bash

print_usage() {
	echo "Genetic variants vizualization."
	echo "run script in bcftools_filtration folder"
}

while getopts '' flag; do
	case "${flag}" in
		*) print_usage
			exit 1 ;;
	esac
done

# paths
assembly=${1}

# script
for file in *.snp.hetero.vcf.gz; do
	echo ${file};
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o ${file%.*.*}.w1mb.s100kb -w 1000000 -s 100000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --hide_track_label --density_thresholds 0,0.1,0.5,1,1.5,2  --rounded --subplots_adjust_left 0.2 --only_count;
done



