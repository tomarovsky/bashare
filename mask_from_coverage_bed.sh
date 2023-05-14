#!/bin/bash

print_usage() {
	echo "Generate mask from coverage BED (mosdepth)" 
	echo ""
}

while getopts '' flag; do
	case "${flag}" in
		*) print_usage
			exit 1 ;;
	esac
done

# script
for file in *.hic.purged.mkdup.mapq20.per-base.bed.gz; do 
	echo $file; 
	python3 $TOOLS/MAVR/scripts/alignment/coverage/generate_mask_from_coverage_bed.py -c ${file} -m $(cat *.mapq20_whole_genome_stats.csv | sed -n 2p | awk '{print $2}') --max_coverage_threshold 2.5 --min_coverage_threshold 0.33 -o ${file%_*_*_*}.max250.min33
done

for file in *.hic.purged.mkdup.per-base.bed.gz; do 
	echo $file; 
	python3 $TOOLS/MAVR/scripts/alignment/coverage/generate_mask_from_coverage_bed.py -c ${file} -m $(cat *_whole_genome_stats.csv | sed -n 4p | awk '{print $2}') --max_coverage_threshold 2.5 --min_coverage_threshold 0.33 -o ${file%_*_*_*}.max250.min33
done

