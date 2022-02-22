#!/bin/bash

print_usage() {
	echo "Creating files for coverage (mosdepth) vizualization." 
	echo ""
}

while getopts '' flag; do
	case "${flag}" in
		*) print_usage
			exit 1 ;;
	esac
done

# script
for file in *_1000000_windows_stats.csv; do 
	echo $file; 
	python3 $TOOLS/MACE/scripts/draw_coverage.py --scaffold_column_name '#scaffold' --window_column_name frame --coverage_column_name median -i ${file} -o ${file%_*_*_*}.1000000.track --subplots_adjust_left 0.35 -l "'Coverage of ${file%_*_*_*}'" -m $(cat *_whole_genome_stats.csv | sed -n 2p | awk '{print $2}') -w 1000000 -n ../../../../assemblies/hic/*.len -a ../../../../assemblies/hic/*.whitelist -z ../../../../assemblies/hic/*.renamelist --scaffold_syn_file ../../../../assemblies/hic/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet;
done

for file in *_100000_windows_stats.csv; do 
	echo $file; 
	python3 $TOOLS/MACE/scripts/draw_coverage.py --scaffold_column_name '#scaffold' --window_column_name frame --coverage_column_name median -i ${file} -o ${file%_*_*_*}.100000.track --subplots_adjust_left 0.35 -l "'Coverage of ${file%_*_*_*}'" -m $(cat *_whole_genome_stats.csv | sed -n 2p | awk '{print $2}') -w 100000 -n ../../../../assemblies/hic/*.len -a ../../../../assemblies/hic/*.whitelist -z ../../../../assemblies/hic/*.renamelist --scaffold_syn_file ../../../../assemblies/hic/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet;
done


