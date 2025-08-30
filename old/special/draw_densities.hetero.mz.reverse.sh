#!/bin/bash

print_usage() {
	echo "Creating files for genome vizualization." 
	echo "run script in bcftools_filtration folder"
}

while getopts '' flag; do
	case "${flag}" in
		*) print_usage
			exit 1 ;;
	esac
done

# paths
assembly="../../../../../assemblies/hic/"

# script
for file in *.SNPs.hetero.vcf; do
	SAMPLE=$(echo $file | cut -d'.' -f 3);
	echo "${file} is being processed ...";
	python3 $TOOLS/forks/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $SAMPLE.snps.100kb.hetero --subplots_adjust_left 0.35 -l "Heterozygous SNPs for $SAMPLE (M.zibellina)" -w 100000 -s 100000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --scaffold_reverse_list ${assembly}/*.reverselist --colormap jet --hide_track_label --density_thresholds 0,0.1,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7;
	python3 $TOOLS/forks/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $SAMPLE.snps.1mb.hetero --subplots_adjust_left 0.35 -l "Heterozygous SNPs for $SAMPLE (M.zibellina)" -w 1000000 -s 1000000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --scaffold_reverse_list ${assembly}/*.reverselist --colormap jet --hide_track_label --density_thresholds 0,0.1,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7;
done


for file in *.indels.hetero.vcf; do
	SAMPLE=$(echo $file | cut -d'.' -f 3);
	echo "${file} is being processed ...";
	python3 $TOOLS/forks/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $SAMPLE.indels.100kb.hetero --subplots_adjust_left 0.35 -l "Heterozygous indels for $SAMPLE (M.zibellina)" -w 100000 -s 100000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --scaffold_reverse_list ${assembly}/*.reverselist --colormap jet --hide_track_label --density_thresholds 0,0.1,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7;
	python3 $TOOLS/forks/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $SAMPLE.indels.1mb.hetero --subplots_adjust_left 0.35 -l "Heterozygous indels for $SAMPLE (M.zibellina)" -w 1000000 -s 1000000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --scaffold_reverse_list ${assembly}/*.reverselist --colormap jet --hide_track_label --density_thresholds 0,0.1,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7;
done


