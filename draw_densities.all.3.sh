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
	echo "${file} is being processed ...";
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 3).snps.100kb.hetero --subplots_adjust_left 0.35 -l "Heterozygous SNPs for $(echo $file | cut -d'.' -f 3)" -w 100000 -s 100000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet --hide_track_label;
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 3).snps.1mb.hetero --subplots_adjust_left 0.35 -l "Heterozygous SNPs for $(echo $file | cut -d'.' -f 3)" -w 1000000 -s 1000000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet --hide_track_label;
done


for file in *.indels.hetero.vcf; do
	echo "${file} is being processed ...";
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 3).indels.100kb.hetero --subplots_adjust_left 0.35 -l "Heterozygous indels for $(echo $file | cut -d'.' -f 3)" -w 100000 -s 100000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet --hide_track_label;
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 3).indels.1mb.hetero --subplots_adjust_left 0.35 -l "Heterozygous indels for $(echo $file | cut -d'.' -f 3)" -w 1000000 -s 1000000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet --hide_track_label;
done

for file in *.SNPs.homo.vcf; do
	echo "${file} is being processed ...";
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 3).snps.100kb.homo --subplots_adjust_left 0.35 -l "Homozygous SNPs for $(echo $file | cut -d'.' -f 3)" -w 100000 -s 100000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet --hide_track_label;
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 3).snps.1mb.homo --subplots_adjust_left 0.35 -l "Homozygous SNPs for $(echo $file | cut -d'.' -f 3)" -w 1000000 -s 1000000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet --hide_track_label;
done


for file in *.indels.homo.vcf; do
	echo "${file} is being processed ...";
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 3).indels.100kb.homo --subplots_adjust_left 0.35 -l "Homozygous indels for $(echo $file | cut -d'.' -f 3)" -w 100000 -s 100000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet --hide_track_label;
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 3).indels.1mb.homo --subplots_adjust_left 0.35 -l "Homozygous indels for $(echo $file | cut -d'.' -f 3)" -w 1000000 -s 1000000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet --hide_track_label;
done


