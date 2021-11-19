#!/bin/bash

print_usage() {
	echo "Creating files for genome vizualization." 
	echo "run script in 'same_filtration_but_using_bcftools' folder"
}

while getopts '' flag; do
	case "${flag}" in
		*) print_usage
			exit 1 ;;
	esac
done


# script
#for file in *.SNPs.hetero.vcf; do
#	echo "${file} is being processed ...";
#	python3 ~/tools/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 1).snps.100kb.hetero --subplots_adjust_left 0.35 -l "'SNPs for ${file%.*.*.*}'" -w 100000 -s 100000  -a ../../*X.whitelist -z ../../*X.orderlist --scaffold_syn_file ../../*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet;
#	python3 ~/tools/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 1).snps.1mb.hetero --subplots_adjust_left 0.35 -l "'SNPs for ${file%.*.*.*}'" -w 1000000 -s 1000000  -a ../../*X.whitelist -z ../../*X.orderlist --scaffold_syn_file ../../*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet;
#done


#for file in *.indels.hetero.vcf; do
#	echo "${file} is being processed ...";
#	python3 ~/tools/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 1).indels.100kb.hetero --subplots_adjust_left 0.35 -l "'indels for ${file%.*.*.*}'" -w 100000 -s 100000  -a ../../*X.whitelist -z ../../*X.orderlist --scaffold_syn_file ../../*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet;
#	python3 ~/tools/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 1).indels.1mb.hetero --subplots_adjust_left 0.35 -l "'indels for ${file%.*.*.*}'" -w 1000000 -s 1000000  -a ../../*X.whitelist -z ../../*X.orderlist --scaffold_syn_file ../../*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet;
#done

# script
for file in *.SNPs.hetero.vcf; do
	echo "${file} is being processed ...";
	python3 ~/tools/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 1).snps.100kb.hetero --subplots_adjust_left 0.35 -l "'SNPs for ${file%.*.*.*}'" -w 100000 -s 100000  -a ../../*r.whitelist -z ../../*r.orderlist --scaffold_syn_file ../../*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet;
	python3 ~/tools/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 1).snps.1mb.hetero --subplots_adjust_left 0.35 -l "'SNPs for ${file%.*.*.*}'" -w 1000000 -s 1000000  -a ../../*r.whitelist -z ../../*r.orderlist --scaffold_syn_file ../../*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet;
done


for file in *.indels.hetero.vcf; do
	echo "${file} is being processed ...";
	python3 ~/tools/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 1).indels.100kb.hetero --subplots_adjust_left 0.35 -l "'indels for ${file%.*.*.*}'" -w 100000 -s 100000  -a ../../*r.whitelist -z ../../*r.orderlist --scaffold_syn_file ../../*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet;
	python3 ~/tools/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 1).indels.1mb.hetero --subplots_adjust_left 0.35 -l "'indels for ${file%.*.*.*}'" -w 1000000 -s 1000000  -a ../../*r.whitelist -z ../../*r.orderlist --scaffold_syn_file ../../*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet;
done


