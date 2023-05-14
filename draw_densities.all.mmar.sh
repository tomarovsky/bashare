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
assembly="../../../../../assemblies/hic.purged/"

# script
for file in *.snp.hetero.vcf; do
	echo ${file};
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o ${file%.*}.100kb --subplots_adjust_left 0.35 -l "Heterozygous SNP for ${file%.*.*.*.*.*.*.*} (M. martes)" -w 100000 -s 100000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --hide_track_label;
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o ${file%.*}.1mb --subplots_adjust_left 0.35 -l "Heterozygous SNP for ${file%.*.*.*.*.*.*.*} (M. martes)" -w 1000000 -s 1000000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --hide_track_label;
done


for file in *.indel.hetero.vcf; do
	echo ${file};
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o ${file%.*}.100kb --subplots_adjust_left 0.35 -l "Heterozygous indels for ${file%.*.*.*.*.*.*.*} (M. martes)" -w 100000 -s 100000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --hide_track_label;
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o ${file%.*}.1mb --subplots_adjust_left 0.35 -l "Heterozygous indels for ${file%.*.*.*.*.*.*.*} (M. martes)" -w 1000000 -s 1000000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --hide_track_label;
done

for file in *.snp.homo.vcf; do
	echo ${file};
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o ${file%.*}.100kb --subplots_adjust_left 0.35 -l "Homozygous SNP for ${file%.*.*.*.*.*.*.*} (M. martes)" -w 100000 -s 100000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --hide_track_label;
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o ${file%.*}.1mb --subplots_adjust_left 0.35 -l "Homozygous SNP for ${file%.*.*.*.*.*.*.*} (M. martes)" -w 1000000 -s 1000000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --hide_track_label;
done


for file in *.indel.homo.vcf; do
	echo ${file};
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o ${file%.*}.100kb --subplots_adjust_left 0.35 -l "Homozygous indels for ${file%.*.*.*.*.*.*.*} (M. martes)" -w 100000 -s 100000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --hide_track_label;
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o ${file%.*}.1mb --subplots_adjust_left 0.35 -l "Homozygous indels for ${file%.*.*.*.*.*.*.*} (M. martes)" -w 1000000 -s 1000000  -a ${assembly}/*.whitelist -z ${assembly}/*.orderedlist --scaffold_syn_file ${assembly}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --hide_track_label;
done


