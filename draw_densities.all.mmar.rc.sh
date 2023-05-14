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
assembly="/mnt/tank/scratch/skliver/common/mustelidae/martes_martes/genome/assemblies/hic.purged/"

# script
for file in *.snp.hetero.vcf.gz; do
	echo ${file};
	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o ${file%.*}.w100kb.s10kb -l "HeteroSNPs densities for ${file%.*.*.*.*.*.*.*.*.*} (Sex: M, Reference: M. martes)" -w 100000 -s 10000 -a ${assembly}/*.revcomp.whitelist -z ${assembly}/*.revcomp.orderedlist --scaffold_syn_file ${assembly}/*.revcomp.syn --syn_file_key_column 0 --syn_file_value_column 1 --hide_track_label --density_thresholds 0,0.1,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6 --rounded --subplots_adjust_left 0.2 --only_count;
done

#for file in *.snp.homo.vcf.gz; do
#	echo ${file};
#	python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o ${file%.*}.w1mb.s100kb -l "HomoSNPs densities for ${file%.*.*.*.*.*.*.*.*.*} (Sex: M, Reference: M. martes)" -w 1000000 -s 100000  -a ${assembly}/*.revcomp.whitelist -z ${assembly}/*.revcomp.orderedlist --scaffold_syn_file ${assembly}/*.revcomp.syn --syn_file_key_column 0 --syn_file_value_column 1 --hide_track_label --density_thresholds 0,0.1,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6 --rounded --subplots_adjust_left 0.2;
#done


