#!/bin/bash

file=''
scafname=''
window=''
whole_genome_stats=''

print_usage() {
	echo "Creating files for genome vizualization." 
	echo "'-s' sex scsffold name"
	echo "'-i' file *_windows_stats name"
	echo "'-w' window size in *_windows_stats file"
	echo "'-g' *_whole_genome_stats file"

}

while getopts 'i:s:w:' flag; do
	case "${flag}" in
		s) scafname="${OPTARG}" ;;
		i) file="${OPTARG}" ;;
		w) window="${OPTARG}" ;;
		g) whole_genome_stats="${OPTARG}" ;;
		*) print_usage
			exit 1 ;;
	esac
done

# script
mkdir PAR/; cd PAR/;

cat ../${file} | awk '{ if ($1 == "'${scafname}'") print $0}' > ${file%_*_*_*}_chrscaf.csv

python3 $TOOLS/Biocrutch/scripts/PAR/pseudoautosomal_region.py -f ${window} -i *_chrscaf.csv -s ${scafname} -m $(cat ${whole_genome_stats} | sed -n 2p | awk '{print $2}') -o ${file%_*_*_*} | tee ${file%_*_*_*}_pseudo.log

cd -
