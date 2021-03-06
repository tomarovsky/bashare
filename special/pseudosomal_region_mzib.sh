#!/bin/bash

file=''
scaf=''

print_usage() {
	echo "Creating files for genome vizualization." 
	echo "'-s' sex scsffold name"
	echo "'-i' file *_10000_windows_stats name"
}

while getopts 'i:s:' flag; do
	case "${flag}" in
		s) scaf="${OPTARG}" ;;
		i) file="${OPTARG}" ;;
		*) print_usage
			exit 1 ;;
	esac
done

# script
mkdir PAR;
cd PAR/;
# mv ../${file} .

cat ../${file} | awk '{ if ($1 == "'${scaf}'") print $0}' > ${file%_*_*_*}_chrscaf.csv

$TOOLS/Biocrutch/scripts/Coverage/pseudoautosomal_region.py -s ${scaf} -f 10000 -r 10 -i ${file%_*_*_*}_chrscaf.csv -m $(cat ../${file%_*_*_*}_whole_genome_stats.csv | sed -n 2p | awk '{print $2}') -o ${file%_*_*_*} | tee ${file%_*_*_*}_pseudo.log

cd ..
