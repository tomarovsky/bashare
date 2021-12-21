#!/bin/bash

tabfile=''

print_usage() {
	echo "Usage: '-i' your coverage.tab.gz file and use -t 'genomecov' or 'mosdepth'"
}

while getopts 'i:t:' flag; do
	case "${flag}" in
		i) tabfile="${OPTARG}" ;;
		t) tool="${OPTARG}" ;;
		*) print_usage
			exit 1 ;;
	esac
done

# script
prefix=$(echo ${tabfile:10} | sed 's/\..*//')

echo "whole genome stats..."
$TOOLS/Biocrutch/scripts/Coverage/coverage_statistics.py -i ${tabfile} --tool-name ${tool} -g -o ${prefix}
echo "scaffolds stats..."
$TOOLS/Biocrutch/scripts/Coverage/coverage_statistics.py -i ${tabfile} --tool-name ${tool} -s -o ${prefix}
echo "nonoverlapping windows stats..."
$TOOLS/Biocrutch/scripts/Coverage/coverage_statistics.py -i ${tabfile} --tool-name ${tool} -n -f 1000000 -o ${prefix}
$TOOLS/Biocrutch/scripts/Coverage/coverage_statistics.py -i ${tabfile} --tool-name ${tool} -n -f 100000 -o ${prefix}
$TOOLS/Biocrutch/scripts/Coverage/coverage_statistics.py -i ${tabfile} --tool-name ${tool} -n -f 10000 -o ${prefix}
