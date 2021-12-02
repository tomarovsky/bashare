#!/bin/bash

genomecov=''

print_usage() {
	echo "Usage: '-i' your genomecov.tab.gz file"
}

while getopts 'i:' flag; do
	case "${flag}" in
		i) genomecov="${OPTARG}" ;;
		*) print_usage
			exit 1 ;;
	esac
done

# script
prefix=$(echo ${genomecov:10} | sed 's/\..*//')

echo "whole genome stats..."
$TOOLS/Biocrutch/scripts/genomecov/coverage_statistics.py -i ${genomecov} -g -o ${prefix}
echo "scaffolds stats..."
$TOOLS/Biocrutch/scripts/genomecov/coverage_statistics.py -i ${genomecov} -s -o ${prefix}
echo "nonoverlapping windows stats..."
$TOOLS/Biocrutch/scripts/genomecov/coverage_statistics.py -i ${genomecov} -n -f 1000000 -o ${prefix}
$TOOLS/Biocrutch/scripts/genomecov/coverage_statistics.py -i ${genomecov} -n -f 100000 -o ${prefix}
$TOOLS/Biocrutch/scripts/genomecov/coverage_statistics.py -i ${genomecov} -n -f 10000 -o ${prefix}  
