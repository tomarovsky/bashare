#!/bin/bash

fasta=''

print_usage() {
	echo "Creating files for genome vizualization." 
	echo "Usage: '-f' your genome.fasta file"
}

while getopts 'i:' flag; do
	case "${flag}" in
		i) fasta="${OPTARG}" ;;
		*) print_usage
			exit 1 ;;
	esac
done

# script
mkdir bwa_index
cd bwa_index
ln -s ../${fasta} ${fasta}
bwa index ${fasta}
cd ..

