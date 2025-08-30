#!/bin/bash

fasta=''

print_usage() {
	echo "Creating bowtie2 index." 
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
mkdir bwtie2_index
cd bowtie2_index
ln -s ../${fasta} ${fasta}
bowtie2-build ${fasta} ${fasta%.*}
cd ..

