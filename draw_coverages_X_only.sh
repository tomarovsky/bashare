#!/bin/bash

fasta=''

print_usage() {
	echo "Creating files for genome vizualization." 
	echo "Usage: '-i' your file"
}

while getopts 'i:' flag; do
	case "${flag}" in
		f) fasta="${OPTARG}" ;;
		*) print_usage
			exit 1 ;;
	esac
done

# script
