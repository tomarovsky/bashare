#!/bin/bash

fastq_pref=''
m='23'
s='30G'
t='10'

print_usage() {
	echo "Creating files for genome vizualization." 
	echo "Usage: '-f' your fastq file prefix (file name without .final_[12].fastq)"
}

while getopts 'i:m:s:t:' flag; do
	case "${flag}" in
		i) fastq_pref="${OPTARG}" ;;
		m) m="${OPTARG}" ;;
		s) s="${OPTARG}" ;;
		t) t="${OPTARG}" ;;
		*) print_usage
			exit 1 ;;
	esac
done

# script
python3 ~/tools/krater/scripts/jf/draw_kmer_distribution_from_fastq.py -i ${fastq_pref}.final_1.fastq,${fastq_pref}.final_2.fastq -o ${fastq_pref} -m ${m} -s ${s} -t ${t} -b
