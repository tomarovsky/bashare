#!/bin/bash

fasta=''
chr_number=''

print_usage() {
	echo "Creating files for genome vizualization." 
	echo "Usage: '-i' your genome.fasta file"
	echo "'-n' the number of chromosomes"
}

while getopts 'i:n:' flag; do
	case "${flag}" in
		i) fasta="${OPTARG}" ;;
		n) chr_number="${OPTARG}" ;;
		*) print_usage
			exit 1 ;;
	esac
done

PREFIX=${fasta%.*}

# script
# samtools index
samtools faidx $fasta

# file.len
cat ${fasta}.fai | awk '{print $1"\t"$2}' | sort -nr -k2 > ${PREFIX}.len
echo "len file is done"
# file.lengths
cat ${fasta}.fai | awk '{print $1"\t"$2}' | sort -nr -k2 | head -n ${chr_number} > ${PREFIX}.lengths
echo "lengths file is done"
# file.whitelist 
cat ${PREFIX}.lengths | awk '{print $1}' > ${PREFIX}.whitelist
echo "whitelist file is done"
# file.orderedlist 
cat ${PREFIX}.lengths | awk '{print $1}' > ${PREFIX}.orderedlist
echo "orderedlist file is done"
# file.syn 
cat ${PREFIX}.whitelist | awk '{print $1"\t"$1}' > ${PREFIX}.syn
echo "syn file is done"
# file.renamelist 
cat ${PREFIX}.syn | awk '{print $2}' > ${PREFIX}.renamelist
echo "rename file is done"

# bwa index
# /home/atomarovsky/bashare/bwa_index.sh -f ${fasta}
