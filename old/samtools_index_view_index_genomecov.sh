#!/bin/bash

bam=''

print_usage() {
	echo "Remove low quality alignments and PCR duplicates" 
	echo "Usage: '-i' your aligment.bam file"
}

while getopts 'i:' flag; do
	case "${flag}" in
		i) bam="${OPTARG}" ;;
		*) print_usage
			exit 1 ;;
	esac
done

# script
echo "samtools index $(bam)"
samtools index ${bam}
echo "samtools view..."
samtools view -b -@ 10 -q 20 -F 1024 -o ${bam%.*.*}.bam ${bam}
rm ${bam}
rm ${bam}.bai
echo "samtools index ${bam%.*.*}.bam"
samtools index ${bam%.*.*}.bam 

echo "genomecov for ${bam%.*.*}.bam"
bedtools genomecov -ibam ${bam%.*.*}.bam -d | gzip > genomecov_${bam%.*.*}.tab.gz
