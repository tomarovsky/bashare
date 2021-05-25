#!/bin/bash

vcf=''

print_usage() {
	echo "Creating files for genome vizualization." 
	echo "Usage: '-i' your VCF file"
}

while getopts 'i:' flag; do
	case "${flag}" in
		i) vcf="${OPTARG}" ;;
		*) print_usage
			exit 1 ;;
	esac
done

# script
# mkdir same_filtration_but_using_bcftools/
# cd same_filtration_but_using_bcftools/
# vcf prefix
pref=${vcf%.*}

echo "bcftool filter..."
bcftools filter  -o ${pref}.filt.vcf --exclude 'QUAL < 20.0 || (FORMAT/SP > 60.0 | FORMAT/DP < 5.0 | FORMAT/GQ < 20.0)' ../${vcf}

fvcf=${pref}.filt.vcf
pref=${fvcf%.*}

touch tmp.txt
echo "bcftools query..."
for SAMPLE in `bcftools query -l ${fvcf}`; do echo "Handling ${SAMPLE}..."; vcf-subset --exclude-ref -c ${SAMPLE} ${fvcf} > ${pref}.${SAMPLE}.vcf; echo ${pref}.${SAMPLE}.vcf >> tmp.txt ; done

for file in $(cat tmp.txt); do
	echo "vcf-subset...";
	vcf-subset -t indels ${file} > ${file%.*}.indels.vcf && vcf-subset -t SNPs ${file} > ${file%.*}.SNPs.vcf;
	echo "get_hyterozygous_variants..."
	python3 ~/tools/MACE/scripts/get_heterozygous_variants.py -i ${file%.*}.indels.vcf -o ${file%.*}.indels && python3 ~/tools/MACE/scripts/get_heterozygous_variants.py -i ${file%.*}.SNPs.vcf -o ${file%.*}.SNPs;
done

rm tmp.txt

# Draw densities of heterozygous files
# /home/atomarovsky/bashare/draw_densities_of_hetero.sh
