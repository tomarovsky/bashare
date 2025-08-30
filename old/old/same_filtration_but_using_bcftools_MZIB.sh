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
mkdir same_filtration_but_using_bcftools/
cd same_filtration_but_using_bcftools/
# vcf prefix
pref=${vcf%.*.*}

#echo "bcftool filter..."
bcftools filter  -o ${pref}.filt.vcf --exclude 'QUAL < 20.0 || (FORMAT/SP > 60.0 | FORMAT/DP < 5.0 | FORMAT/GQ < 20.0)' ../${vcf}

fvcf=${pref}.filt.vcf
pref=${fvcf%.*}

touch tmp.txt
echo "bcftools query..."
for SAMPLE in `bcftools query -l ${fvcf}`; do echo "Handling ${SAMPLE}..."; vcf-subset --exclude-ref -c ${SAMPLE} ${fvcf} > ${pref}.${SAMPLE}.vcf; echo ${pref}.${SAMPLE}.vcf >> tmp.txt ; done

for file in $(cat tmp.txt); do
	echo "not vcf-subset. snps and indels by bcftools filter";
	bcftools  filter -i  'TYPE="indel"' -O v ${file} > ${file%.*}.indels.vcf && bcftools  filter -i  'TYPE="snp"' -O v ${file} > ${file%.*}.SNPs.vcf;
	echo "get_hyterozygous_variants by bcftools filter"
	bcftools filter -i 'FMT/GT = "het"' -O v ${file%.*}.indels.vcf > ${file%.*}.indels.hetero.vcf && bcftools filter -i 'FMT/GT = "het"' -O v ${file%.*}.SNPs.vcf > ${file%.*}.SNPs.hetero.vcf
done

rm tmp.txt

# Draw densities of heterozygous files
# /home/atomarovsky/bashare/draw_densities_of_hetero.sh
