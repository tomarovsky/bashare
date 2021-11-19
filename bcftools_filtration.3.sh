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

# paths
export PATH=${PATH}:/mnt/tank/scratch/atomarovsky/tools/bcftools-1.12/bin:/mnt/tank/scratch/atomarovsky/tools/parallel-20210622/bin:/mnt/tank/scratch/atomarovsky/tools/vcftools-0.1.16/bin

# vcf prefix
pref=${vcf%.*.*}

# script
mkdir bcftools_filtration;
cd  bcftools_filtration/;

echo "bcftool filter..."
bcftools filter -S . -o ${pref}.filt.vcf --exclude 'QUAL < 20.0 || (FORMAT/SP > 60.0 | FORMAT/DP < 5.0 | FORMAT/GQ < 20.0)' ../${vcf}

fvcf=${pref}.filt.vcf
pref=${fvcf%.*}

touch tmp.list
echo "bcftools query..."
for SAMPLE in `bcftools query -l ${fvcf}`; do
	echo "Handling ${SAMPLE}..." ; 
	vcf-subset --exclude-ref -c ${SAMPLE} ${fvcf} > ${pref}.${SAMPLE}.vcf ; 
	echo ${pref}.${SAMPLE}.vcf >> tmp.list ;
done

for file in $(cat tmp.list); do
	echo "bcftools filter...";
	bcftools  filter -i  'TYPE="indel"' -O v ${file} > ${file%.*}.indels.vcf && bcftools  filter -i  'TYPE="snp"' -O v ${file} > ${file%.*}.SNPs.vcf;
	echo "get hyterozygous variants..."
	bcftools filter -i 'FMT/GT = "het"' -O v ${file%.*}.indels.vcf > ${file%.*}.indels.hetero.vcf && bcftools filter -i 'FMT/GT = "het"' -O v ${file%.*}.SNPs.vcf > ${file%.*}.SNPs.hetero.vcf
	echo "get homozygous variants..."
	bcftools filter -i 'FMT/GT = "hom"' -O v ${file%.*}.indels.vcf > ${file%.*}.indels.homo.vcf && bcftools filter -i 'FMT/GT = "hom"' -O v ${file%.*}.SNPs.vcf > ${file%.*}.SNPs.homo.vcf
done

rm tmp.list

# Draw densities of heterozygous files
# /home/atomarovsky/bashare/draw_densities_of_hetero.sh
