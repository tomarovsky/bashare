#!/bin/bash

print_usage() {
	echo "bcftools_filtration with mask specification." 
	echo "'-i' VCF file name"
	echo "'-m' mask file in BED format"
}

while getopts 'i:m:' flag; do
	case "${flag}" in
		i) VCF="${OPTARG}" ;;
		m) MASK="${OPTARG}" ;;
		*) print_usage
			exit 1 ;;
	esac
done

PREFIX=${VCF%.*.*}
export PATH=/mnt/tank/scratch/atomarovsky/tools/bcftools-1.15.1/bin:/mnt/tank/scratch/atomarovsky/tools/bedtools-2.30.0/bin:/mnt/tank/scratch/atomarovsky/tools/vcftools-0.1.16/bin:${PATH}

# ---- script ----
mkdir ${PREFIX}.bcftools_filtration/
cd ${PREFIX}.bcftools_filtration/

echo "VCF filtration"
bcftools filter --threads 4 -S . -O z -o ${PREFIX}.filt.vcf.gz --exclude 'QUAL < 20.0 || (FORMAT/SP > 60.0 || FORMAT/DP < 5.0 || FORMAT/GQ < 20.0)' ../${VCF}

echo "Sample separation"
touch samples.tmp
for SAMPLE in `bcftools query -l ${PREFIX}.filt.vcf.gz`; do
	echo "Sample: ${SAMPLE}";
	bcftools view --threads 4 --min-ac 1 --with-header -s ${SAMPLE} -O z -o ${SAMPLE}.${PREFIX}.filt.vcf.gz ${PREFIX}.filt.vcf.gz;
	bedtools intersect -header -v -a ${SAMPLE}.${PREFIX}.filt.vcf.gz -b ${MASK} > ${SAMPLE}.${PREFIX}.filt.masked.vcf;
	echo ${SAMPLE}.${PREFIX}.filt.masked.vcf >> samples.tmp;
done

for FILE in $(cat samples.tmp); do
	echo "indel and snp";
	bcftools filter --threads 4 -i  'TYPE="indel"' -O z -o ${FILE%.*.*}.indel.vcf.gz ${FILE};
	bcftools filter --threads 4 -i  'TYPE="snp"' -O z -o ${FILE%.*.*}.snp.vcf.gz ${FILE};
	echo "hetero variants"
	bcftools filter --threads 4 -i 'FMT/GT = "het"' -O z -o ${FILE%.*.*}.indel.hetero.vcf.gz ${FILE%.*.*}.indel.vcf.gz;
	bcftools filter --threads 4 -i 'FMT/GT = "het"' -O z -o ${FILE%.*.*}.snp.hetero.vcf.gz ${FILE%.*.*}.snp.vcf.gz;
	echo "homo variants"
	bcftools filter --threads 4 -i 'FMT/GT = "hom"' -O z -o ${FILE%.*.*}.indel.homo.vcf.gz ${FILE%.*.*}.indel.vcf.gz;
	bcftools filter --threads 4 -i 'FMT/GT = "hom"' -O z -o ${FILE%.*.*}.snp.homo.vcf.gz ${FILE%.*.*}.snp.vcf.gz;
done

rm samples.tmp



