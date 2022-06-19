#!/bin/bash

BWA_ALIGNMENT="/mnt/tank/scratch/skliver/common/mustelidae/martes_martes/genome/alignment/hic.purged/bwa_alignment"
VCF="mmar.allsamples.vcf.gz"
PREFIX=${VCF%.*.*}

#---- script ----
mkdir bcftools_filtration/
cd bcftools_filtration/

echo "VCF filtration"
bcftools filter -S . -O v -o ${PREFIX}.filt.vcf --exclude 'QUAL < 20.0 || (FORMAT/SP > 60.0 | FORMAT/DP < 5.0 | FORMAT/GQ < 20.0)' ../${VCF}

echo "Sample separation"
touch samples.tmp
for SAMPLE in `bcftools query -l ${PREFIX}.filt.vcf`; do
	echo "Sample: ${SAMPLE}";
	vcf-subset --exclude-ref -c ${SAMPLE} ${PREFIX}.filt.vcf > ${SAMPLE}.${PREFIX}.filt.vcf;
	bedtools intersect -header -v -a ${SAMPLE}.${PREFIX}.filt.vcf -b ${BWA_ALIGNMENT}/${SAMPLE}/${SAMPLE}.mmar.hic.purged.mkdup.per-base.bed.gz.max250.min33 > ${SAMPLE}.${PREFIX}.filt.masked.vcf
	echo ${SAMPLE}.${PREFIX}.filt.masked.vcf >> samples.tmp;
done

for FILE in $(cat samples.tmp); do
	echo "indel and snp";
	bcftools  filter -i  'TYPE="indel"' -O v ${FILE} > ${FILE%.*}.indel.vcf;
	bcftools  filter -i  'TYPE="snp"' -O v ${FILE} > ${FILE%.*}.snp.vcf;
	echo "hetero variants"
	bcftools filter -i 'FMT/GT = "het"' -O v ${FILE%.*}.indel.vcf > ${FILE%.*}.indel.hetero.vcf;
	bcftools filter -i 'FMT/GT = "het"' -O v ${FILE%.*}.snp.vcf > ${FILE%.*}.snp.hetero.vcf;
	echo "homo variants"
	bcftools filter -i 'FMT/GT = "hom"' -O v ${FILE%.*}.indel.vcf > ${FILE%.*}.indel.homo.vcf;
	bcftools filter -i 'FMT/GT = "hom"' -O v ${file%.*}.snp.vcf > ${file%.*}.snp.homo.vcf;
done

rm samples.tmp