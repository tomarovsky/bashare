#!/bin/bash
# Usage:
# $TOOLS/bashare/localADMIXTURE.sh VCF REGION_FILE "SCAFFOLDS"
# for S in {1..18}; do $TOOLS/bashare/localADMIXTURE.sh mzib.allsamples.filt.masked.auto.snp.plink.0.2.vcf.gz mzib.w1mb.s100kb.features.auto.num.bed $S > /dev/null 2>&1 & done


VCF=$1 # absolute path
REGION_FILE=$2 # absolute path
SCAFFOLDS=$3

source $(conda info --base)/etc/profile.d/conda.sh

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 VCF REGION_FILE 'SCAFFOLDS'"
    exit 1
fi

mkdir -p w1mb.s100kb;
cd w1mb.s100kb;

for CHR in $SCAFFOLDS; do
    mkdir ${CHR};
    cd ${CHR};
    cat ${REGION_FILE} | grep -P "^$CHR\t" | while read REG; do
        echo $REG;
        REGION=$(echo $REG | awk '{print $1":"$2"-"$3}')
        echo $REGION;
        conda activate varcall;
        bcftools view --threads 1 --regions ${REGION} -O v -o ${CHR}.${REGION}.vcf ${VCF};
        conda deactivate && conda activate plink1.9;
        if plink --vcf ${CHR}.${REGION}.vcf --out ${CHR}.${REGION}.plink --make-bed --threads 1; then
            conda deactivate && conda activate admixture;
            admixture -j1 --cv ${CHR}.${REGION}.plink.bed 2 | tee ${CHR}.${REGION}.admixture.log;
            for SAMPLE in {1..33}; do
                ADMIXTURE_OUTLINE=$(sed "${SAMPLE}q;d" ${CHR}.${REGION}.plink.2.Q | awk '{print $1"\t"$2}');
                echo $ADMIXTURE_OUTLINE;
                echo -e "HiC_scaffold_${REG}\t${ADMIXTURE_OUTLINE}" >> ${SAMPLE}.HiC_scaffold_${CHR}.admixture.concat.Q;
            done;
            rm *plink* *.vcf *.admixture.log;
        else
            echo 'WARNING! Empty PLINK results!';
            for SAMPLE in {1..33}; do
                echo -e "HiC_scaffold_${REG}\tEMPTY\tEMPTY" >> ${SAMPLE}.HiC_scaffold_${CHR}.admixture.concat.Q;
            done;
            rm *plink* *.vcf;
        conda deactivate;
        fi;
    done;
    cd ..;
done
