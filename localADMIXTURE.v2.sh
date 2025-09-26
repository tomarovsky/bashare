#!/bin/bash
# Usage:
# $TOOLS/bashare/localADMIXTURE.sh VCF REGION_FILE "SCAFFOLDS"

VCF=$1 # absolute path
REGION_FILE=$2 # absolute path mzib.autosomes_and_PAR.w1mb.s100kb.features.bed
SCAFFOLDS=$3

source $(conda info --base)/etc/profile.d/conda.sh

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 "
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
                echo -e "HiC_scafflod_${REG}\t${ADMIXTURE_OUTLINE}" >> ${SAMPLE}.${CHR}.admixture.concat.Q;
            done;
            cat ${CHR}.${REGION}.admixture.log >> ${SAMPLE}.${CHR}.admixture.concat.log;
            rm *plink* *.vcf *.admixture.log;
        else
            echo 'WARNING! Empty PLINK results!';
            for SAMPLE in {1..33}; do
                echo -e "HiC_scafflod_${REG}\tEMPTY\tEMPTY" >> ${SAMPLE}.${CHR}.admixture.concat.Q;
            done;
            rm *plink* *.vcf;
        conda deactivate;
        fi;
    done;
    cd ..;
done
