#!/bin/bash

set -euo pipefail

# GATK
export PATH=$(conda info --base)/envs/gatk/bin/:${TOOLS}/gatk-4.6.2.0/:${PATH}

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 VCF MASK ASSEMBLY THREADS"
    exit 1
fi

export VCF=$1 # Output from GenotypeGVCFs
export MASK=$2 # BED file
export ASSEMBLY=$3 # Reference FASTA
export THREADS=$4

export PREFIX=${VCF%.*.*}

mkdir -p gatk_filtration/ROH/
cd gatk_filtration/

echo "Step 1: Marking variants with GATK Hard Filters"
# QD < 2.0 - Low quality score normalized by depth. Помогает отсеять ложноположительные варианты, которые имеют высокий QUAL только из-за огромной глубины покрытия, но на самом деле плохи.
# QUAL < 20.0 - Low raw quality score of the variant.
# SOR > 3.0 - (Strand Odds Ratio) Strand bias score.
# FS > 60.0 - Fisher Strand for SNPs - Extreme strand bias specifically for SNPs (Phred-scaled).
# FS > 200.0 - Fisher Strand for Indels - Extreme strand bias specifically for insertions/deletions (Phred-scaled).
# MQ < 40.0 - Mapping quality of the supporting reads.
#
# Genotype-Level Filter is applied to individual sample genotypes within a variant.
# FAIL_GT: DP < 5 || GQ < 20 - Low read depth or low genotype quality for a specific sample.
# Note: Filtered variants and genotypes are marked as FAIL in the output VCF but are not removed.

# --filter-expression "MQ < 40.0" --filter-name "MQ40" \
# --filter-expression "QD < 2.0" --filter-name "QD2" \

gatk --java-options "-Xmx8g" VariantFiltration \
    -R ${ASSEMBLY} \
    -V ../${VCF} \
    -O ${PREFIX}.marked.vcf.gz \
    --filter-expression "QUAL < 20.0" --filter-name "QUAL20" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3" \
    --filter-expression "vc.isSNP() && FS > 60.0" --filter-name "SNP_FS60" \
    --filter-expression "vc.isIndel() && FS > 200.0" --filter-name "INDEL_FS200" \
    --genotype-filter-expression "DP < 5 || GQ < 20" \
    --genotype-filter-name "FAIL_GT"


echo "$(date) | Step 2: Applying filters and Masking"
# --exclude-filtered: Remove sites that failed any filter
# --exclude-intervals: Masking by BED
# --set-filtered-gt-to-nocall: If one sample has DP < 5 (FAIL_GT),
# its genotype will become "./.", but the site will be kept if other samples are fine.

gatk --java-options "-Xmx8g" SelectVariants \
    -R ${ASSEMBLY} \
    -V ${PREFIX}.marked.vcf.gz \
    -O ${PREFIX}.filt.mask.vcf.gz \
    --exclude-intervals ${MASK} \
    --exclude-filtered \
    --set-filtered-gt-to-nocall

echo "$(date) | Step 3: Splitting samples"

SAMPLES=$(bcftools query -l ${PREFIX}.filt.mask.vcf.gz)

export ASSEMBLY PREFIX

process_sample() {
    SAMPLE=$1
    echo "Processing: ${SAMPLE}"

    # 1. Extract Sample & SNP Only
    gatk --java-options '-Xmx4g' SelectVariants \
        -R "${ASSEMBLY}" \
        -V "${PREFIX}.filt.mask.vcf.gz" \
        -O "${SAMPLE}.${PREFIX}.filt.mask.snp.vcf.gz" \
        -sn "${SAMPLE}" \
        -select-type SNP \
        --exclude-non-variants

    # 2. Heterozygous
    gatk --java-options '-Xmx4g' SelectVariants \
        -R "${ASSEMBLY}" \
        -V "${SAMPLE}.${PREFIX}.filt.mask.snp.vcf.gz" \
        -O "${SAMPLE}.${PREFIX}.filt.mask.snp.hetero.vcf.gz" \
        -select 'vc.getGenotype("'${SAMPLE}'").isHet()'

    # 3. Homozygous Alternative
    gatk --java-options '-Xmx4g' SelectVariants \
        -R "${ASSEMBLY}" \
        -V "${SAMPLE}.${PREFIX}.filt.mask.snp.vcf.gz" \
        -O "${SAMPLE}.${PREFIX}.filt.mask.snp.homo.vcf.gz" \
        -select 'vc.getGenotype("'${SAMPLE}'").isHom()'
}
export -f process_sample

echo "$SAMPLES" | parallel -j $THREADS process_sample {}

echo "$(date) | Done."
