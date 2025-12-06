#!/bin/bash
set -euo pipefail

source "$TOOLS/bashare/lib/log_functions.sh"

VCF=$1
OUTPREFIX=$2
LD=$3
THREADS=$4

if [[ $# -ne 4 ]]; then
    log_error "Usage: $0 VCF OUTPREFIX LD THREADS"
    exit 1
fi

# Step 1: geno and maf
log_info "Step 1: geno and maf"
plink --vcf $VCF --out ${OUTPREFIX}.prefiltered \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --geno 0.05 \
    --maf 0.03 \
    --snps-only \
    --make-bed \
    --threads $THREADS |& tee -a ${OUTPREFIX}.1.log

# Step 2: Scaffold renaming
log_info "Step 2: Renaming scaffolds to integers"

cp ${OUTPREFIX}.prefiltered.bim ${OUTPREFIX}.prefiltered.bim.raw

awk -v synfile="${OUTPREFIX}.prefiltered.bim.syn" '
BEGIN {
    OFS="\t";
    count=0;
}
{
    # $1
    if (!($1 in map)) {
        count++;
        map[$1] = count;    # Запоминаем mapping
        print $1, count > synfile
    }

    $1 = map[$1];
    print $0;
}' ${OUTPREFIX}.prefiltered.bim.raw > ${OUTPREFIX}.prefiltered.bim

# Step 3: LD pruning
log_info "Step 3: LD pruning"

plink --bfile ${OUTPREFIX}.prefiltered --out ${OUTPREFIX} \
    --indep-pairwise 50 10 $LD \
    --allow-extra-chr \
    --threads $THREADS |& tee -a ${OUTPREFIX}.2.log

# Step 4: PCA
log_info "Step 4: PCA"

plink --bfile ${OUTPREFIX}.prefiltered --out ${OUTPREFIX} \
    --extract ${OUTPREFIX}.prune.in \
    --allow-extra-chr \
    --pca \
    --make-bed \
    --threads $THREADS |& tee -a ${OUTPREFIX}.3.log

log_info "Done!"
