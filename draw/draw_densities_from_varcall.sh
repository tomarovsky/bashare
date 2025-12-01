#!/bin/bash

set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 FILT_MASK_VCF ASSEMBLY_BASEPATH THREADS"
    exit 1
fi

export FILT_MASK_VCF=$1
export PREFIX=${FILT_MASK_VCF%.filt.mask.vcf.gz}
export ASSEMBLY_BASEPATH=$2
export THREADS=$3


echo "$(date) | Draw variant window densities"

bcftools query -l ${FILT_MASK_VCF} | parallel -j $THREADS '
    SAMPLE={}
    echo "Sample: ${SAMPLE}"

    # Heterozygous SNP density
    $TOOLS/MACE/scripts/draw_variant_window_densities.py \
        -i ${SAMPLE}.${PREFIX}.filt.mask.snp.hetero.vcf.gz \
        -o ${SAMPLE}.${PREFIX}.filt.mask.snp.hetero.w1mb.s100kb \
        --scaffold_ordered_list ${ASSEMBLY_BASEPATH}.orderlist \
        --scaffold_white_list ${ASSEMBLY_BASEPATH}.whitelist \
        --scaffold_length_file ${ASSEMBLY_BASEPATH}.len \
        --scaffold_syn_file ${ASSEMBLY_BASEPATH}.syn \
        -l "HeteroSNP densities for ${SAMPLE}" \
        --density_thresholds 0,0.1,0.5,1,2,3,4,5,6,7 \
        -w 1000000 \
        -s 100000 \
        --subplots_adjust_left 0.2 \
        --hide_track_label \
        --rounded

    # Homozygous SNP density
    $TOOLS/MACE/scripts/draw_variant_window_densities.py \
        -i ${SAMPLE}.${PREFIX}.filt.mask.snp.homo.vcf.gz \
        -o ${SAMPLE}.${PREFIX}.filt.mask.snp.homo.w1mb.s100kb \
        --scaffold_ordered_list ${ASSEMBLY_BASEPATH}.orderlist \
        --scaffold_white_list ${ASSEMBLY_BASEPATH}.whitelist \
        --scaffold_length_file ${ASSEMBLY_BASEPATH}.len \
        --scaffold_syn_file ${ASSEMBLY_BASEPATH}.syn \
        -l "HeteroSNP densities for ${SAMPLE}" \
        --density_thresholds 0,0.1,0.5,1,2,3,4,5,6,7 \
        -w 1000000 \
        -s 100000 \
        --subplots_adjust_left 0.2 \
        --hide_track_label \
        --rounded

    # Heterozygous SNP density (w100kb.s10kb) only counts
    $TOOLS/MACE/scripts/draw_variant_window_densities.py \
        -i ${SAMPLE}.${PREFIX}.filt.mask.snp.homo.vcf.gz \
        -o ${SAMPLE}.${PREFIX}.filt.mask.snp.homo.w100kb.s10kb \
        --scaffold_ordered_list ${ASSEMBLY_BASEPATH}.orderlist \
        --scaffold_white_list ${ASSEMBLY_BASEPATH}.whitelist \
        --scaffold_length_file ${ASSEMBLY_BASEPATH}.len \
        --scaffold_syn_file ${ASSEMBLY_BASEPATH}.syn \
        -w 100000 \
        -s 10000 \
        --only_count

    # ROH (w1mb.s100kb)
    $TOOLS/Biocrutch/scripts/ROH/get_ROH_regions.py \
        -i ${SAMPLE}.${PREFIX}.filt.mask.snp.hetero.w1mb.s100kb.features.bed \
        -o ROH/${SAMPLE}.${PREFIX}.filt.mask.snp.hetero.w1mb.s100kb.features.roh

    # Draw ROH
    $TOOLS/MACE/scripts/draw_features.py \
        -i ROH/${SAMPLE}.${PREFIX}.filt.mask.snp.hetero.w1mb.s100kb.features.roh \
        -o ROH/${SAMPLE}.${PREFIX}.filt.mask.snp.hetero.w1mb.s100kb.features \
        -t bed \
        -l "ROHs for ${SAMPLE}" \
        --scaffold_ordered_list ${ASSEMBLY_BASEPATH}.orderlist \
        --scaffold_length_file ${ASSEMBLY_BASEPATH}.len \
        --scaffold_syn_file ${ASSEMBLY_BASEPATH}.syn \
        --hide_track_label \
        --rounded \
        --subplots_adjust_left 0.2 \
        --figure_width 10 \
        --default_color "tab:blue"

    # ROH (w100kb.s10kb)
    $TOOLS/Biocrutch/scripts/ROH/get_ROH_regions.py \
        -i ${SAMPLE}.${PREFIX}.filt.mask.snp.hetero.w100kb.s10kb.features.bed \
        -o ROH/${SAMPLE}.${PREFIX}.filt.mask.snp.hetero.w100kb.s10kb.features.roh

    $TOOLS/MACE/scripts/draw_features.py \
        -i ROH/${SAMPLE}.${PREFIX}.filt.mask.snp.hetero.w100kb.s10kb.features.roh \
        -o ROH/${SAMPLE}.${PREFIX}.filt.mask.snp.hetero.w100kb.s10kb.features \
        -t bed \
        -l "ROHs for ${SAMPLE}" \
        --scaffold_ordered_list ${ASSEMBLY_BASEPATH}.orderlist \
        --scaffold_length_file ${ASSEMBLY_BASEPATH}.len \
        --scaffold_syn_file ${ASSEMBLY_BASEPATH}.syn \
        --hide_track_label \
        --rounded \
        --subplots_adjust_left 0.2 \
        --figure_width 10 \
        --default_color "tab:blue"
'
