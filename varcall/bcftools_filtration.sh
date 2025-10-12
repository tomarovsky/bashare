#!/bin/bash
# Usage:
# $TOOLS/bashare/bcftools_filtration.sh VCF MASK ASSEMBLY THREADS

set -ex

export VCF=$1 # vcf.gz
export PREFIX=${VCF%.*.*}
export MASK=$2 # common mask from common_mask.sh
export ASSEMBLY=$3 # absolute path to assembly files
export THREADS=$4

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 VCF MASK ASSEMBLY THREADS"
    exit 1
fi

#---- script ----
mkdir -p bcftools_filtration/ROH/
cd bcftools_filtration/


source $(conda info --base)/etc/profile.d/conda.sh
# mamba create -n varcall bcftools bedtools samtools
conda activate varcall

echo "$(date) | VCF filtration"
bcftools filter --threads 30 -S . -O z -o $PREFIX.filt.vcf.gz --exclude 'QUAL < 20.0 || (FORMAT/SP > 60.0 | FORMAT/DP < 5.0 | FORMAT/GQ < 20.0)' ../${VCF}

echo "$(date) | VCF masking"
bedtools intersect -header -v -a ${PREFIX}.filt.vcf.gz -b ${MASK} | bgzip -c > ${PREFIX}.filt.mask.vcf.gz

echo "$(date) | Sample separation"
bcftools query -l $PREFIX.filt.masked.vcf.gz | parallel -j 6 '
    SAMPLE={}
    echo "Sample: ${SAMPLE}"
    bcftools view --threads 10 --min-ac 1 --with-header -s ${SAMPLE} -O z -o ${SAMPLE}.$PREFIX.filt.masked.vcf.gz $PREFIX.filt.masked.vcf.gz
    bcftools filter --threads 10 -i "TYPE=\"snp\"" -O z -o ${SAMPLE}.$PREFIX.filt.masked.snp.vcf.gz ${SAMPLE}.$PREFIX.filt.masked.vcf.gz
    bcftools filter --threads 10 -i "FMT/GT = \"het\"" -O z -o ${SAMPLE}.$PREFIX.filt.masked.snp.hetero.vcf.gz ${SAMPLE}.$PREFIX.filt.masked.snp.vcf.gz
    bcftools filter --threads 10 -i "FMT/GT = \"hom\"" -O z -o ${SAMPLE}.$PREFIX.filt.masked.snp.homo.vcf.gz ${SAMPLE}.$PREFIX.filt.masked.snp.vcf.gz
    # bcftools filter --threads 10 -i "TYPE=\"indel\"" -O z -o ${SAMPLE}.$PREFIX.filt.masked.indel.vcf.gz ${SAMPLE}.$PREFIX.filt.masked.vcf.gz
    # bcftools filter --threads 10 -i "FMT/GT = \"het\"" -O z -o ${SAMPLE}.$PREFIX.filt.masked.indel.hetero.vcf.gz ${SAMPLE}.$PREFIX.filt.masked.indel.vcf.gz
    # bcftools filter --threads 10 -i "FMT/GT = \"hom\"" -O z -o ${SAMPLE}.$PREFIX.filt.masked.indel.homo.vcf.gz ${SAMPLE}.$PREFIX.filt.masked.indel.vcf.gz
'

conda deactivate && conda activate py38;

echo "$(date) | Draw variant window densities"
bcftools query -l $PREFIX.filt.masked.vcf.gz | parallel -j 64 '
    SAMPLE={}
    echo "Sample: ${SAMPLE}"
    $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${SAMPLE}.$PREFIX.filt.masked.snp.hetero.vcf.gz -o ${SAMPLE}.$PREFIX.filt.masked.snp.hetero.w1mb.s100kb -l "HeteroSNP densities for ${SAMPLE}" -w 1000000 -s 100000 --scaffold_white_list $ASSEMBLY/*d.whitelist -z $ASSEMBLY/*d.orderedlist --scaffold_length_file $ASSEMBLY/*d.len --scaffold_syn_file $ASSEMBLY/*d.syn --syn_file_key_column 0 --syn_file_value_column 1 --hide_track_label --density_thresholds 0,0.1,0.5,1,2,3,4,5,6,7 --rounded --subplots_adjust_left 0.2 --centromere_bed $ASSEMBLY/*.centromere.bed

    $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${SAMPLE}.$PREFIX.filt.masked.snp.homo.vcf.gz -o ${SAMPLE}.$PREFIX.filt.masked.snp.homo.w1mb.s100kb -l "HomoSNP densities for ${SAMPLE}" -w 1000000 -s 100000 --scaffold_white_list $ASSEMBLY/*d.whitelist -z $ASSEMBLY/*d.orderedlist --scaffold_length_file $ASSEMBLY/*d.len --scaffold_syn_file $ASSEMBLY/*d.syn --syn_file_key_column 0 --syn_file_value_column 1 --hide_track_label --density_thresholds 0,0.1,0.5,1,2,3,4,5,6,7 --rounded --subplots_adjust_left 0.2 --centromere_bed $ASSEMBLY/*.centromere.bed

    $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${SAMPLE}.$PREFIX.filt.masked.snp.hetero.vcf.gz -o ${SAMPLE}.$PREFIX.filt.masked.snp.hetero.w100kb.s10kb -w 100000 -s 10000 -a $ASSEMBLY/*d.whitelist -z $ASSEMBLY/*d.orderedlist --scaffold_length_file $ASSEMBLY/*d.len --scaffold_syn_file $ASSEMBLY/*d.syn --syn_file_key_column 0 --syn_file_value_column 1 --only_count

    $TOOLS/Biocrutch/scripts/ROH/get_ROH_regions.py -i ${SAMPLE}.$PREFIX.filt.masked.snp.hetero.w1mb.s100kb.features.bed -o ROH/${SAMPLE}.$PREFIX.filt.masked.snp.hetero.w1mb.s100kb.features.roh
    $TOOLS/MACE/scripts/draw_features.py -i ROH/${SAMPLE}.$PREFIX.filt.masked.snp.hetero.w1mb.s100kb.features.roh -o ROH/${SAMPLE}.$PREFIX.filt.masked.snp.hetero.w1mb.s100kb.features -t bed -l "ROHs for ${SAMPLE}" --scaffold_ordered_list $ASSEMBLY/*d.orderedlist --scaffold_length_file $ASSEMBLY/*d.len --scaffold_syn_file $ASSEMBLY/*d.syn --hide_track_label --rounded --subplots_adjust_left 0.2 --figure_width 10 --default_color "tab:blue" --centromere_bed $ASSEMBLY/*.centromere.bed

    $TOOLS/Biocrutch/scripts/ROH/get_ROH_regions.py -i ${SAMPLE}.$PREFIX.filt.masked.snp.hetero.w100kb.s10kb.features.bed -o ROH/${SAMPLE}.$PREFIX.filt.masked.snp.hetero.w100kb.s10kb.features.roh
    $TOOLS/MACE/scripts/draw_features.py -i ROH/${SAMPLE}.$PREFIX.filt.masked.snp.hetero.w100kb.s10kb.features.roh -o ROH/${SAMPLE}.$PREFIX.filt.masked.snp.hetero.w100kb.s10kb.features -t bed -l "ROHs for ${SAMPLE}" --scaffold_ordered_list $ASSEMBLY/*d.orderedlist --scaffold_length_file $ASSEMBLY/*d.len --scaffold_syn_file $ASSEMBLY/*d.syn --hide_track_label --rounded --subplots_adjust_left 0.2 --figure_width 10  --default_color "tab:blue" --centromere_bed $ASSEMBLY/*.centromere.bed
'
