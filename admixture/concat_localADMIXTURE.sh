#!/bin/bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 4 ]]; then
    echo "Usage: $0 VCF_FILE REGION_FILE [THREADS] [SCAFFOLD_PREFIX]"
    echo "Example: $0 mzib.allsamples.filt.masked.auto.snp.vcf.gz mzib.w1mb.s100kb.features.auto.num.bed 18 \"HiC_scaffold_\""
    exit 1
fi

VCF=$1
REGION_FILE=$2
THREADS=${3:-18}
SCAFFOLD_PREFIX=${4:-"HiC_scaffold_"}

source $(conda info --base)/etc/profile.d/conda.sh

mkdir -p figures/

SAMPLES=($(conda run -n varcall bcftools query -l "$VCF"))
SCAFFOLDS=($(awk '{print $1}' "$REGION_FILE" | sort | uniq))

echo "Samples: ${SAMPLES[@]}"
echo "Scaffolds: ${SCAFFOLDS[@]}"

export SCAFFOLDS
export SCAFFOLD_PREFIX

process_sample() {
    local sample="$1"
    local idx="$2"
    local outfile="figures/${sample}.local_admixture.features.bed"
    : > "$outfile"

    echo "Processing $sample (index $idx)"

    for scaf in "${SCAFFOLDS[@]}"; do
        local file="${scaf}/${idx}.${SCAFFOLD_PREFIX}${scaf}.admixture.concat.Q"
        [[ -f "$file" ]] || { echo "Warning: missing $file" >&2; continue; }
        awk -v chr="$scaf" -v sample="$sample" -v prefix="$SCAFFOLD_PREFIX" 'BEGIN{OFS="\t"} {print prefix chr, $1, $2, $3, sample}' "$file" >> "$outfile"
    done
}

export -f process_sample

parallel -j "$THREADS" process_sample {1} {#} ::: "${SAMPLES[@]}"
