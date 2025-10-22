#!/bin/bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
    echo "Usage: $0 <prefix> <source1:samples> <source2:samples> <hybrid:samples> <outgroups:samples> <outfile>"
    echo "Example:"
    echo "  $0 mzib.mfoi.allsamples.filt.mask.auto.snp.plink \\"
    echo "     mzib:10xmzib,S26,T8,T26,T50,T72,T90,T104,T118,T148,T150,T194,china \\"
    echo "     mmar:10xmmar,S44,S46,S49,T149,T24,T76,T77,T82 \\"
    echo "     hybrid:T84,T87 \\"
    echo "     outgroups:Mfoina,Siberian_marten \\"
    echo "     qpAdm_results.txt"
    exit 1
fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate admixtools

PREFIX=$1
SOURCE1=$2
SOURCE2=$3
HYBRID=$4
OUTGROUPS=$5
OUTFILE=$6

GENO_FILE="${PREFIX}.geno"
SNP_FILE="${PREFIX}.snp"
IND_FILE="${PREFIX}.ind"

# --- Temporary directory ---
TMPDIR="${TMPDIR:-/tmp}/admixtools_qpAdm_$(date +%s)_$$"
mkdir -p "$TMPDIR"
trap "rm -rf ${TMPDIR}" EXIT

# --- Parse input groups ---
parse_group() {
    local arg="$1"
    local group_name=$(echo "$arg" | cut -d':' -f1)
    local group_samples=$(echo "$arg" | cut -d':' -f2)
    echo "$group_name" "$group_samples"
}

read SOURCE1_NAME SOURCE1_SAMPLES <<< "$(parse_group "$SOURCE1")"
read SOURCE2_NAME SOURCE2_SAMPLES <<< "$(parse_group "$SOURCE2")"
read HYBRID_NAME HYBRID_SAMPLES <<< "$(parse_group "$HYBRID")"
read -a OUTGROUPS_ARR <<< $(echo "$OUTGROUPS" | tr ',' ' ')

IFS=',' read -ra SOURCE1_ARR <<< "$SOURCE1_SAMPLES"
IFS=',' read -ra SOURCE2_ARR <<< "$SOURCE2_SAMPLES"
IFS=',' read -ra HYBRID_ARR <<< "$HYBRID_SAMPLES"

# --- Read original .ind ---
ORIG_IND=()
while read -r line; do
    sample=$(echo "$line" | awk '{print $1}')
    ORIG_IND+=("$sample")
done < "$IND_FILE"

IND_TEMP="${TMPDIR}/dataset.ind"
POPFILE="${TMPDIR}/poplist.txt"
PARFILE="${TMPDIR}/parfile.par"

# --- Create .ind ---
> "$IND_TEMP"
for sample in "${ORIG_IND[@]}"; do
    if [[ " ${SOURCE1_ARR[*]} " =~ " $sample " ]]; then
        pop="$SOURCE1_NAME"
    elif [[ " ${SOURCE2_ARR[*]} " =~ " $sample " ]]; then
        pop="$SOURCE2_NAME"
    elif [[ " ${HYBRID_ARR[*]} " =~ " $sample " ]]; then
        pop="$HYBRID_NAME"
    elif [[ " ${OUTGROUPS_ARR[*]} " =~ " $sample " ]]; then
        pop="$sample"
    else
        pop="ignore"
    fi
    echo -e "${sample}\tU\t${pop}" >> "$IND_TEMP"
done

# --- Create poplist.txt ---
{
    echo "$HYBRID_NAME"
    echo "$SOURCE1_NAME"
    echo "$SOURCE2_NAME"
    for og in "${OUTGROUPS_ARR[@]}"; do
        echo "$og"
    done
} > "$POPFILE"

# --- Create parfile ---
cat > "$PARFILE" <<EOF
genotypename: $GENO_FILE
snpname:      $SNP_FILE
indivname:    $IND_TEMP
popfilename:  $POPFILE
details:      YES
allsnps:      YES
inbreed:      NO
EOF

# --- Run qpAdm ---
qpAdm -p "$PARFILE" > "$OUTFILE"
