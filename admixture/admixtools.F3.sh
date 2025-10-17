#!/bin/bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: $0 <prefix> <group1:samples> <group2:samples> <hybrid:samples> <outfile>"
    echo "Note: Outgroup not needed for F3 statistic"
    echo "Example:"
    echo "  $0 mzib.mfoi.allsamples.filt.mask.auto.snp.plink \\"
    echo "     mzib:10xmzib,S26,T8,T26,T50,T72,T90,T104,T118,T148,T150,T194,china \\"
    echo "     mmar:10xmmar,S44,S46,S49,T149,T24,T76,T77,T82 \\"
    echo "     hybrid:T84,T87 \\"
    echo "     F3stat_hybrids_vs_mzib_mmar.txt"
    exit 1
fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate admixtools

PREFIX=$1
GROUP1=$2
GROUP2=$3
HYBRID=$4
OUTFILE=$5

GENO_FILE="${PREFIX}.geno"
SNP_FILE="${PREFIX}.snp"
IND_FILE="${PREFIX}.ind"

# --- Temporary directory ---
TMPDIR="${TMPDIR:-/tmp}/admixtools_F3_$(date +%s)_$$"
mkdir -p "$TMPDIR"
trap "rm -rf ${TMPDIR}" EXIT

# --- Parse input groups ---
parse_group() {
    local arg="$1"
    local group_name=$(echo "$arg" | cut -d':' -f1)
    local group_samples=$(echo "$arg" | cut -d':' -f2)
    echo "$group_name" "$group_samples"
}

read GROUP1_NAME GROUP1_SAMPLES <<< "$(parse_group "$GROUP1")"
read GROUP2_NAME GROUP2_SAMPLES <<< "$(parse_group "$GROUP2")"
read HYBRID_NAME HYBRID_SAMPLES <<< "$(parse_group "$HYBRID")"

IFS=',' read -ra GROUP1_ARR <<< "$GROUP1_SAMPLES"
IFS=',' read -ra GROUP2_ARR <<< "$GROUP2_SAMPLES"
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
    if [[ " ${GROUP1_ARR[*]} " =~ " $sample " ]]; then
        pop="$GROUP1_NAME"
    elif [[ " ${GROUP2_ARR[*]} " =~ " $sample " ]]; then
        pop="$GROUP2_NAME"
    elif [[ " ${HYBRID_ARR[*]} " =~ " $sample " ]]; then
        pop="$HYBRID_NAME"
    else
        pop="ignore"
    fi
    echo -e "${sample}\tU\t${pop}" >> "$IND_TEMP"
done

# --- Create poplist.txt ---
echo -e "${GROUP1_NAME}\t${GROUP2_NAME}\t${HYBRID_NAME}\n" > "$POPFILE"

# --- Run F3-statistic ---
cat > "$PARFILE" <<EOF
genotypename: $GENO_FILE
snpname:      $SNP_FILE
indivname:    $IND_TEMP
popfilename:  $POPFILE
inbreed:      NO
EOF

echo "=== [INFO] qp3Pop -> F3-statistic ==="
qp3Pop -p "$PARFILE" > "$OUTFILE"

echo "[INFO] Output file: $OUTFILE"
