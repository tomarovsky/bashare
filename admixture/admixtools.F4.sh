#!/bin/bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
    echo "Usage: $0 <prefix> <W:samples> <X:samples> <Y:samples> <Z:samples> <outfile>"
    echo "Example:"
    echo "  $0 mzib.mfoi.allsamples.filt.mask.auto.snp.plink \\"
    echo "     hybrid:T84,T87 \\"
    echo "     mzib:10xmzib,S26,T8,T26,T50,T72,T90,T104,T118,T148,T150,T194,china \\"
    echo "     mmar:10xmmar,S44,S46,S49,T149,T24,T76,T77,T82 \\"
    echo "     outgroup:10xmfoi \\"
    echo "     Dstat_hybrids_vs_mzib_mmar.txt"
    exit 1
fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate admixtools

PREFIX=$1
W=$2
X=$3
Y=$4
Z=$5
OUTFILE=$6

GENO_FILE="${PREFIX}.geno"
SNP_FILE="${PREFIX}.snp"
IND_FILE="${PREFIX}.ind"

# --- Temporary directory ---
TMPDIR="${TMPDIR:-/tmp}/admixtools_D_$(date +%s)_$$"
mkdir -p "$TMPDIR"
trap "rm -rf ${TMPDIR}" EXIT

# --- Parse input groups ---
parse_group() {
    local arg="$1"
    local group_name=$(echo "$arg" | cut -d':' -f1)
    local group_samples=$(echo "$arg" | cut -d':' -f2)
    echo "$group_name" "$group_samples"
}

read W_GROUP_NAME W_GROUP_SAMPLES <<< "$(parse_group "$W")"
read X_GROUP_NAME X_GROUP_SAMPLES <<< "$(parse_group "$X")"
read Y_GROUP_NAME Y_GROUP_SAMPLES <<< "$(parse_group "$Y")"
read Z_OUTGROUP_NAME Z_OUTGROUP_SAMPLES <<< "$(parse_group "$Z")"

IFS=',' read -ra W_GROUP_ARR <<< "$W_GROUP_SAMPLES"
IFS=',' read -ra X_GROUP_ARR <<< "$X_GROUP_SAMPLES"
IFS=',' read -ra Y_GROUP_ARR <<< "$Y_GROUP_SAMPLES"
IFS=',' read -ra Z_OUTGROUP_ARR <<< "$Z_OUTGROUP_SAMPLES"

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
    if [[ " ${W_GROUP_ARR[*]} " =~ " $sample " ]]; then
        pop="$W_GROUP_NAME"
    elif [[ " ${X_GROUP_ARR[*]} " =~ " $sample " ]]; then
        pop="$X_GROUP_NAME"
    elif [[ " ${Y_GROUP_ARR[*]} " =~ " $sample " ]]; then
        pop="$Y_GROUP_NAME"
    elif [[ " ${Z_OUTGROUP_ARR[*]} " =~ " $sample " ]]; then
        pop="$Z_OUTGROUP_NAME"
    else
        pop="ignore"
    fi
    echo -e "${sample}\tU\t${pop}" >> "$IND_TEMP"
done

# --- Create poplist.txt ---
echo -e "${W_GROUP_NAME}\t${X_GROUP_NAME}\t${Y_GROUP_NAME}\t${Z_OUTGROUP_NAME}\n" > "$POPFILE"

# --- Run D-statistic ---
cat > "$PARFILE" <<EOF
genotypename: $GENO_FILE
snpname:      $SNP_FILE
indivname:    $IND_TEMP
popfilename:  $POPFILE
f4mode:       YES
inbreed:      NO
printsd:      YES
EOF

qpDstat -p "$PARFILE" > "$OUTFILE"

echo "[INFO] Output file: $OUTFILE"

echo -e "W,X,Y,Z(Outgroup),D/F4,stderr,Z-score,ABBA,BABA,SNPs"

for file in "$OUTFILE"; do
    grep "result:" "$file" | while read -r line; do
        line_clean=$(echo "$line" | sed 's/result:[[:space:]]*//')
        echo -e "$(echo "$line_clean" | awk '{$1=$1; print}' | tr -s ' ' ',')"
    done
done
