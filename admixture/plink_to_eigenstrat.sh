#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 PLINK_PREFIX"
    echo "Converts a PLINK dataset to an Eigenstrat dataset: .geno, .snp, .ind"
    exit 1
fi

PLINK_PREFIX="$1"

source $(conda info --base)/etc/profile.d/conda.sh
conda activate admixtools

# Create parameter file for convertf
PARFILE="plink_to_eigenstrat.par"
cat > "$PARFILE" <<EOF
genotypename:    ${PLINK_PREFIX}.bed
snpname:         ${PLINK_PREFIX}.bim
indivname:       ${PLINK_PREFIX}.fam
outputformat:    EIGENSTRAT
genotypeoutname: ${PLINK_PREFIX}.geno
snpoutname:      ${PLINK_PREFIX}.snp
indivoutname:    ${PLINK_PREFIX}.ind
EOF

# Run conversion
echo "[INFO] PLINK → EIGENSTRAT..."
convertf -p "$PARFILE" > "${PLINK_PREFIX}.convertf.log"

# Fix .ind
IND_FILE="${PLINK_PREFIX}.ind"
cp "$IND_FILE" "${IND_FILE}.bak"

# "10xmzib:10xmzib U ???" → "10xmzib U 10xmzib"
awk '{
    # remove possible colons
    gsub(/:.*/, "", $1);
    sample=$1;
    sex=($2=="") ? "U" : $2;
    printf "%s\t%s\t%s\n", sample, sex, sample;
}' "${IND_FILE}.bak" > "$IND_FILE"
