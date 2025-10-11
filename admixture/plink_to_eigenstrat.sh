#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 PLINK_PREFIX"
    exit 1
fi

PLINK_PREFIX="$1"

export PATH=${PATH}:/mnt/tank/scratch/atomarovsky/tools/AdmixTools-7.0.2/bin

# Create parameter file for convertf
PARFILE="plink_to_eigenstrat.par"
cat > "$PARFILE" <<EOF
genotypename:   ${PLINK_PREFIX}.bed
snpname:        ${PLINK_PREFIX}.bim
indivname:      ${PLINK_PREFIX}.fam
outputformat:   EIGENSTRAT
genotypeoutname: ${PLINK_PREFIX}.geno
snpoutname:     ${PLINK_PREFIX}.snp
indivoutname:   ${PLINK_PREFIX}.ind
EOF

# Run conversion
echo "[INFO] PLINK → EIGENSTRAT..."
convertf -p "$PARFILE" > "${OUTDIR}/${BASENAME}.convertf.log"
