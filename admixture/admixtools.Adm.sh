#!/bin/bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: $0 <prefix> <left_groups> <right_groups> <target_group> <outfile>"
    echo ""
    echo "  left_groups  : multiple groups separated by ';', each in format name:sample1,sample2,..."
    echo "  right_groups : multiple groups separated by ';', each in format name:sample1,sample2,..."
    echo "  target_group : single group in format name:sample1,sample2,..."
    echo ""
    echo "Example:"
    echo "  $0 mzib.mfoi.allsamples.filt.mask.auto.snp.plink \\"
    echo "     Martes_martes:S1,S2,S3;Martes_zibellina:T1,T2,T3 \\"
    echo "     Martes_foina:R1,R2 \\"
    echo "     hybrids:T84,T87 \\"
    echo "     qpAdm_results.txt"
    exit 1
fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate admixtools

PREFIX=$1
LEFT=$2
RIGHT=$3
TARGET=$4
OUTFILE=$5

GENO_FILE="${PREFIX}.geno"
SNP_FILE="${PREFIX}.snp"
IND_FILE="${PREFIX}.ind"

# --- Temporary directory ---
TMPDIR="${TMPDIR:-/tmp}/admixtools_qpAdm_$(date +%s)_$$"
mkdir -p "$TMPDIR"
trap "rm -rf ${TMPDIR}" EXIT

# --- Function to parse groups ---
parse_group() {
    local arg="$1"
    local name=$(echo "$arg" | cut -d':' -f1)
    local samples=$(echo "$arg" | cut -d':' -f2)
    IFS=',' read -ra arr <<< "$samples"
    echo "$name" "${arr[*]}"
}

# --- Parse left groups ---
declare -A LEFT_MAP
IFS=';' read -ra LEFT_GROUPS <<< "$LEFT"
for group in "${LEFT_GROUPS[@]}"; do
    read name samples <<< "$(parse_group "$group")"
    LEFT_MAP[$name]="$samples"
done

# --- Parse right groups ---
declare -A RIGHT_MAP
IFS=';' read -ra RIGHT_GROUPS <<< "$RIGHT"
for group in "${RIGHT_GROUPS[@]}"; do
    read name samples <<< "$(parse_group "$group")"
    RIGHT_MAP[$name]="$samples"
done

# --- Parse target group ---
read TARGET_NAME TARGET_SAMPLES <<< "$(parse_group "$TARGET")"
IFS=' ' read -ra TARGET_ARR <<< "$TARGET_SAMPLES"

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
    assigned="ignore"
    # Check left
    for pop in "${!LEFT_MAP[@]}"; do
        if [[ " ${LEFT_MAP[$pop]} " =~ " $sample " ]]; then
            assigned="$pop"
            break
        fi
    done
    # Check right
    if [[ "$assigned" == "ignore" ]]; then
        for pop in "${!RIGHT_MAP[@]}"; do
            if [[ " ${RIGHT_MAP[$pop]} " =~ " $sample " ]]; then
                assigned="$pop"
                break
            fi
        done
    fi
    # Check target
    if [[ "$assigned" == "ignore" ]]; then
        if [[ " ${TARGET_ARR[*]} " =~ " $sample " ]]; then
            assigned="$TARGET_NAME"
        fi
    fi
    echo -e "${sample}\tU\t${assigned}" >> "$IND_TEMP"
done

# --- Create poplist.txt ---
{
    echo "$TARGET_NAME"
    for pop in "${!LEFT_MAP[@]}"; do
        echo "$pop"
    done
    for pop in "${!RIGHT_MAP[@]}"; do
        echo "$pop"
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

# --- Optional: extract admixture proportions ---
grep -E "left pop|admixture proportions|stderr" "$OUTFILE"
