#!/bin/bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: $0 <prefix> <group1_samples> <group2_samples> <hybrid_sample> <outgroup_sample>"
    echo "Example:"
    echo "  $0 mzib.mfoi.allsamples.filt.mask.auto.snp.plink \\"
    echo "     10xmzib,S26,T8,T26,T50,T72,T90,T104,T118,T148,T150,T194,china \\"
    echo "     10xmmar,S44,S46,S49,T149,T24,T76,T77,T82 \\"
    echo "     T87 \\"
    echo "     10xmfoi"
    exit 1
fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate admixtools

PREFIX=$1
GROUP1_SAMPLES=$2
GROUP2_SAMPLES=$3
HYBRID=$4
OUTGROUP=$5

GENO_FILE="${PREFIX}.geno"
SNP_FILE="${PREFIX}.snp"
IND_FILE="${PREFIX}.ind"

IND_TEMP="dataset.ind"
POPFILE="poplist.txt"
PARFILE="qpDstat.par"

# --- Преобразуем строки с запятыми в массивы ---
IFS=',' read -ra GROUP1_ARR <<< "$GROUP1_SAMPLES"
IFS=',' read -ra GROUP2_ARR <<< "$GROUP2_SAMPLES"
IFS=',' read -ra HYBRID_ARR <<< "$HYBRID"
OUTGROUP_ARR=("$OUTGROUP")

# --- Считываем все образцы из оригинального .ind ---
ORIG_IND=()
while read -r line; do
    sample=$(echo "$line" | awk '{print $1}')
    ORIG_IND+=("$sample")
done < "$IND_FILE"

# --- Создаём новый .ind ---
> "$IND_TEMP"

for sample in "${ORIG_IND[@]}"; do
    if [[ " ${GROUP1_ARR[*]} " =~ " $sample " ]]; then
        pop="mzib"
    elif [[ " ${GROUP2_ARR[*]} " =~ " $sample " ]]; then
        pop="mmar"
    elif [[ " ${HYBRID_ARR[*]} " =~ " $sample " ]]; then
        pop="hybrid"
    elif [[ " ${OUTGROUP_ARR[*]} " =~ " $sample " ]]; then
        pop="outgroup"
    else
        pop="ignore"
    fi
    echo -e "${sample}\tU\t${pop}" >> "$IND_TEMP"
done

# --- Создаём poplist.txt ---
echo -e "mzib\tmmar\thybrid\toutgroup\n" > "$POPFILE"

# --- Создаём парфайл для qpDstat ---
cat > "$PARFILE" <<EOF
genotypename: $GENO_FILE
snpname:      $SNP_FILE
indivname:    $IND_TEMP
popfilename:  $POPFILE
f4mode:       NO
EOF

# --- Запуск qpDstat ---
OUTFILE="Dstat_${HYBRID}_vs_mzib_mmar.txt"
echo "=== Расчёт D-статистики для гибрида $HYBRID ==="
qpDstat -p "$PARFILE" > "$OUTFILE"
echo "Результаты сохранены: $OUTFILE"
