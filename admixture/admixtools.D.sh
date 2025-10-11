#!/bin/bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: $0 <prefix> <group1_samples> <group2_samples> <hybrid_sample> <outgroup_sample>"
    echo "Example:"
    echo "  $0 mzib.mfoi.allsamples.filt.namefix.mask.aut_and_PAR.snp.plink \\"
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

# --- Создаём временную директорию ---
TMPDIR=$(mktemp -d)
IND_TEMP="${TMPDIR}/dataset.ind"
POPFILE="${TMPDIR}/poplist.txt"

# --- Генерация .ind файла ---
echo "=== Генерация .ind файла ==="
> "$IND_TEMP"

# Функция для добавления группы
add_group() {
    local samples=$1
    local popname=$2
    IFS=',' read -ra arr <<< "$samples"
    for s in "${arr[@]}"; do
        echo -e "${s}\tU\t${popname}" >> "$IND_TEMP"
    done
}

add_group "$GROUP1_SAMPLES" "mzib"
add_group "$GROUP2_SAMPLES" "mmar"
add_group "$HYBRID" "hybrid"
add_group "$OUTGROUP" "outgroup"

# --- Генерация poplist.txt ---
echo "mzib   mmar   hybrid   outgroup" > "$POPFILE"

# --- Конфигурация qpDstat ---
PARFILE="${TMPDIR}/qpDstat.par"
cat > "$PARFILE" <<EOF
genotypename: $GENO_FILE
snpname:      $SNP_FILE
indivname:    $IND_TEMP
popfilename:  $POPFILE
f4mode:       YES
EOF

# --- Запуск qpDstat ---
echo "=== Расчёт D-статистики для гибрида $HYBRID ==="
OUTFILE="Dstat_${HYBRID}_vs_mzib_mmar.txt"

# временно заменяем "hybrid" на имя гибрида в poplist
sed "s/hybrid/$HYBRID/" "$POPFILE" > "${TMPDIR}/popfile_temp.txt"
qpDstat -p <(sed "s/popfilename:.*/popfilename: ${TMPDIR}\/popfile_temp.txt/" "$PARFILE") > "$OUTFILE"

echo "Результаты сохранены: $OUTFILE"

# --- Уборка ---
rm -r "$TMPDIR"
