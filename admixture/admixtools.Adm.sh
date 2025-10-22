#!/bin/bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: $0 <plink_prefix> <left_groups> <right_groups> <target_samples> <outfile>"
    echo "Example:"
    echo "  $0 mzib.mfoi.allsamples.filt.mask.auto.snp.plink \\"
    echo "     \"sables:10xmzib,S26,T8; pines:10xmmar,S44,S46\" \\"
    echo "     \"foina:10xmfoi; out:1344\" \\"
    echo "     \"T18,T87\" \\"
    echo "     qpAdm_results.txt"
    exit 1
fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate admixtools

PREFIX=$1
LEFT_INPUT=$2
RIGHT_INPUT=$3
TARGET_SAMPLES_INPUT=$4
OUTFILE=$5

GENO_FILE="${PREFIX}.geno"
SNP_FILE="${PREFIX}.snp"
IND_FILE="${PREFIX}.ind"

TMPDIR="${TMPDIR:-/tmp}/admixtools_qpAdm_$(date +%s)_$$"
mkdir -p "$TMPDIR"
trap "rm -rf ${TMPDIR}" EXIT

# --- Функция для парсинга групп ---
parse_groups() {
    local arg="$1"
    declare -A grp_samples
    IFS=';' read -ra parts <<< "$arg"
    for part in "${parts[@]}"; do
        part=$(echo "$part" | xargs)  # убираем лишние пробелы
        name=$(echo "$part" | cut -d':' -f1)
        samples=$(echo "$part" | cut -d':' -f2)
        IFS=',' read -ra arr <<< "$samples"
        grp_samples["$name"]="${arr[*]}"
    done
    echo "$(declare -p grp_samples)"
}

# --- Парсим left и right ---
eval "$(parse_groups "$LEFT_INPUT")"
LEFT_NAMES=(${!grp_samples[@]})
LEFT_SAMPLES=()
for n in "${LEFT_NAMES[@]}"; do
    IFS=' ' read -ra tmp <<< "${grp_samples[$n]}"
    LEFT_SAMPLES+=("${tmp[@]}")
done

eval "$(parse_groups "$RIGHT_INPUT")"
RIGHT_NAMES=(${!grp_samples[@]})
RIGHT_SAMPLES=()
for n in "${RIGHT_NAMES[@]}"; do
    IFS=' ' read -ra tmp <<< "${grp_samples[$n]}"
    RIGHT_SAMPLES+=("${tmp[@]}")
done

# --- Target ---
IFS=',' read -ra TARGET_SAMPLES <<< "$TARGET_SAMPLES_INPUT"
TARGET_NAME="target"

# --- Read original .ind ---
ORIG_IND=()
while read -r line; do
    sample=$(echo "$line" | awk '{print $1}')
    ORIG_IND+=("$sample")
done < "$IND_FILE"

IND_TEMP="${TMPDIR}/dataset.ind"
POPFILE="${TMPDIR}/poplist.txt"
PARFILE="${TMPDIR}/parfile.par"

# --- Создаём .ind ---
> "$IND_TEMP"
for sample in "${ORIG_IND[@]}"; do
    pop="ignore"
    if [[ " ${TARGET_SAMPLES[*]} " =~ " $sample " ]]; then
        pop="$TARGET_NAME"
    else
        for name in "${LEFT_NAMES[@]}"; do
            IFS=' ' read -ra tmp <<< "${grp_samples[$name]}"
            if [[ " ${tmp[*]} " =~ " $sample " ]]; then
                pop="$name"
                break
            fi
        done
        if [[ "$pop" == "ignore" ]]; then
            for name in "${RIGHT_NAMES[@]}"; do
                IFS=' ' read -ra tmp <<< "${grp_samples[$name]}"
                if [[ " ${tmp[*]} " =~ " $sample " ]]; then
                    pop="$name"
                    break
                fi
            done
        fi
    fi
    echo -e "${sample}\tU\t${pop}" >> "$IND_TEMP"
done

# --- Создаём poplist.txt ---
{
    echo "target: $TARGET_NAME"
    echo -n "left: "
    echo "${LEFT_NAMES[*]}" | tr ' ' ','
    echo -n "right: "
    echo "${RIGHT_NAMES[*]}" | tr ' ' ','
} > "$POPFILE"

# --- Создаём parfile ---
cat > "$PARFILE" <<EOF
genotypename: $GENO_FILE
snpname:      $SNP_FILE
indivname:    $IND_TEMP
popfilename:  $POPFILE
allsnps:      YES
inbreed:      NO
EOF

# --- Запуск qpAdm ---
qpAdm -p "$PARFILE" > "$OUTFILE"

# Optional: краткий вывод
grep -E "estimated proportions|std errors" "$OUTFILE"
