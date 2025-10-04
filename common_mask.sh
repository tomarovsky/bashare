#!/bin/bash
# common_mask.sh BED_FILES OUTPREFIX THREADS

BED_FILES=$1
OUTPREFIX=$2
THREADS=$3

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 BED_FILES OUTPREFIX THREADS"
    exit 1
fi

# Сортируем файлы через parallel
echo "$BED_FILES" | tr ' ' '\n' | parallel -j "$THREADS" '
    echo 'Sorting {}'
    sort -k1,1 -k2,2n {} > {}.sorted.bed
'

# Получаем список отсортированных файлов
sorted_bed_files=()
for bed in $BED_FILES; do
    sorted_bed_files+=("${bed}.sorted.bed")
done

tmpdir="./tmp_intersect"
mkdir -p "$tmpdir"

for i in "${!sorted_bed_files[@]}"; do
    for j in "${sorted_bed_files[@]:$((i+1))}"; do
        echo "${sorted_bed_files[i]} $j"
    done
done | parallel -j "$THREADS" '
    a={1}
    b={2}
    tmpfile=${tmpdir}/${a%%.*}_${b%%.*}.intersect
    bedtools intersect \
        -a "${a}" \
        -b "${b}" \
        -sorted > "$tmpfile"
'

bedtools merge -i <(cat "$tmpdir"/*.intersect | sort -k1,1 -k2,2n) > ${OUTPREFIX}.merge_all.intersect_2.mapq10.bed

# rm -r "$tmpdir"


