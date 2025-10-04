#!/bin/bash
# common_mask.sh BED_FILES OUTPREFIX THREADS

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 BED_FILE1 [BED_FILE2 ...] OUTPREFIX THREADS"
    exit 1
fi

ARGS=("$@")
OUTPREFIX="${ARGS[-2]}"
THREADS="${ARGS[-1]}"
BED_FILES=("${ARGS[@]:0:${#ARGS[@]}-2}")

# sort all input BED files
printf "%s\n" "${BED_FILES[@]}" | parallel -j "$THREADS" '
    echo "Sorting {}"
    sort -k1,1 -k2,2n {} > {}.sorted.bed
'

# array of sorted BED files
sorted_bed_files=()
for bed in "${BED_FILES[@]}"; do
    sorted_bed_files+=("${bed}.sorted.bed")
done

tmpdir="./tmp_intersect"
mkdir -p "$tmpdir"

# all pairwise combinations and compute intersections in parallel
for i in "${!sorted_bed_files[@]}"; do
    for j in "${sorted_bed_files[@]:$((i+1))}"; do
        echo "${sorted_bed_files[i]} $j"
    done
done | parallel -j "$THREADS" --colsep ' ' '
    a={1}
    b={2}
    a_sample=$(basename "${a}" | cut -d"." -f1)
    b_sample=$(basename "${b}" | cut -d"." -f1)
    tmpfile="'"$tmpdir"'/${a_sample}_${b_sample}.intersect"
    echo "Processing ${a_sample} vs ${b_sample} -> $tmpfile"
    bedtools intersect \
        -a "${a}" \
        -b "${b}" \
        -sorted > "$tmpfile"
'

echo "Sorting..."
sort -k1,1 -k2,2n "$tmpdir"/*.intersect > "$tmpdir/all_intersect.sorted.bed"

echo "Merging..."
bedtools merge -i "$tmpdir/all_intersect.sorted.bed" > "${OUTPREFIX}.merge_all.intersect_2.mapq10.bed"

# cleanup temporary files
rm -r "$tmpdir"

