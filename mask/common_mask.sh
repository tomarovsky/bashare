#!/bin/bash
# Usage: 200 GB RAM
# common_mask_universal.sh BED_FILES N_COMB OUTPREFIX THREADS
# Example: ./common_mask_universal.sh *.bed 2 result_prefix 64

set -euo pipefail

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 BED_FILE1 [BED_FILE2 ...] N_COMB OUTPREFIX THREADS"
    echo "  N_COMB: Number of files for intersection (e.g., 2, 3, 4...)"
    echo "  Tools: bedtools, python3"
    exit 1
fi

ARGS=("$@")
THREADS="${ARGS[-1]}"
OUTPREFIX="${ARGS[-2]}"
N_COMB="${ARGS[-3]}"
BED_FILES=("${ARGS[@]:0:${#ARGS[@]}-3}")

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

tmpdir="./tmp_intersect_${OUTPREFIX}"
mkdir -p "${tmpdir}/tmp_sort/"

# Export tmpdir for parallel
export tmpdir

# 2. Generate combinations and intersect
python3 -c "
import sys, itertools
n = int(sys.argv[1])
files = sys.argv[2:]
for combo in itertools.combinations(files, n):
    print(' '.join(combo))
" "$N_COMB" "${sorted_bed_files[@]}" | parallel -j "$THREADS" '
    read -ra files <<< "{}"
    outfile="$tmpdir/job_{#}.intersect"

    cmd="cat ${files[0]}"

    for ((i=1; i<${#files[@]}; i++)); do
        f="${files[$i]}"
        cmd="$cmd | bedtools intersect -a stdin -b $f -sorted"
    done

    echo "Processing combination #{#} (Size ${#files[@]}) -> $outfile"

    eval "$cmd" > "$outfile"
'

echo "Sorting merged results..."
sort -S 200G --parallel="$THREADS" -T "${tmpdir}/tmp_sort" --merge -k1,1 -k2,2n "$tmpdir"/*.intersect > "$tmpdir/all_intersect.sorted.bed"

echo "Merging final mask..."
bedtools merge -i "$tmpdir/all_intersect.sorted.bed" > "${OUTPREFIX}.merge_all.intersect_N${N_COMB}.bed"

echo "Total masked:"
awk '{sum += $3 - $2} END {print sum}' "${OUTPREFIX}.merge_all.intersect_N${N_COMB}.bed"

# cleanup temporary files
echo "Cleaning up intersection tmp dir..."
rm -r "$tmpdir"

echo "Cleaning up sorted input files..."
for f in "${sorted_bed_files[@]}"; do
    rm "$f"
done
