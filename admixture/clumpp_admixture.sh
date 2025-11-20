#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 INDFILE Q1.Q [Q2.Q ...]"
    exit 1
fi

# parse arguments
INDFILE="$1" # sample per line (VCF order)
shift
QFILES=("$@") # Q-files from ADMIXTURE (identical K values)

R=${#QFILES[@]}                                 # number of runs
C=$(wc -l < "${QFILES[0]}")                     # number of individuals == number of lines in Q-file
K=$(head -1 "${QFILES[0]}" | awk '{print NF}')  # number of clusters == number of columns

# Output files
CLUMPP_INFILE="clumpp.indfile"
CLUMPP_OUTFILE="clumpp.output"
CLUMPP_MISCFILE="clumpp.log"
CLUMPP_PARAMFILE="clumpp.paramfile"
CLUMPP_PERMUTED_DATAFILE="clumpp.permutted.R"
CLUMPP_OUTFILE_TSV="clumpp.output.tsv"

# 1. Generate clumpp.indfile
> "$CLUMPP_INFILE"
line_nr=1
for ((r=0; r<R; r++)); do
    Q="${QFILES[$r]}"
    indnr=1
    while read -r line; do
        echo -e "${line_nr}\t${indnr}\t(1)\t${indnr} : $line" >> "$CLUMPP_INFILE"
        ((indnr++))
        ((line_nr++))
    done < "$Q"
done

# 2. Generate clumpp.paramfile
cat > "$CLUMPP_PARAMFILE" <<EOF
DATATYPE 0
INDFILE $CLUMPP_INFILE
OUTFILE $CLUMPP_OUTFILE
MISCFILE $CLUMPP_MISCFILE
K $K
C $C
R $R
M 1
W 0
S 2
GREEDY_OPTION 2
REPEATS 100
PRINT_PERMUTED_DATA 2
PERMUTED_DATAFILE $CLUMPP_PERMUTED_DATAFILE
PRINT_EVERY_PERM 0
PRINT_RANDOM_INPUTORDER 0
OVERRIDE_WARNINGS 0
ORDER_BY_RUN 1
EOF

# 3. Run CLUMPP
CLUMPP "$CLUMPP_PARAMFILE" > /dev/null
echo $K
echo "Done! K=$K"
grep -P "^The highest value" $CLUMPP_MISCFILE

# 4. CLUMPP output to TSV
python3 - <<EOF
import sys

with open('$INDFILE', 'r') as f:
    sample_ids = [line.strip() for line in f]

output_lines = []
with open('$CLUMPP_OUTFILE', 'r') as f:
    for line in f:
        clumpp_clusters = line.split(':')[1].strip()
        clumpp_clusters = '\\t'.join(clumpp_clusters.split())
        output_lines.append(clumpp_clusters)

if len(sample_ids) != len(output_lines):
    print(f"Warning: number of samples ({len(sample_ids)}) not equal to the number of lines in the output file ({len(output_lines)})", file=sys.stderr)

with open('$CLUMPP_OUTFILE_TSV', 'w') as f:
    for sample_id, q_values in zip(sample_ids, output_lines):
        f.write(f"{sample_id}\\t{q_values}\\n")
EOF
