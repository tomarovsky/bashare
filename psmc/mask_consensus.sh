#!/bin/bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input.fastq.gz> <mask.bed> <output.fastq.gz>"
    exit 1
fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate py38

CONSENSUS_FASTQ="$1"
MASK_BED="$2"
OUTPUT_FASTQ_GZ="$3"

# temporary files
FASTA_TMP=${CONSENSUS_FASTQ%.*}.tmp.fasta
QUAL_TMP=${CONSENSUS_FASTQ%.*}.tmp.qual
FASTA_MASKED_TMP=${CONSENSUS_FASTQ%.*}.tmp.masked.fasta
QUAL_MASKED_TMP=${CONSENSUS_FASTQ%.*}.tmp.masked.qual

# cleanup function
# cleanup() {
#     rm -f "$FASTA_TMP" "$QUAL_TMP" "$FASTA_MASKED_TMP" "$QUAL_MASKED_TMP"
# }
# trap cleanup EXIT # Execute cleanup on exit


echo "$(date) | Split $CONSENSUS_FASTQ into FASTA and QUAL..."

# We pass the Python script and file names as arguments
python3 - "$CONSENSUS_FASTQ" "$FASTA_TMP" "$QUAL_TMP" <<'EOF'
import sys
import gzip

def split_consensus_fastq(fastq_gz_path, fasta_path, qual_path):
    state = None # 'seq' or 'qual'
    current_header = None
    current_seq = []
    current_qual = []

    with gzip.open(fastq_gz_path, 'rt') as f_in, \
         open(fasta_path, 'w') as f_out, \
         open(qual_path, 'w') as q_out:

        def write_record():
            if current_header:
                f_out.write(f">{current_header}\n")
                f_out.write("".join(current_seq) + "\n")
                q_out.write(f">{current_header}\n")
                q_out.write("".join(current_qual) + "\n")

        for line in f_in:
            line = line.strip()
            if not line:
                continue

            if line.startswith('@'):
                write_record()
                current_header = line[1:]
                current_seq = []
                current_qual = []
                state = 'seq'

            elif line == '+' and state == 'seq':
                state = 'qual'

            elif state == 'seq':
                current_seq.append(line)

            elif state == 'qual':
                current_qual.append(line)

        write_record()

if __name__ == "__main__":
    # sys.argv[0] - '-' (script name)
    in_fastq = sys.argv[1]
    out_fasta = sys.argv[2]
    out_qual = sys.argv[3]
    split_consensus_fastq(in_fastq, out_fasta, out_qual)
EOF


conda deactivate && conda activate varcall

echo "$(date) | Mask sequences (replace with 'n') and quality (replace with '!')..."
bedtools maskfasta -fi "$FASTA_TMP" -bed "$MASK_BED" -fo "$FASTA_MASKED_TMP" -mc 'n' &
bedtools maskfasta -fi "$QUAL_TMP" -bed "$MASK_BED" -fo "$QUAL_MASKED_TMP" -mc '!' &

wait
conda deactivate && conda activate py38

echo "$(date) | Merge masked FASTA and QUAL files into ${OUTPUT_FASTQ_GZ}..."
python3 - "$FASTA_MASKED_TMP" "$QUAL_MASKED_TMP" <<'EOF'
import sys
import textwrap
import gzip

def fasta_parser(filename):
    header = None
    seq_parts = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if header:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header:
            yield header, "".join(seq_parts)

def merge_to_fastq(fasta_file, qual_file):
    seq_gen = fasta_parser(fasta_file)
    qual_gen = fasta_parser(qual_file)

    wrapper = textwrap.TextWrapper(width=60, break_long_words=True, replace_whitespace=False)

    for (seq_header, seq), (qual_header, qual) in zip(seq_gen, qual_gen):
        if seq_header != qual_header:
            print(f"Critical error: Mismatched headers: {seq_header} != {qual_header}", file=sys.stderr)
            sys.exit(1)

        # Print to stdout
        print(f"@{seq_header}")

        # Print with line wrapping
        if seq:
            print("\n".join(wrapper.wrap(seq)))
        else:
            print() # For empty sequence

        print("+")

        if qual:
            print("\n".join(wrapper.wrap(qual)))
        else:
            print()

if __name__ == "__main__":
    fasta_in = sys.argv[1]
    qual_in = sys.argv[2]
    merge_to_fastq(fasta_in, qual_in)
EOF
 | gzip -c > "$OUTPUT_FASTQ_GZ"

echo "$(date) | DONE | ${CONSENSUS_FASTQ} -> ${OUTPUT_FASTQ_GZ}"
