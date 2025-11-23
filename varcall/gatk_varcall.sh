#!/bin/bash

set -euo pipefail

SCATTER_COUNT=64
JAVA_MEM="4g"

# GATK
export PATH=$(conda info --base)/envs/gatk/bin/:${TOOLS}/gatk-4.6.2.0/:${PATH}

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 genome.fasta <BAM_LIST_FILE> <comma-separated male IDs> haploid.bed outprefix"
    echo "  <BAM_LIST_FILE>: Path to a text file containing one BAM path per line."
    exit 1
fi

REF=$1
BAM_LIST_FILE=$2
MALES_LIST_RAW=$3
HAPLOID_BED=$4
OUTPREFIX=$5

if [ ! -f "$REF" ]; then
    echo "[ERROR] Reference file not found: $REF"
    exit 1
fi

if [ ! -f "$BAM_LIST_FILE" ]; then
    echo "[ERROR] BAM list file not found: $BAM_LIST_FILE"
    exit 1
fi

if [ ! -f "$HAPLOID_BED" ]; then
    echo "[ERROR] Haploid BED file not found: $HAPLOID_BED"
    exit 1
fi

mapfile -t BAM_ARRAY < <(grep -v '^[[:space:]]*$' "$BAM_LIST_FILE")
IFS=',' read -r -a MALES_ARRAY <<< "$MALES_LIST_RAW"

if [[ ! -f "${REF%.*}.dict" ]]; then
    picard CreateSequenceDictionary -R "$REF"
fi

# output directories
mkdir -p gatk_tmp/intervals
mkdir -p gatk_tmp/chunks
mkdir -p gatk_tmp/logs
mkdir -p gatk_tmp/gvcfs

echo "[1/4] SplitIntervals"
if [ -z "$(ls -A gatk_tmp/intervals)" ]; then
    gatk --java-options "-Xmx8g" SplitIntervals \
        -R "$REF" \
        --scatter-count "$SCATTER_COUNT" \
        -O gatk_tmp/intervals
else
    echo "[INFO] Intervals already exist, using them."
fi

INTERVAL_FILES=(gatk_tmp/intervals/*.interval_list)

ALL_SAMPLE_GVCFS=()

echo "[2/4] BAM processing..."
for BAM in "${BAM_ARRAY[@]}"; do
    # Sample ID: path/to/SAMPLE.alignment.bam -> SAMPLE
    NAME=$(basename "$BAM" | cut -d. -f1)
    FINAL_GVCF="gatk_tmp/gvcfs/${NAME}.g.vcf.gz"

    # Check sex
    IS_MALE=false
    for M in "${MALES_ARRAY[@]}"; do
        if [[ "$M" == "$NAME" ]]; then IS_MALE=true; break; fi
    done

    if [ "$IS_MALE" = true ]; then
        echo "[INFO] ${NAME}: MALE"
    else
        echo "[INFO] ${NAME}: FEMALE"
    fi

    CHUNK_FILES_ARGS=""
    for i in "${!INTERVAL_FILES[@]}"; do
        INTERVAL="${INTERVAL_FILES[$i]}"
        ID=$(printf "%04d" $i) # 0001, 0002...

        # Files for this chunk
        CHUNK_OUT="gatk_tmp/chunks/${NAME}.${ID}.g.vcf.gz"
        LOG="gatk_tmp/logs/${NAME}.${ID}.log"

        # Collect arguments for future merging (-I file1 -I file2...)
        CHUNK_FILES_ARGS="${CHUNK_FILES_ARGS} -I $CHUNK_OUT"

        (
            # Skip if chunk is already ready
            if [ -f "$CHUNK_OUT" ] && [ -f "${CHUNK_OUT}.tbi" ]; then
                exit 0
            fi

            if [ "$IS_MALE" = true ]; then
                # === Logic for MALES ===
                P1="gatk_tmp/chunks/${NAME}.${ID}.p1.g.vcf.gz"
                P2="gatk_tmp/chunks/${NAME}.${ID}.p2.g.vcf.gz"

                # 1. Ploidy 1 (intersection of interval with haploid_bed)
                gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
                    -R "$REF" -I "$BAM" -O "$P1" -ERC GVCF \
                    --sample-ploidy 1 \
                    -L "$INTERVAL" -L "$HAPLOID_BED" \
                    |& tee -a "$LOG" || true

                # 2. Ploidy 2 (interval WITHOUT haploid_bed)
                gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
                    -R "$REF" -I "$BAM" -O "$P2" -ERC GVCF \
                    --sample-ploidy 2 \
                    -L "$INTERVAL" -XL "$HAPLOID_BED" \
                    |& tee -a "$LOG" || true

                # 3. Merge p1 and p2
                # Check if files exist (GATK does not create a file if the region is empty)
                MERGE_LIST=""
                if [ -f "$P1" ]; then MERGE_LIST="$MERGE_LIST -I $P1"; fi
                if [ -f "$P2" ]; then MERGE_LIST="$MERGE_LIST -I $P2"; fi

                gatk --java-options "-Xmx2g" MergeVcfs $MERGE_LIST -O "$CHUNK_OUT" > "$LOG" 2>&1

                # Remove p1 and p2
                rm -f "$P1"* "$P2"*

            else
                # === Logic for FEMALES ===
                gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
                    -R "$REF" -I "$BAM" -O "$CHUNK_OUT" -ERC GVCF \
                    --sample-ploidy 2 \
                    -L "$INTERVAL" \
                    |& tee -a "$LOG"
            fi
        ) &
    done

    echo "[INFO] ${NAME}: waiting for all chunks ..."
    wait

    echo "[INFO] ${NAME}: merging chunks ..."
    gatk --java-options "-Xmx16g" MergeVcfs \
        $CHUNK_FILES_ARGS \
        -O "$FINAL_GVCF"

    # Remove chunk files
    rm gatk_tmp/chunks/${NAME}.*.g.vcf.gz*

    ALL_SAMPLE_GVCFS+=("-V $FINAL_GVCF")
    echo "[INFO] $NAME Done."
done

echo "[3/4] CombineGVCFs"
gatk --java-options "-Xmx64g" CombineGVCFs \
    -R "$REF" \
    "${ALL_SAMPLE_GVCFS[@]}" \
    -O ${OUTPREFIX}.g.vcf.gz

echo "[4/4] GenotypeGVCFs"
gatk --java-options "-Xmx64g" GenotypeGVCFs \
    -R "$REF" \
    -V ${OUTPREFIX}.g.vcf.gz \
    -O ${OUTPREFIX}.vcf.gz

echo "[INFO] Done! Result: ${OUTPREFIX}.vcf.gz"
