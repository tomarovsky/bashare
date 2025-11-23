#!/bin/bash

set -euo pipefail

SCATTER_COUNT=64
JAVA_MEM="4g"

# GATK
export PATH=$(conda info --base)/envs/gatk/bin/:${TOOLS}/gatk-4.6.2.0/:${PATH}

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 genome.fasta <BAM_LIST_FILE> <MALES_LIST_FILE> haploid.bed outprefix"
    echo "  <BAM_LIST_FILE>: Path to a text file containing one BAM path per line."
    echo "  <MALES_LIST_FILE>: Path to a text file containing one male sample ID path per line."
    exit 1
fi

REF=$1
BAM_LIST_FILE=$2
MALES_LIST_FILE=$3
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
mapfile -t MALES_ARRAY < <(grep -v '^[[:space:]]*$' "$MALES_LIST_FILE")

if [[ ! -f "${REF%.*}.dict" ]]; then
    picard CreateSequenceDictionary -R "$REF"
fi

# output directories
mkdir -p gatk_tmp/intervals
mkdir -p gatk_tmp/chunks
mkdir -p gatk_tmp/logs
mkdir -p gatk_tmp/gvcfs

echo "[1/4] SplitIntervals"
if [ -z "$(ls -A gatk_tmp/intervals 2>/dev/null)" ]; then
    gatk --java-options "-Xmx8g" SplitIntervals \
        -R "$REF" \
        --scatter-count "$SCATTER_COUNT" \
        -O gatk_tmp/intervals
    echo "[INFO] Intervals created."
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

    # Skip if sample GVCF already exists
    if [ -f "$FINAL_GVCF" ] && [ -f "${FINAL_GVCF}.tbi" ]; then
        echo "[INFO] ${NAME}: GVCF already exists, skipping."
        ALL_SAMPLE_GVCFS+=("-V" "$FINAL_GVCF")
        continue
    fi

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
    CHUNKS_TO_PROCESS=0
    for i in "${!INTERVAL_FILES[@]}"; do
        INTERVAL="${INTERVAL_FILES[$i]}"
        ID=$(printf "%04d" $i) # 0001, 0002...

        # Files for this chunk
        CHUNK_OUT="gatk_tmp/chunks/${NAME}.${ID}.g.vcf.gz"
        LOG="gatk_tmp/logs/${NAME}.${ID}.log"

        # Collect arguments for future merging (-I file1 -I file2...)
        CHUNK_FILES_ARGS="${CHUNK_FILES_ARGS} -I $CHUNK_OUT"

        # Check if chunk already exists
        if [ -f "$CHUNK_OUT" ] && [ -f "${CHUNK_OUT}.tbi" ]; then
            echo "[INFO] ${NAME}.${ID}: Chunk already exists, skipping."
            continue
        fi

        ((CHUNKS_TO_PROCESS++))
        (
            if [ "$IS_MALE" = true ]; then
                # === Logic for MALES ===
                P1="gatk_tmp/chunks/${NAME}.${ID}.p1.g.vcf.gz"
                P2="gatk_tmp/chunks/${NAME}.${ID}.p2.g.vcf.gz"

                # 1. Ploidy 1 (intersection of interval with haploid_bed)
                if [ ! -f "$P1" ] || [ ! -f "${P1}.tbi" ]; then
                    gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
                        -R "$REF" -I "$BAM" -O "$P1" -ERC GVCF \
                        --sample-ploidy 1 \
                        -L "$INTERVAL" -L "$HAPLOID_BED" \
                        |& tee -a "$LOG" || true
                else
                    echo "[INFO] ${NAME}.${ID}.p1: Already exists, skipping."
                fi

                # 2. Ploidy 2 (interval WITHOUT haploid_bed)
                if [ ! -f "$P2" ] || [ ! -f "${P2}.tbi" ]; then
                    gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
                        -R "$REF" -I "$BAM" -O "$P2" -ERC GVCF \
                        --sample-ploidy 2 \
                        -L "$INTERVAL" -XL "$HAPLOID_BED" \
                        |& tee -a "$LOG" || true
                else
                    echo "[INFO] ${NAME}.${ID}.p2: Already exists, skipping."
                fi

                # 3. Merge p1 and p2
                # Check if files exist (GATK does not create a file if the region is empty)
                MERGE_LIST=""
                if [ -f "$P1" ]; then MERGE_LIST="$MERGE_LIST -I $P1"; fi
                if [ -f "$P2" ]; then MERGE_LIST="$MERGE_LIST -I $P2"; fi

                if [ -n "$MERGE_LIST" ] && ([ ! -f "$CHUNK_OUT" ] || [ ! -f "${CHUNK_OUT}.tbi" ]); then
                    gatk --java-options "-Xmx2g" MergeVcfs $MERGE_LIST -O "$CHUNK_OUT" >> "$LOG" 2>&1
                elif [ -f "$CHUNK_OUT" ] && [ -f "${CHUNK_OUT}.tbi" ]; then
                    echo "[INFO] ${NAME}.${ID}: Merged chunk already exists, skipping merge."
                fi

                # Remove p1 and p2 (only if we created them and chunk merge was successful)
                if [ -f "$CHUNK_OUT" ] && [ -f "${CHUNK_OUT}.tbi" ]; then
                    rm -f "$P1"* "$P2"*
                fi

            else
                # === Logic for FEMALES ===
                if [ ! -f "$CHUNK_OUT" ] || [ ! -f "${CHUNK_OUT}.tbi" ]; then
                    gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
                        -R "$REF" -I "$BAM" -O "$CHUNK_OUT" -ERC GVCF \
                        --sample-ploidy 2 \
                        -L "$INTERVAL" \
                        |& tee -a "$LOG"
                else
                    echo "[INFO] ${NAME}.${ID}: Chunk already exists, skipping."
                fi
            fi
        ) &
    done

    if [ $CHUNKS_TO_PROCESS -gt 0 ]; then
        echo "[INFO] ${NAME}: waiting for $CHUNKS_TO_PROCESS chunks to process..."
        wait
        echo "[INFO] ${NAME}: all chunks processed."
    else
        echo "[INFO] ${NAME}: all chunks already exist."
    fi

    echo "[INFO] ${NAME}: merging chunks ..."
    if [ ! -f "$FINAL_GVCF" ] || [ ! -f "${FINAL_GVCF}.tbi" ]; then
        gatk --java-options "-Xmx16g" MergeVcfs \
            $CHUNK_FILES_ARGS \
            -O "$FINAL_GVCF"
        echo "[INFO] ${NAME}: chunks merged into final GVCF."
    else
        echo "[INFO] ${NAME}: final GVCF already exists, skipping merge."
    fi

    # Remove chunk files only if final GVCF was successfully created
    if [ -f "$FINAL_GVCF" ] && [ -f "${FINAL_GVCF}.tbi" ]; then
        rm -f gatk_tmp/chunks/${NAME}.*.g.vcf.gz*
        echo "[INFO] ${NAME}: chunk files cleaned up."
    fi

    ALL_SAMPLE_GVCFS+=("-V $FINAL_GVCF")
    echo "[INFO] $NAME Done."
done

# All intervals (needed for GenomicsDBImport and GenotypeGVCFs)
ALL_INTERVALS=()
for INTERVAL_FILE in gatk_tmp/intervals/*.interval_list; do
    ALL_INTERVALS+=("-L" "$INTERVAL_FILE")
done

echo "[3/4] GenomicsDBImport"
DB_WORKSPACE="gatk_tmp/genomics_db"
# GenomicsDBImport do not work if outdir already exists
if [ ! -d "$DB_WORKSPACE" ]; then
    echo "[INFO] Creating new GenomicsDB..."
    # use --genomicsdb-update-workspace-path "$DB_WORKSPACE" \ to update existing db
    gatk --java-options "-Xmx64g" GenomicsDBImport \
        "${ALL_SAMPLE_GVCFS[@]}" \
        "${ALL_INTERVALS[@]}" \
        --genomicsdb-workspace-path "$DB_WORKSPACE" \
        --tmp-dir gatk_tmp \
        --merge-input-intervals \
        --reader-threads 1 \
        --batch-size 50
    echo "[INFO] GenomicsDBImport completed."

else
    echo "[INFO] GenomicsDB directory ($DB_WORKSPACE) already exists, skipping import."
fi

echo "[4/4] GenotypeGVCFs"
FINAL_VCF="${OUTPREFIX}.vcf.gz"
if [ ! -f "$FINAL_VCF" ] || [ ! -f "${FINAL_VCF}.tbi" ]; then
    gatk --java-options "-Xmx64g" GenotypeGVCFs \
        -R "$REF" \
        -V gendb://"$DB_WORKSPACE" \
        -O "$FINAL_VCF" \
        -G StandardAnnotation \
        "${INTERVAL_ARGS[@]}" \
        --only-output-calls-starting-in-intervals \
        --merge-input-intervals

    echo "[INFO] GenotypeGVCFs completed."
else
    echo "[INFO] Final VCF already exists, skipping."
fi

echo "[INFO] Done! Result: $FINAL_VCF"
