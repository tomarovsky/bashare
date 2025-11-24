#!/bin/bash

set -euo pipefail

SCATTER_COUNT=64
JAVA_MEM="4g"

# GATK
export PATH=$(conda info --base)/envs/gatk/bin/:${TOOLS}/gatk-4.6.2.0/:${PATH}

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 genome.fasta <BAM_LIST_FILE> <MALES_LIST_FILE> haploid.bed outprefix"
    echo "  <BAM_LIST_FILE>: Path to a file containing one BAM path per line."
    echo "  <MALES_LIST_FILE>: Path to a file containing one male sample ID path per line."
    exit 1
fi

REF=$1
BAM_LIST_FILE=$2
MALES_LIST_FILE=$3
HAPLOID_BED=$4
OUTPREFIX=$5

if [ ! -f "$REF" ]; then echo "[ERROR] Reference file not found: $REF"; exit 1; fi
if [ ! -f "$BAM_LIST_FILE" ]; then echo "[ERROR] BAM list file not found: $BAM_LIST_FILE"; exit 1; fi
if [ ! -f "$HAPLOID_BED" ]; then echo "[ERROR] Haploid BED file not found: $HAPLOID_BED"; exit 1; fi

mapfile -t BAM_ARRAY < <(grep -v '^[[:space:]]*$' "$BAM_LIST_FILE")
mapfile -t MALES_ARRAY < <(grep -v '^[[:space:]]*$' "$MALES_LIST_FILE")

if [[ ! -f "${REF%.*}.dict" ]]; then
    picard CreateSequenceDictionary -R "$REF"
fi

# ---- Directories ----
mkdir -p gatk_tmp/intervals gatk_tmp/chunks gatk_tmp/logs gatk_tmp/gvcfs

# ---- SplitIntervals ----
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

# ---- BAM processing ----
echo "[2/4] BAM processing..."
for BAM in "${BAM_ARRAY[@]}"; do
    NAME=$(basename "$BAM" | cut -d. -f1)
    FINAL_GVCF="gatk_tmp/gvcfs/${NAME}.g.vcf.gz"

    # Skip if sample GVCF already exists
    if [ -f "$FINAL_GVCF" ] && [ -f "${FINAL_GVCF}.tbi" ]; then
        echo "[INFO] ${NAME}: GVCF already exists."
        ALL_SAMPLE_GVCFS+=("-V" "$FINAL_GVCF")
        continue
    fi

    # Check sex
    IS_MALE=false
    for M in "${MALES_ARRAY[@]}"; do
        if [[ "$M" == "$NAME" ]]; then IS_MALE=true; break; fi
    done

    if [ "$IS_MALE" = true ]; then echo "[INFO] ${NAME}: MALE"; else echo "[INFO] ${NAME}: FEMALE"; fi

    CHUNK_FILES_ARGS=""

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
            echo "[INFO] ${NAME}.${ID}: Chunk already exists."
            continue
        fi

        (
            if [ "$IS_MALE" = true ]; then
                # === Logic for MALES ===
                P1="gatk_tmp/chunks/${NAME}.${ID}.p1.g.vcf.gz"
                P2="gatk_tmp/chunks/${NAME}.${ID}.p2.g.vcf.gz"

                # 1. Ploidy 1
                # Check if INTERVAL is fully within HAPLOID_BED
                HAS_HAPLOID=$(bedtools intersect -a <(cat "$INTERVAL" | grep -v "^@") -b "$HAPLOID_BED")

                if [ -z "$HAS_HAPLOID" ]; then
                    echo "[INFO] ${NAME}.${ID}: INTERVAL fully without HAPLOID, skipping haploid varcall." | tee -a "$LOG"
                else
                    if [ ! -f "$P1" ] || [ ! -f "${P1}.tbi" ]; then
                        gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
                            -R "$REF" -I "$BAM" -O "$P1" -ERC GVCF \
                            --sample-ploidy 1 \
                            -L "$INTERVAL" -L "$HAPLOID_BED" \
                            >> "$LOG" 2>&1
                    fi
                fi

                # 2. Ploidy 2
                # Check if INTERVAL is fully without HAPLOID_BED
                HAS_DIPLOID=$(bedtools subtract -a <(cat "$INTERVAL" | grep -v "^@") -b "$HAPLOID_BED")

                if [ -z "$HAS_DIPLOID" ]; then
                    echo "[INFO] ${NAME}.${ID}: INTERVAL fully within HAPLOID, skipping diploid varcall." | tee -a "$LOG"
                else
                    if [ ! -f "$P2" ] || [ ! -f "${P2}.tbi" ]; then
                        gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
                            -R "$REF" -I "$BAM" -O "$P2" -ERC GVCF \
                            --sample-ploidy 2 \
                            -L "$INTERVAL" -XL "$HAPLOID_BED" \
                            >> "$LOG" 2>&1
                    fi
                fi

                # 3. Combine p1 and p2
                if [ ! -f "$CHUNK_OUT" ] || [ ! -f "${CHUNK_OUT}.tbi" ]; then
                    if [ -f "$P1" ] && [ -f "$P2" ]; then
                        gatk --java-options "-Xmx4g" CombineGVCFs \
                            -R "$REF" \
                            -V "$P1" -V "$P2" \
                            -O "$CHUNK_OUT" \
                            >> "$LOG" 2>&1
                    elif [ -f "$P1" ]; then
                        cp "$P1" "$CHUNK_OUT"
                        cp "${P1}.tbi" "${CHUNK_OUT}.tbi"
                    elif [ -f "$P2" ]; then
                        cp "$P2" "$CHUNK_OUT"
                        cp "${P2}.tbi" "${CHUNK_OUT}.tbi"
                    else
                        echo "[ERROR] No P1 or P2 produced for ${NAME}.${ID}" | tee -a "$LOG"
                        exit 1
                    fi
                fi

                # Cleanup
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
                        >> "$LOG" 2>&1
                fi
            fi
        ) &
    done

    wait

    echo "[INFO] ${NAME}: merging chunks ..."
    if [ ! -f "$FINAL_GVCF" ] || [ ! -f "${FINAL_GVCF}.tbi" ]; then
        gatk --java-options "-Xmx16g" MergeVcfs \
            $CHUNK_FILES_ARGS \
            -O "$FINAL_GVCF" \
            >> "$LOG" 2>&1
        echo "[INFO] ${NAME}: chunks merged into final GVCF."
    else
        echo "[INFO] ${NAME}: final GVCF already exists."
    fi

    # Cleanup chunks
    if [ -f "$FINAL_GVCF" ] && [ -f "${FINAL_GVCF}.tbi" ]; then
        rm -f gatk_tmp/chunks/${NAME}.*.g.vcf.gz*
        echo "[INFO] ${NAME}: chunk files cleaned up."
    fi

    ALL_SAMPLE_GVCFS+=("-V" "$FINAL_GVCF")
    echo "[INFO] $NAME Done."
done

echo "[3/4] CombineGVCFs (Merging all samples)"
COMBINED_GVCF="${OUTPREFIX}.combined.g.vcf.gz"

if [ ! -f "$COMBINED_GVCF" ] || [ ! -f "${COMBINED_GVCF}.tbi" ]; then
    echo "[INFO] Merging sample GVCFs into one..."

    gatk --java-options "-Xmx64g" CombineGVCFs \
        -R "$REF" \
        "${ALL_SAMPLE_GVCFS[@]}" \
        -O "$COMBINED_GVCF"

    echo "[INFO] CombineGVCFs completed."
else
    echo "[INFO] Combined GVCF already exists: $COMBINED_GVCF"
fi

echo "[4/4] GenotypeGVCFs"
FINAL_VCF="${OUTPREFIX}.vcf.gz"

if [ ! -f "$FINAL_VCF" ] || [ ! -f "${FINAL_VCF}.tbi" ]; then
    echo "[INFO] Genotyping..."

    gatk --java-options "-Xmx64g" GenotypeGVCFs \
        -R "$REF" \
        -V "$COMBINED_GVCF" \
        -O "$FINAL_VCF" \
        -G StandardAnnotation

    echo "[INFO] GenotypeGVCFs completed."
else
    echo "[INFO] Final VCF already exists."
fi

echo "[INFO] Done! Result: $FINAL_VCF"
