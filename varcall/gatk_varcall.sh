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

# Check dependencies
command -v bedtools >/dev/null 2>&1 || { echo >&2 "[ERROR] bedtools not found in PATH."; exit 1; }
command -v gatk >/dev/null 2>&1 || { echo >&2 "[ERROR] gatk not found in PATH."; exit 1; }

if [ ! -f "$REF" ]; then echo "[ERROR] Reference file not found: $REF"; exit 1; fi
if [ ! -f "$BAM_LIST_FILE" ]; then echo "[ERROR] BAM list file not found: $BAM_LIST_FILE"; exit 1; fi
if [ ! -f "$HAPLOID_BED" ]; then echo "[ERROR] Haploid BED file not found: $HAPLOID_BED"; exit 1; fi

mapfile -t BAM_ARRAY < <(grep -v '^[[:space:]]*$' "$BAM_LIST_FILE")
mapfile -t MALES_ARRAY < <(grep -v '^[[:space:]]*$' "$MALES_LIST_FILE")

if [[ ! -f "${REF%.*}.dict" ]]; then
    picard CreateSequenceDictionary -R "$REF"
fi

# ---- Directories ----
mkdir -p gatk_tmp/intervals gatk_tmp/chunks gatk_tmp/joint_chunks gatk_tmp/gvcfs gatk_tmp/logs


# ---- SplitIntervals ----
echo "[1/4] | $(date) | SplitIntervals"
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

# Pre-convert Intervals to BED
echo "[INFO] Pre-converting intervals to BED..."
for INTERVAL in "${INTERVAL_FILES[@]}"; do
    if [ ! -f "${INTERVAL}.bed" ]; then
        # Picard (1-based) to BED (0-based start)
        awk -v OFS='\t' '!/^@/ {print $1, $2-1, $3}' "$INTERVAL" > "${INTERVAL}.bed"
    fi
done


# ---- BAM processing ----
echo "[2/4] | $(date) | BAM processing..."
for BAM in "${BAM_ARRAY[@]}"; do
    NAME=$(basename "$BAM" | cut -d. -f1)
    FINAL_GVCF="gatk_tmp/gvcfs/${NAME}.g.vcf.gz"
    SAMPLE_LOG="gatk_tmp/logs/${NAME}.main.log"

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
        INTERVAL_BED="${INTERVAL}.bed"
        ID=$(printf "%04d" $i) # 0001, 0002...

        # Files for this chunk
        CHUNK_OUT="gatk_tmp/chunks/${NAME}.${ID}.g.vcf.gz"
        CHUNK_LOG="gatk_tmp/logs/${NAME}.${ID}.log"

        # Collect arguments for future merging
        CHUNK_FILES_ARGS="${CHUNK_FILES_ARGS} -I $CHUNK_OUT"

        if [ -f "$CHUNK_OUT" ] && [ -f "${CHUNK_OUT}.tbi" ]; then
            echo "[INFO] ${NAME}.${ID}: Chunk already exists."
            continue
        fi

        (
            if [ "$IS_MALE" = true ]; then
                # === Logic for MALES ===
                P1="gatk_tmp/chunks/${NAME}.${ID}.p1.g.vcf.gz"
                P2="gatk_tmp/chunks/${NAME}.${ID}.p2.g.vcf.gz"

                # Temp interval files for this chunk
                P1_BED="gatk_tmp/intervals/${NAME}.${ID}.p1.bed"
                P2_BED="gatk_tmp/intervals/${NAME}.${ID}.p2.bed"

                # 1. Ploidy 1 (Intersection with Haploid regions)
                bedtools intersect -a "$INTERVAL_BED" -b "$HAPLOID_BED" > "$P1_BED"
                if [ -s "$P1_BED" ]; then
                    if [ ! -f "$P1" ] || [ ! -f "${P1}.tbi" ]; then
                        echo "${NAME}.${ID} | Ploidy 1 processing..."
                        gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
                            -R "$REF" -I "$BAM" -O "$P1" -ERC GVCF \
                            --sample-ploidy 1 -L "$P1_BED" >> "$CHUNK_LOG" 2>&1
                    fi
                fi

                # 2. Ploidy 2 (Subtraction of Haploid regions = Diploid/PAR)
                bedtools subtract -a "$INTERVAL_BED" -b "$HAPLOID_BED" > "$P2_BED"
                if [ -s "$P2_BED" ]; then
                    if [ ! -f "$P2" ] || [ ! -f "${P2}.tbi" ]; then
                        gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
                            -R "$REF" -I "$BAM" -O "$P2" -ERC GVCF \
                            --sample-ploidy 2 -L "$P2_BED" >> "$CHUNK_LOG" 2>&1
                    fi
                fi

                # 3. Combine p1 and p2
                if [ ! -f "$CHUNK_OUT" ] || [ ! -f "${CHUNK_OUT}.tbi" ]; then
                    if [ -f "$P1" ] && [ -f "$P2" ]; then
                        gatk --java-options "-Xmx4g" CombineGVCFs \
                            -R "$REF" -V "$P1" -V "$P2" -O "$CHUNK_OUT" >> "$CHUNK_LOG" 2>&1
                    elif [ -f "$P1" ]; then
                        cp "$P1" "$CHUNK_OUT" && cp "${P1}.tbi" "${CHUNK_OUT}.tbi"
                    elif [ -f "$P2" ]; then
                        cp "$P2" "$CHUNK_OUT" && cp "${P2}.tbi" "${CHUNK_OUT}.tbi"
                    else
                        echo "[ERROR] No P1 or P2 produced for ${NAME}.${ID}" >> "$CHUNK_LOG"
                        exit 1
                    fi
                fi

                # Cleanup temps
                if [ -f "$CHUNK_OUT" ] && [ -f "${CHUNK_OUT}.tbi" ]; then
                    rm -f "$P1_BED" "$P2_BED" "$P1" "$P1.tbi" "$P2" "$P2.tbi"
                fi

            else
                # === Logic for FEMALES ===
                if [ ! -f "$CHUNK_OUT" ] || [ ! -f "${CHUNK_OUT}.tbi" ]; then
                    gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
                        -R "$REF" -I "$BAM" -O "$CHUNK_OUT" -ERC GVCF \
                        --sample-ploidy 2 -L "$INTERVAL" >> "$CHUNK_LOG" 2>&1
                fi
            fi
        ) &
    done

    wait

    echo "[INFO] ${NAME}: gathering chunks..."
    if [ ! -f "$FINAL_GVCF" ] || [ ! -f "${FINAL_GVCF}.tbi" ]; then
        gatk --java-options "-Xmx16g" GatherVcfs \
            $CHUNK_FILES_ARGS \
            -O "$FINAL_GVCF" \
            >> "$SAMPLE_LOG" 2>&1
        gatk IndexFeatureFile -I "$FINAL_GVCF" >> "$SAMPLE_LOG" 2>&1
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


echo "[3/4] | $(date) | Combine and Genotype..."
GATHER_ARGS=""
for i in "${!INTERVAL_FILES[@]}"; do
    INTERVAL="${INTERVAL_FILES[$i]}"
    ID=$(printf "%04d" $i) # 0001, 0002...

    CHUNK_COMBINED="gatk_tmp/joint_chunks/${OUTPREFIX}.${ID}.combined.g.vcf.gz"
    CHUNK_FINAL="gatk_tmp/joint_chunks/${OUTPREFIX}.${ID}.vcf.gz"
    CHUNK_LOG="gatk_tmp/logs/joint_call.${ID}.log"

    # Save name of file for final GatherVcfs
    GATHER_ARGS="${GATHER_ARGS} -I ${CHUNK_FINAL}"

    (
        # 1. CombineGVCFs (only for this interval)
        if [ ! -f "$CHUNK_COMBINED" ] || [ ! -f "${CHUNK_COMBINED}.tbi" ]; then
            gatk --java-options "-Xmx8g" CombineGVCFs \
                -R "$REF" \
                "${ALL_SAMPLE_GVCFS[@]}" \
                -L "$INTERVAL" \
                -O "$CHUNK_COMBINED" >> "$CHUNK_LOG" 2>&1
        fi

        # 2. GenotypeGVCFs (also only for the interval)
        if [ ! -f "$CHUNK_FINAL" ] || [ ! -f "${CHUNK_FINAL}.tbi" ]; then
            gatk --java-options "-Xmx8g" GenotypeGVCFs \
                -R "$REF" \
                -V "$CHUNK_COMBINED" \
                -L "$INTERVAL" \
                -O "$CHUNK_FINAL" \
                -G StandardAnnotation >> "$CHUNK_LOG" 2>&1
        fi

        # rm -f "$CHUNK_COMBINED" "${CHUNK_COMBINED}.tbi"

    ) &

done

wait

echo "[4/4] | $(date) | Merging all chunks into final VCF..."
FINAL_VCF="${OUTPREFIX}.vcf.gz"

if [ ! -f "$FINAL_VCF" ] || [ ! -f "${FINAL_VCF}.tbi" ]; then
    gatk --java-options "-Xmx8g" GatherVcfs \
        $GATHER_ARGS \
        -O "$FINAL_VCF"

    gatk IndexFeatureFile -I "$FINAL_VCF"
else
    echo "[INFO] Final VCF already exists."
fi

echo "[INFO] Done! Result: $FINAL_VCF"
