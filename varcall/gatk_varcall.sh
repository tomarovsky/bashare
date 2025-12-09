#!/bin/bash
set -euo pipefail

source "$TOOLS/bashare/lib/log_functions.sh"

SCATTER_COUNT=64
JAVA_MEM="32g"

# GATK path
export PATH=${TOOLS}/gatk-4.6.2.0/:${PATH}

if [ "$#" -ne 5 ]; then
    log_error "Usage: $0 genome.fasta <BAM_LIST_FILE> <MALES_LIST_FILE> haploid.bed outprefix"
    log_error "  <BAM_LIST_FILE>: Path to a file containing one BAM path per line."
    log_error "  <MALES_LIST_FILE>: Path to a file containing one male sample ID path per line."
    exit 1
fi

REF=$1
BAM_LIST_FILE=$2
MALES_LIST_FILE=$3
HAPLOID_BED=$4
OUTPREFIX=$5

# Check dependencies
command -v bedtools >/dev/null 2>&1 || { log_error "bedtools not found."; exit 1; }
command -v gatk >/dev/null 2>&1 || { log_error "gatk not found."; exit 1; }
command -v bcftools >/dev/null 2>&1 || { log_error "bcftools not found."; exit 1; }
command -v samtools >/dev/null 2>&1 || { log_error "samtools not found."; exit 1; }

if [ ! -f "$REF" ]; then log_error "Reference file not found: $REF"; exit 1; fi
if [ ! -f "$BAM_LIST_FILE" ]; then log_error "BAM list file not found: $BAM_LIST_FILE"; exit 1; fi
if [ ! -f "$HAPLOID_BED" ]; then log_error "Haploid BED file not found: $HAPLOID_BED"; exit 1; fi

if [[ ! -f "${REF%.*}.dict" ]]; then
    gatk CreateSequenceDictionary -R "$REF"
fi

# ---- 0. Prepare BAM file lists ----
log_info "Categorizing BAM files..."

# Males list into array
declare -A MALES_MAP
while read -r line; do
    [[ -z "$line" ]] && continue
    MALES_MAP["$line"]=1
done < "$MALES_LIST_FILE"

# Declare arrays
ALL_BAMS_ARGS=()
MALE_BAMS_ARGS=()
FEMALE_BAMS_ARGS=()

mapfile -t BAM_ARRAY < <(grep -v '^[[:space:]]*$' "$BAM_LIST_FILE")

for BAM in "${BAM_ARRAY[@]}"; do
    NAME=$(basename "$BAM" | cut -d. -f1)

    ALL_BAMS_ARGS+=("-I" "$BAM")

    # Check sex
    if [[ -n "${MALES_MAP[$NAME]-}" ]]; then
        MALE_BAMS_ARGS+=("-I" "$BAM")
    else
        FEMALE_BAMS_ARGS+=("-I" "$BAM")
    fi
done

log_info "Total BAMs: $((${#ALL_BAMS_ARGS[@]}/2))"
log_info "Males: $((${#MALE_BAMS_ARGS[@]}/2))"
log_info "Females: $((${#FEMALE_BAMS_ARGS[@]}/2))"

# ---- Extract sample names from BAMs ----
SAMPLE_ORDER=()
for BAM in "${BAM_ARRAY[@]}"; do
    SAMPLE=$(samtools view -H "$BAM" | grep "^@RG" | grep -o "SM:[^[:space:]]*" | cut -d: -f2 | head -n1)
    SAMPLE_ORDER+=("$SAMPLE")
done

# ---- Directories ----
mkdir -p gatk_tmp/intervals gatk_tmp/chunks gatk_tmp/logs

# ---- 1. SplitIntervals ----
log_info "[1/3] SplitIntervals"
if [ -z "$(ls -A gatk_tmp/intervals 2>/dev/null)" ]; then
    gatk --java-options "-Xmx4g" SplitIntervals \
        -R "$REF" \
        --scatter-count "$SCATTER_COUNT" \
        -O gatk_tmp/intervals
else
    log_info "Intervals already exist."
fi

INTERVAL_FILES=(gatk_tmp/intervals/*.interval_list)

# ---- 2. Processing Chunks ----
log_info "[2/3] Processing chunks with Multi-sample HaplotypeCaller..."

CHUNK_GATHER_ARGS=""

for i in "${!INTERVAL_FILES[@]}"; do
    INTERVAL="${INTERVAL_FILES[$i]}"
    INTERVAL_BED="${INTERVAL}.bed"
    if [ ! -f "$INTERVAL_BED" ]; then
        awk -v OFS='\t' '!/^@/ {print $1, $2-1, $3}' "$INTERVAL" > "$INTERVAL_BED"
    fi

    ID=$(printf "%04d" $i) # 0000, 0001...
    CHUNK_FINAL="gatk_tmp/chunks/${OUTPREFIX}.${ID}.vcf.gz"
    CHUNK_LOG="gatk_tmp/logs/chunk.${ID}.log"

    CHUNK_GATHER_ARGS="${CHUNK_GATHER_ARGS} -I $CHUNK_FINAL"

    # Skip if chunk already exists
    if [ -f "$CHUNK_FINAL" ] && [ -f "${CHUNK_FINAL}.tbi" ]; then
        log_info "Chunk $ID exists."
        continue
    fi

    (
        P1_BED="gatk_tmp/intervals/${ID}.p1.bed"
        P2_BED="gatk_tmp/intervals/${ID}.p2.bed"

        P1_VCF="gatk_tmp/chunks/${ID}.p1.vcf.gz"
        P1_M_VCF="gatk_tmp/chunks/${ID}.p1.m.vcf.gz"
        P1_F_VCF="gatk_tmp/chunks/${ID}.p1.f.vcf.gz"

        P2_VCF="gatk_tmp/chunks/${ID}.p2.vcf.gz"

        # P1 (Ploidy 1 regions) = Intersection with Haploid BED
        bedtools intersect -a "$INTERVAL_BED" -b "$HAPLOID_BED" > "$P1_BED"
        # P2 (Ploidy 2 regions) = Subtraction of Haploid BED
        bedtools subtract -a "$INTERVAL_BED" -b "$HAPLOID_BED" > "$P2_BED"

        LISTS_TO_MERGE=()

        # --- Ploidy 1 ---
        if [ -s "$P1_BED" ]; then
            # For P1 we need to run Males (Ploidy 1) and Females (Ploidy 2) separately

            # A. Males (Ploidy 1)
            if [ ${#MALE_BAMS_ARGS[@]} -gt 0 ]; then
                gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
                    -R "$REF" "${MALE_BAMS_ARGS[@]}" -O "$P1_M_VCF" \
                    --sample-ploidy 1 -L "$P1_BED" >> "$CHUNK_LOG" 2>&1
            fi

            # B. Females (Ploidy 2)
            if [ ${#FEMALE_BAMS_ARGS[@]} -gt 0 ]; then
                gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
                    -R "$REF" "${FEMALE_BAMS_ARGS[@]}" -O "$P1_F_VCF" \
                    --sample-ploidy 2 -L "$P1_BED" >> "$CHUNK_LOG" 2>&1
            fi

            # C. Merge Males and Females for P1
            if [ -f "$P1_M_VCF" ] && [ -f "$P1_F_VCF" ]; then
                # Merge VCFs
                bcftools merge -O z -o "${P1_VCF%.vcf.gz}.unsorted.vcf.gz" "$P1_M_VCF" "$P1_F_VCF" >> "$CHUNK_LOG" 2>&1
                # Sort samples
                bcftools view -s <(printf "%s," "${SAMPLE_ORDER[@]}") -O z -o "$P1_VCF" "${P1_VCF%.vcf.gz}.unsorted.vcf.gz"
                gatk IndexFeatureFile -I "$P1_VCF"
            elif [ -f "$P1_M_VCF" ]; then
                mv "$P1_M_VCF" "$P1_VCF" && mv "${P1_M_VCF}.tbi" "${P1_VCF}.tbi"
            elif [ -f "$P1_F_VCF" ]; then
                mv "$P1_F_VCF" "$P1_VCF" && mv "${P1_F_VCF}.tbi" "${P1_VCF}.tbi"
            fi

            if [ -f "$P1_VCF" ]; then
                LISTS_TO_MERGE+=("-I" "$P1_VCF")
            fi

            # rm -f "$P1_M_VCF" "$P1_F_VCF" "${P1_M_VCF}.tbi" "${P1_F_VCF}.tbi"
        fi

        # --- Ploidy 2 ---
        if [ -s "$P2_BED" ]; then
            gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
                -R "$REF" "${ALL_BAMS_ARGS[@]}" -O "$P2_VCF" \
                --sample-ploidy 2 -L "$P2_BED" >> "$CHUNK_LOG" 2>&1

            if [ -f "$P2_VCF" ]; then
                LISTS_TO_MERGE+=("-I" "$P2_VCF")
            fi
        fi

        gatk --java-options "-Xmx4g" MergeVcfs \
            "${LISTS_TO_MERGE[@]}" \
            -O "$CHUNK_FINAL" >> "$CHUNK_LOG" 2>&1
        gatk IndexFeatureFile -I "$CHUNK_FINAL" >> "$CHUNK_LOG" 2>&1

        # rm -f "$P1_BED" "$P2_BED" "$P1_VCF" "${P1_VCF}.tbi" "$P2_VCF" "${P2_VCF}.tbi"

    ) &

done

wait

# ---- 3. Gather Final VCF ----
log_info "[3/3] Gather all chunks..."
FINAL_VCF="${OUTPREFIX}.vcf.gz"

if [ ! -f "$FINAL_VCF" ]; then
    gatk --java-options "-Xmx8g" GatherVcfs \
        $CHUNK_GATHER_ARGS \
        -O "$FINAL_VCF"
    gatk IndexFeatureFile -I "$FINAL_VCF"
    log_info "Done! Result: $FINAL_VCF"
else
    log_info "Final VCF already exists."
fi
