#!/bin/bash
# based on https://github.com/atotickov/BerryTart/bash_scr/PSMS.sh
# mamba create -n PSMC -c bioconda -c conda-forge -c mahajrod routoolpa mace mavr parallel gnuplot samtools=0.1.19
# export PATH="${TOOLS}/psmc:${TOOLS}/psmc/utils:${TOOLS}/Biocrutch/scripts/psmc/:${PATH}"

source "$TOOLS/bashare/lib/log_functions.sh"

ASSEMBLY=""
BAM_FILE=""
VCF_FILE=""
UNMASKED_FQ=""
STATS_FILE=""
MASK_FILE=""
GEN_TIME=""
MU_RATE=""
CHRX_ID=""

THREADS=32
N_PSMC=25
T_PSMC=15
R_PSMC=5
PATTERN="4+25*2+4+6"
ROUNDS=false
STOP_AFTER_VCF=false

print_usage() {
    echo "Required options:"
    echo " -f    Genome assembly file in .fasta format (with .fasta file there should be a .fasta.fai file)"
    echo " -b    Full path to the alignment .bam file (a .bam.bai file must exist) (required if -v or -q are not used)"
    echo " -v    Full path to the input SAMPLE.vcf.gz file (required if -b or -q are not used)"
    echo " -q    Full path to the unmasked consensus .fq.gz file (required if -b or -v are not used)"
    echo " -w    Full path to the sample's whole genome stats file (e.g., _whole_genome_stats.csv)"
    echo " -g    Generation time"
    echo " -u    Mutation rate"
    echo " -c    ChrX scaffold name"
    echo ""
    echo "Optional:"
    echo " -m    Full path to the BED mask file (optional, if not provided no masking will be performed)"
    echo " -j    Number of threads for parallel processes (Default: ${THREADS})"
    echo " -N    PSMC parameter N (Max number of iterations) (Default: ${N_PSMC})"
    echo " -t    PSMC parameter t (Max number of hidden states) (Default: ${T_PSMC})"
    echo " -r    PSMC parameter r (Recombination rate parameter) (Default: ${R_PSMC})"
    echo " -p    PSMC time segment pattern (e.g., '4+25*2+4+6') (Default: \"${PATTERN}\")"
    echo " -z    (Optional) Run 100 bootstrap analysis. Default: ${ROUNDS}"
    echo " -x    (Optional) Stop after generating VCF file. Default: ${STOP_AFTER_VCF}"
    echo ""
    echo "Usage: $0 -f assembly.fasta [-b sample.bam | -v sample.vcf.gz | -q sample.unmasked.fq.gz] -w whole_genome_stats.csv -g 5 -u 4.64e-9 -c scaff_19 [-m mask.bed] [-x]"
}

# --- Process command-line options ---
while getopts 'f:b:v:q:w:g:u:c:m:j:n:t:r:p:zx' flag; do
    case "${flag}" in
        f) ASSEMBLY="${OPTARG}" ;;
        b) BAM_FILE="${OPTARG}" ;;
        v) VCF_FILE="${OPTARG}" ;;
        q) UNMASKED_FQ="${OPTARG}" ;;
        w) STATS_FILE="${OPTARG}" ;;
        g) GEN_TIME="${OPTARG}" ;;
        u) MU_RATE="${OPTARG}" ;;
        c) CHRX_ID="${OPTARG}" ;;
        m) MASK_FILE="${OPTARG}" ;;
        j) THREADS="${OPTARG}" ;;
        n) N_PSMC="${OPTARG}" ;;
        t) T_PSMC="${OPTARG}" ;;
        r) R_PSMC="${OPTARG}" ;;
        p) PATTERN="${OPTARG}" ;;
        z) ROUNDS=true ;;
        x) STOP_AFTER_VCF=true ;;
        *) print_usage
           exit 1 ;;
    esac
done

# --- Check for required arguments ---
if [ -z "$BAM_FILE" ] && [ -z "$VCF_FILE" ] && [ -z "$UNMASKED_FQ" ]; then
    log_error "Missing required argument. You must provide either -b (BAM), -v (VCF) or -q (unmasked FQ)."
    print_usage
    exit 1
fi

# Check for conflicting arguments
input_count=0
[ -n "$BAM_FILE" ] && input_count=$((input_count + 1))
[ -n "$VCF_FILE" ] && input_count=$((input_count + 1))
[ -n "$UNMASKED_FQ" ] && input_count=$((input_count + 1))

if [ $input_count -gt 1 ]; then
    log_error "Conflicting arguments. You can only provide one of -b (BAM), -v (VCF) or -q (unmasked FQ)."
    print_usage
    exit 1
fi

if [ -z "$ASSEMBLY" ] || [ -z "$STATS_FILE" ] || [ -z "$GEN_TIME" ] || [ -z "$MU_RATE" ] || [ -z "$CHRX_ID" ]; then
    log_error "Missing one or more required arguments (-f, -w, -g, -u, -c)."
    print_usage
    exit 1
fi

# --- Global script variables ---
WORKDIR=$(pwd)

if [ -n "$BAM_FILE" ]; then
    SAMPLE=$(basename "$BAM_FILE" | cut -d. -f1)
elif [ -n "$VCF_FILE" ]; then
    SAMPLE=$(basename "$VCF_FILE" | cut -d. -f1)
elif [ -n "$UNMASKED_FQ" ]; then
    SAMPLE=$(basename "$UNMASKED_FQ" | cut -d. -f1)
fi

# --- Setup directories ---
ALL_CHR_DIR="${WORKDIR}/all_Chr/${SAMPLE}"
NO_CHRX_DIR="${WORKDIR}/no_ChrX/${SAMPLE}"

mkdir -p "${ALL_CHR_DIR}"
mkdir -p "${NO_CHRX_DIR}"

log_info "Working directory: ${WORKDIR}"
log_info "Sample: ${SAMPLE}"
log_info "Mask file: ${MASK_FILE:-Not provided, skipping masking}"

if [ -n "$BAM_FILE" ]; then
    # Variant Calling and VCF/FQ creation (all chromosomes)
    log_info "${SAMPLE} | Variant calling from BAM"
    mkdir -p "${ALL_CHR_DIR}/split/bcf"

    prepare_region_list.py -r "${ASSEMBLY}.fai" -s -m 1500000 -n 1 -g samtools -x 1000 2>/dev/null | \
        parallel -j "${THREADS}" "samtools mpileup -C50 -uf ${ASSEMBLY} -r {} ${BAM_FILE} 2>/dev/null | bcftools view -b -c - > ${ALL_CHR_DIR}/split/bcf/tmp.{#}.bcf"

    log_info "${SAMPLE} | BCF -> VCF"
    bcftools cat $(ls ${ALL_CHR_DIR}/split/bcf/tmp.*.bcf | sort -V) | bcftools view - | gzip > "${ALL_CHR_DIR}/${SAMPLE}.vcf.gz"
    rm -r "${ALL_CHR_DIR}/split" # Remove split dir after use

    if $STOP_AFTER_VCF; then
        log_info "${SAMPLE} | VCF file created: ${ALL_CHR_DIR}/${SAMPLE}.vcf.gz"
        log_info "STOPPING as requested by -x flag"
        exit 0
    fi

    # -D and -d parameters from whole genome stats file
    MIN_DEPTH=$(awk 'NR==2{printf "%.0f", $2/3}' "${STATS_FILE}")   # Minimum (33% of median)
    MAX_DEPTH=$(awk 'NR==2{printf "%.0f", $2*2.5}' "${STATS_FILE}") # Maximum (250% of median)
    log_info "Consensus file | -d:${MIN_DEPTH} -D:${MAX_DEPTH}"
    zcat "${ALL_CHR_DIR}/${SAMPLE}.vcf.gz" | vcfutils.pl vcf2fq -d "${MIN_DEPTH}" -D "${MAX_DEPTH}" | gzip > "${ALL_CHR_DIR}/${SAMPLE}.fq.gz"

elif [ -n "$VCF_FILE" ]; then
    log_info "${SAMPLE} | Using provided VCF file: ${VCF_FILE}"

    ln -s "$(realpath ${VCF_FILE})" "${ALL_CHR_DIR}/${SAMPLE}.vcf.gz"

    log_info "${SAMPLE} | Created symlink: $(realpath ${VCF_FILE}) -> ${ALL_CHR_DIR}/${SAMPLE}.vcf.gz"

    # -D and -d parameters from whole genome stats file
    MIN_DEPTH=$(awk 'NR==2{printf "%.0f", $2/3}' "${STATS_FILE}")   # Minimum (33% of median)
    MAX_DEPTH=$(awk 'NR==2{printf "%.0f", $2*2.5}' "${STATS_FILE}") # Maximum (250% of median)
    log_info "Consensus file | -d:${MIN_DEPTH} -D:${MAX_DEPTH}"
    zcat "${ALL_CHR_DIR}/${SAMPLE}.vcf.gz" | vcfutils.pl vcf2fq -d "${MIN_DEPTH}" -D "${MAX_DEPTH}" | gzip > "${ALL_CHR_DIR}/${SAMPLE}.fq.gz"

elif [ -n "$UNMASKED_FQ" ]; then
    log_info "${SAMPLE} | Using provided unmasked FQ file: ${UNMASKED_FQ}"

    # Create symlink to the provided unmasked FQ file
    ln -s "$(realpath ${UNMASKED_FQ})" "${ALL_CHR_DIR}/${SAMPLE}.fq.gz"

    log_info "${SAMPLE} | Created symlink: $(realpath ${UNMASKED_FQ}) -> ${ALL_CHR_DIR}/${SAMPLE}.fq.gz"
fi

# Masking (optional)
if [ -n "$MASK_FILE" ]; then
    log_info "${SAMPLE} | Consensus masking"
    mask_consensus.py -i "${ALL_CHR_DIR}/${SAMPLE}.fq.gz" -m "${MASK_FILE}" -o "${ALL_CHR_DIR}/${SAMPLE}.masked.fq.gz"
    FQ_FOR_PSMC="${ALL_CHR_DIR}/${SAMPLE}.masked.fq.gz"
    log_info "${SAMPLE} | Masked FQ file prepared: ${FQ_FOR_PSMC}"
else
    FQ_FOR_PSMC="${ALL_CHR_DIR}/${SAMPLE}.fq.gz"
    log_info "${SAMPLE} | No mask file provided, using unmasked consensus: ${FQ_FOR_PSMC}"
fi

# PSMC (all chromosomes)
log_info "${SAMPLE} | Fasta-like consensus file preparation";
fq2psmcfa -q20 "${FQ_FOR_PSMC}" > "${ALL_CHR_DIR}/${SAMPLE}.diploid.psmcfa"

log_info "${SAMPLE} | PSMC calculation";
psmc -N"${N_PSMC}" -t"${T_PSMC}" -r"${R_PSMC}" -p "${PATTERN}" -o "${ALL_CHR_DIR}/${SAMPLE}.diploid.psmc" "${ALL_CHR_DIR}/${SAMPLE}.diploid.psmcfa"

log_info "${SAMPLE} | PSMC plot";
psmc_plot.pl -u "${MU_RATE}" -g "${GEN_TIME}" -R "${ALL_CHR_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.diploid" "${ALL_CHR_DIR}/${SAMPLE}.diploid.psmc"
mv ${ALL_CHR_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.diploid.0.txt ${ALL_CHR_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.diploid.txt
rm ${ALL_CHR_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.diploid.{eps,gp,par}

# Bootstrapping (optional)
if $ROUNDS; then
    log_info "${SAMPLE} | Fasta-like consensus file preparation for bootstrapping";
    splitfa "${ALL_CHR_DIR}/${SAMPLE}.diploid.psmcfa" > "${ALL_CHR_DIR}/${SAMPLE}.diploid.split.psmcfa"
    log_info "${SAMPLE} | PSMC 100 bootstrapping"
    seq 100 | xargs -I{} echo "psmc -N${N_PSMC} -t${T_PSMC} -r${R_PSMC} -b -p \"${PATTERN}\" -o ${ALL_CHR_DIR}/round-{}.psmc ${ALL_CHR_DIR}/${SAMPLE}.diploid.split.psmcfa" > "${ALL_CHR_DIR}/seq100.txt"
    parallel -j "${THREADS}" < "${ALL_CHR_DIR}/seq100.txt"
    cat ${ALL_CHR_DIR}/round-*.psmc > "${ALL_CHR_DIR}/${SAMPLE}.round.psmc"
    rm ${ALL_CHR_DIR}/round-*.psmc
    rm ${ALL_CHR_DIR}/seq100.txt

    log_info "${SAMPLE} | PSMC bootstrap plot";
    psmc_plot.pl -u "${MU_RATE}" -g "${GEN_TIME}" -R "${ALL_CHR_DIR}/round" "${ALL_CHR_DIR}/${SAMPLE}.round.psmc"
    cat ${ALL_CHR_DIR}/round-*.txt > "${ALL_CHR_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.round.txt"
    rm ${ALL_CHR_DIR}/round-*.txt
fi

# PSMC excluding ChrX
log_info "${SAMPLE} | Removing ${CHRX_ID}"
if [ -n "$MASK_FILE" ]; then
    # Use masked file if masking was performed
    INPUT_FQ="${ALL_CHR_DIR}/${SAMPLE}.masked.fq.gz"
    log_info "Use masked FQ: ${INPUT_FQ}"
else
    # Use unmasked file if no masking
    INPUT_FQ="${ALL_CHR_DIR}/${SAMPLE}.fq.gz"
    log_info "Use UNmasked FQ: ${INPUT_FQ}"
fi

consensus_remove_scaffold.py -i "${INPUT_FQ}" -o "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.fq.gz" -s "${CHRX_ID}"

log_info "${SAMPLE} | Fasta-like consensus file preparation"
fq2psmcfa -q20 "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.fq.gz" > "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.diploid.psmcfa"

log_info "${SAMPLE} | PSMC calculation"
psmc -N"${N_PSMC}" -t"${T_PSMC}" -r"${R_PSMC}" -p "${PATTERN}" -o "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.diploid.psmc" "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.diploid.psmcfa"

log_info "${SAMPLE} | PSMC plot"
psmc_plot.pl -u "${MU_RATE}" -g "${GEN_TIME}" -R "${NO_CHRX_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.no_ChrX.diploid" "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.diploid.psmc"
mv ${NO_CHRX_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.no_ChrX.diploid.0.txt ${NO_CHRX_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.no_ChrX.diploid.txt
rm ${NO_CHRX_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.no_ChrX.diploid.{eps,gp,par}

# Bootstrapping (optional)
if $ROUNDS; then
    log_info "${SAMPLE} | Fasta-like consensus file preparation for bootstrapping";
    splitfa "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.diploid.psmcfa" > "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.diploid.split.psmcfa"
    log_info "${SAMPLE} | PSMC 100 bootstrapping"
    seq 100 | xargs -I{} echo "psmc -N${N_PSMC} -t${T_PSMC} -r${R_PSMC} -b -p '${PATTERN}' -o ${NO_CHRX_DIR}/round-{}.psmc ${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.diploid.split.psmcfa" > "${NO_CHRX_DIR}/seq100.txt"
    parallel -j "${THREADS}" < "${NO_CHRX_DIR}/seq100.txt"
    cat ${NO_CHRX_DIR}/round-*.psmc > "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.round.psmc"
    rm ${NO_CHRX_DIR}/round-*.psmc
    rm ${NO_CHRX_DIR}/seq100.txt

    log_info "${SAMPLE} | PSMC bootstrap plot";
    psmc_plot.pl -u "${MU_RATE}" -g "${GEN_TIME}" -R "${NO_CHRX_DIR}/round" "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.round.psmc"
    cat ${NO_CHRX_DIR}/round.*.txt > "${NO_CHRX_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.no_ChrX.round.txt"
    rm ${NO_CHRX_DIR}/round.*.txt
fi

log_info "Done! | ${SAMPLE} | Gen=${GEN_TIME} | Mu=${MU_RATE} | Masking=$([ -n "$MASK_FILE" ] && echo Yes || echo No) | Rounds=${ROUNDS}"
