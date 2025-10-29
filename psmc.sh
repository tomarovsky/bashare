#!/bin/bash
# based on https://github.com/atotickov/BerryTart/bash_scr/PSMS.sh
# mamba create -n PSMC -c bioconda -c conda-forge -c mahajrod routoolpa mace mavr parallel gnuplot samtools=0.1.19
# export PATH="${TOOLS}/psmc:${TOOLS}/psmc/utils:${TOOLS}/bedtools-2.31.1/bin:${TOOLS}/htslib-1.22.1/bin:${PATH}"

ASSEMBLY=""
BAM_FILE=""
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

print_usage() {
    echo "Required options:"
    echo " -f   Genome assembly file in .fasta format (with .fasta file there should be a .fasta.fai file)"
    echo " -b   Full path to the alignment .bam file (a .bam.bai file must exist)"
    echo " -w   Full path to the sample's whole genome stats file (e.g., _whole_genome_stats.csv)"
    echo " -m   Full path to the BED mask file"
    echo " -g   Generation time"
    echo " -u   Mutation rate"
    echo " -c   ChrX scaffold name"
    echo ""
    echo "Optional:"
    echo " -j   Number of threads for parallel processes (Default: ${THREADS})"
    echo " -N   PSMC parameter N (Max number of iterations) (Default: ${N_PSMC})"
    echo " -t   PSMC parameter t (Max number of hidden states) (Default: ${T_PSMC})"
    echo " -r   PSMC parameter r (Recombination rate parameter) (Default: ${R_PSMC})"
    echo " -p   PSMC time segment pattern (e.g., '4+25*2+4+6') (Default: \"${PATTERN}\")"
    echo " -z   (Optional) Run 100 bootstrap analysis. Default: ${ROUNDS}"
    echo ""
    echo "Usage: $0 -f assembly.fasta -b sample.bam -w whole_genome_stats.csv -m mask.bed -g 5 -u 4.64e-9 -c scaff_19"
}

# --- Process command-line options ---
while getopts 'f:b:w:m:g:u:c:j:n:t:r:p:Z' flag; do
    case "${flag}" in
        f) ASSEMBLY="${OPTARG}" ;;
        b) BAM_FILE="${OPTARG}" ;;
        w) STATS_FILE="${OPTARG}" ;;
        m) MASK_FILE="${OPTARG}" ;;
        g) GEN_TIME="${OPTARG}" ;;
        u) MU_RATE="${OPTARG}" ;;
        c) CHRX_ID="${OPTARG}" ;;
        j) THREADS="${OPTARG}" ;;
        n) N_PSMC="${OPTARG}" ;;
        t) T_PSMC="${OPTARG}" ;;
        r) R_PSMC="${OPTARG}" ;;
        p) PATTERN="${OPTARG}" ;;
        Z) ROUNDS=true ;;
        *) print_usage
           exit 1 ;;
    esac
done

# --- Check for required arguments ---
if [ -z "$ASSEMBLY" ] || [ -z "$BAM_FILE" ] || [ -z "$STATS_FILE" ] || [ -z "$GEN_TIME" ] || [ -z "$MU_RATE" ] || [ -z "$CHRX_ID" ]; then
    echo "Error: Missing required arguments."
    print_usage
    exit 1
fi

# --- Global script variables ---
WORKDIR=$(pwd)
SAMPLE=$(basename "$BAM_FILE" | cut -d. -f1)

ALL_CHR_DIR="${WORKDIR}/all_Chr/${SAMPLE}"
NO_CHRX_DIR="${WORKDIR}/no_ChrX/${SAMPLE}"

# --- Setup directories ---
mkdir -p "${ALL_CHR_DIR}/split/bcf"
mkdir -p "${NO_CHRX_DIR}"

echo "$(date) | Sample: ${SAMPLE}"

# Variant Calling and VCF/FQ creation (all chromosomes)
echo "$(date) | ${SAMPLE} | Variant calling"
prepare_region_list.py -r "${ASSEMBLY}.fai" -s -m 1500000 -n 1 -g samtools -x 1000 2>/dev/null | \
    parallel -j "${THREADS}" "samtools mpileup -C50 -uf ${ASSEMBLY} -r {} ${BAM_FILE} | bcftools view -b -c - > ${ALL_CHR_DIR}/split/bcf/tmp.{#}.bcf"

echo "$(date) | ${SAMPLE} | BCF -> VCF"
bcftools cat $(ls ${ALL_CHR_DIR}/split/bcf/tmp.*.bcf | sort -V) | bcftools view - | gzip > "${ALL_CHR_DIR}/${SAMPLE}.vcf.gz"
rm -r "${ALL_CHR_DIR}/split"

echo "$(date) | ${SAMPLE} | VCF masking"
bedtools intersect -header -v -a "${ALL_CHR_DIR}/${SAMPLE}.vcf.gz" -b ${MASK_FILE} | bgzip -c >"${ALL_CHR_DIR}/${SAMPLE}.masked.vcf.gz"

# -D and -d parameters from whole genome stats file
MIN_DEPTH=$(awk 'NR==2{printf "%.0f", $2/3}' "${STATS_FILE}")
MAX_DEPTH=$(awk 'NR==2{printf "%.0f", $2*2.5}' "${STATS_FILE}")
echo "$(date) | Consensus file | -d:${MIN_DEPTH} -D:${MAX_DEPTH}"
zcat "${ALL_CHR_DIR}/${SAMPLE}.masked.vcf.gz" | vcfutils.pl vcf2fq -d "${MIN_DEPTH}" -D "${MAX_DEPTH}" | gzip > "${ALL_CHR_DIR}/${SAMPLE}.fq.gz"

# PSMC (all chromosomes)
echo "$(date) | ${SAMPLE} | Fasta-like consensus file preparation";
fq2psmcfa -q20 "${ALL_CHR_DIR}/${SAMPLE}.fq.gz" > "${ALL_CHR_DIR}/${SAMPLE}.diploid.psmcfa"

echo "$(date) | ${SAMPLE} | Fasta-like consensus file preparation for bootstrapping";
splitfa "${ALL_CHR_DIR}/${SAMPLE}.diploid.psmcfa" > "${ALL_CHR_DIR}/${SAMPLE}.diploid.split.psmcfa"

echo "$(date) | ${SAMPLE} | PSMC calculation";
psmc -N"${N_PSMC}" -t"${T_PSMC}" -r"${R_PSMC}" -p "${PATTERN}" -o "${ALL_CHR_DIR}/${SAMPLE}.diploid.psmc" "${ALL_CHR_DIR}/${SAMPLE}.diploid.psmcfa"

echo "$(date) | ${SAMPLE} | PSMC plot";
psmc_plot.pl -u "${MU_RATE}" -g "${GEN_TIME}" -R "${ALL_CHR_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.diploid" "${ALL_CHR_DIR}/${SAMPLE}.diploid.psmc"
mv ${ALL_CHR_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.diploid.0.txt ${ALL_CHR_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.diploid.txt
rm ${ALL_CHR_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.diploid.{eps,gp,par}

# Bootstrapping (optional -r)
if $ROUNDS; then
    echo "$(date) | ${SAMPLE} | PSMC 100 bootstrapping"
    seq 100 | xargs -I{} echo "psmc -N${N_PSMC} -t${T_PSMC} -r${R_PSMC} -b -p \"${PATTERN}\" -o ${ALL_CHR_DIR}/round-{}.psmc ${ALL_CHR_DIR}/${SAMPLE}.diploid.split.psmcfa" > "${ALL_CHR_DIR}/seq100.txt"
    parallel -j "${THREADS}" < "${ALL_CHR_DIR}/seq100.txt"
    cat ${ALL_CHR_DIR}/round-*.psmc > "${ALL_CHR_DIR}/${SAMPLE}.round.psmc"
    rm ${ALL_CHR_DIR}/round-*.psmc
    rm ${ALL_CHR_DIR}/seq100.txt

    echo "$(date) | ${SAMPLE} | PSMC bootstrap plot";
    psmc_plot.pl -u "${MU_RATE}" -g "${GEN_TIME}" -R "${ALL_CHR_DIR}/round" "${ALL_CHR_DIR}/${SAMPLE}.round.psmc"
    cat ${ALL_CHR_DIR}/round-*.txt > "${ALL_CHR_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.round.txt"
    rm ${ALL_CHR_DIR}/round-*.txt
fi

# PSMC excluding ChrX
echo "$(date) | ${SAMPLE} | Removing ${CHRX_ID}"
zcat "${ALL_CHR_DIR}/${SAMPLE}.fq.gz" | \
    awk -v ID_TO_REMOVE="${CHRX_ID}" 'BEGIN {RS="@"; ORS="@"} /^$/ {ORS=""; next} { if ($0 !~ "^" ID_TO_REMOVE "[ \t\n]") { print $0 }}' | gzip > "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.fq.gz"

echo "$(date) | ${SAMPLE} | Fasta-like consensus file preparation";
fq2psmcfa -q20 "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.fq.gz" > "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.diploid.psmcfa"

echo "$(date) | ${SAMPLE} | Fasta-like consensus file preparation for bootstrapping";
splitfa "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.diploid.psmcfa" > "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.diploid.split.psmcfa"

echo "$(date) | ${SAMPLE} | PSMC calculation";
psmc -N"${N_PSMC}" -t"${T_PSMC}" -r"${R_PSMC}" -p "${PATTERN}" -o "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.diploid.psmc" "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.diploid.psmcfa"

echo "$(date) | ${SAMPLE} | PSMC plot";
psmc_plot.pl -u "${MU_RATE}" -g "${GEN_TIME}" -R "${NO_CHRX_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.no_ChrX.diploid" "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.diploid.psmc"
mv ${NO_CHRX_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.no_ChrX.diploid.0.txt ${NO_CHRX_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.no_ChrX.diploid.txt
rm ${NO_CHRX_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.no_ChrX.diploid.{eps,gp,par}

# Bootstrapping (optional -r)
if $ROUNDS; then
    echo "$(date) | ${SAMPLE} | PSMC 100 bootstrapping"
    seq 100 | xargs -I{} echo "psmc -N${N_PSMC} -t${T_PSMC} -r${R_PSMC} -b -p '${PATTERN}' -o ${NO_CHRX_DIR}/round-{}.psmc ${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.diploid.split.psmcfa" > "${NO_CHRX_DIR}/seq100.txt"
    parallel -j "${THREADS}" < "${NO_CHRX_DIR}/seq100.txt"
    cat ${NO_CHRX_DIR}/round-*.psmc > "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.round.psmc"
    rm ${NO_CHRX_DIR}/round-*.psmc
    rm ${NO_CHRX_DIR}/seq100.txt

    echo "$(date) | ${SAMPLE} | PSMC bootstrap plot";
    psmc_plot.pl -u "${MU_RATE}" -g "${GEN_TIME}" -R "${NO_CHRX_DIR}/round" "${NO_CHRX_DIR}/${SAMPLE}.no_ChrX.round.psmc"
    cat ${NO_CHRX_DIR}/round.*.txt > "${NO_CHRX_DIR}/${SAMPLE}.G${GEN_TIME}_U${MU_RATE}.no_ChrX.round.txt"
    rm ${NO_CHRX_DIR}/round.*.txt
fi

echo "$(date) | DONE | ${SAMPLE} | Gen=${GEN_TIME} | Mu=${MU_RATE} | Rounds=${ROUNDS}"
