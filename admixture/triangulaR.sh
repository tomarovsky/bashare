#!/bin/bash
# Usage: ./triangulaR.sh VCF_FILE POP_FILE DIFF
# Example: ./triangulaR.sh mzib.vcf mzib.pop 0.75

set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 VCF_FILE POP_FILE DIFF"
    echo "Example: $0 mzib.allsamples.filt.masked.auto.snp.plink.0.2.vcf mzib.allsamples.filt.masked.auto.snp.plink.0.2.pop 0.75"
    exit 1
fi

VCF_FILE=$1
POP_FILE=$2
DIFF=$3

# mamba create -n triangular -c conda-forge r-base r-vcfr r-devtools; devtools::install_github("omys-omics/triangulaR")
source $(conda info --base)/etc/profile.d/conda.sh
conda activate triangulaR

Rscript --vanilla - <<EOF
library(triangulaR)
library(vcfR)

data <- read.vcfR("$VCF_FILE", verbose = FALSE)
pops <- read.table("$POP_FILE", header = FALSE, stringsAsFactors = FALSE)
names(pops) <- c("id", "pop")

data.diff <- alleleFreqDiff(vcfR = data, pm = pops,
                            p1 = "M.martes", p2 = "M.zibellina",
                            difference = $DIFF)
hi.het <- hybridIndex(vcfR = data.diff, pm = pops,
                      p1 = "M.martes", p2 = "M.zibellina")

output_file <- paste0("triangulaR.diff_", $DIFF, ".tsv")
write.table(hi.het, file = output_file, quote = FALSE, row.names = FALSE, col.names = TRUE)
cat("Saved hybridIndex for difference =", $DIFF, "to", output_file, "\n")
EOF
