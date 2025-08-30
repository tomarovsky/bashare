#!/usr/bin/bash

set -u -e

export PATH=${PATH}:/mnt/tank/scratch/atomarovsky/tools/AdmixTools-7.0.2/bin

INDIR="/mnt/tank/scratch/skliver/common/mustelidae/martes_zibellina/genome/varcall/hic.purged/bcftools/mzib_and_mfoi/plink/admixtr/"
OUTDIR="/mnt/tank/scratch/skliver/common/mustelidae/martes_zibellina/genome/varcall/hic.purged/bcftools/mzib_and_mfoi/plink/admixtr/"
OUTGROUP="mfoi"
for TARGET in ${1}; do
	GENO=${INDIR}/mzib.mfoi.allsamples.filt.namefix.mask.aut_and_PAR.snp.plink.geno
	SNP=${INDIR}/mzib.mfoi.allsamples.filt.namefix.mask.aut_and_PAR.snp.plink.snp
	IND=${INDIR}/mzib.mfoi.allsamples.filt.namefix.mask.aut_and_PAR.snp.plink.ind
	POPLIST=${OUTDIR}/${TARGET}.f4.poplist.txt
	PARAMSFILE=$OUTDIR/$TARGET.f4.params.txt
	OUT=$OUTDIR/$TARGET.groups.f4.out
	for SAMPLE_1 in 10xmmar S44 S46 S49 T149 T24 T76 T77 T82; do
		printf "" > $POPLIST
		for SAMPLE_2 in 10xmzib S26 T8 T26 T50 T72 T90 T104 T118 T148 T150 T194 china; do
			printf "${TARGET}\t${SAMPLE_1}\t${SAMPLE_2}\t${OUTGROUP}\n" >> $POPLIST
		done
		printf "genotypename:\t$GENO\n" > $PARAMSFILE
		printf "snpname:\t$SNP\n" >> $PARAMSFILE
		printf "indivname:\t$IND\n" >> $PARAMSFILE
		printf "popfilename:\t$POPLIST\n" >> $PARAMSFILE
		printf "f4mode: YES\n" >> $PARAMSFILE
		qpDstat -p $PARAMSFILE >> $OUT
	done
done




