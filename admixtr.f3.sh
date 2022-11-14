#!/usr/bin/bash

set -u -e

export PATH=${PATH}:/mnt/tank/scratch/atomarovsky/tools/AdmixTools-7.0.2/bin

INDIR="/mnt/tank/scratch/skliver/common/mustelidae/martes_martes/genome/admixture/autosomes_and_PAR/admixtr"
OUTDIR="/mnt/tank/scratch/skliver/common/mustelidae/martes_martes/genome/admixture/autosomes_and_PAR/admixtr"
for TARGET in ${1}; do
	GENO=${INDIR}/mmar.allsamples.filt.mapq10.max250.min33.intersect.2.merge.all.snp.autosomes_and_PAR.admixtr.geno
	SNP=${INDIR}/mmar.allsamples.filt.mapq10.max250.min33.intersect.2.merge.all.snp.autosomes_and_PAR.admixtr.snp
	IND=${INDIR}/mmar.allsamples.filt.mapq10.max250.min33.intersect.2.merge.all.snp.autosomes_and_PAR.admixtr.ind
	POPLIST=${OUTDIR}/${TARGET}.f3stats.poplist.txt
	PARAMSFILE=$OUTDIR/$TARGET.f3stats.qp3Pop.params.txt
	OUT=$OUTDIR/$TARGET.qp3Pop.out
	for SAMPLE_1 in 10xmmar S44 S46 S49 S50 T149 T151 10xmzib china S26 T104 T118 T148 T150 T18 T194 T26 T50 T72 T8 T90 T24 T76 T77 T78 T79 T81 T82 T83 T84 T85 T86 T87; do
		printf "" > $POPLIST
		for SAMPLE_2 in 10xmmar S44 S46 S49 S50 T149 T151 10xmzib china S26 T104 T118 T148 T150 T18 T194 T26 T50 T72 T8 T90 T24 T76 T77 T78 T79 T81 T82 T83 T84 T85 T86 T87; do
			printf "${SAMPLE_1}\t${SAMPLE_2}\t${TARGET}\n" >> $POPLIST
		done
		printf "genotypename:\t$GENO\n" > $PARAMSFILE
		printf "snpname:\t$SNP\n" >> $PARAMSFILE
		printf "indivname:\t$IND\n" >> $PARAMSFILE
		printf "popfilename:\t$POPLIST\n" >> $PARAMSFILE
		qp3Pop -p $PARAMSFILE >> $OUT
	done
done




