#!/bin/bash

if [ -f ${RESULTS_DIR}/Polysolver/winners.hla.txt ];then
	echo "winners.hla.txt already exists."
    sleep 1m
	exit 0
fi

module load samtools/1.6
module load perl/5.20.2

mkdir ${RESULTS_DIR}/Polysolver
cd ${RESULTS_DIR}/Polysolver

PSHOME=/icgc/dkfzlsdf/analysis/G200/immuno/tools/polysolver/polysolver

#source config file
source ${PSHOME}/scripts/config.bash

#genotyping:polysolver
#run on  data (bam-race-includeFreq-build-format-insertCalc-outDir)
${PSHOME}/scripts/shell_call_hla_type $PATH_TO_BAM Unknown 1 hg19 STDFQ 0 $RESULTS_DIR/Polysolver

rm ${RESULTS_DIR}/Polysolver/*.bam
rm ${RESULTS_DIR}/Polysolver/hla_* 
rm ${RESULTS_DIR}/Polysolver/merged.*.fastq 
rm ${RESULTS_DIR}/Polysolver/*.sam
rm ${RESULTS_DIR}/Polysolver/ids_*
rm ${RESULTS_DIR}/Polysolver/tag.*.fastq
rm ${RESULTS_DIR}/Polysolver/nv.complete.chr6region.R0k6.csorted.bam.bai