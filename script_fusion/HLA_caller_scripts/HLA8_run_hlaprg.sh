#!/bin/bash

module load boost/1.64.0
module load bwa/0.7.15
module load samtools/1.6
module load picard/1.61
module load gcc/7.2.0

export PATH=$PATH:/icgc/dkfzlsdf/analysis/G200/pfitzerl/software/bin/bin

HLAPRG_DIR=${SOFTWARE_DIR}/HLA-PRG-LA/src

mkdir ${RESULTS_DIR}/HLA_PRG
cd ${RESULTS_DIR}/HLA_PRG

${HLAPRG_DIR}/HLA-PRG-LA.pl --BAM ${PATH_TO_BAM} --graph PRG_MHC_GRCh38_withIMGT --sampleID ${SAMPLE} --maxThreads 7 --workingDir ${RESULTS_DIR}/HLA_PRG

rm ${RESULTS_DIR}/HLA_PRG/*/*bam
rm ${RESULTS_DIR}/HLA_PRG/*/*bai