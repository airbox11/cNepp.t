#!/bin/bash

module load bwa/0.7.15
module load java/1.8.0_131
module load perl/5.20.2 

mkdir ${RESULTS_DIR}/HLAVBSeq
cd ${RESULTS_DIR}/HLAVBSeq
HLAVBSeq_DIR=${PIPELINE_DIR}/external_scripts/HLAVBSeq

#FASTQ1=${RESULTS_DIR}/*total_extracted_1.fastq
#FASTQ2=${RESULTS_DIR}/*total_extracted_2.fastq

FASTQ1="/icgc/dkfzlsdf/analysis/G200/pfitzerl/test/K28A-1QNKZZ/buffy/K28A-1QNKZZ_total_extracted_1.fastq"
FASTQ2="/icgc/dkfzlsdf/analysis/G200/pfitzerl/test/K28A-1QNKZZ/buffy/K28A-1QNKZZ_total_extracted_2.fastq"



echo $FASTQ1
echo $FASTQ2

# Alignment by BWA-MEM allowing multiple alignments for each read

echo "bwa index ${HLAVBSeq_DIR}/hla_gen.fasta"
bwa index ${HLAVBSeq_DIR}/hla_gen.fasta
echo "bwa mem -t 10 -P -L 10000 -a ${HLAVBSeq_DIR}/hla_gen.fasta $FASTQ1 $FASTQ2 > ${RESULTS_DIR}/HLAVBSeq/${NAME}_HLAVBSeq.sam"
bwa mem -t 10 -P -L 10000 -a ${HLAVBSeq_DIR}/hla_gen.fasta $FASTQ1 $FASTQ2 > ${RESULTS_DIR}/HLAVBSeq/${NAME}_HLAVBSeq.sam

# Estimation of HLA types by HLA-VBSeq
# For paired-end read data:

echo "java -jar ${HLAVBSeq_DIR}/HLAVBSeq.jar ${HLAVBSeq_DIR}/hla_gen.fasta ${RESULTS_DIR}/HLAVBSeq/${NAME}_HLAVBSeq.sam ${RESULTS_DIR}/HLAVBSeq/${NAME}_HLAVBSeq_result.txt --alpha_zero 0.01 --is_paired"
java -jar ${HLAVBSeq_DIR}/HLAVBSeq.jar ${HLAVBSeq_DIR}/hla_gen.fasta ${RESULTS_DIR}/HLAVBSeq/${NAME}_HLAVBSeq.sam ${RESULTS_DIR}/HLAVBSeq/${NAME}_HLAVBSeq_result.txt --alpha_zero 0.01 --is_paired

# Here, alpha_zero is a hyperparameter as described in the paper and we recommend to use 0.01.


echo "perl ${HLAVBSeq_DIR}/parse_result.pl ${HLAVBSeq_DIR}/Allelelist.txt ${RESULTS_DIR}/HLAVBSeq/${NAME}_HLAVBSeq_result.txt | sort -k2 -n -r > ${RESULTS_DIR}/HLAVBSeq/${NAME}_HLAVBSeq_result_parsed.txt"
perl ${HLAVBSeq_DIR}/parse_result.pl ${HLAVBSeq_DIR}/Allelelist.txt ${RESULTS_DIR}/HLAVBSeq/${NAME}_HLAVBSeq_result.txt | sort -k2 -n -r > ${RESULTS_DIR}/HLAVBSeq/${NAME}_HLAVBSeq_result_parsed.txt


rm ${RESULTS_DIR}/HLAVBSeq/*.sam
#rm ${RESULTS_DIR}/HLAVBSeq/${NAME}_HLAVBSeq_result.txt