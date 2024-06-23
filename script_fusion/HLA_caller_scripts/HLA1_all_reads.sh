#!/bin/bash

module load samtools/1.6
module load bedtools/2.24.0


# Convert bam to fastq

samtools sort -n ${PATH_TO_BAM} -o ${RESULTS_DIR}/${NAME}.aln.qsort.bam -m 10G -@ 10

bedtools bamtofastq -i ${RESULTS_DIR}/${NAME}.aln.qsort.bam -fq ${RESULTS_DIR}/${NAME}_total_extracted_1.fastq -fq2 ${RESULTS_DIR}/${NAME}_total_extracted_2.fastq


rm ${RESULTS_DIR}/*.bam