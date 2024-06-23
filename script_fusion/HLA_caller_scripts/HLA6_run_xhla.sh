#!/bin/bash

mkdir ${RESULTS_DIR}/xHLA
cd ${RESULTS_DIR}/xHLA

module load samtools/1.3.1
module load bwa/0.7.15 
module load R/3.4.0
module load bedtools/2.24.0
module load sambamba/0.6.6
module load python/2.7.9

export PYTHONPATH=/icgc/dkfzlsdf/analysis/G200/pfitzerl/software/CrossMap-0.2.8/lib:$PYTHONPATH

pip install --user --upgrade pysam==0.10.0
## works with pysam version 0.10.0 (pip install --user --upgrade pysam==0.10.0)

# add chr prefix to bam file of extracted reads for CrossMap.py 
echo xHLA_test_runtime

samtools view -H ${RESULTS_DIR}/${PID}_intermediate.bam | sed -e 's/SN:\([0-9XY]*\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/'  > ${RESULTS_DIR}/xHLA/${NAME}_mapped_extracted_chr.sam

samtools view ${RESULTS_DIR}/${NAME}_mapped_extracted.bam  | awk 'BEGIN {OFS="\t"} { print $1,$2,"chr"$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15 }' >> ${RESULTS_DIR}/xHLA/${NAME}_mapped_extracted_chr.sam

samtools view -b ${RESULTS_DIR}/xHLA/${NAME}_mapped_extracted_chr.sam > ${RESULTS_DIR}/xHLA/${NAME}_mapped_extracted_chr.bam

# use CrossMap.py to remap the reads from hg19 to hg38

CrossMap.py bam ${SOFTWARE_DIR}/CrossMap-0.2.8/data/UCSC_chain/hg19ToHg38.over.chain.gz ${RESULTS_DIR}/xHLA/${NAME}_mapped_extracted_chr.bam ${RESULTS_DIR}/xHLA/${NAME}_mapped_extracted_chr_hg38

#run xHLA

${SOFTWARE_DIR}/HLA/bin/typer.sh ${RESULTS_DIR}/xHLA/${NAME}_mapped_extracted_chr_hg38.sorted.bam ${NAME}_mapped_extracted_chr_hg38

#rm ${RESULTS_DIR}/xHLA/${NAME}_mapped_extracted_chr.sam
#rm ${RESULTS_DIR}/xHLA/${NAME}_mapped_extracted_chr.bam