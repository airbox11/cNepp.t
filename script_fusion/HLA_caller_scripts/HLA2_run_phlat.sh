#!/bin/bash

if [ -f ${RESULTS_DIR}/PHLAT/PHLAT_HLA.sum ];then
	echo "PHLAT_HLA.sum already exists."
	exit 0
fi

module load python/2.7.9

#changes user's pysam version!
pip install --user --upgrade pysam==0.8.3

PHLAT_DIR=/icgc/dkfzlsdf/analysis/G200/immuno/tools/phlat-release
INDEX_DIR=/icgc/dkfzlsdf/analysis/G200/immuno/tools/phlat-release/b2folder
B2URL=/ibios/tbi_cluster/13.1/x86_64/bowtie2/bowtie2-2.2.1/bin/bowtie2

FASTQ1=${RESULTS_DIR}/*total_extracted_1.fastq
FASTQ2=${RESULTS_DIR}/*total_extracted_2.fastq

mkdir ${RESULTS_DIR}/PHLAT

echo $RESULTS_DIR/PHLAT
echo $FASTQ1
echo ${NAME}_PHLAT

python -O ${PHLAT_DIR}/dist/PHLAT.py -1 $FASTQ1 -2 $FASTQ2 -index $INDEX_DIR -b2url $B2URL -tag PHLAT -e $PHLAT_DIR -o $RESULTS_DIR/PHLAT

### Note: the length of the -tag inut variable seems to be limited and can lead to an 'index out of range' error