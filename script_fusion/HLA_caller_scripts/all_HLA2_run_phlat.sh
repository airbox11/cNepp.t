#!/bin/bash

module load python/2.7.9

#changes user's pysam version!
pip install --user --upgrade pysam==0.8.3

PHLAT_DIR=/icgc/dkfzlsdf/analysis/G200/immuno/tools/phlat-release
INDEX_DIR=/icgc/dkfzlsdf/analysis/G200/immuno/tools/phlat-release/b2folder
B2URL=/ibios/tbi_cluster/13.1/x86_64/bowtie2/bowtie2-2.2.1/bin/bowtie2

FASTQ1=
FASTQ2=

mkdir ${RESULTS_DIR}/PHLAT


BAM_DIR=`dirname ${PATH_TO_BAM}`

READ1=`ls ${BAM_DIR}/../run*/sequence/*R1.fastq.gz`
READ2=`ls ${BAM_DIR}/../run*/sequence/*R2.fastq.gz`

FASTQ1=/icgc/dkfzlsdf/analysis/G200/pfitzerl/data/HLA/K28A-1CNM9W/buffy/K28A-1CNM9W_buffy_total_extracted_1.fastq
FASTQ2=/icgc/dkfzlsdf/analysis/G200/pfitzerl/data/HLA/K28A-1CNM9W/buffy/K28A-1CNM9W_buffy_total_extracted_2.fastq

#python -O ${PHLAT_DIR}/dist/PHLAT.py -1 ${FASTQ1} -2 ${FASTQ2} -index $INDEX_DIR -b2url $B2URL -tag ${NAME}_PHLAT -e $PHLAT_DIR -o ${RESULTS_DIR}/PHLAT

if [ ! -f "${RESULTS_DIR}/PHLAT/*HLA.sum" ];then
    #### create a named pipe for decompress fastq files

    mkfifo ${PBS_SCRATCH_DIR}/${PBS_JOBID}/tmp_fifo_READ1.fastq
    mkfifo ${PBS_SCRATCH_DIR}/${PBS_JOBID}/tmp_fifo_READ2.fastq

    zcat ${READ1} >  ${PBS_SCRATCH_DIR}/${PBS_JOBID}/tmp_fifo_READ1.fastq &
    zcat ${READ2} >  ${PBS_SCRATCH_DIR}/${PBS_JOBID}/tmp_fifo_READ2.fastq &  

    #### to run software PHLAT for hla typing
    
    python -O ${PHLAT_DIR}/dist/PHLAT.py -1 ${PBS_SCRATCH_DIR}/${PBS_JOBID}/tmp_fifo_READ1.fastq -2 ${PBS_SCRATCH_DIR}/${PBS_JOBID}/tmp_fifo_READ2.fastq -index $INDEX_DIR -b2url $B2URL -tag PHLAT -e $PHLAT_DIR -o ${RESULTS_DIR}/PHLAT
    
    
    #### remove temporary files
    rm ${PBS_SCRATCH_DIR}/${PBS_JOBID}/tmp_fifo_READ1.fastq
    rm ${PBS_SCRATCH_DIR}/${PBS_JOBID}/tmp_fifo_READ2.fastq
else
    echo "${RESULTS_DIR}/PHLAT/*HLA.sum is already available, running of Phlat stops"
fi
