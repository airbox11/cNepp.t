#!/bin/bash

mkdir ${RESULTS_DIR}/Optitype
cd $RESULTS_DIR/Optitype

OPTITYPE_HOME=/icgc/dkfzlsdf/analysis/G200/immuno/tools/optitype/OptiType-master

#FASTQ1=${RESULTS_DIR}/*total_extracted_1.fastq
#FASTQ2=${RESULTS_DIR}/*total_extracted_2.fastq

BAM_DIR=`dirname ${PATH_TO_BAM}`

READ1=`ls ${BAM_DIR}/../run*/sequence/*R1.fastq.gz`
READ2=`ls ${BAM_DIR}/../run*/sequence/*R2.fastq.gz`

if [ ! -f "${RESULTS_DIR}/Optitype/20*/*result.tsv" ];then
    #### create a named pipe for decompress fastq files

    mkfifo $SCRATCHDIR/$LSB_JOBID/tmp_fifo_READ1.fastq
    mkfifo $SCRATCHDIR/$LSB_JOBID/tmp_fifo_READ2.fastq

    zcat ${READ1} >  $SCRATCHDIR/$LSB_JOBID/tmp_fifo_READ1.fastq &
    zcat ${READ2} >  $SCRATCHDIR/$LSB_JOBID/tmp_fifo_READ2.fastq &    
    
    singularity run \
        -B "$SCRATCHDIR/$LSB_JOBID:$SCRATCHDIR/$LSB_JOBID:ro" -B "$RESULTS_DIR:$RESULTS_DIR:rw" \
        "${OPTITYPE_SIF}" -i ${FASTQ1} ${FASTQ2} -d -o ${RESULTS_DIR}/Optitype

    #### remove temporary files
    rm $SCRATCHDIR/$LSB_JOBID/tmp_fifo_READ1.fastq
    rm $SCRATCHDIR/$LSB_JOBID/tmp_fifo_READ2.fastq
else
    echo "${RESULTS_DIR}/Optitype/20*/*result.tsv is already available, running of Optitype stops"
fi

