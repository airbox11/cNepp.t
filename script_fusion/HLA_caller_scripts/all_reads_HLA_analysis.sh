#!/bin/bash

#run script: ./all_reads_HLA_analysis.sh /path/to/bam /path/to/ouput/directory project_ID

#check number of input variables
if [ $# -ne 3 ]; then
    echo 'Number of given arguments does not equal 3!'
    echo 'Usage: ./HLA_analysis_pipeline.sh /path/to/bam /path/to/ouput/directory project_ID'
    exit 1
fi

#declare input variables
PATH_TO_BAM=$1
PATH_TO_OUT_DIR=$2
PROJECT_ID=$3


# declare variables
PIPELINE_DIR=/icgc/dkfzlsdf/analysis/G200/pfitzerl/HLA-typing
SOFTWARE_DIR=/icgc/dkfzlsdf/analysis/G200/pfitzerl/software

# create folder structure
BASENAME=`basename $PATH_TO_BAM .bam`

cd $PATH_TO_OUT_DIR
PID=`echo $BASENAME | sed -n -e "s/^.*\(${PROJECT_ID}-[A-Z0-9]*\).*/\1/p"`
mkdir $PID
cd $PID
SAMPLE=`echo $BASENAME | sed -n -e "s/^\(.*\)_${PROJECT_ID}.*/\1/p"`
mkdir $SAMPLE
cd $SAMPLE


RESULTS_DIR=${PATH_TO_OUT_DIR}/${PID}/${SAMPLE}

# Extract reads mapped to HLA loci and unmapped reads and convert to fastq files
 
#HLA1_all_reads=`qsub -N HLA1_all_reads_${PID}_${SAMPLE} -l walltime=5:00:00,mem=11G,nodes=01:ppn=10 -j oe -o ${RESULTS_DIR} -v PATH_TO_BAM=${PATH_TO_BAM},NAME=${PID}_${SAMPLE},RESULTS_DIR=${RESULTS_DIR} ${PIPELINE_DIR}/HLA1_all_reads.sh`

#echo $HLA1_all_reads

# run phlat on extracted reads fastq
###########################-W depend=afterok:${HLA1_all_reads}
HLA2_run_phlat=`qsub -N HLA2_run_phlat_${PID}_${SAMPLE} -l walltime=12:00:00,mem=4G,nodes=01:ppn=9 -j oe -o ${RESULTS_DIR} -v NAME=${PID}_${SAMPLE},RESULTS_DIR=${RESULTS_DIR},PATH_TO_BAM=${PATH_TO_BAM} ${PIPELINE_DIR}/all_HLA2_run_phlat.sh`

#-W depend=afterok:${HLA1_all_reads}

echo $HLA2_run_phlat

# run HLA-VBSeq on extracted reads

#HLA3_run_HLAVBSeq=`qsub -N HLA3_run_HLAVBSeq_${PID}_${SAMPLE} -l walltime=12:00:00,mem=20G,nodes=01:ppn=10 -j oe -o ${RESULTS_DIR} -W depend=afterok:${HLA1_all_reads} -v NAME=${PID}_${SAMPLE},RESULTS_DIR=${RESULTS_DIR},PIPELINE_DIR=${PIPELINE_DIR} ${PIPELINE_DIR}/HLA3_run_HLAVBSeq.sh`

#-W depend=afterok:${HLA1_extract_reads}

#echo $HLA3_run_HLAVBSeq

# run Polysolver

#HLA4_run_polysolver=`qsub -N HLA4_run_polysolver_${PID}_${SAMPLE} -l walltime=10:00:00,mem=20G,nodes=01:ppn=8 -j oe -o ${RESULTS_DIR} -v PATH_TO_BAM=${PATH_TO_BAM},RESULTS_DIR=${RESULTS_DIR} ${PIPELINE_DIR}/HLA4_run_polysolver.sh`

#echo $HLA4_run_polysolver

# run Optitype on extracted reads
#########################-W depend=afterok:${HLA1_all_reads}
HLA5_run_optitype=`qsub -N HLA5_run_optitype_${PID}_${SAMPLE} -l walltime=12:00:00,mem=30G,nodes=01:ppn=16 -j oe -o ${RESULTS_DIR} -v RESULTS_DIR=${RESULTS_DIR},PIPELINE_DIR=${PIPELINE_DIR},PATH_TO_BAM=${PATH_TO_BAM} ${PIPELINE_DIR}/all_HLA5_run_optitype.sh`

echo $HLA5_run_optitype

# run xHLA on extracted reads
#############################
#HLA6_run_xhla=`qsub -N HLA6_run_xhla_${PID}_${SAMPLE} -l walltime=2:00:00,mem=30G,nodes=01:ppn=16 -j oe -o ${RESULTS_DIR} -W depend=afterok:${HLA1_extract_reads} -v RESULTS_DIR=${RESULTS_DIR},SOFTWARE_DIR=${SOFTWARE_DIR},PATH_TO_BAM=${PATH_TO_BAM},NAME=${PID}_${SAMPLE} ${PIPELINE_DIR}/HLA6_run_xhla.sh`

#echo $HLA6_run_xhla

# run kourami on extracted reads
############################

#HLA7_run_kourami=`qsub -N HLA7_run_kourami_${PID}_${SAMPLE} -l walltime=2:00:00,mem=30G,nodes=01:ppn=16 -j oe -o ${RESULTS_DIR} -W depend=afterok:${HLA6_run_xhla} -v RESULTS_DIR=${RESULTS_DIR},SOFTWARE_DIR=${SOFTWARE_DIR},NAME=${PID}_${SAMPLE} ${PIPELINE_DIR}/HLA7_run_kourami.sh`

#echo $HLA7_run_kourami

