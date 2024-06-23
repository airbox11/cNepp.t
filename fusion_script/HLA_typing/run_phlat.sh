#! /bin/bash
##Author: Zhiqin Huang
##To run software PHLAT for HLA typing

# Predefine parameters in main scripts
#PID=$1
#OUTPUT_DIR=$2
#READ1=$3	# if multiple read1 files, then use ";" to combine them.
#READ2=$4
#threads=$5	# the number CPU to use

# Those are defined in main scripts, if not, define them here
#phlatRelease=/icgc/dkfzlsdf/analysis/G200/immuno/tools/phlat-release
#indexdir=$phlatRelease/b2folder
#bowtie2=/odcf/cluster/13.1/x86_64/bin/bowtie2 # update the path in main script if changed

if [ ! -f "${OUTPUT_DIR}/phlat_${PID}_HLA.sum" ];then
	#### create a named pipe for decompress fastq files
	if [ -e ${OUTPUT_DIR}/tmp_fifo_READ1 ];then
		rm ${OUTPUT_DIR}/tmp_fifo_READ1
	fi

	if [ -e ${OUTPUT_DIR}/tmp_fifo_READ2 ];then
		rm ${OUTPUT_DIR}/tmp_fifo_READ2
	fi

	mkfifo ${OUTPUT_DIR}/tmp_fifo_READ1
	mkfifo ${OUTPUT_DIR}/tmp_fifo_READ2
	
	cREAD1=`echo "${READ1}" | awk '{gsub(";"," ",$0);print $0}'`
	cREAD2=`echo "${READ2}" | awk '{gsub(";"," ",$0);print $0}'`

	zcat ${cREAD1} >  ${OUTPUT_DIR}/tmp_fifo_READ1 &
	zcat ${cREAD2} >  ${OUTPUT_DIR}/tmp_fifo_READ2 &

	#### to run software PHLAT for hla typing
	touch ${OUTPUT_DIR}/running_phlat

	echo "python -O $phlatRelease/dist/PHLAT.py -1 ${OUTPUT_DIR}/tmp_fifo_READ1  -2 ${OUTPUT_DIR}/tmp_fifo_READ2 -index $phlatRelease/b2folder -b2url ${bowtie2} -tag "phlat_$PID" -e $phlatRelease -o ${OUTPUT_DIR} -p $threads"
	python -O $phlatRelease/dist/PHLAT.py -1 ${OUTPUT_DIR}/tmp_fifo_READ1  -2 ${OUTPUT_DIR}/tmp_fifo_READ2 -index $phlatRelease/b2folder -b2url ${bowtie2} -tag "phlat_$PID" -e $phlatRelease -o ${OUTPUT_DIR} -p ${threads}

	#### remove temporary files
	rm ${OUTPUT_DIR}/running_phlat
	rm ${OUTPUT_DIR}/tmp_fifo_READ1
	rm ${OUTPUT_DIR}/tmp_fifo_READ2
else
	echo "${OUTPUT_DIR}/phlat_$PID_HLA.sum is already available, running of Phlat stops"
fi

