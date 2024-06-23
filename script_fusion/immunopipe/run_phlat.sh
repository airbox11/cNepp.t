#! /bin/bash
##Author: Zhiqin Huang

#OUTPUT_DIR=$1
#READ1=$2
#READ2=$3
#PID=$4

#phlatPY=/icgc/dkfzlsdf/analysis/D120/kosalogl/tools/phlat-release/dist/PHLAT.py
#phlatRelease=/icgc/dkfzlsdf/analysis/D120/kosalogl/tools/phlat-release
#indexdir=$phlatRelease/b2folder


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

	zcat ${READ1} >  ${OUTPUT_DIR}/tmp_fifo_READ1 &
	zcat ${READ2} >  ${OUTPUT_DIR}/tmp_fifo_READ2 &

	#### to run software PHLAT for hla typing
	touch ${OUTPUT_DIR}/running_phlat

	echo "python -O $phlatRelease/dist/PHLAT.py -1 ${OUTPUT_DIR}/tmp_fifo_READ1  -2 ${OUTPUT_DIR}/tmp_fifo_READ2 -index $phlatRelease/b2folder -b2url bowtie2 -tag "phlat_$PID" -e $phlatRelease -o ${OUTPUT_DIR}"
	python -O $phlatRelease/dist/PHLAT.py -1 ${OUTPUT_DIR}/tmp_fifo_READ1  -2 ${OUTPUT_DIR}/tmp_fifo_READ2 -index $phlatRelease/b2folder -b2url bowtie2 -tag "phlat_$PID" -e $phlatRelease -o ${OUTPUT_DIR}


	#### remove temporary files
	rm ${OUTPUT_DIR}/running_phlat
	rm ${OUTPUT_DIR}/tmp_fifo_READ1
	rm ${OUTPUT_DIR}/tmp_fifo_READ2
else
	echo "${OUTPUT_DIR}/phlat_$PID_HLA.sum is already available, running of Phlat stops"
fi

