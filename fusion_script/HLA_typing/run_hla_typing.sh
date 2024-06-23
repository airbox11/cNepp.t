#! /bin/bash

usage() 
{ echo "Usage: $0 -i PID -1 READ1 -2 READ2 -o OUTPUT_DIR
	-i: sample ID
	-1: reads 1, if multiple fastq files, then quote symbol is needed. e.g. \"fq1_1;fq2_1;fq3_1\" 
	-2: reads 2, like read 1.
	-o: Output directory (will be created if not available)
>>> Note: Absolute path needed " 1>&2; exit 1;}

while getopts ":i:r1:r2:o:" opt; do
    case "${opt}" in
        i)
            PID=${OPTARG}
            ;;
        1)
            READ1=${OPTARG}
            ;;
        2)
            READ2=${OPTARG}
            ;;
        o)
            OUTPUT_DIR=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${PID}" ] || [ -z "${READ1}" ] || [ -z "${READ2}" ] || [ -z "${OUTPUT_DIR}" ] ; then
	usage
fi

# Parameters, check the following path before you use it. Make sure they are available.
phlatRelease=/icgc/dkfzlsdf/analysis/G200/immuno/tools/phlat-release
#bowtie2=/odcf/cluster/13.1/x86_64/bin/bowtie2 # ODCF changed the path in 2021
bowtie2=/software/bowtie2/2.3.5.1/bowtie2 
PIPELINE_DIR=/icgc/dkfzlsdf/analysis/D120/scripts/HLA_typing
threads=8 # the number CPU to use, setting as 8


###########
# create folder for PHLAT output
phlatDir=${OUTPUT_DIR}/phlat
if [ ! -d $phlatDir ];then
	mkdir -p $phlatDir
	echo "mkdir -p $phlatDir"
fi

module load perl/5.24.1
module load R/3.4.0
module load python/2.7.9
module load samtools/1.5

bsub -r -J phlat_${PID} -o ${OUTPUT_DIR}/log_run_phlat  -W 250:00 -R "rusage[mem=10G]"  -n $threads -env "all, OUTPUT_DIR=${phlatDir},PID=${PID},READ1=${READ1},READ2=${READ2},phlatRelease=${phlatRelease},bowtie2=${bowtie2},threads=${threads}" ${PIPELINE_DIR}/run_phlat.sh

#echo "bsub -r -J phlat_${PID} -o ${OUTPUT_DIR}/log_run_phlat  -W 250:00 -R \"rusage[mem=10G]\"  -n $threads -env \"all, OUTPUT_DIR=${phlatDir},PID=${PID},READ1=${READ1},READ2=${READ2},phlatRelease=${phlatRelease},bowtie2=${bowtie2},threads=${threads}\" ${PIPELINE_DIR}/run_phlat.sh"

