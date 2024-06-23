#! /bin/bash

usage() 
{ echo "Usage: $0 -i PID 
    -i: sample ID 

Result folder is in /icgc/dkfzlsdf/analysis/D120/project/fusion_liposarcoma/results_per_pid

The generated file will be saved like:
/icgc/dkfzlsdf/analysis/D120/project/fusion_liposarcoma/results_per_pid/H021-1MRCLC/fastq/exon_sequencing/blood/phlat/execute_log.sh " 1>&2; exit 1;}

while getopts ":i:" opt; do
    case "${opt}" in
        i)
            PID=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${PID}" ] ; then
    usage
fi

#PID=$1

#controlType=(control buffy_coat blood buffy_coat02)
VIEW_PER_PID=/icgc/dkfzlsdf/project/hipo/hipo_02*/sequencing/*/view-by-pid
RESULTS_PER_PID=/icgc/dkfzlsdf/analysis/D120/project/fusion_liposarcoma/results_per_pid

phlatRelease=/icgc/dkfzlsdf/analysis/G200/immuno/tools/phlat-release
bowtie2=/odcf/cluster/13.1/x86_64/bin/bowtie2

module load python/2.7.9
module load R/3.4.0
module load perl/5.20.2
module load biobambam/0.0.191

for bamFile in `find $VIEW_PER_PID/$PID -name *.bam | grep -v ".merging\|chimeric\|results"`
do
	seqType=`echo $bamFile | awk '{split($1,a,"/");print a[8]}'`
	PID=`echo $bamFile | awk '{split($1,a,"/");print a[10]}'`
	sampleType=`echo $bamFile | awk '{split($1,a,"/");print a[11]}'`
	fqDir=$RESULTS_PER_PID/$PID/fastq/$seqType/$sampleType
	if [ ! -d "$fqDir/phlat" ];then
		mkdir -p $fqDir/phlat
	fi
	echo "module load python/2.7.9
module load R/3.4.0
module load perl/5.24.1
module load biobambam/0.0.191
bamtofastq collate=1  F=$fqDir/read1.fastq F2=$fqDir/read2.fastq filename=$bamFile ranges=6 S=$fqDir/single.fastq O=$fqDir/unmatch.1.fastq O2=$fqDir/unmatch.2.fastq
python -O /icgc/dkfzlsdf/analysis/G200/immuno/tools/phlat-release/dist/PHLAT.py -1 $fqDir/read1.fastq -2 $fqDir/read2.fastq  -index ${phlatRelease}/b2folder -b2url ${bowtie2}  -tag phlat_$PID  -e ${phlatRelease} -o $fqDir/phlat " &> $fqDir/phlat/execute_log.sh

done

