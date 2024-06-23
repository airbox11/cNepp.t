#! /bin/bash


# usage message for users
usage() { 
    echo "Usage: $0 -i PID -1 read1 -2 read2 -o outputDir
    -1: Path to read1 file, if multiple fastq files, then quote symbol is needed. e.g. \"fq1_1;fq2_1;fq3_1\"
    -2: Path to read2 file
    -i: Sample ID
    -o: Path to result folder
>>> Note: Absolute path needed" 1>&2; exit 1; 
}

while getopts ":i:1:2:d:o:f" opt; do
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
        d)
            s2Dir=${OPTARG}
            ;;
        o)
            OUTPUT_DIR=${OPTARG}
            ;;
        f)
            defuse=${OPTARG}
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


deFuseOutputDir=$OUTPUT_DIR/$s2Dir
mkdir -p $deFuseOutputDir

## unzip ====
if [ ! -f $OUTPUT_DIR/reads.1.fastq ];then
    echo 'unzip fastq.gz'
    zcat ${read1} > $OUTPUT_DIR/reads.1.fastq 
fi
if [ ! -f $OUTPUT_DIR/reads.2.fastq ];then
    echo 'unzip fastq.gz'
    zcat ${read2} > $OUTPUT_DIR/reads.2.fastq 
fi


if [ $defuse == 'yanhong' ];then
    echo $defuse

    config=/home/lyuya/miniconda3/opt/defuse/scripts/config.txt
    dataset_directory=/icgc/dkfzlsdf/analysis/D120/yanhong/git/deFuse_reference

    script=/home/lyuya/miniconda3/opt/defuse/scripts/defuse_run.pl
    $script -c $config -d $dataset_directory -1 $OUTPUT_DIR/reads.1.fastq -2 $OUTPUT_DIR/reads.2.fastq -o $deFuseOutputDir -p 8 -n libraryName

fi


## run zhiqin's defuse====
if [ $defuse == 'zh' ];then
    script=/lsdf/dkfz/analysis/D120/tools/deFuse/defuse-0.6.1/scripts/defuse.pl
    config=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline/fusion_script/fusion/config.txt
    perl $script -c $config -1 $OUTPUT_DIR/reads.1.fastq -2 $OUTPUT_DIR/reads.2.fastq -o $deFuseOutputDir -p 8 -n libraryName
fi
# deFuse_sh=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline/fusion_script/fusion/run_defuse.sh
# ${deFuse_sh} ${READ1} ${READ2} ${OUTPUT_DIR}
# echo "defuse_run.pl -c $config -1 $OUTPUT_DIR/reads.1.fastq -2 $OUTPUT_DIR/reads.2.fastq -o $deFuseOutputDir -p 8"

## run fusion tool deFuse
# bsub -R "rusage[mem=200G]" -u lvyhwind@hotmail.com -r -o log_${PID}_defuse -J ${PID}_deFuse -W 200:00 -n $threads -env "all, threads=$threads,LeftReads=${READ1},RightReads=${READ2},smID_path=${OUTPUT_DIR},libraryName=${PID}" ${deFuse_sh}


## run yanhong's deFuse ====


## end of the run 
echo ">>>> running defuse : done"

