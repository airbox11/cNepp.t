#! /bin/bash
# update on 2019-05
# the deFuse scritps folder
deFuseScriptsFoder=/lsdf/dkfz/analysis/D120/tools/deFuse/defuse-0.6.1/scripts
config=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline/fusion_script/fusion/config.txt

# leftRead=$LeftReads
# rightRead=$RightReads
# outputDir=${smID_path}

read1=$1
read2=$2
outputDir=$3

threads=8

# set path
deFuseOutputDir=$outputDir/deFuse
# create deFuse output folder is not exsited.
if [ ! -d $deFuseOutputDir ];then
    mkdir $deFuseOutputDir
fi

    # zcat ${read1} |head 


if [ ! -f $deFuseOutputDir/reads.1.fastq ];then
    zcat ${read1} > $deFuseOutputDir/reads.1.fastq 
fi
if [ ! -f $deFuseOutputDir/reads.2.fastq ];then
    zcat ${read2} > $deFuseOutputDir/reads.2.fastq 
fi

# deFuseFile=$deFuseOutputDir/results.filtered.tsv
# if [ ! -f "$deFuseFile" ];then
    # set buffer path for FIFO files 
    #cd $deFuseOutputDir
    #makes a FIFO special file

        # mkfifo $deFuseOutputDir/tmp_fifo_file1
        # mkfifo $deFuseOutputDir/tmp_fifo_file2

    #r1=`echo $leftRead | sed 's\%%\ \g'`
    #r2=`echo $rightRead | sed 's\%%\ \g'`
    #zcat $r1 > $deFuseOutputDir/tmp_fifo_file1 & # & , 2016-02-15
    #zcat $r2 > $deFuseOutputDir/tmp_fifo_file2 & # & , 2016-02-15

    # 2018-08. Merge multiple fastq files

        # cREAD1=`echo "${leftRead}" | awk '{gsub(";"," ",$0);print $0}'`
        # cREAD2=`echo "${rightRead}" | awk '{gsub(";"," ",$0);print $0}'`

        # echo "zcat ${cREAD1} >  ${phlatOutDir}/tmp_fifo_READ1 &"
        # echo "zcat ${cREAD2} >  ${phlatOutDir}/tmp_fifo_READ2 &"



    # touch $deFuseOutputDir/running_defuse

    # echo "perl $deFuseScriptsFoder/defuse.pl -c $deFuseScriptsFoder/config.txt -1 $deFuseOutputDir/tmp_fifo_file1 -2 $deFuseOutputDir/tmp_fifo_file2 -o $deFuseOutputDir -p $threads -n ${libraryName}"
## my code

# perl $deFuseScriptsFoder/defuse.pl -c $deFuseScriptsFoder/config.txt -1 $deFuseOutputDir/reads.1.fastq -2 $deFuseOutputDir/reads.2.fastq -o $deFuseOutputDir -p $threads -n libraryName
perl $deFuseScriptsFoder/defuse.pl -c $config -1 $deFuseOutputDir/reads.1.fastq -2 $deFuseOutputDir/reads.2.fastq -o $deFuseOutputDir -p $threads -n libraryName

    # rm $deFuseOutputDir/tmp_fifo_file1
    # rm $deFuseOutputDir/tmp_fifo_file2
    # rm $deFuseOutputDir/running_defuse
# else
#     echo "$deFuseFile is already there! Please check!"
#     exit 0
# fi

