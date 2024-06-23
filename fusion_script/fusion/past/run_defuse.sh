#! /bin/bash
# update on 2019-05
# the deFuse scritps folder
deFuseScriptsFoder=/lsdf/dkfz/analysis/D120/tools/deFuse/defuse-0.6.1/scripts

leftRead=$LeftReads
rightRead=$RightReads
outputDir=${smID_path}

threads=$threads

# set path
deFuseOutputDir=$outputDir/deFuse
# create deFuse output folder is not exsited.
if [ ! -d $deFuseOutputDir ];then
    mkdir $deFuseOutputDir
fi

deFuseFile=$deFuseOutputDir/results.filtered.tsv
if [ ! -f "$deFuseFile" ];then
    # set buffer path for FIFO files 
    #cd $deFuseOutputDir
    #makes a FIFO special file
    mkfifo $deFuseOutputDir/tmp_fifo_file1
    mkfifo $deFuseOutputDir/tmp_fifo_file2
    #r1=`echo $leftRead | sed 's\%%\ \g'`
    #r2=`echo $rightRead | sed 's\%%\ \g'`
    #zcat $r1 > $deFuseOutputDir/tmp_fifo_file1 & # & , 2016-02-15
    #zcat $r2 > $deFuseOutputDir/tmp_fifo_file2 & # & , 2016-02-15

    # 2018-08. Merge multiple fastq files
    cREAD1=`echo "${leftRead}" | awk '{gsub(";"," ",$0);print $0}'`
    cREAD2=`echo "${rightRead}" | awk '{gsub(";"," ",$0);print $0}'`

    echo "zcat ${cREAD1} >  ${phlatOutDir}/tmp_fifo_READ1 &"
    echo "zcat ${cREAD2} >  ${phlatOutDir}/tmp_fifo_READ2 &"

    zcat ${cREAD1} > $deFuseOutputDir/tmp_fifo_file1 & # & , 2016-02-15
    zcat ${cREAD2} > $deFuseOutputDir/tmp_fifo_file2 & # & , 2016-02-15

    touch $deFuseOutputDir/running_defuse

    echo "perl $deFuseScriptsFoder/defuse.pl -c $deFuseScriptsFoder/config.txt -1 $deFuseOutputDir/tmp_fifo_file1 -2 $deFuseOutputDir/tmp_fifo_file2 -o $deFuseOutputDir -p $threads -n ${libraryName}"
    perl $deFuseScriptsFoder/defuse.pl -c $deFuseScriptsFoder/config.txt -1 $deFuseOutputDir/tmp_fifo_file1 -2 $deFuseOutputDir/tmp_fifo_file2 -o $deFuseOutputDir -p $threads -n ${libraryName}

    rm $deFuseOutputDir/tmp_fifo_file1
    rm $deFuseOutputDir/tmp_fifo_file2
    rm $deFuseOutputDir/running_defuse
else
    echo "$deFuseFile is already there! Please check!"
    exit 0
fi


if [ -f "$deFuseFile" ];then
    touch $deFuseOutputDir/finish_defuse
fi
    
