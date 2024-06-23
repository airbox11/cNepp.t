set -e 
PATH=~/miniconda3/bin:$PATH
export runID=$1
export s2Dir=$2
export steps=$3
export defuse=$4

module load perl/5.24.1
module load R/3.4.0
module load python/3.6.1
module load samtools/1.5
module load blat/3.6

## parameter sets
workDir=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline_collection/mhc4.1/${runID}/9_Fusion_gene_based_neoepitope_identification
scriptsDir=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline/fusion_script



function s0_clearData () {
    rm -rf ${workDir}
}

function s1_prepareDir () {
    mkdir -p ${workDir}
    mkdir -p ${workDir}/1_prepare_fastq
    mkdir -p ${workDir}/2_mainRun
    mkdir -p ${workDir}/3_rerun
    mkdir -p ${workDir}/4_confuse
    mkdir -p ${workDir}/5_prediction

    fastqDir=/icgc/dkfzlsdf/analysis/D120/yanhong/data_collection_fastq.gz
    cp -P ${fastqDir}/*${runID}* ${workDir}/1_prepare_fastq
}

function s2_mainRun () {
    OUTPUT_DIR=$workDir/2_mainRun
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
        script=/home/lyuya/miniconda3/opt/defuse/scripts/defuse_run.pl
        config=/home/lyuya/miniconda3/opt/defuse/scripts/config.txt
        dataset_directory=/icgc/dkfzlsdf/analysis/D120/yanhong/git/deFuse_reference
        
        $script -c $config -d $dataset_directory -1 $OUTPUT_DIR/reads.1.fastq -2 $OUTPUT_DIR/reads.2.fastq -o $deFuseOutputDir -p 8 -n libraryName
    fi

    if [ $defuse == 'zh' ];then
        script=/lsdf/dkfz/analysis/D120/tools/deFuse/defuse-0.6.1/scripts/defuse.pl
        config=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline/fusion_script/fusion/config.txt
        echo "perl $script -c $config -1 $OUTUPT_DIR/reads.1.fastq -2 $OUTPUT_DIR/reads.2.fastq -o $deFuseOutputDir -p 8 -n libraryName"
        perl $script -c $config -1 $OUTPUT_DIR/reads.1.fastq -2 $OUTPUT_DIR/reads.2.fastq -o $deFuseOutputDir -p 8 -n libraryName
    fi
}


function s3_reRun () {
    inputDir=${workDir}/2_mainRun/$s2Dir 
    outputDir=${workDir}/3_rerun/$s2Dir
    sh ${scriptsDir}/fusion/2_rerun_defuse.sh \
    -i ${runID} \
    -o $inputDir
    # -o ${workDir}/2_mainRun/deFuse
    mkdir -p $outputDir
    mv $inputDir/results* $inputDir/annotations $inputDir/log $outputDir
}


function s4_confuse () {
    python /icgc/dkfzlsdf/analysis/D120/scripts/fusion/conffuse/confFuse.py \
    ${workDir}/3_rerun/results.filtered.tsv \
    ${workDir}/4_confuse/confFuse.tsv \
    /icgc/dkfzlsdf/analysis/D120/scripts/fusion/conffuse/artefact_list.tab \
    /icgc/dkfzlsdf/analysis/D120/reference/gencode.v19.chr_patch_hapl_scaff.annotation.gtf
}

function s5_prediction () {
    script=/icgc/dkfzlsdf/analysis/D120/scripts/fusion/epitopePrediction/main_run_fusion_neoantigen_prediction.sh
    sh ${script} -f /icgc/dkfzlsdf/analysis/D120/scripts/fusion/epitopePrediction/test_fusion/fasta_myb_nfib.fa \
    -n ${runID} \
    -m /icgc/dkfzlsdf/analysis/D120/scripts/fusion/epitopePrediction/test_fusion/format_hla_phlat_H021-RFEX31 \
    -o ${workDir}/5_prediction
}

if echo $steps | grep -o s0; then
    s0_clearData
fi

if echo $steps | grep -o s1; then
    s1_prepareDir
fi
if echo $steps | grep -o s2; then
    s2_mainRun
fi
if echo $steps | grep -o s3; then
    s3_reRun
fi
if echo $steps | grep -o s4; then
    s4_confuse
fi
if echo $steps | grep -o s5; then
    s5_prediction
fi

