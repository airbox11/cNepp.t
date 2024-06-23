workDir=/icgc/dkfzlsdf/analysis/D120/yanhong/test_fusion
scriptsDir=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline/fusion_script
mkdir -p ${workDir}

sampleID=H021-RFEX31-metastasis07
testDir=deFuse



sampleID=libraryNam
testDir=test2


function step1 () {
    sh ${scriptsDir}/fusion/main_run_fusion.sh \
    -i ${sampleID} \
    -1 ${workDir}/AS-356663-LR-44266_R1.fastq.gz \
    -2 ${workDir}/AS-356663-LR-44266_R2.fastq.gz \
    -o ${workDir}/${testDir}
}


function step2 () {
    sh ${scriptsDir}/fusion/2_rerun_defuse.sh \
    -i ${sampleID} \
    -o ${workDir}/${testDir}
}


function step3_confuse () {
    python /icgc/dkfzlsdf/analysis/D120/scripts/fusion/conffuse/confFuse.py \
    ${workDir}/${testDir}/results.filtered.tsv \
    ${workDir}/${testDir}/confFuse.tsv \
    /icgc/dkfzlsdf/analysis/D120/scripts/fusion/conffuse/artefact_list.tab \
    /icgc/dkfzlsdf/analysis/D120/reference/gencode.v19.chr_patch_hapl_scaff.annotation.gtf
}

function step4_prediction () {
    script=/icgc/dkfzlsdf/analysis/D120/scripts/fusion/epitopePrediction/main_run_fusion_neoantigen_prediction.sh
    sh ${script} -f /icgc/dkfzlsdf/analysis/D120/scripts/fusion/epitopePrediction/test_fusion/fasta_myb_nfib.fa \
    -n NYB_NFIB \
    -m /icgc/dkfzlsdf/analysis/D120/scripts/fusion/epitopePrediction/test_fusion/format_hla_phlat_H021-RFEX31 \
    -o /icgc/dkfzlsdf/analysis/D120/yanhong/test_fusion
}

# step1
# step2
# step3_confuse
step4_prediction