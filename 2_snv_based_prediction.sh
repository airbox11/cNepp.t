# code=/icgc/dkfzlsdf/analysis/D120/scripts/neoepitope_SNVs/main_run_neoepitope_snvs.sh -c 
# CONFIG_FILE -i PID -v vcf -m hla_type -o outputDir
# -c: Path to config file
# -i: Sample ID
# -v: Path to VCF file
# -m: Path to HLA type result file
# -o: Path to result folder



sampleID=$1
patientID=$2
workDir=$3
mhcIversion=$4


if [[ $mhcIversion == 4.0 ]]
then
    script=/icgc/dkfzlsdf/analysis/D120/scripts/neoepitope_SNVs/main_run_neoepitope_snvs.sh
    config=/icgc/dkfzlsdf/analysis/D120/scripts/neoepitope_SNVs/config.sh
elif [[ $mhcIversion == 4.1 ]]
then
    script=/icgc/dkfzlsdf/analysis/D120/yanhong/update_netMHCpan/netMHCpan_4.1/main_run_neoepitope_snvs.sh
    config=/icgc/dkfzlsdf/analysis/D120/yanhong/update_netMHCpan/netMHCpan_4.1/config.sh
fi


if [ ! -e $workDir/2_SNVs_based_neoepitope_prediction/*vcf ]; then
    vcf1=/icgc/dkfzlsdf/project/OE0422/immuno_patients_nct/sequencing/exon_sequencing/view-by-pid/IPNCT_${patientID}/snv_results/paired/*/results*/*somatic_functional_snvs_conf_8_to_10.vcf
    vcf2=${workDir}/2_SNVs_based_neoepitope_prediction/snvs_IPNCT_${patientID}_somatic_functional_snvs_conf_8_to_10.vcf
    ln -sf $vcf1 $vcf2
    if [ ! -e $vcf2 ]; then
        echo "file not exist!! >>>>>>"
        echo $vcf2
        exit 1
    fi

fi


mkdir -p ${workDir}/2_SNVs_based_neoepitope_prediction/result


sh $script \
-c $config \
-i ${sampleID} \
-v ${workDir}/2_SNVs_based_neoepitope_prediction/*.vcf \
-m ${workDir}/1.2_format_hla_result/format_hla_phlat \
-o ${workDir}/2_SNVs_based_neoepitope_prediction/result


# -v ${workDir}/2_SNVs_based_neoepitope_prediction/snvs_IPNCT_${patientID}_somatic_functional_snvs_conf_8_to_10.vcf \