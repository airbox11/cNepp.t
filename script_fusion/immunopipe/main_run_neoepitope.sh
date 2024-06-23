#! /bin/bash

# to submit this script: run_SNV_mpileup_pipeline.sh [-v] -c <CONFIG_FILE> PID1 PID2 PID3

#processing arguments
usage() { echo "Usage: $0 -c CONFIG_FILE -i PID
        >>>> Note: Absolute path of config file" 1>&2; exit 1; }

while getopts ":c:i:" o; do
    case "${o}" in
        c)
            CONFIG_FILE=${OPTARG}
            ;;
        i)
            PID=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${CONFIG_FILE}" ] || [ ! -f "${CONFIG_FILE}" ] || [ -z "${PID}" ]; then
    usage
fi

source ${CONFIG_FILE}

# identifiy path of raw fastq files of control samples
for i in "${controlType[@]}"
do
        controlFolder=$VIEW_PER_PID/${PID}/$i/
        if [ -d $controlFolder ];then
                READ1=`find $controlFolder -name *.fastq.gz | sort | sed '1q;d'`
                READ2=`find $controlFolder -name *.fastq.gz | sort | sed '2q;d'`
        fi
done
echo ">>>>>Input control Read1 is: " $READ1
echo ">>>>>Input control Read2 is: " $READ2

if [ ! -f $READ1 ];then
	echo $READ1 " is not found. This run stops. Pleas recheck it."
	exit
fi

# identify the path of SNV file
VCF_ORIGINAL=${RESULTS_PER_PID}/${PID}/mpileup_indels/snvs_${PID}_somatic_functional_snvs_conf_7_to_10.vcf

if [ ! -f ${VCF_ORIGINAL} ];then
	echo "${VCF_ORIGINAL} is not available. Running stops and exits";
	exit 1
fi

# create a new folder 
OUTPUT_DIR=${RESULTS_PER_PID}/${PID}/immuno
if [ ! -d "${OUTPUT_DIR}" ];then
	mkdir -p ${OUTPUT_DIR}
fi
cd ${OUTPUT_DIR}

# make a soft link of mutation calling VCF
if [ ! -f "snvs_${PID}_somatic_functional_snvs_conf_7_to_10.vcf" ];then
	ln -s $VCF_ORIGINAL snvs_${PID}_somatic_functional_snvs_conf_7_to_10.vcf
	#continue
fi
VCF=${OUTPUT_DIR}/snvs_${PID}_somatic_functional_snvs_conf_7_to_10.vcf

# running Phlat
run_phlat=`qsub -N phlat_${PID} -o ${OUTPUT_DIR}/log_run_phlat  -l walltime=15:00:00,mem=5g,nodes=1:ppn=8 -v OUTPUT_DIR=${OUTPUT_DIR},PID=${PID},READ1=${READ1},READ2=${READ2},phlatRelease=${phlatRelease} ${PIPELINE_DIR}/run_phlat.sh`

# predict mutated protein and reference protein
run_predict_protein=`qsub -N predictProtein_${PID} -o ${OUTPUT_DIR}/log_predictProtein -l walltime=1:00:00,mem=1g,nodes=1:ppn=1 -v ANNOVAR=${ANNOVAR},VCF=${VCF},predictedProtein=${OUTPUT_DIR}/${PID}_predictedProtein.fa,PIPELINE_DIR=${PIPELINE_DIR} ${PIPELINE_DIR}/run_predict_protein_fa.sh`

# extract mutated and reference peptides from predicted protein fasta (e.g. 29aa)
run_separate_ref_mut_peptide=`qsub -W depend=afterok:${run_predict_protein} -N separate_ref_mut_peptide_${PID} -o ${OUTPUT_DIR}/log_separatePeptide -l walltime=1:00:00,mem=1g,nodes=1:ppn=1 -v fastaFile=${OUTPUT_DIR}/${PID}_predictedProtein.fa,outputMut=${OUTPUT_DIR}/${PID}_mut.fa,outputRef=${OUTPUT_DIR}/${PID}_ref.fa,extendLen=${extendLen},PIPELINE_DIR=${PIPELINE_DIR} ${PIPELINE_DIR}/run_separate_ref_mut_peptide.sh`

# run netMHC for mutated peptides
run_netMHC_mut=`qsub -W depend=afterok:${run_phlat}:${run_separate_ref_mut_peptide} -N run_netMHC_mutated_${PID} -o ${OUTPUT_DIR}/log_netMHC_mutated -l walltime=10:00:00,mem=4g,nodes=1:ppn=8 -v phlatResultFile=${OUTPUT_DIR}/phlat_${PID}_HLA.sum,peptideFasta=${OUTPUT_DIR}/${PID}_mut.fa,OUTPUT_DIR=${OUTPUT_DIR},netMHCpan_SIF=${netMHCpan_SIF},netMHCIIpan_SIF=${netMHCIIpan_SIF},pepType=mut $PIPELINE_DIR/run_netMHC.sh`

# run netMHC for reference peptides
run_netMHC_ref=`qsub -W depend=afterok:${run_phlat}:${run_separate_ref_mut_peptide} -N run_netMHC_ref_${PID} -o ${OUTPUT_DIR}/log_netMHC_ref -l walltime=10:00:00,mem=4g,nodes=1:ppn=8 -v phlatResultFile=${OUTPUT_DIR}/phlat_${PID}_HLA.sum,peptideFasta=${OUTPUT_DIR}/${PID}_ref.fa,OUTPUT_DIR=${OUTPUT_DIR},netMHCpan_SIF=${netMHCpan_SIF},netMHCIIpan_SIF=${netMHCIIpan_SIF},pepType=ref $PIPELINE_DIR/run_netMHC.sh `

# Summarize results for MHC-I and MHC-II neoepitopes
run_summary=`qsub -W depend=afterok:${run_netMHC_ref}:${run_netMHC_mut} -o ${OUTPUT_DIR}/log_summary -N run_summary_results_${PID} -l walltime=1:00:00,mem=1g,nodes=1:ppn=1 -v OUTPUT_DIR=${OUTPUT_DIR},VCF=${VCF},PID=${PID},PIPELINE_DIR=${PIPELINE_DIR},extendLen=${extendLen}  ${PIPELINE_DIR}/run_summary.sh`

			
