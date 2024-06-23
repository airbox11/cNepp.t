#! /bin/bash

## This pipeline script aims to predict neoantigen based on VCF and HLA-type files

usage() { echo "Usage: $0 -c CONFIG_FILE -i PID -v vcf -m hla_type -o outputDir
	-c: Path to config file
	-i: Sample ID
	-v: Path to VCF file
	-m: Path to HLA type result file
	-o: Path to result folder
>>> Note: Absolute path needed" 1>&2; exit 1; }

while getopts ":c:i:v:m:o:" o; do
    case "${o}" in
        c)
            CONFIG_FILE=${OPTARG}
            ;;
        i)
            PID=${OPTARG}
            ;;
		v)
			VCF=${OPTARG}
			;;
		m)
			HLA=${OPTARG}
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

if [ -z "${CONFIG_FILE}" ] || [ ! -f "${CONFIG_FILE}" ] || [ ! -f "${VCF}" ] || [ ! -f "${HLA}" ]; then
    usage
fi

source ${CONFIG_FILE}

# create a new folder 
if [ ! -d "${OUTPUT_DIR}" ];then
	mkdir -p ${OUTPUT_DIR}
fi
cd ${OUTPUT_DIR}

module load R/3.4.0
module load python/2.7.9
module load perl/5.20.2

# predict mutated protein and reference protein
echo "run peptide prediction"
run_predict_protein=`bsub -r -J predictProtein_${PID} -o ${OUTPUT_DIR}/log_predictProtein -W 01:00 -R "rusage[mem=2G]" -n 1 -env "all, ANNOVAR=${ANNOVAR},VCF=${VCF},predictedProtein=${OUTPUT_DIR}/${PID}_predictedProtein.fa,PIPELINE_DIR=${PIPELINE_DIR},OUTPUT_DIR=${OUTPUT_DIR}" ${PIPELINE_DIR}/run_predict_protein_fa.sh | cut -d'<' -f2 | cut -d'>' -f1 `

# extract mutated and reference peptides from predicted protein fasta (e.g. 29aa)
echo "run seperating mut/ref peptides"
run_separate_ref_mut_peptide=`bsub -r -w "done(${run_predict_protein})" -J separate_ref_mut_peptide_${PID} -o ${OUTPUT_DIR}/log_separatePeptide -W 01:00 -R "rusage[mem=2G]" -n 1 -env "all, fastaFile=${OUTPUT_DIR}/${PID}_predictedProtein.fa,outputMut=${OUTPUT_DIR}/${PID}_mut.fa,outputRef=${OUTPUT_DIR}/${PID}_ref.fa,extendLen=${extendLen},PIPELINE_DIR=${PIPELINE_DIR}" ${PIPELINE_DIR}/run_separate_ref_mut_peptide.sh | cut -d'<' -f2 | cut -d'>' -f1 `

echo "run netMHC"
# run netMHC for mutated peptides
run_netMHC_mut=`bsub -r -w "done(${run_separate_ref_mut_peptide})" -J run_netMHC_mutated_${PID} -o ${OUTPUT_DIR}/log_netMHC_mutated -W 20:00 -R "rusage[mem=4G]" -n 8 -env "all, HLA=${HLA},peptideFasta=${OUTPUT_DIR}/${PID}_mut.fa,OUTPUT_DIR=${OUTPUT_DIR},netMHCpan=${netMHCpan},netMHCIIpan=${netMHCIIpan},pepType=mut" $PIPELINE_DIR/run_netMHC.sh | cut -d'<' -f2 | cut -d'>' -f1 `

# run netMHC for reference peptides
run_netMHC_ref=`bsub -r -w "done(${run_separate_ref_mut_peptide})" -J run_netMHC_ref_${PID} -o ${OUTPUT_DIR}/log_netMHC_ref -W 20:00 -R "rusage[mem=4G]" -env "all, HLA=${HLA},peptideFasta=${OUTPUT_DIR}/${PID}_ref.fa,OUTPUT_DIR=${OUTPUT_DIR},netMHCpan=${netMHCpan},netMHCIIpan=${netMHCIIpan},pepType=ref" $PIPELINE_DIR/run_netMHC.sh | cut -d'<' -f2 | cut -d'>' -f1`

echo "run filtering"
# Summarize results for MHC-I and MHC-II neoepitopes
# to update run_summary.sh because "key" in python scripts not consistant
####  Use map_NM_geneID file to solve the problem
####
run_summary=`bsub -r -w "done(${run_netMHC_ref}) && done(${run_netMHC_mut})" -o ${OUTPUT_DIR}/log_summary -J run_summary_results_${PID} -W 01:00 -R "rusage[mem=2G]" -n 1  -env "all, OUTPUT_DIR=${OUTPUT_DIR},VCF=${VCF},PID=${PID},PIPELINE_DIR=${PIPELINE_DIR},extendLen=${extendLen}"  ${PIPELINE_DIR}/run_summary.sh | cut -d'<' -f2 | cut -d'>' -f1`

