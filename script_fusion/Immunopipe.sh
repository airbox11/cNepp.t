#!/bin/bash

set -e 
########################################################
## Immunopipe runscript                               ##                                     
########################################################

function print_usage {
    echo '
    Usage: ./Immunopipe.sh  
		--dna-bam /path/to/bam [ 
		--snvs /path/to/snv.vcf] [ 
		--indels /path/to/indel.vcf] [ 
		--fusions /path/to/arriba_fusions.tsv] [ 
		--fusion-confidence "medium high"]  
		--output /path/to/output/directory  
		--pid PID
    ' 1>&2
    exit 1
}

# parse command-line arguments
FUSION_CONFIDENCE="medium high"
while [ $# -gt 0 ]; do
	PARAMETER="$1"
	shift
	case "$PARAMETER" in
		--dna-bam) PATH_TO_BAM="$1";;
		--hla) HLA="$1";;
		--snvs)
			SNV="$1"
			if [ ! -f "${SNV}" ];then
				echo "SNV-VCF file is not available at ${SNV}. Running stops and exits";
				exit 1
			fi
		;;
		--indels)
			INDEL="$1"
			if [ ! -f "${INDEL}" ];then
				echo "INDEL-VCF file is not available at ${INDEL}. Running stops and exits";
				exit 1
			fi
		;;
		--fusions)
			FUSION="$1"
			if [ ! -f "${FUSION}" ];then
				echo "Arriba fusion file is not available at ${FUSION}. Running stops and exits";
				exit 1
			fi
		;;
		--fusion-confidence) FUSION_CONFIDENCE="$1";;
		--output) PATH_TO_OUT_DIR="$1";;
		--OUT_FUSION) OUT_FUSION="$1";;
		--pid) PID="$1";;
		*) echo "Unknown parameter: $PARAMETER" 1>&2; print_usage;;
	esac
	shift
done
if [ -z "$PATH_TO_BAM" -o -z "$PATH_TO_OUT_DIR" -o -z "$PID" ]; then
	print_usage
fi

PIPELINE_DIR=$(readlink -f "$(dirname "$0")")
	# /omics/groups/OE0422/internal/yanhong/git/neoepitope_prediction-master1

source ${PIPELINE_DIR}/immunopipe/config.sh
### Yanhong code: =====
# auto-detect installation directory



# netMHCpan=/icgc/dkfzlsdf/analysis/D120/yanhong/update_netMHCpan/netMHCpan_4.1/netMHCpan
# netMHCIIpan=/icgc/dkfzlsdf/analysis/D120/yanhong/update_netMHCpan/netMHCIIpan-4.0/netMHCIIpan

netMHCpan=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/update_netMHCpan/netMHCpan_4.1/netMHCpan
netMHCIIpan=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/update_netMHCpan/netMHCIIpan-4.0/netMHCIIpan



########################### HLA-Typing ##################################
if [ -n "$HLA" ]; then
	RESULTS_DIR=${PATH_TO_OUT_DIR}
	# rm -rf $RESULTS_DIR
	mkdir -p $RESULTS_DIR

	## Extract reads mapped to HLA loci and unmapped reads and convert to fastq files ====
	HLA1_extract_reads=`bsub -J HLA1_extract_reads_${PID} -W 24:00 -n 2 -R "rusage[mem=21504]" -o ${RESULTS_DIR}/stdout_extract_reads -env "all,PATH_TO_BAM=${PATH_TO_BAM},NAME=${PID},RESULTS_DIR=${RESULTS_DIR}" < ${PIPELINE_DIR}/HLA_caller_scripts/HLA1_extract_reads.sh`

	## run Polysolver ====
	HLA4_run_polysolver=`bsub -J HLA4_run_polysolver_${PID} -W 24:00 -n 8 -R "rusage[mem=21504]" -o ${RESULTS_DIR}/stdout_run_polysolver -w "done(HLA1_extract_reads_${PID})" -env "all,PATH_TO_BAM=${PATH_TO_BAM},RESULTS_DIR=${RESULTS_DIR}" < ${PIPELINE_DIR}/HLA_caller_scripts/HLA4_run_polysolver.sh`


	## run Optitype on extracted reads ====
	rm -rf ${RESULTS_DIR}/stdout_run_optitype

	HLA5_run_optitype=`bsub -J HLA5_run_optitype_${PID} -W 240:00 -n 16 -R "rusage[mem=200G]" -o ${RESULTS_DIR}/stdout_run_optitype -w "done(HLA1_extract_reads_${PID})" -env "all,RESULTS_DIR=${RESULTS_DIR},OPTITYPE_SIF=${OPTITYPE_SIF}" < ${PIPELINE_DIR}/HLA_caller_scripts/HLA5_run_optitype.sh`


	# HLA5_run_optitype=`bsub -J HLA5_run_optitype_${PID} -W 240:00 -n 16 -R "rusage[mem=50G]" -env "all,RESULTS_DIR=${RESULTS_DIR},OPTITYPE_SIF=${OPTITYPE_SIF}" < ${PIPELINE_DIR}/HLA_caller_scripts/HLA5_run_optitype.sh`

	## run kourami on extracted reads ====

	HLA7_run_kourami=`bsub -J HLA7_run_kourami_${PID} -W 24:00 -n 5 -R "rusage[mem=5000]" -o ${RESULTS_DIR}/stdout_run_kourami_%J -w "done(HLA1_extract_reads_${PID})" -env "all,RESULTS_DIR=${RESULTS_DIR},KOURAMI_SIF=${KOURAMI_SIF},NAME=${PID}" < ${PIPELINE_DIR}/HLA_caller_scripts/HLA7_run_kourami.sh`

	## compute consensus ====
	compute_consensus=`bsub -J compute_consensus_${PID} -W 00:05 -R "rusage[mem=100]" -o ${RESULTS_DIR}/stdout_compute_consensus_%J -w "done(HLA4_run_polysolver_${PID}) && done(HLA5_run_optitype_${PID})" -env "all,RESULTS_DIR=${RESULTS_DIR},PIPELINE_DIR=${PIPELINE_DIR}/HLA_caller_scripts" < ${PIPELINE_DIR}/HLA_caller_scripts/run_compute_consensus.sh`

	echo $compute_consensus
fi


########################### Neoepitope-Prediction ########################
##################### FUSION

if [ -n "$FUSION" ]; then

	## create a new folders
	mkdir -p ${OUT_FUSION}
	echo $OUT_FUSION

# export OUT=${OUT_FUSION}/peptides_for_binding_prediction.fa
# export FUSION=${FUSION}
# export FUSION_CONFIDENCE=${FUSION_CONFIDENCE}
# ${PIPELINE_DIR}/immunopipe_fusion/run_extract_protein_fa.sh

# export peptideFasta=${OUT_FUSION}/peptides_for_binding_prediction.fa
# export sourceDir=${PATH_TO_OUT_DIR}
# export OUTPUT_DIR=${OUT_FUSION}
# export netMHCpan=${netMHCpan}
# export netMHCIIpan=${netMHCIIpan}
# export pepType=fusion
# $PIPELINE_DIR/immunopipe/run_netMHC_on_consensus.sh



# export OUTPUT_DIR=${OUT_FUSION}
# export FUSION=${FUSION}
# export PID=${PID}
# ${PIPELINE_DIR}/immunopipe_fusion/run_summary.sh



	## FUSION: extract predicted fusion peptides from Arriba output
	rm -f ${OUT_FUSION}/stdout_extractProtein_fusion
	bsub -K -J extractProtein_fusion_${PID} -W 00:10 -R "rusage[mem=2000]" -o ${OUT_FUSION}/stdout_extractProtein_fusion -env "all,OUT=${OUT_FUSION}/peptides_for_binding_prediction.fa,FUSION=${FUSION},FUSION_CONFIDENCE=${FUSION_CONFIDENCE}" < ${PIPELINE_DIR}/immunopipe_fusion/run_extract_protein_fa.sh

	## FUSION: run netMHC for fusion peptides
	rm -f ${OUT_FUSION}/stdout_netMHC_fusion
	bsub -K -J run_netMHC_fusion_${PID} -W 10:00 -n 8 -R "rusage[mem=4000]" -o ${OUT_FUSION}/stdout_netMHC_fusion -env "all,peptideFasta=${OUT_FUSION}/peptides_for_binding_prediction.fa,sourceDir=${PATH_TO_OUT_DIR},OUTPUT_DIR=${OUT_FUSION},netMHCpan=${netMHCpan},netMHCIIpan=${netMHCIIpan},pepType=fusion" < $PIPELINE_DIR/immunopipe/run_netMHC_on_consensus.sh



	## FUSION: Summarize results for MHC-I and MHC-II neoepitopes
	rm -f ${OUT_FUSION}/stdout_summary
	bsub -K -J run_summary_results_fusion_${PID} -W 00:10 -R "rusage[mem=1000]" -o ${OUT_FUSION}/stdout_summary -env "all,OUTPUT_DIR=${OUT_FUSION},FUSION=${FUSION},PID=${PID}" < ${PIPELINE_DIR}/immunopipe_fusion/run_summary.sh


fi # fusion file given

########################### SNV 

if [ -n "$SNV" ]; then

# create a new folders
mkdir ${RESULTS_DIR}/neoepitopes_snv
OUT_SNV=${RESULTS_DIR}/neoepitopes_snv

# SNV: predict mutated protein and reference protein sequences
run_predict_protein_snv=`bsub -J predictProtein_snv_${PID} -W 01:00 -R "rusage[mem=10000]" -o ${OUT_SNV}/stdout_predictProtein_%J -env "all,OUT_DIR=${OUT_SNV},ANNOVAR=${ANNOVAR},SNV=${SNV},predictedProtein=${PID}_predictedProtein.fa,PIPELINE_DIR=${PIPELINE_DIR}/immunopipe" < ${PIPELINE_DIR}/immunopipe/run_predict_protein_fa.sh`

echo $run_predict_protein_snv

# SNV: extract mutated and reference peptides from predicted protein fasta (e.g. 29aa)
run_separate_ref_mut_peptide_snv=`bsub -J separate_ref_mut_peptide_snv_${PID} -W 01:00 -R "rusage[mem=10000]" -o ${OUT_SNV}/stdout_separatePeptide_%J -w "done(predictProtein_snv_${PID})" -env "all,fastaFile=${OUT_SNV}/${PID}_predictedProtein.fa,outputMut=${OUT_SNV}/${PID}_mut.fa,outputRef=${OUT_SNV}/${PID}_ref.fa,extendLen=${extendLen},PIPELINE_DIR=${PIPELINE_DIR}/immunopipe" < ${PIPELINE_DIR}/immunopipe/run_separate_ref_mut_peptide.sh`

echo $run_separate_ref_mut_peptide_snv

# SNV: run netMHC for mutated peptides
run_netMHC_mut_snv=`bsub -J run_netMHC_mutated_snv_${PID} -W 10:00 -n 8 -R "rusage[mem=4000]" -o ${OUT_SNV}/stdout_netMHC_mutated_%J -w "done(compute_consensus_${PID}) && done(separate_ref_mut_peptide_snv_${PID})" -env "all,peptideFasta=${OUT_SNV}/${PID}_mut.fa,OUTPUT_DIR=${OUT_SNV},netMHCpan_SIF=${netMHCpan_SIF},netMHCIIpan_SIF=${netMHCIIpan_SIF},pepType=mut" < $PIPELINE_DIR/immunopipe/run_netMHC_on_consensus.sh`

echo $run_netMHC_mut_snv

# SNV: run netMHC for reference peptides
run_netMHC_ref_snv=`bsub -J run_netMHC_ref_snv_${PID} -W 10:00 -n 8 -R "rusage[mem=4000]" -o ${OUT_SNV}/stdout_netMHC_ref_%J -w "done(compute_consensus_${PID}) && done(separate_ref_mut_peptide_snv_${PID})" -env "all,peptideFasta=${OUT_SNV}/${PID}_ref.fa,OUTPUT_DIR=${OUT_SNV},netMHCpan_SIF=${netMHCpan_SIF},netMHCIIpan_SIF=${netMHCIIpan_SIF},pepType=ref" < $PIPELINE_DIR/immunopipe/run_netMHC_on_consensus.sh `

echo $run_netMHC_ref_snv

# SNV: Summarize results for MHC-I and MHC-II neoepitopes
run_summary_snv=`bsub -J run_summary_results_snv_${PID} -W 01:00 -R "rusage[mem=1000]" -o ${OUT_SNV}/stdout_summary_%J -w "done(run_netMHC_ref_snv_${PID}) && done(run_netMHC_mutated_snv_${PID})" -env "all,OUTPUT_DIR=${OUT_SNV},VCF=${SNV},PID=${PID},PIPELINE_DIR=${PIPELINE_DIR}/immunopipe,extendLen=${extendLen}" < ${PIPELINE_DIR}/immunopipe/run_summary.sh`

echo $run_summary_snv

fi # SNV file given

##################### INDEL 

if [ -n "$INDEL" ]; then

# create a new folders
mkdir ${RESULTS_DIR}/neoepitopes_indel
OUT_INDEL=${RESULTS_DIR}/neoepitopes_indel

# INDEL: predict mutated protein and reference protein sequences
run_predict_protein_indel=`bsub -J predictProtein_indel_${PID} -W 01:00 -R "rusage[mem=10000]" -o ${OUT_INDEL}/stdout_predictProtein_%J -env "all,OUT_DIR=${OUT_INDEL},ANNOVAR=${ANNOVAR},INDEL=${INDEL},predictedProtein=${PID}_predictedProtein.fa,PIPELINE_DIR=${PIPELINE_DIR}/immunopipe_indel" < ${PIPELINE_DIR}/immunopipe_indel/run_predict_protein_indel.sh`

echo $run_predict_protein_indel

# INDEL: extract mutated and reference peptides from predicted protein fasta (e.g. 29aa)
run_separate_ref_mut_peptide_indel=`bsub -J separate_ref_mut_peptide_indel_${PID} -W 01:00 -R "rusage[mem=10000]" -o ${OUT_INDEL}/stdout_separatePeptide_%J -w "done(predictProtein_indel_${PID})" -env "all,fastaFile=${OUT_INDEL}/${PID}_predictedProtein.fa,outputMut=${OUT_INDEL}/${PID}_mut.fa,outputRef=${OUT_INDEL}/${PID}_ref.fa,extendLen=${extendLen},PIPELINE_DIR=${PIPELINE_DIR}/immunopipe_indel" < ${PIPELINE_DIR}/immunopipe_indel/run_separate_ref_mut_peptide_indel.sh`

echo $run_separate_ref_mut_peptide_indel

# INDEL: run netMHC for mutated peptides
run_netMHC_mut_indel=`bsub -J run_netMHC_mutated_indel_${PID} -W 10:00 -n 8 -R "rusage[mem=4000]" -o ${OUT_INDEL}/stdout_netMHC_mutated_%J -w "done(compute_consensus_${PID}) && done(separate_ref_mut_peptide_indel_${PID})" -env "all,peptideFasta=${OUT_INDEL}/${PID}_mut.fa,OUTPUT_DIR=${OUT_INDEL},netMHCpan_SIF=${netMHCpan_SIF},netMHCIIpan_SIF=${netMHCIIpan_SIF},pepType=mut" < $PIPELINE_DIR/immunopipe_indel/run_netMHC_on_consensus.sh`

echo $run_netMHC_mut_indel

# INDEL: run netMHC for reference peptides
run_netMHC_ref_indel=`bsub -J run_netMHC_ref_indel_${PID} -W 10:00 -n 8 -R "rusage[mem=4000]" -o ${OUT_INDEL}/stdout_netMHC_ref_%J -w "done(compute_consensus_${PID}) && done(separate_ref_mut_peptide_indel_${PID})" -env "all,peptideFasta=${OUT_INDEL}/${PID}_ref.fa,OUTPUT_DIR=${OUT_INDEL},netMHCpan_SIF=${netMHCpan_SIF},netMHCIIpan_SIF=${netMHCIIpan_SIF},pepType=ref" < $PIPELINE_DIR/immunopipe_indel/run_netMHC_on_consensus.sh `

echo $run_netMHC_ref_indel

# INDEL: Summarize results for MHC-I and MHC-II neoepitopes
run_summary_indel=`bsub -J run_summary_results_indel_${PID} -W 01:00 -R "rusage[mem=1000]" -o ${OUT_INDEL}/stdout_summary_%J -w "done(run_netMHC_ref_indel_${PID}) && done(run_netMHC_mutated_indel_${PID})" -env "all,OUTPUT_DIR=${OUT_INDEL},VCF=${INDEL},PID=${PID},PIPELINE_DIR=${PIPELINE_DIR}/immunopipe_indel,extendLen=${extendLen}" < ${PIPELINE_DIR}/immunopipe_indel/run_summary.sh`

echo $run_summary_indel

fi # indel file given

