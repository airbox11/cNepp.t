#!/usr/bin/env bash
set -e

module load perl/5.24.1
module load R/3.6.2
module load samtools/1.5

STAR=/home/lyuya/miniconda3/bin/STAR
arriba=/software/arriba/2.1.0/bin/arriba
python=/home/lyuya/miniconda3/bin/python
export PATH="/home/lyuya/miniconda3/bin/:$PATH"
export PATH="/software/r/3.6.2/lib64/R/bin/:$PATH"



## prepare data_manually ====
# hla fastq
# snv vcf, [bam, tsv, TCGA]
# indel vcf

## input parameters: ====
dir_work=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1

if [ $hlaID = promise ]; then

	tumorID=`echo $runID | cut -d _ -f 2`
	wegs=`echo $runID | cut -d _ -f 3`
	tcga=`echo $runID | cut -d _ -f 4-5 --output-delimiter='_'`
	runID=`echo $runID | cut -d _ -f 1`

	if [ $tumorID == 'tumor1' -o $tumorID == 'tumor11' ];then
		stage=T1T2
	elif [ $tumorID == 'tumor12' ];then
		stage=T12
	elif [ $tumorID == 'tumor31' -o $tumorID == 'tumor3' ];then
		stage=T2ERW
	elif [ $tumorID == 'tumor41' ];then
		stage=T3
	elif [ $tumorID == 'tumor51' ];then
		stage=T3ERW
	elif [ $tumorID == 'tumor61' ];then
		stage=T4
	fi

	dir_input=/omics/odcf/project/OE0422/promise/sequencing
	workDir=$dir_work/promise/batch_process_20230302/result/${runID}_${stage}_${tumorID}
	outputDir=/omics/odcf/analysis/OE0422_projects/promise/results_per_pid_2/$runID

	if [ $wegs == 'wes' ]; then
		dir_input_ws=exon_sequencing
	else
		dir_input_ws=whole_genome_sequencing
	fi


else
	dir_input=/omics/odcf/project/OE0422/immuno_patients_nct/sequencing
	workDir=$dir_work/${runID}

	wegs=exon_sequencing
	# wegs=whole_genome_sequencing
	outputDir=/omics/odcf/analysis/OE0422_projects/Immuno-Patients-NCT/sequencing/${wegs}/results_per_pid/$runID 
fi



workDir9=$workDir/9_Fusion_gene_based_neoepitope_identification
workDir92=$workDir9/2_arriba_result_${dataType}
workDir93=$workDir9/3_neoPrediction_${dataType}
scriptsDir=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline


## functions ================================

## SNV BASED PREIDCTION ====
## 0) prepare folder ====

function prepare_folder () {
	echo '=== ===> prepare_folder: start...'

	mkdir -p $workDir
	mkdir -p $workDir/1_hla_type
	mkdir -p $workDir/2_SNVs_based_neoepitope_prediction
	mkdir -p $workDir/3_add_expression
	mkdir -p $workDir/4_indel_based_prediction
	mkdir -p $workDir/8_chose_neoepitode 
	mkdir -p $workDir/9_Fusion_gene_based_neoepitope_identification
	mkdir -p $workDir92

	mkdir -p $outputDir/Mutation_analysis/snv 
	mkdir -p $outputDir/Mutation_analysis/CGI 
	mkdir -p $outputDir/Mutation_analysis/indel 
	mkdir -p $outputDir/Epitope_prediction/snv_based 
	mkdir -p $outputDir/Epitope_prediction/indel_based 
	mkdir -p $outputDir/Epitope_prediction/fusion_genes 
	mkdir -p $outputDir/HLA/DKMS 
	mkdir -p $outputDir/HLA/In_silico 
	mkdir -p $outputDir/HLA/LOH
	mkdir -p $outputDir/Gene_Expression


	## creat link to mhc/runID: ====

	## rna input:

	if [ $hlaID = promise ]; then

		echo ${runID}_${tumorID}

		## get tumorID (complete folder name ) ======== ======== ========

		tumor_sample_number=$(find $dir_input/$dir_input_ws/view-by-pid/$runID -maxdepth 1 -type d -name "*$tumorID*" | wc -l)
		if [ $tumor_sample_number == 1 ]; then
			tumorID_2=$(find $dir_input/$dir_input_ws/view-by-pid/$runID -maxdepth 1 -type d -name "*$tumorID*" | xargs -I {} basename {})
		else
			echo 'Error: multple tumor samples with same tumorID!'
			exit 1
		fi
		# tumorID_2=${tumorID}-01

		## prepare RNA bam file ======== ======== ========


		# if [[ $tcga == 'RNAseq' ]]; then
		if [[ `echo $tcga | grep -o 'RNAseq' | wc -l` == 1 ]]; then
			dir_RNA=$dir_input/rna_sequencing/view-by-pid/${runID}/${tumorID_2}/paired/merged-alignment/.merging_0
			if [ ! -d $dir_RNA ];then
				tumorIDcount=$(ls -d $dir_input/rna_sequencing/view-by-pid/${runID}/${tumorID}* | wc -l)
				if [ $tumorIDcount == 1 ]; then
					dir_RNA=$dir_input/rna_sequencing/view-by-pid/${runID}/${tumorID}*/paired/merged-alignment/.merging_0
				else
					echo 'Error: there is no or multiple tumor RNA data, please check the details!'
					exit 1
				fi
			fi

			file1=`find $dir_RNA -name "*${runID}_merged.mdup.bam"| xargs -I {}  realpath {}`
			file2=`find $dir_RNA -name '*fpkm_tpm.featureCounts.tsv'| xargs -I {}  realpath {}`
			ln -sf $file1 $workDir/3_add_expression
			ln -sf $file2 $workDir/3_add_expression
		fi

		dir_wgs=$dir_input/${dir_input_ws}/view-by-pid/${runID}

		## prepare WGS/WES bam file for HLA-prediction ======== ======== ========
		tumorID_index=`echo $tumorID | sed 's/tumor//'`

		buffcoatID_1=buffy-coat${tumorID_index}-01
		buffcoatID_2=buffy-coat11-01


		if [ -d "${dir_wgs}/${buffcoatID_1}" ];then
			dir_wgs2=$dir_input/${dir_input_ws}/view-by-pid/${runID}/$buffcoatID_1/paired/merged-alignment
		elif [ -d "${dir_wgs}/$buffcoatID_2" ];then
			dir_wgs2=$dir_input/${dir_input_ws}/view-by-pid/${runID}/$buffcoatID_2/paired/merged-alignment
			buffcoatID_1=$buffcoatID_2
		else
			echo 'Error: no control bam file found for HLA typping prediction!' 
			echo $workDir/error_report
			echo 'Error: no control bam file found for HLA typping prediction!' >> $workDir/error_report
			ls -ld $dir_input/${dir_input_ws}/view-by-pid/${runID} >> $workDir/error_report
			exit 1
		fi


		file1=`find $dir_wgs2 -maxdepth 1 -name '*merged.mdup.bam'`
		ln -sf $file1 $workDir/1_hla_type

		## snv, indel ======== ======== ========
		dir_snv=$dir_wgs/snv_results/paired/${tumorID_2}_${buffcoatID_1}
		file1=`find ${dir_snv} -name '*somatic_functional_snvs_conf_8_to_10.vcf'`
		ln -sf $file1 $workDir/2_SNVs_based_neoepitope_prediction

		dir_indel=$dir_wgs/indel_results/paired/${tumorID_2}_${buffcoatID_1}
		file1=`find ${dir_indel} -type f -name '*somatic_functional_indels_conf_8_to_10.vcf'`
		ln -sf $file1 $workDir/4_indel_based_prediction

		echo '>>>>>>>>>>>>>>> s0: prepare_folder: done'

	else
		runID_for_links=$runID

		if [ $(echo $runID_for_links | grep -i _tumor | wc -l) -gt 0 ];then
			runID_for_links=$(echo $runID_for_links | sed 's/_tumor.*//g')
		fi


		## links for RNAseq data ==== 
		dir_RNA=$dir_input/rna_sequencing/view-by-pid/IPNCT_${runID_for_links}
		if [ -d $dir_RNA ]; then
			cd $dir_RNA
			if [ $(ls $dir_RNA | wc -l) -eq 1 ];then
				input1=`ls`
			elif [ $(ls $dir_RNA | wc -l) -gt 1 ];then
				ls $dir_RNA
				read -p "Input Selection for RNAseq tumor ID:" input1
			fi
			cd $input1/paired/merged-alignment/.merging_0
			file1=`find . -name "*${runID_for_links}_merged.mdup.bam"| xargs -I {}  realpath {}`
			file2=`find . -name '*fpkm_tpm.featureCounts.tsv'| xargs -I {}  realpath {}`

			ln -sf $file1 $workDir/3_add_expression
			ln -sf $file2 $workDir/3_add_expression
		fi

		## bam for HLA-prediction ==== 
		cd $dir_input/${wegs}/view-by-pid/IPNCT_${runID_for_links}
		input_control=control
		if [ $(ls ./control* | wc -l) -gt 1 ];then
			ls 
			read -p "Input Selection for exon/wes control ID:" input_control
		fi
		cd $input_control/paired/merged-alignment/.merging_0

		file1=`find . -name '*merged.mdup.bam'| xargs -I {}  realpath {}`
		ln -sf $file1 $workDir/1_hla_type

		## snv, indel ====
		cd $dir_input/${wegs}/view-by-pid/IPNCT_${runID_for_links}

		input_tumor=tumor
		if [ $(ls ./tumor* | wc -l) -gt 1 ];then
			ls -d tumor*
			read -p "Input Selection for wes/wgs tumor ID:" input_tumor
		fi

		files=`find . -name '*_functional_snvs_conf_8_to_10.vcf'| grep "${input_tumor}_${input_control}/" |xargs -I {}  realpath {}`
		for file in ${files[@]}; do
			ln -sf $file $workDir/2_SNVs_based_neoepitope_prediction
		done


		files=`find . -name '*_functional_indels_conf_8_to_10.vcf'|grep "${input_tumor}_${input_control}/" | xargs -I {}  realpath {}`
		for file in ${files[@]}; do
			ln -sf $file $workDir/4_indel_based_prediction
		done

	fi
	ls -l --color=auto $workDir/*

	rm -f $file_run_status
	touch $file_run_status

}



function rename_zip {

	echo '=== ===> rename_zip: start...'
	zipDir=$outputDir/zip
	# rm -fr ${runID}.zip
	rm -fr $outputDir/*.zip
	rm -fr $zipDir
	rsync -a $outputDir $zipDir 2>&1 > /dev/null


	IFS=$'\n'
	files=($(find $zipDir -type f))
	unset IFS
	for i in ${files[*]}
	do
		mv $i `echo $i | sed "s;^\(.*/\)\([^/]\+\)$;\1${runID}_\2;"` 
	done

	cd $zipDir/$runID
	zip -r $outputDir/${runID}.zip *
	cd -

	rm -fr $zipDir
	echo $outputDir
	ls -l --color $outputDir
	find $outputDir -type f
}

## 1) HLA typing prediction ====
function s1a_HLA () {
	echo '=== ===> s1a_HLA: start...'
	phlat_dir=${workDir}/1_hla_type/hla_${hlaID}
	rm -rf $phlat_dir
	mkdir -p $phlat_dir


	if [ $hlaID = 'promise' ]; then
		sourceDir=/omics/odcf/project/OE0422/promise/sequencing/${dir_input_ws}/view-by-pid/$runID/buffy-coat*/paired/run*/sequence/

	elif [ $hlaID = 'phlat_all' ]; then
		sourceDir=$workDir/1_hla_type
	fi
	fq1=${workDir}/1_hla_type/fq1.fastq.gz
	fq2=${workDir}/1_hla_type/fq2.fastq.gz

	if [ ! -f $fq1 ] || [ ! -f $fq2 ]; then
		f1=`find $sourceDir -name '*R1.fastq.gz' | sort | sed 's:\n: :g'`
		f2=`find $sourceDir -name '*R2.fastq.gz' | sort | sed 's:\n: :g'`

		if [ -z $f1 ] || [ -z $f2 ]; then
			echo '>>>>>>>>>>> >>>>>>>>>>> no fastq files found to dump together'
			exit 1
		fi

		cat $f1 > $fq1 &
		cat $f2 > $fq2 &
		wait 
	fi

	phlatRelease=/icgc/dkfzlsdf/analysis/G200/immuno/tools/phlat-release
	bowtie2=/software/bowtie2/2.3.5.1/bowtie2 
	threads=16
	$python -O $phlatRelease/dist/PHLAT.py -1 $fq1  -2 $fq2 -index $phlatRelease/b2folder -b2url ${bowtie2} -tag "phlat_$hlaID" -e $phlatRelease -o ${phlat_dir} -p ${threads}

	# rm -f ${workDir}/1_hla_type/fq*.fastq.gz

	## format hla ====
	input=${workDir}/1_hla_type/phlat_$hlaID/*_HLA.sum
	output=${workDir}/1_hla_type/format_hla
	script=/omics/groups/OE0422/internal/scripts/HLA_typing/run_format_phlat.sh
	sh $script -p $input -o $output


}


function s1b_hla_sab () {
	echo '=== ===> s1b_hla_sab: start...'
	PID=hla_sab_$runID


	workDirHla=$workDir/1_hla_type
	scriptDir=/omics/groups/OE0422/internal/yanhong/git/neoepitope_prediction-master_gitlab/neoepitope_prediction
	script=$scriptDir/Immunopipe.sh
	referenceDir=$scriptDir/data


	if [ $hlaID = 'promise' ]; then
		sourceDir=/omics/odcf/project/OE0422/promise/sequencing/${dir_input_ws}/view-by-pid/$runID/buffy-coat*/paired/merged-alignment
	else
		sourceDir=$workDir/1_hla_type
	fi

	if [ -f ${workDirHla}/format_hla ] && [ $(wc -l < ${workDirHla}/format_hla) -gt 1 ];then
		cat $workDir/1_hla_type/format_hla
		return 0
	else 
		if [ -d $outputDir/HLA ]; then
			format_hla=`find ${outputDir}/HLA -type f -name 'format_hla'`
		fi

		if [[ -f $format_hla ]] && [ $(wc -l < $format_hla) -gt 1 ]; then
			ln -fs $format_hla $workDir/1_hla_type
		else
			hlaFile=`find $workDir/1_hla_type -name '*.result_formatted.tsv'`
			if [ -z $hlaFile ];then
				bam=`ls $sourceDir | grep -P '^(?!.*chimeric).*mdup.bam$'| xargs -I {} echo $sourceDir/{}`
				bam=`realpath $bam`
				$script --control-dna-bam $bam --references $referenceDir --output $workDirHla --pid $PID --hla hla
			fi

			if [ -f $hlaFile ] && [ $(wc -l < $hlaFile) -ge 1 ]; then
				less $hlaFile | sed -E 's/^([ABC].*)\t(.*)/HLA-\1\2/g' | sed -E 's/^(D.*)\t(.*):(.*)/\1_\2\3/g' > $workDir/1_hla_type/format_hla
				if [ $hlaID != 'promise' ]; then
					cp $workDir/1_hla_type/format_hla $outputDir/HLA/In_silico
				fi
			else
				echo '==== ==== Error: no format_hla file '
				exit 1
			fi
		fi
	fi

}




## 2) neoepitode predition with snv calling result ====


function f_vcfOnly () {
	rm -rf 1.vcf
	cat /omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/vcf_header.txt >> 1.vcf
	cat *org | grep -v '^#' >> 1.vcf
	script=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/prepare_vcf.r
	Rscript $script $workDir/2_SNVs_based_neoepitope_prediction
}


function s2_snv () {
	echo '=== ===> s2_snv: start...'
	# if [ $hlaID = 'promise' ]; then
		# vcf=`find /omics/odcf/project/OE0422/promise/sequencing/${dir_input_ws}/view-by-pid/$runID -name "snv*somatic_functional_snvs_conf_8_to_10.vcf" | grep ${tumorID}-`

		# format_hla=`find $workDir/1_hla_type -name "format_hla"`
		# ln -fs $vcf $workDir/2_SNVs_based_neoepitope_prediction
		# ln -fs $vcf $outputDir/Mutation_analysis/snv/$stage_$tumorID
		# if [ ! -f ${workDir}/1_hla_type/format_hla ]; then
		# 	ln -fs $format_hla ${workDir}/1_hla_type/format_hla
		# fi
	# fi


	cd ${workDir}/2_SNVs_based_neoepitope_prediction
	if [ -f *org ] & [ $vcfOnly == 'Yes' ] || [ $vcfOnly == 'yes' ] ;then
		f_vcfOnly
		vcf=${workDir}/2_SNVs_based_neoepitope_prediction/2.vcf

	else
		vcf=${workDir}/2_SNVs_based_neoepitope_prediction/*somatic*.vcf
	fi
	cd -

	# module load python/2.7.9
	script=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/update_netMHCpan/netMHCpan_4.1/main_run_neoepitope_snvs.sh
	format_hla=${workDir}/1_hla_type/format_hla


	sh $script \
	-i ${netMHCpanID} \
	-v $vcf \
	-m $format_hla \
	-o ${workDir}/2_SNVs_based_neoepitope_prediction/${netMHCpanID}


}


## 3) If RNA-seq exist, excute 3.2, otherwise excute 3.1 ====
function s3_rna () {
	echo '=== ===> s3_rna: start...'
	inputMHCI=${workDir}/2_SNVs_based_neoepitope_prediction/$netMHCpanID/results_MHCI_epitopes.tab_splitGenes
	inputMHCII=${workDir}/2_SNVs_based_neoepitope_prediction/$netMHCpanID/results_MHCII_epitopes.tab_splitGenes

	rm -rf ${workDir}/3_add_expression/*tab 


	tab_RNA_tmp=0

	# if [[ $tcga == 'RNAseq' ]];then
	if [[ `echo $tcga | grep -o 'RNAseq' | wc -l` == 1 ]]; then
		outputMHCI=${workDir}/3_add_expression/MHCI_epitopes_RNAseq_$netMHCpanID.tab
		outputMHCII=${workDir}/3_add_expression/MHCII_epitopes_RNAseq.tab

		script=/omics/groups/OE0422/internal/scripts/add_expression/main_add_RNA.sh
		pipelineDir=/omics/groups/OE0422/internal/scripts/add_expression
		vcf=${workDir}/2_SNVs_based_neoepitope_prediction/*somatic*.vcf

		if [ $hlaID = 'promise' ]; then

			rnaDir=/omics/odcf/project/OE0422/promise/sequencing/rna_sequencing/view-by-pid/$runID/$tumorID-*
			RNA_bam=`find $rnaDir -name '*merged.mdup.bam'  | grep -P '^(?!.*chimeric).*merged.mdup.bam$'`
			expression=`find $rnaDir -name '*.fpkm_tpm.featureCounts.tsv'`

			ln -fs $expression $workDir/3_add_expression
			ln -fs $RNA_bam $workDir/3_add_expression
		fi


		RNA_bam=`ls ${workDir}/3_add_expression | grep -P '^(?!.*chimeric).*mdup.bam$'| xargs -I {} echo ${workDir}/3_add_expression/{}`
		expression=${workDir}/3_add_expression/*fpkm_tpm.featureCounts.tsv

		if [ -f $inputMHCI ]; then
			echo 
			sh $script -p $pipelineDir -v $vcf -b $RNA_bam -e $expression -n $inputMHCI -o $outputMHCI
		else
			echo '>>>> >>>> WARNING: no file exist:'
			echo $inputMHCI
		fi

		if [ -f $inputMHCII ]; then
			echo 
			sh $script -p $pipelineDir -v $vcf -b $RNA_bam -e $expression -n $inputMHCII -o $outputMHCII 
		else
			echo '>>>> >>>> WARNING: no file exist:'
			echo $inputMHCII
		fi

		script=$scriptsDir/check_wish_list_genes.r
		outputFile=$workDir/8_chose_neoepitode/wish_list_genes_expression.csv
		Rscript $script $expression $outputFile

		tab_RNA_tmp=1



	fi
	if [[ `echo $tcga | grep -oi 'tcga' | wc -l` == 1 ]]; then

		tcgaID=`echo $tcga | sed 's/\(TCGA-[^_]\+\).*/\1/g'`
		code=/omics/groups/OE0422/internal/scripts/add_expression/run_add_refExpression.R
		tcgaFile=/omics/groups/OE0422/internal/tcga_fpkm/${tcgaID}/${tcgaID}_expression_addName.tab

		if [[ -f ${workDir}/3_add_expression/MHCI_epitopes_RNAseq_$netMHCpanID.tab ]]; then
			inputMHCI=${workDir}/3_add_expression/MHCI_epitopes_RNAseq_$netMHCpanID.tab
		fi
		outputMHCI=${workDir}/3_add_expression/MHCI_epitopes_${tcga}_$netMHCpanID.tab
		Rscript $code $tcgaFile $inputMHCI $outputMHCI 



		if [[ -f ${workDir}/3_add_expression/MHCII_epitopes_RNAseq.tab ]]; then
			inputMHCII=${workDir}/3_add_expression/MHCII_epitopes_RNAseq.tab
		fi
		outputMHCII=${workDir}/3_add_expression/MHCII_epitopes_${tcga}.tab
		if [ -f $inputMHCII ];then
			Rscript $code $tcgaFile $inputMHCII $outputMHCII 
		else
			echo ">>>> >>>> no file: $inputMHCII"
		fi
		if [ $tab_RNA_tmp -eq 1 ];then
			mkdir -p  ${workDir}/3_add_expression/past
			mv $inputMHCI $inputMHCII $workDir/3_add_expression/past
		fi
	fi


}

function s8a_filter () {
	rm -rf $workDir/8_chose_neoepitode/MHC*
	echo '=== ===> s8a_filter: start...'
	script=$scriptsDir/8a_filter_neoepitode.r
	Rscript $script $netMHCpanID $workDir $tcga

	## rename columns ==== ==== 
	cd $workDir/8_chose_neoepitode
	IFS=$'\n'
	files=($(find . -maxdepth 1 -name '*tab'))
	unset IFS
	for i in ${files[*]}
	do
		awk 'BEGIN{FS='\t'; OFS='\t'}{if (NR == 1) {gsub(" +","_",$0)}; print $0}' $i > ${i}_renameCol
	done
	cd -
}


function s8b_xlsx_to_public (){
	echo '=== ===> s8b_xlsx_to_public: start...'
	if [ $hlaID = 'promise' ]; then
		outputDir_snv=${outputDir}/Epitope_prediction/snv_based/${stage}_${tumorID}
	else
		outputDir_snv=${outputDir}/Epitope_prediction/snv_based
	fi
	mkdir -p $outputDir_snv
	/tbi/software/x86_64/R/R-3.4.0/el7/bin/Rscript ${scriptsDir}/convert_to_xlsx.r $workDir $outputDir_snv snv

	## convert wishList to xlsx file
	if [ `echo $tcga | grep 'RNAseq' | wc -l` == 1 ];then
		inputFile=$workDir/8_chose_neoepitode/wish_list_genes_expression.csv
		/tbi/software/x86_64/R/R-3.4.0/el7/bin/Rscript /omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/convert_to_xlsx.r $inputFile $outputDir/Gene_Expression wishList
	fi


}


## INDEL BASED PREIDCTION ====

function i4a_indel_predict () {
	echo '=== ===> i4a_indel_predict: start...'

	if [ $hlaID = 'promise' ]; then
		vcf=`find $indel_source_dir -name 'indel*somatic_functional_indels_conf_8_to_10.vcf'| xargs realpath`
		indel_vcf_dir=$outputDir/Mutation_analysis/indel/$stage_$tumorID
		mkdir -p $indel_vcf_dir
		ln -fs $vcf $indel_vcf_dir
	fi

	vcf=${workDir}/4_indel_based_prediction/*somatic*.vcf
	if [ ! -e ${workDir}/4_indel_based_prediction/*vcf ]; then
		echo ">>>> >>>> indel vcf file not exist!! "
		exit 1
	fi

	lineCount=`wc -l ${workDir}/4_indel_based_prediction/*vcf | cut -d ' ' -f1`
	if (( $lineCount <= 1 )); then
		echo ">>>> >>>> linecount = ${lineCount};  no indel predicted"
		echo 'indel_stop=1' >> $file_run_status
		return
	fi

	indel_result=$workDir/4_indel_based_prediction/result
	mkdir -p $indel_result

	script=/omics/groups/OE0422/internal/scripts/neoepitope_indels/main_indels.sh
	hla=${workDir}/1_hla_type/format_hla

	if [ `less $vcf | wc -l ` -gt 1 ]; then
		sh $script -f $vcf -l 21 -a $hla -o $indel_result
	else
		echo 'indel_stop=1' >> $file_run_status
	fi
}


function i4b_indel_tsv () {
	source $file_run_status
	if [[ $indel_stop == 1 ]]; then
		return
	fi
	echo '=== ===> i4b_indel_tsv: start...'
	script=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/4b_indel_table.r
	Rscript $script $workDir

}


function i4c_xlsx_to_public () {
	echo '=== ===> i4c_xlsx_to_public: start...'

	if [ $indel_stop == 1 ]; then
		return
	fi

	if [ $hlaID = 'promise' ]; then
		outputDir_indel=${outputDir}/Epitope_prediction/indel_based/${stage}_${tumorID}
	else
		outputDir_indel=${outputDir}/Epitope_prediction/indel_based
	fi

	mkdir -p $outputDir_indel

	/tbi/software/x86_64/R/R-3.4.0/el7/bin/Rscript ${scriptsDir}/convert_to_xlsx.r $workDir $outputDir_indel indel
}

## fusion BASED PREIDCTION ====
## step 1 ====

function f1a_run_star_arriba () {
	echo '=== ===> f1a_run_star_arriba: start...'
	bam=$workDir/3.2_add_local_rna/*merged.mdup.bam
	dirSource=`realpath  $bam | xargs dirname `


	fq1n=`find $dirSource -name '*1.fastq.gz'| wc -l `
	fq2n=`find $dirSource -name '*2.fastq.gz'| wc -l `
	if [ $fq1n -ne 1 ] || [ $fq2n -ne 1 ]; then
		echo '>>>> >>>> 1/2-fastq.gz for STAR is not unique!!'
		exit 0
	fi


	READ1=`find $dirSource -name '*1.fastq.gz' | head -1`
	READ2=`find $dirSource -name '*2.fastq.gz' | head -1`


	STAR_INDEX_DIR=/omics/groups/OE0422/internal/yanhong/git/arriba/STAR_index_GRCh37_ENSEMBL87
	annotation=/omics/groups/OE0422/internal/yanhong/git/arriba/ENSEMBL87.gtf
	assembly=/omics/groups/OE0422/internal/yanhong/git/arriba/GRCh37.fa


	cd $workDir92

	# READ1=`find ../3.2* -name '*fastq.gz' | xargs -I {} echo $workDir92/{}| grep -P  '^(?!.*?RNA).*fastq.gz' | grep -P 'R1' | tr '\n' ',' | sed 's/,$/\n/'`
	# READ2=`find ../1_prepare_fastq -name '*fastq.gz' | xargs -I {} echo $workDir92/{}| grep -P  '^(?!.*?RNA).*fastq.gz' | grep -P 'R2' | tr '\n' ',' | sed 's/,$/\n/'`
	$STAR \
		--runThreadN 8 \
		--genomeDir "$STAR_INDEX_DIR" --genomeLoad NoSharedMemory \
		--readFilesIn "$READ1" "$READ2" \
		--readFilesCommand zcat \
		--outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
		--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
		--chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 |
	tee Aligned.out.bam |
	$samtools view -h |
	$arriba -x /dev/stdin  \
	-f blacklist \
	-g $annotation -a $assembly \
	-o $workDir92/fusions.tsv -O $workDir92/fusions.discarded.tsv 

}


function f1b_run_arriba () {
	echo '=== ===> f1b_run_arriba: start...'

	if [ $hlaID = 'promise' ]; then
		rnaDir=/omics/odcf/project/OE0422/promise/sequencing/rna_sequencing/view-by-pid/$runID/$tumorID-*/paired/merged-alignment/.merging_0
		RNA_bam=`find $rnaDir -name '*merged.mdup.bam'  | grep -P '^(?!.*chimeric).*merged.mdup.bam$'`
		RNA_chimeric_bam=`find $rnaDir -name '*chimeric_merged.mdup.bam'`
	else
		bam=`ls ${workDir}/3_add_expression | grep -P '^(?!.*chimeric).*mdup.bam$'| xargs -I {} echo ${workDir}/3_add_expression/{}`
		if [[ $(echo $bam | wc -w) == 0 ]]; then
			echo '>>>> >>>> no bam file for fusion-based epitiope predicton'
			return 0
		fi
		dirSource=`realpath  $bam | xargs dirname `
		RNA_chimeric_bam=`find $dirSource -name '*chimeric_merged.mdup.bam'`
		RNA_bam=`find $dirSource -maxdepth 1 -name '*merged.mdup.bam'  | grep -P '^(?!.*chimeric).*merged.mdup.bam$'`
	fi

	sourceDir=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/fusion_arriba/arriba_reference_hg19
	blacklist=$sourceDir/blacklist_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz
	knowFusion=$sourceDir/known_fusions_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz
	proteinDomain=$sourceDir/protein_domains_hg19_hs37d5_GRCh37_v2.1.0.gff3
	gtf=$sourceDir/gencode.v19.annotation_plain.gtf
	fa=$sourceDir/hs37d5_PhiX.fa

	if [ -z $RNA_chimeric_bam ]; then
		echo '>>>> >>>> RNA_chimeric_bam needed.'
		return 0
	fi

	if [ -f $RNA_bam ]; then
		ls -l  $RNA_chimeric_bam
		ls -l  $RNA_bam
		$arriba \
			-c $RNA_chimeric_bam \
			-x $RNA_bam \
			-b $blacklist \
			-k $knowFusion \
			-t $knowFusion \
			-p $proteinDomain \
			-g $gtf \
			-a $fa \
			-o $workDir92/fusions.tsv \
			-O $workDir92/fusions.discarded.tsv 
	fi
}



## step 2: ====
function f2_prepare_HLA () {
	echo '=== ===> f2_prepare_HLA: start...'
	script=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/fusion_arriba/prepare_HLA.r
	Rscript $script $workDir $workDir92
}

## step 3: ====

function f3_neo_prediction () {
	echo '=== ===> f3_neo_prediction: start...'
	if [ ! -f $workDir92/fusions.tsv ];then
		return
	fi
	## re-format fusion.tsv
	script=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/fusion_arriba/splitRow.py
	python $script $workDir92

	## predict
	PID=${runID}_neoPrediction
	OUT_FUSION=$workDir9/3_neoPrediction_$dataType

	script=/omics/groups/OE0422/internal/yanhong/git/neoepitope_prediction-master1/Immunopipe.sh
	bam=$workDir92/Aligned.out.bam # not neccessary for fusion
	fusions_tsv=$workDir92/fusions.tsv_splitGenes


	$script --dna-bam $bam --fusions $fusions_tsv --fusion-confidence "medium high" --output $workDir92 --OUT_FUSION $OUT_FUSION --pid $PID
}

## step 4:
function f4_mer21 () {
	echo '=== ===> f4_mer21: start...'
	if [ ! -d $workDir93 ];then
		echo ">>>> >>>> folder not exist: $workDir93 "
		exit 0
	fi
	script=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/fusion_arriba/mer21.r
	Rscript $script $workDir93
}


## step 5: ====
function f5_to_xlsx () {
	echo '=== ===> f5_to_xlsx: start...'
	if [[ `echo $tcga |grep 'RNAseq' | wc -l` != 1 ]];then
		return
	fi
	script=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/fusion_arriba/convert_to_xlsx_arribaFusion.r

	if [ $hlaID = promise ]; then
		outputDir_fusion=$outputDir/Epitope_prediction/fusion_genes/${stage}_${tumorID}/
	else
		outputDir_fusion=$outputDir/Epitope_prediction/fusion_genes
	fi

	/tbi/software/x86_64/R/R-3.4.0/el7/bin/Rscript $script $workDir $outputDir_fusion $dataType

	## convert fusionTsv to xlsx
	fusions_tsv=$workDir92/fusions.tsv_splitGenes
	/tbi/software/x86_64/R/R-3.4.0/el7/bin/Rscript ${scriptsDir}/convert_to_xlsx.r ${fusions_tsv} $outputDir_fusion fusionsTsv
}

## loh ==== ====

function l1_lohAnalysis () {
	echo '=== ===> l1_lohAnalysis: start...'

	if [[ $vcfOnly == Yes ]]; then
		return
	fi

	## run LOH ==== ==== 
	set +e
	bash ${scriptsDir}/loh_hla.sh $workDir $runID $tumorID
	if [[ $? -ne 0 ]] ; then
		loh_stop=1
		echo 'loh_stop=1' >> $file_run_status
	fi

	set -e
}

function l2_addLoh_cpFile () {
	echo '=== ===> l2_addLoh_cpFile: start...'
	source $file_run_status

	if [ ! -f $workDir/5_LOHHLA/example-out/example.10.DNA.HLAlossPrediction_CI.xls ]; then
		echo no loh l1 result exist.
		return
	fi

	if [ -n $loh_stop ] & [ $loh_stop == 1 ] || [ $vcfOnly == 'Yes' ] ;then
		return
	fi
	## merge with mhc results ==== ====
	if [ ! $(find $workDir/8_chose_neoepitode -name 'MHCI_*renameCol' | wc -l) -gt 0 ];then
		echo '>>>> >>>> no file found:MHCI_*renameCol : '
		ls -l --color $workDir/8_chose_neoepitode
		return
	fi
	mhc_input=`find $workDir/8_chose_neoepitode -name 'MHCI_*renameCol' | xargs basename`
	Rscript ${scriptsDir}/add_loh.r $workDir $mhc_input
	
	## copy files to final /HLA/LOH ==== ====
	if [ $hlaID = 'promise' ]; then
		outputDir_loh=$outputDir/HLA/LOH/${stage}_${tumorID}
	else
		outputDir_loh=$outputDir/HLA/LOH
	fi
	mkdir -p $outputDir_loh
	if [ $(find $workDir/5_LOHHLA/example-out -type f -name '*xls' | wc -l) -gt 0 ];then
		# ls  $workDir/5_LOHHLA/example-out/*xls
		cp $workDir/5_LOHHLA/example-out/*xls $workDir/5_LOHHLA/example-out/Figures/*pdf $outputDir_loh
		echo $outputDir_loh
		ls -l --color $outputDir_loh
	fi

	echo '=== ===> l2_addLoh_cpFile: done'

}


## cgi download and merge
function c1 () {
	echo '=== ===> c1_cgi_onlineLaunch: start...'
	## prepare folder

	cd $workDir/1_hla_type
	mkdir -p cgi
	rm -rf cgi*log
	rm -rf $workDir/1_hla_type/mutation.csv $workDir/1_hla_type/mutation_germline.csv

	if [ $vcfOnly == 'Yes' ]; then
		script=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/vcfOnly_snv_indel_seperate.py
		python $script $workDir
	fi

	find $workDir/2* -name '*somatic*vcf' | xargs -I {} awk '(NR>1) { print $1, $2, $4, $5, $17}' {} >> mutation.csv
	find $workDir/4* -name '*somatic*vcf' | xargs -I {} awk '(NR>1) { print $1, $2, $4, $5, $18}' {} >> mutation.csv
	find $workDir/2* -name '*germline*vcf' | xargs -I {} awk '(NR>1) { print $1, $2, $4, $5, $17}' {} >> mutation_germline.csv
	find $workDir/4* -name '*germline*vcf' | xargs -I {} awk '(NR>1) { print $1, $2, $4, $5, $18}' {} >> mutation_germline.csv

	echo 'chr pos ref alt gene' > $workDir/1_hla_type/tmp.tsv
	cat $workDir/1_hla_type/mutation.csv  >> tmp.tsv
	sed 's/[[:space:]]\{1,\}/\t/g' tmp.tsv > mutation.tsv
	rm -rf tmp.tsv

	## run cgi and download data
	script=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/cgi/cgi_api.py
	python $script $workDir/1_hla_type $cgiTumorType
	cat *log
	cd -

	##
	cp $workDir/1_hla_type/cgi/* $outputDir/Mutation_analysis/CGI
	cp $workDir/1_hla_type/*csv $outputDir/Mutation_analysis
}

function c2 () {
	echo '=== ===> c2_cgi_merge2file: start...'
	## merge to snv/indel/fusion results
	script=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/cgi/merge_cgi_wishList.r
	Rscript $script $workDir
}

## EXCUTION (step order sensitive) =======================================================================================

OLDIFS=$IFS
IFS='-' read -r -a array1 <<< $steps

for step in "${array1[@]}"
do
	## prepare_folder ==== ====
	export file_run_status=$workDir/tmp_run_status

	if [ $step == 's0' ]; then
		prepare_folder
	fi
	
	if [ ! -f $file_run_status ]; then
		touch $file_run_status
	fi

	## format_hla  ==== ====

	if [ $step == 's1a' ]; then
		s1a_HLA
	fi

	if [ $step == 's1b' ]; then
		s1b_hla_sab
	fi

	## snv ==== ====

	if [ $step == 's2' ]; then
		if [ $netMHCpanID = 'both' ]; then
			netMHCpanID=netMHCpan4_1
			s2_snv
			netMHCpanID=netMHCstabpan
			s2_snv
		else
			s2_snv
		fi
	fi

	if [ $step == 's3' ]; then
		s3_rna
	fi

	if [ $step == 's8a' ]; then
		s8a_filter
	fi

	if [ $step == 's8b' ]; then
		s8b_xlsx_to_public
	fi

	### indel ==== ====

	if [ $step == 'i4a' ]; then
		i4a_indel_predict
	fi

	if [ $step == 'i4b' ]; then
		i4b_indel_tsv
	fi

	if [ $step == 'i4c' ]; then
		i4c_xlsx_to_public
	fi



	## fusion ==== ====

	if [ $step == 'f1a' ]; then
		f1a_run_star_arriba
	fi

	if [ $step == 'f1b' ]; then
		f1b_run_arriba
	fi

	if [ $step == 'f2' ]; then
		f2_prepare_HLA
	fi

	if [ $step == 'f3' ]; then
		f3_neo_prediction
	fi

	if [ $step == 'f4' ]; then
		f4_mer21
	fi

	if [ $step == 'f5' ]; then
		f5_to_xlsx
	fi

	## snv, indel, fusion  ==== ====

	if [ $step == 'snv' ]; then
		if [ $netMHCpanID = 'both' ]; then
			netMHCpanID=netMHCpan4_1
			s2_snv
			netMHCpanID=netMHCstabpan
			s2_snv
		else
			s2_snv
		fi
		s3_rna
		s8a_filter
	fi

	if [ $step == 'indel' ]; then
		i4a_indel_predict
		i4b_indel_tsv
	fi

	if [ $step == 'fusion' ]; then
		f1b_run_arriba
		f2_prepare_HLA
		f3_neo_prediction
		f4_mer21
	fi


	## lohhla ==== ====

	if [ $step == 'l1' ]; then
		l1_lohAnalysis
	fi

	if [ $step == 'l2' ]; then
		l2_addLoh_cpFile
	fi

	if [ $step == 'loh' ]; then
		l1_lohAnalysis
		l2_addLoh_cpFile
	fi

	## cgi ==== ====

	if [ $step == 'c1' ]; then
		c1
	fi

	if [ $step == 'c2' ]; then
		c2
	fi

	if [ $step == 'cgi' ]; then
		c1
		c2
	fi

	## blast score ==== ==== 

	if [ $step == 'blst' ]; then
		mhc_input=`find $workDir/8_chose_neoepitode -name 'MHCI_*loh' | xargs basename`
		bash $scriptsDir/blast_score.sh $workDir $mhc_input
		Rscript $scriptsDir/blast_score_mergeTable.r $workDir $mhc_input
	fi

	## align BLOSUM62 score ==== ==== 

	if [ $step == 'blhv' ]; then
		mhc_input=`find $workDir/8_chose_neoepitode -name 'MHCI_*loh' | xargs basename`
		Rscript $scriptsDir/align_BL62_score.r $workDir $mhc_input
	fi

	if [ $step == 'blav' ]; then
		$scriptsDir/virus_all_parralel.sh $workDir MHCI_epitopes_${tcga}_netMHCpan4_1.tab
	fi


	## change name , zip ==== ====
	if [ $step == 'nz' ]; then
		rename_zip
	fi

	## xlsx ==== ====
	if [ $step == 'xlsx' ]; then
		l1_lohAnalysis
		l2_addLoh_cpFile
		c1 
		c2

		s8b_xlsx_to_public
		i4c_xlsx_to_public
		f5_to_xlsx

		rename_zip

	fi

done

IFS=$OLDIFS
