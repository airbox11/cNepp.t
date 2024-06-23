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
export vcfOnly



## input parameters: ====
dir_root=/omics/odcf/analysis/OE0422_projects/Immuno-Patients-NCT
dir_result=$dir_root/yanhong/all_in_one_pipeline_result
export dir_pipeline=$dir_root/yanhong/all_in_one_pipeline

script_hla=$dir_pipeline/script_hla
script_fusion=$dir_pipeline/script_fusion


dir_tcga=/omics/groups/OE0422/internal/tcga_fpkm

if [[ $hlaID == promise ]]; then
	runID_tumorID=$runID
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
	workDir=$dir_result/promise/batch_process_20230302/result/${runID}_${stage}_${tumorID}

	outputDir=/omics/odcf/analysis/OE0422_projects/promise/results_per_pid_2/$runID
	outputDir=/omics/odcf/analysis/OE0422_projects/promise/promise_hg38/$runID

	if [ $wegs == 'wes' ]; then
		dir_input_ws=exon_sequencing
	else
		dir_input_ws=whole_genome_sequencing
	fi


else
	runID_tumorID=$runID
	runID=`echo $runID_tumorID | cut -d _ -f 1`
	tumorID=`echo $runID_tumorID | cut -d _ -f 2`

	dir_input=/omics/odcf/project/OE0422/immuno_patients_nct/sequencing
	# dir_input=/omics/odcf/project/OE0317/hpv/sequencing
	# dir_input=/omics/odcf/project/hipo/hipo_021/sequencing/ #debug;master


	workDir=$dir_result/${runID_tumorID}

	wegs=exon_sequencing
	if [[ `echo $wg | grep -oi 'wgs' | wc -l` == 1 ]];then
		wegs=whole_genome_sequencing
	fi

	outputDir=$dir_root/sequencing/${wegs}/results_per_pid/$runID
	outputDir=$dir_root/sequencing/${wegs}/results_per_pid/$runID_tumorID
fi

export file_run_status=$workDir/tmp_run_status
mkdir -p $workDir

if [ ! -f $file_run_status ]; then
	touch $file_run_status
else
	rm -rf $file_run_status 
	touch $file_run_status
fi


workDir9=$workDir/9_Fusion_gene_based_neoepitope_identification
workDir92=$workDir9/2_arriba_result_${dataType}
workDir93=$workDir9/3_neoPrediction_${dataType}
scriptsDir=$dir_root/yanhong/all_in_one_pipeline


## functions ================================

## SNV BASED PREIDCTION ====
## 0) prepare folder ====

mkdir -p $workDir/1_hla_type
mkdir -p $workDir/2_SNVs_based_neoepitope_prediction
mkdir -p $workDir/3_add_expression
mkdir -p $workDir/4_indel_based_prediction
mkdir -p $workDir/8_chose_neoepitode 
mkdir -p $workDir/9_Fusion_gene_based_neoepitope_identification
mkdir -p $workDir92

mkdir -p $outputDir/Mutation_analysis/CGI 

mkdir -p $outputDir/Mutation_analysis/snv 
mkdir -p $outputDir/Mutation_analysis/indel 
mkdir -p $outputDir/Mutation_analysis/fusion

mkdir -p $outputDir/Epitope_prediction/snv_based 
mkdir -p $outputDir/Epitope_prediction/indel_based 
mkdir -p $outputDir/Epitope_prediction/fusion_genes 

mkdir -p $outputDir/HLA/DKMS 
mkdir -p $outputDir/HLA/In_silico 
mkdir -p $outputDir/HLA/LOH

mkdir -p $outputDir/Gene_Expression




function prepare_folder () {

	echo '=== ===> prepare_folder: start...'


	## creat link to mhc/runID: ====

	## rna input:

	if [ $hlaID = promise ]; then
		## get tumorID (complete folder name ) ======== ======== ========

		tumor_sample_number=$(find $dir_input/$dir_input_ws/view-by-pid/$runID -maxdepth 1 -type d -name "*$tumorID*" | wc -l)
		if [ $tumor_sample_number == 1 ]; then
			tumorID_2=$(find $dir_input/$dir_input_ws/view-by-pid/$runID -maxdepth 1 -type d -name "*$tumorID*" | xargs -I {} basename {})
		else
			echo 'Error: multple tumor samples with same tumorID!'
			exit 1
		fi

		## prepare RNA bam file ======== ======== ========


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

			file1=`find $dir_RNA -name "*${runID}_merged.mdup.bam"| xargs -I {}  realpath -s {}`
			file2=`find $dir_RNA -name '*fpkm_tpm.featureCounts.tsv'| xargs -I {}  realpath -s {}`
			ln -sf $file1 $workDir/3_add_expression
			ln -sf $file2 $workDir/3_add_expression
		fi

		dir_wgs=$dir_input/${dir_input_ws}/view-by-pid/${runID}

		## prepare WGS/WES bam file for HLA-prediction ======== ======== ========
		tumorID_index=`echo $tumorID | sed 's/tumor//'`

		buffcoatID_1=buffy-coat${tumorID_index}-01
		buffcoatID_2=buffy-coat*-01

		if [ -d "${dir_wgs}/${buffcoatID_1}" ];then
			dir_wgs2=$dir_input/${dir_input_ws}/view-by-pid/${runID}/$buffcoatID_1/paired/merged-alignment
		elif [ `find ${dir_wgs} -maxdepth 1 -type d -name "$buffcoatID_2" | wc -l ` -eq 1 ];then
			dir_wgs2=$dir_input/${dir_input_ws}/view-by-pid/${runID}/$buffcoatID_2/paired/merged-alignment
			buffcoatID_1=$buffcoatID_2
		else
			ls -l ${dir_wgs}
			set +e
			read -t 15 -p "Input Selection for control ID:" RESP
			if [ -z $RESP ] ; then
				echo "Timeout"
				echo 'Error: no control bam file found for HLA typping prediction!'
				echo 'Error: no control bam file found for HLA typping prediction!' >> $workDir/error_report
				ls -ld $dir_input/${dir_input_ws}/view-by-pid/${runID} >> $workDir/error_report
			else
				buffcoatID_1=$RESP
			fi
			set -e
		fi


		file1=`find $dir_wgs2 -maxdepth 1 -name '*merged.mdup.bam'`
		ln -sf $file1 $workDir/1_hla_type

		## snv, indel ======== ======== ========
		dir_snv=$dir_wgs/snv_results/paired/${tumorID_2}_${buffcoatID_1}
		file1=`find ${dir_snv} -name '*somatic_functional_snvs_conf_8_to_10.vcf'`
		ln -sf $file1 $workDir/2_SNVs_based_neoepitope_prediction
		file1=`find ${dir_snv} -name '*germline_functional_snvs_conf_8_to_10.vcf'`
		ln -sf $file1 $workDir/2_SNVs_based_neoepitope_prediction

		dir_indel=$dir_wgs/indel_results/paired/${tumorID_2}_${buffcoatID_1}
		file1=`find ${dir_indel} -name '*somatic_functional_indels_conf_8_to_10.vcf'`
		ln -sf $file1 $workDir/4_indel_based_prediction
		file1=`find ${dir_indel} -name '*germline_functional_indels_conf_8_to_10.vcf'`
		ln -sf $file1 $workDir/4_indel_based_prediction

		echo '>>>>>>>>>>>>>>> s0: prepare_folder: done'

	else


		## links for RNAseq data ==== 
		dir_RNA=$dir_input/rna_sequencing/view-by-pid/IPNCT_${runID}/$tumorID/paired/merged-alignment
		# dir_RNA=$dir_input/rna_sequencing/view-by-pid/${runID}/$tumorID/paired/merged-alignment #debug;master
		if [ -d $dir_RNA ]; then
			cd $dir_RNA
			file1=`find -L .  -name "*${runID}_merged.mdup.bam"| xargs -I {}  realpath -s {}`
			file2=`find -L .  -name '*fpkm_tpm.featureCounts.tsv'| xargs -I {}  realpath -s {}`

			ln -sf $file1 $workDir/3_add_expression
			ln -sf $file2 $workDir/3_add_expression
		fi

		cd $dir_input/${wegs}/view-by-pid/IPNCT_${runID}
		# cd $dir_input/${wegs}/view-by-pid/${runID} #debug;master
		## bam for HLA-prediction ==== 
		input_control=control01
		# input_control=blood #debug;master

		ls -l $input_control/paired/merged-alignment

		file1=`find $input_control/paired/merged-alignment -name '*merged.mdup.bam'| xargs realpath -s`
		ln -sf $file1 $workDir/1_hla_type


		## snv, indel ====
		files=`find . -name '*_functional_snvs_conf_8_to_10.vcf'| grep "${tumorID}_${input_control}/" |xargs -I {}  realpath -s {}`
		for file in ${files[@]}; do
			ln -sf $file $workDir/2_SNVs_based_neoepitope_prediction
		done

		files=`find . -name '*_functional_indels_conf_8_to_10.vcf'|grep "${tumorID}_${input_control}/" | xargs -I {}  realpath -s {}`
		for file in ${files[@]}; do
			ln -sf $file $workDir/4_indel_based_prediction
		done

	fi
	ls -l --color=auto $workDir/*

}



function rename_zip {

	echo '=== ===> rename and zip : start...'
	find $outputDir -type f | xargs -I {} chmod 660 {}

	# find $outputDir/Epitope_prediction/snv_based/${stage}_${tumorID} -type f 

	## copy file to omics
	if [ $hlaID == 'promise' ]; then
		set +e
		rm -rf $outputDir/DKMS
		dir_omics_cgi=$outputDir/Mutation_analysis/CGI/${stage}_${tumorID}
		dir_omics_snv=$outputDir/Mutation_analysis/snv/${stage}_${tumorID}
		dir_omics_indel=$outputDir/Mutation_analysis/indel/${stage}_${tumorID}
		dir_omics_fusion=$outputDir/Mutation_analysis/fusion/${stage}_${tumorID}
		dir_omics_geneExpression=$outputDir/Gene_Expression/${stage}_$tumorID


		find $outputDir/Mutation_analysis/CGI -maxdepth 1 -type f | xargs -I {} rm {}
		find $outputDir/Mutation_analysis -maxdepth 1 -type f | xargs -I {} rm {}
		rm -rf $outputDir/Mutation_analysis/*tumor*
		rm -rf $dir_omics_cgi
		rm -rf $dir_omics_snv
		rm -rf $dir_omics_indel
		rm -rf $dir_omics_geneExpression

		mkdir -p $dir_omics_cgi
		mkdir -p $dir_omics_snv
		mkdir -p $dir_omics_indel
		mkdir -p $dir_omics_fusion
		mkdir -p $dir_omics_geneExpression

		cp $workDir/1_hla_type/cgi/* $dir_omics_cgi
		if [ ${vcfOnly,,} == 'promise' ] ;then
			rm -rf $dir_omics_snv/*.vcf
			cp $workDir/2_SNVs_based_neoepitope_prediction/2.vcf $dir_omics_snv/snv.vcf
			cp $workDir/4_indel_based_prediction/result/1.vcf $dir_omics_indel/indel.vcf
		else
			cp $workDir/2_SNVs_based_neoepitope_prediction/*vcf $dir_omics_snv
			cp $workDir/4_indel_based_prediction/*vcf $dir_omics_indel
		fi


		file_fusion_mutation=$workDir/9_Fusion_gene_based_neoepitope_identification/2_arriba_result_RNAbamOpt/fusions.tsv_splitGenes
		if [ -f $file_fusion_mutation ];then
			cp $file_fusion_mutation $dir_omics_fusion
		fi

		if [ -f $workDir/3_add_expression/*featureCounts.tsv ]; then
			rm -rf $outputDir/Gene_Expression/*featureCounts.tsv
			cp $workDir/3_add_expression/*tsv $dir_omics_geneExpression
		fi

		cd $outputDir


		## snv vcf
		mv ./Mutation_analysis/snv/${stage}_${tumorID}/*snvs_*_germline_functional_snvs_conf_8_to_10.vcf \
			./Mutation_analysis/snv/${stage}_${tumorID}/${runID}_${stage}_snvs_germline_functional.vcf
		mv ./Mutation_analysis/snv/${stage}_${tumorID}/*snvs_*_somatic_functional_snvs_conf_8_to_10.vcf \
			./Mutation_analysis/snv/${stage}_${tumorID}/${runID}_${stage}_snvs_somatic_functional.vcf


		## indel vcf
		mv ./Mutation_analysis/indel/${stage}_${tumorID}/*indel_*_somatic_functional_indels_conf_8_to_10.vcf \
			./Mutation_analysis/indel/${stage}_${tumorID}/${runID}_${stage}_indel_somatic_functional.vcf
		mv ./Mutation_analysis/indel/${stage}_${tumorID}/*indel_*_germline_functional_indels_conf_8_to_10.vcf \
			./Mutation_analysis/indel/${stage}_${tumorID}/${runID}_${stage}_indel_germline_functional.vcf


		## fusion vcf
		mv ./Mutation_analysis/fusion/${stage}_${tumorID}/fusions.tsv_splitGenes \
			./Mutation_analysis/fusion/${stage}_${tumorID}/${runID}_${stage}_fusions.tsv


		## cgi
		mv ./Mutation_analysis/CGI/${stage}_${tumorID}/*biomarkers.tsv \
			./Mutation_analysis/CGI/${stage}_${tumorID}/${runID}_${stage}_biomarkers.tsv
		mv ./Mutation_analysis/CGI/${stage}_${tumorID}/*summary.txt \
			./Mutation_analysis/CGI/${stage}_${tumorID}/${runID}_${stage}_summary.txt
		mv ./Mutation_analysis/CGI/${stage}_${tumorID}/report.txt \
			./Mutation_analysis/CGI/${stage}_${tumorID}/${runID}_${stage}_report.txt
		mv ./Mutation_analysis/CGI/${stage}_${tumorID}/*alterations.tsv \
			./Mutation_analysis/CGI/${stage}_${tumorID}/${runID}_${stage}_alterations.tsv
		mv ./Mutation_analysis/CGI/${stage}_${tumorID}/*input01.tsv \
			./Mutation_analysis/CGI/${stage}_${tumorID}/${runID}_${stage}_input01.tsv

		## snv xlsx
		mv ./Epitope_prediction/snv_based/${stage}_${tumorID}/*MHCI_*.xlsx \
			./Epitope_prediction/snv_based/${stage}_${tumorID}/${runID}_${stage}_snv_MHCI_${tcga}.xlsx
		mv ./Epitope_prediction/snv_based/${stage}_${tumorID}/*MHCII_*.xlsx \
			./Epitope_prediction/snv_based/${stage}_${tumorID}/${runID}_${stage}_snv_MHCII_${tcga}.xlsx
		## snv wishlist
		mv ./Epitope_prediction/snv_based/${stage}_${tumorID}/*wish_list_genes_expression.xlsx \
			./Epitope_prediction/snv_based/${stage}_${tumorID}/${runID}_${stage}_wish_list_genes_expression.xlsx

		## indel xlsx
		mv ./Epitope_prediction/indel_based/${stage}_${tumorID}/*indel_long_peptides.xlsx \
			./Epitope_prediction/indel_based/${stage}_${tumorID}/${runID}_${stage}_indel_long_peptides.xlsx
		mv ./Epitope_prediction/indel_based/${stage}_${tumorID}/*indel_wildType_MHCI.xlsx \
			./Epitope_prediction/indel_based/${stage}_${tumorID}/${runID}_${stage}_indel_wildType_MHCI.xlsx
		mv ./Epitope_prediction/indel_based/${stage}_${tumorID}/*indel_mutant_MHCI.xlsx \
			./Epitope_prediction/indel_based/${stage}_${tumorID}/${runID}_${stage}_indel_mutant_MHCI.xlsx
		mv ./Epitope_prediction/indel_based/${stage}_${tumorID}/*indel_wildType_MHCII.xlsx \
			./Epitope_prediction/indel_based/${stage}_${tumorID}/${runID}_${stage}_indel_wildType_MHCII.xlsx
		mv ./Epitope_prediction/indel_based/${stage}_${tumorID}/*indel_mutant_MHCII.xlsx \
			./Epitope_prediction/indel_based/${stage}_${tumorID}/${runID}_${stage}_indel_mutant_MHCII.xlsx

		## fusion xlsx
		mv ./Epitope_prediction/fusion_genes/${stage}_${tumorID}/*MHCI_RNAbamOpt.xlsx \
			./Epitope_prediction/fusion_genes/${stage}_${tumorID}/${runID}_${stage}_MHCI_RNAbamOpt.xlsx
		mv ./Epitope_prediction/fusion_genes/${stage}_${tumorID}/*MHCII_RNAbamOpt.xlsx \
			./Epitope_prediction/fusion_genes/${stage}_${tumorID}/${runID}_${stage}_MHCII_RNAbamOpt.xlsx
		mv ./Epitope_prediction/fusion_genes/${stage}_${tumorID}/*fusionsTsv.xlsx \
			./Epitope_prediction/fusion_genes/${stage}_${tumorID}/${runID}_${stage}_fusionsTsv.xlsx

		## loh
		mv ./HLA/LOH/${stage}_${tumorID}/*example.10.DNA.HLAlossPrediction_CI.xls \
			./HLA/LOH/${stage}_${tumorID}/${runID}_${stage}_HLAlossPrediction_CI.xls
		mv ./HLA/LOH/${stage}_${tumorID}/*example_tumor_sorted.minCoverage_10.HLA.pdf \
			./HLA/LOH/${stage}_${tumorID}/${runID}_${stage}_minCoverage_10.HLA.pdf
		mv ./HLA/LOH/${stage}_${tumorID}/*example.10.DNA.IntegerCPN_CI.xls \
			./HLA/LOH/${stage}_${tumorID}/${runID}_${stage}_IntegerCPN_CI.xls

		## rna expression profile
		mv ./Gene_Expression/${stage}_${tumorID}/*tsv \
			./Gene_Expression/${stage}_${tumorID}/${runID}_${stage}_featureCounts.tsv
		# mv ./Gene_Expression/${stage}_${tumorID}/*wish_list_genes_expression.csv \
		# 	./Gene_Expression/${stage}_${tumorID}/${runID}_${stage}_wish_list_genes_expression.csv

		## wish-list 
		mv ./Gene_Expression/*wish*ist* ./Gene_Expression/${stage}_${tumorID}
		mv ./Gene_Expression/${stage}_${tumorID}/*wishList.xlsx \
			./Gene_Expression/${stage}_${tumorID}/${runID}_${stage}_wishList.xlsx
		mv ./Gene_Expression/${stage}_${tumorID}/*wish_list_genes_expression.csv \
			./Gene_Expression/${stage}_${tumorID}/${runID}_${stage}_wish_list_genes_expression.csv


		cd -
		set -e

		## update and zip files
		zipDir=/omics/odcf/analysis/OE0422_projects/promise/results_per_pid_2/zip.tmp/${runID_tumorID}.tmp
		rm -rf $zipDir
		rm -fr $outputDir/*.zip
		find $outputDir -type d -empty -print | xargs -I {} rm -rf {}


		find $outputDir -type d | sed "s,$outputDir,$zipDir,g" | xargs -I {} mkdir -p {}
		find $outputDir -type f -mtime -10 | xargs -I  {} basename {}
		find $outputDir -type f -mtime -10 | while read x; do
			desti=`echo $x | sed "s,$outputDir,$zipDir,g"`
			cp $x $desti
		done
		rm -rf $zipDir/HLA/In_silico/format_hla
		mkdir -p $zipDir/HLA/In_silico
		cat  $workDir/1_hla_type/format_hla > $zipDir/HLA/In_silico/${runID}_${stage}_hla
		# cp $outputDir/HLA/In_silico/format_hla $zipDir/HLA/In_silico/${runID}_${stage}_hla

		cd $zipDir
		zip -r $outputDir/${runID}.zip *
		cd -

		rm -fr $zipDir
		echo $outputDir
		ls -l --color $outputDir
		zipinfo $outputDir/${runID}.zip


		exit 0

	fi


	dir_omics_cgi=$outputDir/Mutation_analysis/CGI
	dir_omics_snv=$outputDir/Mutation_analysis/snv
	dir_omics_indel=$outputDir/Mutation_analysis/indel
	dir_omics_fusion=$outputDir/Mutation_analysis/fusion


	cp $workDir/1_hla_type/cgi/* $dir_omics_cgi
	cp $workDir/2_SNVs_based_neoepitope_prediction/*vcf $dir_omics_snv
	cp $workDir/4_indel_based_prediction/*vcf $dir_omics_indel
	rm -rf $dir_omics_snv/1.vcf
	rm -rf $dir_omics_snv/2.vcf
	rm -rf $dir_omics_indel/*tmp*vcf
	set +e
	cp $workDir/1_hla_type/format_hla $outputDir/HLA/In_silico/format_hla
	set -e

	file_fusion_mutation=$workDir/9_Fusion_gene_based_neoepitope_identification/2_arriba_result_RNAbamOpt/fusions.tsv_splitGenes
	if [ -f $file_fusion_mutation ];then
		cp $file_fusion_mutation $dir_omics_fusion/fusions_splitGenes.tsv
	fi

	echo '=== ===> rename_zip: start...'
	zipDir=$outputDir/zip
	# rm -fr ${runID}.zip
	rm -fr $outputDir/*.zip
	rm -fr $zipDir
	find $outputDir -type d -empty -print | xargs -I {} rm -rf {}
	rsync -a $outputDir $zipDir 2>&1 > /dev/null

	IFS=$'\n'
	files=($(find $zipDir -type f))
	unset IFS
	for i in ${files[*]}
	do
		mv $i `echo $i | sed "s;^\(.*/\)\([^/]\+\)$;\1${runID}_\2;"` 
	done

	cd $zipDir
	zip -r $outputDir/${runID_tumorID}.zip *
	cd -

	rm -fr $zipDir
	echo $outputDir
	ls -l --color $outputDir
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
	script=$dir_pipeline/scripts/HLA_typing/run_format_phlat.sh
	sh $script -p $input -o $output


}


function s1b_hla_sab () {
	echo '=== ===> s1b_hla_sab: start...'
	PID=hla_sab_$runID


	workDirHla=$workDir/1_hla_type
	script=$script_hla/Immunopipe.sh
	referenceDir=$script_hla/data


	if [ $hlaID = 'promise' ]; then
		sourceDir=/omics/odcf/project/OE0422/promise/sequencing/${dir_input_ws}/view-by-pid/$runID/buffy-coat*/paired/merged-alignment
	else
		sourceDir=$workDir/1_hla_type
	fi

	if [ -f ${workDirHla}/format_hla ] && [ $(wc -l < ${workDirHla}/format_hla) -gt 1 ];then
		cat $workDir/1_hla_type/format_hla
		cp $workDir/1_hla_type/format_hla $outputDir/HLA/In_silico/format_hla
		return 0
	else 
		if [ -d $outputDir/HLA ]; then
			format_hla=`find ${outputDir}/HLA -name 'format_hla'`
		fi

		if [[ -f $format_hla ]] && [ $(wc -l < $format_hla) -gt 1 ]; then
			cp $format_hla $workDir/1_hla_type
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

	cat $workDir/1_hla_type/format_hla

}




## 2) neoepitode predition with snv calling result ====


function f_vcfOnly_pathology () {
	rm -rf 1.vcf
	cat $dir_pipeline/vcf_header.txt >> 1.vcf
	cat *org | grep -v '^#' >> 1.vcf

	## input: 1.vcf; output: 2.vcf
	script=$dir_pipeline/prepare_vcf.r
	Rscript $script $workDir/2_SNVs_based_neoepitope_prediction $vcfOnly

	## output indel_tmp.vcf snv_somatic.vcf
	script=$dir_pipeline/vcfOnly_snv_indel_seperate.py
	python $script $workDir

	## reformat indel vcf: input: indel_tmp.vcf; output: indel_somatic.vcf
	script=$dir_pipeline/vcfOnly_indel_reformat.r
	Rscript $script $workDir/4_indel_based_prediction

}


function find_promise_hg38_org () {
	findx=`find /omics/odcf/analysis/OE0422_projects/promise/whole_genome_sequencing/_ALL_SNV_INDEL_2024_03 -name "*${runID}*${stage}*Impact*SNV*txt" | wc -l`
	if [ $findx == 1 ]; then
		org_snv=`find /omics/odcf/analysis/OE0422_projects/promise/whole_genome_sequencing/_ALL_SNV_INDEL_2024_03 -name "*${runID}*${stage}*Impact*SNV*txt"`
		export org_snv
	else
		echo ==== ERROR: no snv vcf found.
	fi

	findx=`find /omics/odcf/analysis/OE0422_projects/promise/whole_genome_sequencing/_ALL_SNV_INDEL_2024_03 -name "*${runID}*${stage}*Impact*INDEL*txt" | wc -l`
	if [ $findx == 1 ]; then
		org_ind=`find /omics/odcf/analysis/OE0422_projects/promise/whole_genome_sequencing/_ALL_SNV_INDEL_2024_03 -name "*${runID}*${stage}*Impact*INDEL*txt"`
		export org_ind
	else
		echo ==== ERROR: no indel vcf found.
	fi
}

if [[ ${vcfOnly,,} == 'promise' ]] ;then
	find_promise_hg38_org
fi

function f_vcfOnly_promise () {
	rm -rf ./netMHCpan4_1/*

	rm -rf 1.vcf
	cat $dir_pipeline/vcf_header.txt > 1.vcf
	cat $org_snv | awk -v OFS='\t' 'NR > 1 {print $1, $2, "NA", $3, $4}' >> 1.vcf

	## input: 1.vcf; output: 2.vcf
	script=$dir_pipeline/prepare_vcf.r
	Rscript $script $workDir/2_SNVs_based_neoepitope_prediction $vcfOnly

}


function f_vcfOnly_origin () {
	findx=`find ${workDir}/2_SNVs_based_neoepitope_prediction  -maxdepth 1 -name "*somatic*.vcf" | wc -l`
	if [ $findx == 1 ]; then
		vcf=${workDir}/2_SNVs_based_neoepitope_prediction/*somatic*.vcf
	else
		echo ==== ERROR: no snv file found.
	fi
	rm -rf 1.vcf
	cat $dir_pipeline/vcf_header.txt > $workDir/2_SNVs_based_neoepitope_prediction/1.vcf
	cat $vcf | awk 'NR > 1 {print $0}' >> 1.vcf

	## input: 1.vcf; output: 2.vcf
	script=$dir_pipeline/prepare_vcf.r
	Rscript $script $workDir/2_SNVs_based_neoepitope_prediction $vcfOnly

}




function s2_snv () {
	echo '=== ===> s2_snv: start...'

	cd ${workDir}/2_SNVs_based_neoepitope_prediction
	# rm -rf netMHCpan4_1 #debug;transcriptIDs
	if [ -f *org ] && [[ ${vcfOnly,,} == 'pathology' ]];then
		f_vcfOnly_pathology
		vcf=${workDir}/2_SNVs_based_neoepitope_prediction/snv_somatic.vcf
	fi

	if [[ ${vcfOnly,,} == 'promise' ]] ;then
		f_vcfOnly_promise
		vcf=${workDir}/2_SNVs_based_neoepitope_prediction/2.vcf
	fi

	if [[ ${vcfOnly,,} == 'origin' ]] ;then
		f_vcfOnly_origin
		vcf=${workDir}/2_SNVs_based_neoepitope_prediction/2.vcf
	fi

	if [ `cat $vcf | wc -l` -lt 2 ]; then
		echo 'Error: no valid lines in snv vcf file'
		echo 'snv_stop=1' >> $file_run_status 
		return
	fi
	cd -
	script=$dir_pipeline/update_netMHCpan/netMHCpan_4.1/main_run_neoepitope_snvs.sh
	format_hla=${workDir}/1_hla_type/format_hla

	if [ $merlength == 21 ];then
		extendLen=10 #21mer
	else
		extendLen=13 #27mer
	fi
	sh $script \
	-i ${netMHCpanID} \
	-v $vcf \
	-m $format_hla \
	-o ${workDir}/2_SNVs_based_neoepitope_prediction/${netMHCpanID} \
	-r $extendLen

}


## 3) If RNA-seq exist, excute 3.2, otherwise excute 3.1 ====

function add_RNAseq_vcfOnly_pathology () {
	input_sf=$workDir/2_SNVs_based_neoepitope_prediction/*sf
	expression=$workDir/3_add_expression/fpkm_tpm.featureCounts.tsv
	script_vcfOnly=$dir_pipeline/vcfOnly/reformat_sf.r

	input_s2=$workDir/2_SNVs_based_neoepitope_prediction/$netMHCpanID/results_MHCI_epitopes.tab_splitGenes
	output_s3=$workDir/3_add_expression/MHCI_epitopes_RNAseq_$netMHCpanID.tab
	Rscript $script_vcfOnly $workDir/3_add_expression $input_sf $expression $input_s2 $output_s3 


	input_s2=$workDir/2_SNVs_based_neoepitope_prediction/$netMHCpanID/results_MHCII_epitopes.tab_splitGenes
	output_s3=$workDir/3_add_expression/MHCII_epitopes_RNAseq.tab
	if [ -f $input_s2 ];then
		Rscript $script_vcfOnly $workDir/3_add_expression $input_sf $expression $input_s2 $output_s3 
	fi
}


function addRNAseq_origin () {
	RNA_bam=`ls ${workDir}/3_add_expression | grep -P '^(?!.*chimeric).*mdup.bam$'| xargs -I {} echo ${workDir}/3_add_expression/{}`
	expression=${workDir}/3_add_expression/*fpkm_tpm.featureCounts.tsv

	if [[ ${vcfOnly,,} == 'promise' ]];then
		vcf=$workDir/2_SNVs_based_neoepitope_prediction/2.vcf
		script_py=$dir_pipeline/hg38_to_hg19_vcf.py
		python $script_py $vcf $workDir/3_add_expression/3_19 $workDir/3_add_expression/3_38
		vcf=$workDir/3_add_expression/3
	fi

	if [ -f $inputMHCI ]; then
		sh $script -p $pipelineDir -v $vcf -b $RNA_bam -e $expression -n $inputMHCI -o $outputMHCI
		:
	else
		echo '>>>> >>>> WARNING: no file exist:'
		echo $inputMHCI
	fi

	if [ -f $inputMHCII ]; then
		sh $script -p $pipelineDir -v $vcf -b $RNA_bam -e $expression -n $inputMHCII -o $outputMHCII
		:
	else
		echo '>>>> >>>> WARNING: no file exist:'
		echo $inputMHCII
	fi

	if [ ! $hlaID = 'promise' ]; then
		## calculate colA3 exon6 FPKM ==== 
		script_col6a3=$dir_pipeline/exon6_col6a3_FPKM/FPKM_calculation.sh
		result=fpkm_COL6A3_exon6.txt

		bash $script_col6a3 $runID $result $workDir
		cp $workDir/3_add_expression/$result $outputDir/Gene_Expression
	fi
}


function s3_rna () {
	echo '=== ===> s3_rna: start...'
	if [[ $snv_stop == 1 ]]; then
		return
	fi

	inputMHCI=${workDir}/2_SNVs_based_neoepitope_prediction/$netMHCpanID/results_MHCI_epitopes.tab.hydro_splitGenes
	inputMHCII=${workDir}/2_SNVs_based_neoepitope_prediction/$netMHCpanID/results_MHCII_epitopes.tab.hydro_splitGenes

	rm -rf ${workDir}/3_add_expression/*tab
	rm -rf ${workDir}/3_add_expression/add*
	tab_RNA_tmp=0

	if [[ `echo $tcga | grep -o 'RNAseq' | wc -l` == 1 ]]; then
		outputMHCI=${workDir}/3_add_expression/MHCI_epitopes_RNAseq_$netMHCpanID.tab
		outputMHCII=${workDir}/3_add_expression/MHCII_epitopes_RNAseq.tab

		script=$dir_pipeline/scripts/add_expression/main_add_RNA.sh
		pipelineDir=$dir_pipeline/scripts/add_expression
		tab_RNA_tmp=1


		if [ $hlaID == 'promise' ]; then
			rnaDir=/omics/odcf/project/OE0422/promise/sequencing/rna_sequencing/view-by-pid/$runID/$tumorID-*
			RNA_bam=`find $rnaDir -name '*merged.mdup.bam'  | grep -P '^(?!.*chimeric).*merged.mdup.bam$'`
			expression=`find $rnaDir -name '*.fpkm_tpm.featureCounts.tsv'`

			ln -fs $expression $workDir/3_add_expression
			ln -fs $RNA_bam $workDir/3_add_expression
		fi

		if [[ ${vcfOnly,,} == 'pathology' ]];then
			add_RNAseq_vcfOnly_pathology
			:
		else
			vcf=${workDir}/2_SNVs_based_neoepitope_prediction/*somatic*.vcf
			vcf=${workDir}/2_SNVs_based_neoepitope_prediction/2.vcf
			addRNAseq_origin
			:
		fi
		script=$scriptsDir/check_wish_list_genes.r
		outputFile=$workDir/8_chose_neoepitode/wish_list_genes_expression.csv
		Rscript $script $expression $outputFile

	fi

	if [[ `echo $tcga | grep -oi 'tcga' | wc -l` == 1 ]] || [[ `echo $tcga | grep -oi 'target' | wc -l` == 1 ]]; then

		if [[ `echo $tcga | grep -oi 'tcga' | wc -l` == 1 ]];then
			tcgaID=`echo $tcga | sed 's/\(TCGA-[^_]\+\).*/\1/g'`
		else
			tcgaID=`echo $tcga | sed 's/\(TARGET-[^_]\+\).*/\1/g'`
		fi

		code=$dir_pipeline/scripts/add_expression/run_add_refExpression.R
		tcgaFile=$dir_tcga/${tcgaID}/${tcgaID}_expression_addName.tab

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
			:
		else
			echo ">>>> >>>> no file: $inputMHCII"
		fi

		if [ $tab_RNA_tmp -eq 1 ];then
			set +e
			mkdir -p  ${workDir}/3_add_expression/past
			mv ${workDir}/3_add_expression/MHCI*_epitopes_RNAseq*.tab $workDir/3_add_expression/past
			set -e
		fi
	fi


}

function s8a_filter () {
	if [[ $snv_stop == 1 ]]; then
		return
	fi

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
	if [[ $snv_stop == 1 ]]; then
		return
	fi

	echo '=== ===> s8b_xlsx_to_public: start...'
	if [ $hlaID = 'promise' ]; then
		outputDir_snv=${outputDir}/Epitope_prediction/snv_based/${stage}_${tumorID}
	else
		outputDir_snv=${outputDir}/Epitope_prediction/snv_based
	fi
	rm -rf $outputDir_snv
	mkdir -p $outputDir_snv
	/tbi/software/x86_64/R/R-3.4.0/el7/bin/Rscript ${scriptsDir}/convert_to_xlsx.r $workDir $outputDir_snv snv

	## convert wishList to xlsx file
	if [ `echo $tcga | grep 'RNAseq' | wc -l` == 1 ];then
		rm -rf $outputDir/Gene_Expression/wishList.xlsx
		inputFile=$workDir/8_chose_neoepitode/wish_list_genes_expression.csv
		/tbi/software/x86_64/R/R-3.4.0/el7/bin/Rscript $dir_pipeline/convert_to_xlsx.r $inputFile $outputDir/Gene_Expression wishList
		inputFile=$workDir/3_add_expression/*fpkm_tpm.featureCounts.tsv
		/tbi/software/x86_64/R/R-3.4.0/el7/bin/Rscript $dir_pipeline/convert_to_xlsx.r $inputFile $outputDir/Gene_Expression wishList
	fi
}


## INDEL BASED PREIDCTION ====

function f_vcfOnly_indel_promise () {
	# findx=`find /omics/odcf/analysis/OE0422_projects/promise/whole_genome_sequencing/_ALL_SNV_INDEL_2024_03 -name "*${runID}*${stage}*Impact*INDEL*txt" | wc -l`
	# if [ $findx == 1 ]; then
	# 	vcf_1=`find /omics/odcf/analysis/OE0422_projects/promise/whole_genome_sequencing/_ALL_SNV_INDEL_2024_03 -name "*${runID}*${stage}*Impact*INDEL*txt"`
	# else
	# 	echo ==== ERROR: no snv file found.
	# fi
	if [ `less $org_ind | wc -l ` -lt 2 ]; then
		touch $workDir/4_indel_based_prediction/result/1.vcf
		return
	fi

	vcf=$workDir/4_indel_based_prediction/result/1.vcf

	awk '{
		header="#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	CONTROL	TUMOR	DBSNP	1K_GENOMES	ExAC	EVS	LocalControlAF	ANNOVAR_FUNCTION	GENE	EXONIC_CLASSIFICATION	ANNOVAR_TRANSCRIPTS	SEGDUP	CYTOBAND	REPEAT_MASKER	DAC_BLACKLIST	DUKE_EXCLUDED	HISEQDEPTH	SELFCHAIN	MAPABILITY	SIMPLE_TANDEMREPEATS	CLASSIFICATION	CONFIDENCE	REGION_CONFIDENCE	Enhancers	CpGislands	TFBScons	ENCODE_DNASE	miRNAs_snoRNAs	miRBase18	COSMIC	miRNAtargets	CgiMountains	phastConsElem20bp	ENCODE_TFBS"

		if (NR==1) {
			print header
		} else {
			for(i=1; i<=43; ++i) {
				if (i<=2) {
					printf ("%s\t", $i)
				} else if (i==3) {
					printf ("%s\t", "NA")
				} else if (i==4) {
					printf ("%s\t", $3)
				} else if (i==5) {
					printf ("%s\t", $4)
				} else if (i==18) {
					printf ("NA\t")
				} else if (i >5 && i < 43) {
					printf ("\t")
				} else if (i == 43 ) {
					printf ("\n")
				} 
			};
		}

	}' ${org_ind} > $vcf

	## convert vcf to protein sequence by ANNOVAR
	script=$dir_pipeline/update_netMHCpan/netMHCpan_4.1/run_predict_protein__indel_from_fake_vcf.sh
	predictedProtein=$workDir/4_indel_based_prediction/predictedProtein.fa
	bash $script $vcf $predictedProtein

	## generate gene.tab for pathology patients
	script=$dir_pipeline/scripts/neoepitope_indels/pathology_get_geneTab.py
	python $script $indel_result $workDir/4_indel_based_prediction/predictedProtein.fa

}

function f_vcfOnly_indel_pathology () {
	## change vcf header
	vcf_1=`find ${workDir}/4_indel_based_prediction -maxdepth 1 -name 'indel_somatic.vcf'`
	if [ `less $vcf_1 | wc -l ` -lt 2 ]; then
		return
	fi
	vcf=$workDir/4_indel_based_prediction/result/1.vcf
	awk '{
		header="#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	CONTROL	TUMOR	DBSNP	1K_GENOMES	ExAC	EVS	LocalControlAF	ANNOVAR_FUNCTION	GENE	EXONIC_CLASSIFICATION	ANNOVAR_TRANSCRIPTS	SEGDUP	CYTOBAND	REPEAT_MASKER	DAC_BLACKLIST	DUKE_EXCLUDED	HISEQDEPTH	SELFCHAIN	MAPABILITY	SIMPLE_TANDEMREPEATS	CLASSIFICATION	CONFIDENCE	REGION_CONFIDENCE	Enhancers	CpGislands	TFBScons	ENCODE_DNASE	miRNAs_snoRNAs	miRBase18	COSMIC	miRNAtargets	CgiMountains	phastConsElem20bp	ENCODE_TFBS"

		if (NR==1) {
			print header
		} else {
			for(i=1; i<=18; ++i) {
				if (i<=11) {
					printf ("%s\t", $i)
				} else if (i >11 && i < 18) {
					printf ("\t")
				} else if (i == 18 ) {
					printf ("%s\n", $12)
				} 
			};
		}

	}' ${vcf_1} > $vcf
	## convert vcf to protein sequence by ANNOVAR
	script=$dir_pipeline/update_netMHCpan/netMHCpan_4.1/run_predict_protein__indel_from_fake_vcf.sh
	predictedProtein=$workDir/4_indel_based_prediction/predictedProtein.fa
	bash $script $vcf $predictedProtein

	## generate gene.tab for pathology patients
	script=$dir_pipeline/scripts/neoepitope_indels/pathology_get_geneTab.py
	python $script $indel_result $workDir/4_indel_based_prediction/predictedProtein.fa
}

function i4a_indel_predict () {
	echo '=== ===> i4a_indel_predict: start...'
	indel_result=$workDir/4_indel_based_prediction/result
	mkdir -p $indel_result

	## promise , hg19
	if [ $hlaID == 'promise' ] && [ $vcfOnly != 'promise' ]; then
		vcf=`find $workDir/4_indel_based_prediction -name 'indel*somatic_functional_indels_conf_8_to_10.vcf'| xargs realpath`
		indel_vcf_dir=$outputDir/Mutation_analysis/indel/$stage_$tumorID
		mkdir -p $indel_vcf_dir
		ln -fs $vcf $indel_vcf_dir
	fi

	## promise , hg38
	if [ ${vcfOnly,,} == 'promise' ] ;then
		find $workDir/4_indel_based_prediction -type f -delete
		f_vcfOnly_indel_promise
		vcf=$workDir/4_indel_based_prediction/result/1.vcf
	fi

	## pathology
	if [ ${vcfOnly,,} == 'pathology' ] ;then
		f_vcfOnly_indel_pathology
		vcf=$workDir/4_indel_based_prediction/result/1.vcf
	fi

	## origin
	if [ ${vcfOnly,,} == 'origin' ] ;then
		if [ ! `find ${workDir}/4_indel_based_prediction -maxdepth 1  -name '*somatic*vcf' | wc -l` -eq 1 ]; then
			echo 'no suitable indel vcf file found!'
			exit 1
		fi
		vcf=`find ${workDir}/4_indel_based_prediction -maxdepth 1 -name '*somatic*vcf'`
	fi

	## predict epitopes ====
	script=$dir_pipeline/scripts/neoepitope_indels/main_indels.sh
	hla=${workDir}/1_hla_type/format_hla

	if [ `less $vcf | wc -l ` -gt 1 ]; then
		if [ $merlength == 21 ];then
			sh $script -f $vcf -l 21 -a $hla -o $indel_result -v $vcfOnly
		else
			sh $script -f $vcf -l 27 -a $hla -o $indel_result -v $vcfOnly #27mer
		fi
	else
		echo 'no indel mutations.'
		echo 'indel_stop=1' >> $file_run_status
	fi
}


function i4b_indel_tsv () {
	source $file_run_status
	if [[ $indel_stop == 1 ]]; then
		return
	fi
	echo '=== ===> i4b_indel_tsv: start...'
	script=$dir_pipeline/4b_indel_table.r
	Rscript $script $workDir


	script=$dir_pipeline/update_netMHCpan/netMHCpan_4.1/hydro.py
	cd $workDir/4_indel_based_prediction
	IFS=$'\n'
	files=($(find . -maxdepth 1 -name 'indel*MHC*tsv'))
	unset IFS
	for i in ${files[*]}
	do
		python $script $i
		mv ${i}.hydro $i
	done
	cd -

}


function i4c_xlsx_to_public () {
	echo '=== ===> i4c_xlsx_to_public: start...'

	if [[ $indel_stop == 1 ]]; then
		return
		:
	fi

	if [ $hlaID = 'promise' ]; then
		outputDir_indel=${outputDir}/Epitope_prediction/indel_based/${stage}_${tumorID}
	else
		outputDir_indel=${outputDir}/Epitope_prediction/indel_based
	fi

	rm -rf $outputDir_indel
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


	STAR_INDEX_DIR=$dir_root/yanhong/git/arriba/STAR_index_GRCh37_ENSEMBL87
	annotation=$dir_root/yanhong/git/arriba/ENSEMBL87.gtf
	assembly=$dir_root/yanhong/git/arriba/GRCh37.fa


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
	if [ `echo $tcga |grep 'RNAseq' | wc -l` -eq 0 ] || [[ $vcfOnly == 'Yes' ]];then
		echo 'fusion_stop=1' >> $file_run_status
		return
	fi

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

	sourceDir=$dir_pipeline/fusion_arriba/arriba_reference_hg19
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
			-u \
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
	if [[ $fusion_stop == 1 ]]; then
		return
	fi
	echo '=== ===> f2_prepare_HLA: start...'
	script=$dir_pipeline/fusion_arriba/prepare_HLA.r
	Rscript $script $workDir $workDir92
}

## step 3: ====

function f3_neo_prediction () {
	if [[ $fusion_stop == 1 ]]; then
		return
	fi

	echo '=== ===> f3_neo_prediction: start...'
	if [ ! -f $workDir92/fusions.tsv ];then
		return
	fi
	## re-format fusion.tsv
	script=$dir_pipeline/fusion_arriba/splitRow.py
	python $script $workDir92

	## predict
	PID=${runID}_neoPrediction
	OUT_FUSION=$workDir9/3_neoPrediction_$dataType

	script=$script_fusion/Immunopipe.sh
	bam=$workDir92/Aligned.out.bam # not neccessary for fusion
	fusions_tsv=$workDir92/fusions.tsv_splitGenes


	$script --dna-bam $bam --fusions $fusions_tsv --fusion-confidence "medium high" --output $workDir92 --OUT_FUSION $OUT_FUSION --pid $PID
}

## step 4:
function f4_mer21 () {

	if [[ $fusion_stop == 1 ]]; then
		return
	fi

	echo '=== ===> f4_mer21: start...'
	if [ ! -d $workDir93 ];then
		echo ">>>> >>>> folder not exist: $workDir93 "
		return
	fi
	script=$dir_pipeline/fusion_arriba/mer21.r
	Rscript $script $workDir93
}


## step 5: ====
function f5_to_xlsx () {
	if [[ $fusion_stop == 1 ]] ; then
		return
	fi

	echo '=== ===> f5_to_xlsx: start...'
	if [[ `echo $tcga |grep 'RNAseq' | wc -l` != 1 ]];then
		return
	fi
	script=$dir_pipeline/fusion_arriba/convert_to_xlsx_arribaFusion.r

	if [ $hlaID = promise ]; then
		outputDir_fusion=$outputDir/Epitope_prediction/fusion_genes/${stage}_${tumorID}/
	else
		outputDir_fusion=$outputDir/Epitope_prediction/fusion_genes
	fi

	rm -rf $outputDir_fusion
	mkdir -p $outputDir_fusion
	/tbi/software/x86_64/R/R-3.4.0/el7/bin/Rscript $script $workDir $outputDir_fusion $dataType


	# set +e
	# # mkdir -p  ${dir_fp}/${runID}_$tumorID
	# # cp $outputDir_fusion/*xlsx ${dir_fp}/${runID}_$tumorID
	# dir_fp=$dir_root/yanhong/all_in_one_pipeline_collection/mhc4.1/promise/batch_3_update/batch.3.1/fusion_promise.mer
	# # cp $workDir93/*MHCI_*mer  ${dir_fp}/${runID}_${tumorID}_MHCI.tsv
	# # cp $workDir93/*MHCII_*mer ${dir_fp}/${runID}_${tumorID}_MHCII.tsv
	# cp $workDir92/fusions.tsv ${dir_fp}/${runID}_${tumorID}_fusionGene.tsv
	# set -e

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


	if [ ! -z $loh_stop ];then
		if [ $loh_stop == 1 ];then
			return
		fi
	fi
	if [ $vcfOnly == 'Yes' ];then
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
	if [ $(find $workDir/5_LOHHLA/example-out -name '*xls' | wc -l) -gt 0 ];then
		cp $workDir/5_LOHHLA/example-out/*xls $workDir/5_LOHHLA/example-out/Figures/*pdf $outputDir_loh
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
	rm -rf $workDir/1_hla_type/mutation.csv $workDir/1_hla_type/mutation_germline.csv

	if [ ${vcfOnly,,} == 'promise' ] ;then
		# less $org_snv | awk -F '\t' '(NR>1) { print $1, $2, $3, $4, $10} NR == 10 {exit}' >> mutation.csv
		# less $org_ind | awk -F '\t' '(NR>1) { print $1, $2, $3, $4, $10} NR == 10 {exit}' >> mutation.csv
		less $org_snv | awk -F '\t' '(NR>1) { print $1, $2, $3, $4, $10}' >> mutation.csv
		less $org_ind | awk -F '\t' '(NR>1) { print $1, $2, $3, $4, $10}' >> mutation.csv
	else 
		find $workDir/2* -name '*somatic*vcf' | xargs -I {} awk -F '\t' '(NR>1) { print $1, $2, $4, $5, $17}' {} >> mutation.csv
		find $workDir/4* -name '*somatic*vcf' | xargs -I {} awk -F '\t' '(NR>1) { print $1, $2, $4, $5, $18}' {} >> mutation.csv
		find $workDir/2* -name '*germline*vcf' | xargs -I {} awk -F '\t' '(NR>1) { print $1, $2, $4, $5, $17}' {} >> mutation_germline.csv
		find $workDir/4* -name '*germline*vcf' | xargs -I {} awk -F '\t' '(NR>1) { print $1, $2, $4, $5, $18}' {} >> mutation_germline.csv
	fi

	echo 'chr pos ref alt gene' > $workDir/1_hla_type/tmp.tsv
	cat $workDir/1_hla_type/mutation.csv  >> tmp.tsv
	sed 's/[[:space:]]\{1,\}/\t/g' tmp.tsv > mutation.tsv
	rm -rf tmp.tsv

	## run cgi and download data
	script=$dir_pipeline/cgi/cgi_api.py
	cgi_try=0
	while [[ `find $workDir/1_hla_type -name 'cgi*log' | wc -l ` == 0 ]] || \
		  [ `cat $workDir/1_hla_type/cgi*log | grep 'Analysis done' | wc -l` -lt 1 ]; do
		rm -rf cgi*log

		cgi_try=$((cgi_try+1))
		if [ $cgi_try -gt 1 ]; then
			echo cgi_try: $cgi_try
			sleep 100
		fi
		if [[ $cgi_try == 5 ]]; then
			echo 'Error >>> >>>> cgi-max;  too much trying on cgi'
			exit 1
		fi

		python $script $workDir/1_hla_type $cgiTumorType
	done
	cd -
}



function c2 () {
	## merge to snv/indel/fusion results
	echo '=== ===> c2_cgi_merge2file: start...'
	script=$dir_pipeline/cgi/merge_cgi_wishList.r
	echo $org_snv
	echo $org_ind
	Rscript $script $workDir 
}

## EXCUTION (step order sensitive) =======================================================================================

OLDIFS=$IFS
IFS='-' read -r -a array1 <<< $steps

for step in "${array1[@]}"
do
	## prepare_folder ==== ====
	source $file_run_status

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
