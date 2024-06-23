#!/usr/bin/env bash
set -e

lohDir=/omics/groups/OE0422/internal/yanhong/git/20220306_hlaloh/lohhla
# workDir=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/JS53
workDir=$1
workDir5=$workDir/5_LOHHLA

runID=$2
tumorID=$3

mkdir -p $workDir5/data
mkdir -p $workDir5/example-file/bam
mkdir -p $workDir5/example-out






function link_bam {

	promise_wgs=/omics/odcf/project/OE0422/promise/sequencing/whole_genome_sequencing/view-by-pid/${runID}
	if [[ ! -z $tumorID ]];then
		# tumorID_index=`echo $tumorID | sed 's/tumor//g'`

		# buffcoatID_1=buffy-coat${tumorID_index}-01
		# buffcoatID_2=buffy-coat11-01

		# if [ ! -d "${promise_wgs}/${buffcoatID_1}" ];then
		# 	if [ -d "${promise_wgs}/$buffcoatID_2" ];then
		# 		buffcoatID_1=$buffcoatID_2
		# 	else
		# 		ls -l $promise_wgs
		# 		read -p "Input Selection for buffy-coat for ${runID}_${tumorID}:" buffcoatID_1
		# 	fi
		# fi


		bam=`find ${promise_wgs}/buffy-coat*/paired/merged-alignment -maxdepth 1 -name '*bam'`
		ln -fs ${bam} $workDir5/example-file/bam/example_BS_GL_sorted.bam
		ln -fs ${bam}.bai $workDir5/example-file/bam/example_BS_GL_sorted.bam.bai

		bam_tumor=`find ${promise_wgs}/${tumorID}*/paired/merged-alignment -maxdepth 1 -name '*bam'`
		ln -fs ${bam_tumor} $workDir5/example-file/bam/example_tumor_sorted.bam
		ln -fs ${bam_tumor}.bai $workDir5/example-file/bam/example_tumor_sorted.bam.bai

	else
		bam=`realpath $workDir/1_hla_type/*bam`
		ln -fs ${bam} $workDir5/example-file/bam/example_BS_GL_sorted.bam
		ln -fs ${bam}.bai $workDir5/example-file/bam/example_BS_GL_sorted.bam.bai

		dir_tumor=`echo $bam | sed 's;control;tumor;' | sed 's;\(.*\)/[^/]\+;\1;g'`

		cd $dir_tumor
		lnqk *bam $workDir5/example-file/bam/example_tumor_sorted.bam
		lnqk *bam.bai $workDir5/example-file/bam/example_tumor_sorted.bam.bai
		cd -
	fi

	if [ ! -f $workDir5/example-file/bam/example_BS_GL_sorted.bam ] | [ ! -f $workDir5/example-file/bam/example_tumor_sorted.bam ]; then
		echo '>>>> >>>> Error: no file exist for loh analysis.'
		echo 'loh_stop=1' >> $file_run_status
		exit 0
	fi

}


function hla_data {
	echo '====================== get hlas'

	hlas=$workDir5/data/hlas
	HLAs=$workDir5/data/HLAs
	cat $workDir/1_hla_type/format_hla | grep HLA |
		sed 's/[A-Z]/\L&/g' |
		sed 's/-/_/g' |
		sed 's/\(hla_.\)\(..\):\(..\)/\1_\2_\3/g' |
		tee $hlas

	cat $workDir/1_hla_type/format_hla | grep HLA |
		sed 's/\(HLA-.\)\(.*\)/\1\*\2/g' |
		tee $HLAs

	echo '====================== '
	fa=$workDir5/data/example.patient.hlaFasta.fa
	rm -f $fa
	hla_dat=${lohDir}/data/hla.dat

	count=0
	while IFS= read -r HLA
	do
		((count+=1))
		hla=`sed -n "${count}p" < $hlas`
		echo -e ">${hla}" | tee -a $fa
		echo $HLA

		hla2=`echo  "${HLA}" |sed 's/\*/\\\*/g'`
		n1=`grep -nP "Human MHC; HLA; Class I; HLA.*Allele; ${hla2}" $hla_dat | head -1 | cut -f1 -d:`

		sed -n "${n1},$ p" $hla_dat | while IFS= read -r line;
		do
			printf '%s\n' "$line" | grep -P '^\s{5}.*' | sed 's/ //g' | sed 's/[0-9]//g' >> $fa
			regex_end='^//'
			if [[ $line =~ $regex_end ]]; then
				break
			fi
		done

		if (($count==10));then
			exit 0
		fi
	done < $HLAs
}

function run_loh {
	script=${lohDir}/example.sh
	bash $script $workDir5 $lohDir
}


## MAIN ===================================
echo '=================================== start loh '

# n1=`ls  --color $workDir/5_LOHHLA/example-out/flagstat | wc -l`
# if (($n1!=2));then
# 	echo $workDir
# fi

# link_bam
hla_data
run_loh

# xls=`find $workDir5/example-out -maxdepth 1 -name '*xls'`
# if [[ -z $xls ]]; then
# fi

echo '=================================== loh Done ]]]]]]]]]]]]]]]]]]]]'
