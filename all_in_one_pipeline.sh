#!/usr/bin/bash
set -e

function usage() { 
	echo "
	Usage: $0 -r runID -s steps [e, p, t, d, h, c, b]

	) -r: runID: 
	] -e: netMHCpanID:
		(netMHCpan4_1), netMHCstabpan, both
	] -p: patientID:
		($ runID)
	] -t: tcga: 
		(RNAseq), TCGA-READ
	) -s: steps: 
		s0: prepare folder
		s1a_HLA
		s1b_HLA

		s2_snv
		s3_add_expression
		s8a_filter
		s8b_xlsx_to_public (snv)
			format names
			wish_list
		snv: (s2,s3,s8a)

		i4a_indel_predict
		i4b_indel_tsv
		i4c_xlsx_to_public
		indel: (i4a, i4b)

		f1a: run_star_arriba , (f1b: run_arriba)
		f2: prepare_HLA
		f3: neo_prediction
		f4: mer21
		f5: to_xlsx
		fusion: (f1b, f2, f3, f4)

		l1_run
		l2_cp_file
		loh: (l1, l2)

		c1
		c2
			add redundant columns
		cgi :(c1, c2)
		(snv, indel, fusion)

		xlsx: convert as xlsx and copy to final folder
			(l1, l2, c1, c2, s8b, i4c, f5, nz)

		nz: change file name
			add format_hla

	] -h: hlaID
		(NCT_IP): NCT internal patient, promise
	] -l: mer length
		(21), 27
	] -c: (CANCER)
		cgi tumor type
	] -b: only vcf processing
		(origin), pathology, promise
	] -w: wgs or wes data
		(wes), wgs

	] -d: debug: 
		debug, empty

	# dependency:
	s1b - l1
	s8a - l2 - c2
	i4a - c2

	main steps:
		s0
		s1b
		snv/indel/fusion
		xlsx


	" 1>&2
	exit 1
	}

while getopts :r:e:p:t:s:d:h:a:c:b:w:o:l: opt; do
	case $opt in
		r)
			runID=${OPTARG}
			;;
		e)
			netMHCpanID=${OPTARG}
			;;
		p)
			patientID=${OPTARG}
			;;
		t)
			tcga=${OPTARG}
			;;
		s)
			steps=${OPTARG}
			;;
		d)
			debug=${OPTARG}
			;;
		h)
			hlaID=${OPTARG}
			;;
		a)
			dataType=${OPTARG}
			;;
		c)
			cgiTumorType=${OPTARG}
			;;
		b)
			vcfOnly=${OPTARG}
			;;
		w)
			wg=${OPTARG}
			;;
		o)
			logDir=${OPTARG}
			;;
		l)
			merlength=${OPTARG}
			;;
		*)
			usage
			;;
	esac
done
shift $((OPTIND-1))

if [ -z "${runID}" ]; then
	usage
fi

## default value 
if [ -z $netMHCpanID ]; then
	netMHCpanID=netMHCpan4_1
fi

if [ -z $patientID ]; then
	patientID=$runID
fi

if [ -z $tcga ]; then
	tcga=RNAseq
fi

if [ -z $hlaID ]; then
	hlaID=NCT_IP
fi

if [ -z $dataType ]; then
	dataType=RNAbamOpt
fi

if [ -z $cgiTumorType ]; then
	cgiTumorType=CANCER
fi

if [ -z $merlength ]; then
	merlength=21
fi

if [ -z $vcfOnly ]; then
	vcfOnly=origin
fi

if [ -z $wg ]; then
	wg=wes
fi



# script=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/all_in_one_pipeline_0.sh
script=/omics/odcf/analysis/OE0422_projects/Immuno-Patients-NCT/yanhong/all_in_one_pipeline/all_in_one_pipeline_0_update.sh

mkdir -p ~/log/tmp_run

if [ -z $debug ]; then
	PID="${runID}-${hlaID}-${netMHCpanID}-${tcga}-${steps}"
	if [ -z $logDir ];then
		logDir=/omics/odcf/analysis/OE0422_projects/Immuno-Patients-NCT/yanhong/log
	fi
	mkdir -p $logDir
	outputReport=$logDir/${PID}.log
	rm -f $outputReport
	env="\
		all,\
		runID=$runID,\
		netMHCpanID=$netMHCpanID,\
		patientID=$patientID,\
		tcga=$tcga,\
		steps=$steps,\
		hlaID=$hlaID,\
		dataType=$dataType,\
		cgiTumorType=$cgiTumorType,\
		vcfOnly=$vcfOnly,\
		wg=$wg\
		merlength=$merlength\
		"
	script_tmp=$logDir/tmp_run/$(date +"%Y_%m_%d__%H_%M_%S")_${PID}
	echo $script_tmp
	rm -rf $logD/tmp_run/*$runID*
	cp $script $script_tmp

	# bsub -q short -u y.lyu@dkfz-heidelberg.de -o $outputReport -r -R "rusage[mem=30G]" -J $PID -env "$env" -n 12 $script_tmp
	bsub -u y.lyu@dkfz-heidelberg.de -o $outputReport -r -R "rusage[mem=30G]" -J $PID -W 50:00 -env "$env" -n 12 $script_tmp
else
	echo ">>>>>>>> DEBUG >>>>>>>> $runID"

	export runID
	export netMHCpanID
	export patientID
	export tcga
	export steps
	export debug
	export hlaID
	export dataType
	export cgiTumorType
	export vcfOnly
	export wg
	export merlength

	bash $script 
fi

