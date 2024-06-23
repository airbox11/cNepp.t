#! /bin/bash

# This script is to do indels-based neoepitope prediction

usage() { echo "Usage: $0 -f indelsFile -l lenAA -a format_hla_phlat -o outputDir
	-f: Path to indels file generating from ODCF pipeline
	-l: Length of predicted neopeptides # e.g. 29 
	-a: Path to file containing HLA class I and DRB allele
	-o: Path to output folder # must be empty
>>>> Note: Absolute path be required" 1>&2; exit 1; }

while getopts ":f:l:a:o:" o; do
	case "${o}" in
	f)
		indelsFile=${OPTARG}
		;;
	l)
		lenAA=${OPTARG}
        ;;
	a)
	    format_hla_phlat=${OPTARG}
	    ;;
	o)
	    outputDir=${OPTARG}
	    ;;
	*)
    	usage
		;;
    esac
done
shift $((OPTIND-1))

if [ ! -f "${indelsFile}" ] || [ ! -f "${format_hla_phlat}" ] || [ -z "${lenAA}" ] || [ -z "${outputDir}" ] ; then
    usage
fi

## create folder if output folder not available
if [ ! -d $outputDir ];then
	mkdir -p $outputDir
else
	if [ ! -z "$(ls -A $outputDir)" ];then
		echo "$outputDir must be empty"
		exit
	fi
fi

pipelineDir=/icgc/dkfzlsdf/analysis/D120/scripts/neoepitope_indels/

runIndels=${pipelineDir}/predict_peptides_indels.R
netMHCbin=${pipelineDir}/indels_netmhcpan.sh
visualBin=${pipelineDir}/run_visualization.R
genHTML=${pipelineDir}/function_generate_html.Rmd
runHTML=${pipelineDir}/run_generate_html.R

module load R/3.4.0
module load samtools/1.5

# to to working directory
cd $outputDir

#predict indel peptides
## Find ANNOVAR_TRANSCRIPTS column position
myColPos=`head -n 1 $indelsFile | awk '{for(i=1;i<=NF;i++){if ($i=="ANNOVAR_TRANSCRIPTS") print i}}'`
# Be noted, it may filter unexpectly because annotion contain "synonymous"
candidates=`grep -v "CHROM\|splicing\|stopgain\|synonymous" $indelsFile | wc -l` 
if [ "$candidates" != 0 ]
then
	grep -v "CHROM\|splicing\|stopgain\|synonymous" $indelsFile | awk -v myColPos=$myColPos -v aaLen=$lenAA -v run_predict_indels=$runIndels 'BEGIN {FS="\t"; OFS="\t"} {split($myColPos,a,":");print "Rscript " run_predict_indels " " $myColPos " " aaLen " " a[1] ".tab" }' | xargs -I {} sh -c {}
else
	echo ">>>> Please check your script in filtering indel mutation step"
	grep -v "CHROM\|splicing\|stopgain" $indelsFile | awk -v myColPos=$myColPos -v aaLen=$lenAA -v run_predict_indels=$runIndels 'BEGIN {FS="\t"; OFS="\t"} {split($myColPos,a,":");print "Rscript " run_predict_indels " " $myColPos " " aaLen " " a[1] ".tab" }' | xargs -I {} sh -c {}
fi

# for each gene
for file in `find ${outputDir}/ -type f -name "*tab"`
#for file in `ls $outputDir/*tab`;
do
	# if file is empty or contain only spaces, TO optimize
	if [[ -s $file && $(wc -l < $file) -ge 2 ]];then

		# run netMHCpan binding affinity prediction
		sh $netMHCbin $format_hla_phlat $file $outputDir
		# extract gene name  from file path
		geneName=`echo $file | awk '{split($1,a,"/");print a[length(a)]}' | sed s/.tab//g`

		# for new version
		grep "<=" $outputDir/${geneName}_neoepitope_indel_mutation > $outputDir/${geneName}_binding

		# if no candidates found, _binding file empty, exit
		if [ ! -s $outputDir/${geneName}_binding ];then
			echo ">>>> $outputDir/${geneName}_binding is empty"
			continue
		fi
		Rscript $visualBin $file $outputDir/${geneName}_binding $outputDir/${geneName}_visualization

		# for old version
		grep "<=" $outputDir/version1/${geneName}_neoepitope_indel_mutation_v1 > $outputDir/version1/${geneName}_binding_v1
		Rscript $visualBin $file $outputDir/version1/${geneName}_binding_v1 $outputDir/version1/${geneName}_visualization_v1
	else
		echo "$file is empty"		
fi
done

# sumarize results
cat  $outputDir/*_visualization > $outputDir/indels_sum_visualization
cat  $outputDir/version1/*_visualization_v1 > $outputDir/version1/indels_sum_visualization_v1

# Generate HTML results
Rscript $runHTML $genHTML $outputDir/indels_sum_visualization $outputDir/indels_sum_visualization.html
Rscript $runHTML $genHTML $outputDir/version1/indels_sum_visualization_v1  $outputDir/version1/indels_sum_visualization_v1.html




