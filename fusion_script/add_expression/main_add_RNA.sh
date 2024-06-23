#! /bin/bash

## This pipeline is to add mutant expression based on RNA-seq profilling

usage() { echo "Usage: $0 -p pipelineFolder -v vcf -b BAM -e expressionProfiling -n neoepitopeResult -o outputDir
    -p: Path to pipeline folder
    -v: Path to VCF file (snvs)
    -b: Path to RNA bam file
    -e: Path to RNA-seq expression file (TPM, FPKM)
    -n: Path to neoepitope prediction result file
    -o: Path to your output file

>>> Note: Absolute path needed" 1>&2; exit 1; }

while getopts ":p:v:b:e:n:o:" o; do
    case "${o}" in
        p)
            pipelineDir=${OPTARG}
            ;;
        v)
            VCF=${OPTARG}
            ;;
        b)
            bam_RNA=${OPTARG}
            ;;
        e)
            expressionProfiling=${OPTARG}
            ;;
        n)
            neoepitope=${OPTARG}
            ;;
		o)
			outFile=${OPTARG}
			;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${pipelineDir}" ] || [ ! -f "${VCF}" ] || [ ! -f "${bam_RNA}" ] || [ ! -f "${expressionProfiling}" ] || [ ! -f "${neoepitope}" ] ; then
    usage
fi

module load samtools/1.5
module load R/3.4.0
module load perl/5.24.1

# predefine parameter
reference_input=/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/1KGRef/hs37d5.fa

## temporary folder for inter-files
## default is in /icgc/dkfzlsdf/analysis/D120/tmp/
## Update it if needed
tmpfile1=$(mktemp /icgc/dkfzlsdf/analysis/D120/tmp/addRNA1-script.XXXXXX)
tmpfile2=$(mktemp /icgc/dkfzlsdf/analysis/D120/tmp/addRNA2-script.XXXXXX)
tmpfile3=$(mktemp /icgc/dkfzlsdf/analysis/D120/tmp/addRNA3-script.XXXXXX)

sh $pipelineDir/1_run_allele_count.sh $VCF ${bam_RNA} $tmpfile1 $reference_input
echo "sh $pipelineDir/1_run_allele_count.sh $VCF ${bam_RNA} $tmpfile1 $reference_input"
Rscript $pipelineDir/2_merge_exom_inRNA.R  $tmpfile1 $VCF $tmpfile2
echo "Rscript $pipelineDir/2_merge_exom_inRNA.R  $tmpfile1 $VCF $tmpfile2"
sh $pipelineDir/3_run_allele_frequence.sh $tmpfile2 > $tmpfile3
echo "sh $pipelineDir/3_run_allele_frequence.sh $tmpfile2 > $tmpfile3"
Rscript $pipelineDir/4_bind_neoepitope_mutExp_hipo.R $tmpfile3 $neoepitope $expressionProfiling $outFile
echo "Rscript $pipelineDir/4_bind_neoepitope_mutExp_hipo.R $tmpfile3 $neoepitope $expressionProfiling $outFile"

rm $tmpfile1
rm $tmpfile2
rm $tmpfile3

##### noted there are same gene symbol in DKFZ pipeline, ENSG00000182378.8 from chromosome X and Y, to be explained
# GENE    TPM   FPKM
# 54851 PLCXD1 2.6394 1.7162
# 57243 PLCXD1      0      0

