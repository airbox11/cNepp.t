#! /bin/bash

usage() { echo "Usage: $0 -f fastaFile -n fusionName -m hla_type -o outputDir
    -f: Path to fusion fasta file 
    -n: Fusion name, e.g. NYB_NFIB
    -m: Path to HLA type result file
    -o: Path to result folder
>>> Note: Absolute path needed" 1>&2; exit 1; }

while getopts ":f:n:m:o:" o; do
    case "${o}" in
        f)
            fastaFile=${OPTARG}
            ;;
        n)
            fusionName=${OPTARG}
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

if [ ! -f "${fastaFile}" ] || [ -z "${fusionName}" ] || [ ! -f "${HLA}" ] || [ -z "${OUTPUT_DIR}" ]; then
    usage
fi

# create a new folder 
if [ ! -d "${OUTPUT_DIR}" ];then
    mkdir -p ${OUTPUT_DIR}
fi
cd ${OUTPUT_DIR}


netMHCRun=/icgc/dkfzlsdf/analysis/D120/scripts/fusion/epitopePrediction/run_netMHCpan_for_single_peptides.sh
visRun=/icgc/dkfzlsdf/analysis/D120/scripts/fusion/epitopePrediction/run_single_visualization.R

sh $netMHCRun $fastaFile $fusionName ${OUTPUT_DIR} ${HLA}


cat ${OUTPUT_DIR}/*${fusionName} > ${OUTPUT_DIR}/fusion_netMHC_${fusionName}.results

Rscript $visRun $fastaFile ${OUTPUT_DIR}/fusion_netMHC_${fusionName}.results ${OUTPUT_DIR}/fusion_${fusionName}_visualization
