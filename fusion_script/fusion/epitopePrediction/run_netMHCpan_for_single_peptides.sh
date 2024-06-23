#! /bin/bash

fastaFile=$1
fusionName=$2  # suffix for output files
OUTPUT_DIR=$3
HLA_phlat=$4

cd ${OUTPUT_DIR}

HLA_phlat=${HLA_phlat}

netMHCIIpan=/icgc/dkfzlsdf/analysis/G200/immuno/tools/netMHCIIpan-3.1/netMHCIIpan
netMHCpan=/icgc/dkfzlsdf/analysis/G200/immuno/tools/netMHCpan-3.0/netMHCpan

# MHC-I epitope
awk -v fa=$fastaFile -v OUTPUT_DIR=$OUTPUT_DIR -v netMHCI=$netMHCpan -v fusionName=$fusionName '{if (substr($1,1,3)=="HLA") print netMHCI " -a " $1 " -l 8,9,10,11 -f " fa " > " OUTPUT_DIR "/netMHCI_" $1 "_" fusionName}' ${HLA_phlat} | xargs -P 8 -n 1 -I {} sh -c {}

# running MHCII -> DRB 
cat ${HLA_phlat} | grep "DRB" |  awk -v fa=$fastaFile -v OUTPUT_DIR=$OUTPUT_DIR -v netMHCII=$netMHCIIpan -v fusionName=$fusionName '{print netMHCII " -a " $1 " -f " fa " > " OUTPUT_DIR "/netMHCII_"$1 "_" fusionName }' |  xargs -P 8 -n 1 -I {} sh -c {}

