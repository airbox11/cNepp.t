#! /bin/bash

# Run netMHC for HLA and DRB 
# HLA format
<<comment1
DQA1_0201
DQA1_0201
DQB1_0202
DQB1_0202
DRB1_0403
DRB1_0701
HLA-A02:01
HLA-A02:01
HLA-B44:02
HLA-B44:02
HLA-C05:01
HLA-C05:01
comment1

# check if HLA file available
if [ ! -f ${HLA} ];then
        printf "File %s is not found. Exiting...\n" $HLA
	exit 2
fi

# running MHCI
awk -v fa=$peptideFasta -v OUTPUT_DIR=$OUTPUT_DIR -v netMHCI=$netMHCpan -v pepType=$pepType '{if (substr($1,1,3)=="HLA") print netMHCI " -a " $1 " -BA -l 8,9,10,11 -f " fa " > " OUTPUT_DIR "/netMHCI_" $1 "_" pepType}' ${HLA} | sort | uniq | xargs -P 8 -n 1 -I {} sh -c {}

# running MHCII, only DRB alleles
cat ${HLA} | grep "DRB" |  awk -v fa=$peptideFasta -v OUTPUT_DIR=$OUTPUT_DIR -v netMHCII=$netMHCIIpan -v pepType=$pepType '{print netMHCII " -a " $1 " -f " fa " > " OUTPUT_DIR "/netMHCII_"$1 "_" pepType}' | sort | uniq | xargs -P 8 -n 1 -I {} sh -c {}


# print command lines, without execution 
awk -v fa=$peptideFasta -v OUTPUT_DIR=$OUTPUT_DIR -v netMHCI=$netMHCpan -v pepType=$pepType '{if (substr($1,1,3)=="HLA") print netMHCI " -a " $1 " -BA -l 8,9,10,11 -f " fa " > " OUTPUT_DIR "/netMHCI_" $1 "_" pepType}' ${HLA} | sort | uniq 

cat ${HLA} | grep "DRB" |  awk -v fa=$peptideFasta -v OUTPUT_DIR=$OUTPUT_DIR -v netMHCII=$netMHCIIpan -v pepType=$pepType '{print netMHCII " -a " $1 " -f " fa " > " OUTPUT_DIR "/netMHCII_"$1 "_" pepType}'  | sort | uniq 
