#! /bin/bash

# Run netMHC for HLA and DRB 

format_hla_phlat=$1
peptideFile=$2
OUTPUT_DIR=$3


ID=`grep -v "geneName" $peptideFile | cut -f 1 | sort | uniq`
grep -v "geneName" $peptideFile | cut -f 1,4 | sort  | uniq | awk '{if (NR==1) print ">"$1 "\n" $2; if (NR!=1) print ">" $1 "_" NR "\n" $2;}' > $OUTPUT_DIR/${ID}_targeted.fa
grep -v "geneName" $peptideFile | cut -f 1,5 | sort  | uniq | awk '{if (NR==1) print ">"$1 "\n" $2; if (NR!=1) print ">" $1 "_" NR "\n" $2;}' > $OUTPUT_DIR/${ID}_refPeptide.fa


ANNOVAR=/icgc/ngs_share/annovar/annovar_Feb2016
netMHCIIpan=/icgc/dkfzlsdf/analysis/G200/immuno/tools/netMHCIIpan-3.1/netMHCIIpan
netMHCpan=/icgc/dkfzlsdf/analysis/G200/immuno/tools/netMHCpan-3.0/netMHCpan
phlatRelease=/icgc/dkfzlsdf/analysis/G200/immuno/tools/phlat-release
netMHCIIpan32=/lsdf/dkfz/analysis/G200/immuno/tools/netMHCIIpan-3.2/netMHCIIpan
netMHCpan40=/lsdf/dkfz/analysis/G200/immuno/tools/netMHCpan-4.0/netMHCpan

## Folder for first version results, netMHCpan-3.0, netMHCIIpan-3.1
OUTPUT_DIR_V1=$OUTPUT_DIR/version1

if [ ! -d "$OUTPUT_DIR_V1" ];then
	mkdir $OUTPUT_DIR_V1
fi

if [ ! -f ${format_hla_phlat} ];then
	printf "File %s is not found. Exiting...\n" $format_hla_phlat
	exit 2
fi

targeted=$OUTPUT_DIR/${ID}_targeted.fa
refPeptide=$OUTPUT_DIR/${ID}_refPeptide.fa

### Old version netMHCpan-3.0, netMHCIIpan-3.1
# running MHCI
awk -v fa=$targeted -v OUTPUT_DIR_V1=$OUTPUT_DIR_V1 -v netMHCI=$netMHCpan -v pepType="mutated" '{if (substr($1,1,3)=="HLA") print netMHCI " -a " $1 " -l 8,9,10,11 -f " fa " > " OUTPUT_DIR_V1 "/netMHCI_" $1 "_" pepType}' ${format_hla_phlat} | xargs -P 8 -n 1 -I {} sh -c {}

awk -v fa=$refPeptide -v OUTPUT_DIR_V1=$OUTPUT_DIR_V1 -v netMHCI=$netMHCpan -v pepType="ref" '{if (substr($1,1,3)=="HLA") print netMHCI " -a " $1 " -l 8,9,10,11 -f " fa " > " OUTPUT_DIR_V1 "/netMHCI_" $1 "_" pepType}' ${format_hla_phlat} | xargs -P 8 -n 1 -I {} sh -c {}

# running MHCII -> DRB 
cat ${format_hla_phlat} | grep "DRB" |  awk -v fa=$targeted -v OUTPUT_DIR_V1=$OUTPUT_DIR_V1 -v netMHCII=$netMHCIIpan -v pepType="mutated" '{print netMHCII " -a " $1 " -f " fa " > " OUTPUT_DIR_V1 "/netMHCII_"$1 "_" pepType}' |  xargs -P 8 -n 1 -I {} sh -c {}
cat ${format_hla_phlat} | grep "DRB" |  awk -v fa=$refPeptide -v OUTPUT_DIR_V1=$OUTPUT_DIR_V1 -v netMHCII=$netMHCIIpan -v pepType="ref" '{print netMHCII " -a " $1 " -f " fa " > " OUTPUT_DIR_V1 "/netMHCII_"$1 "_" pepType}' |  xargs -P 8 -n 1 -I {} sh -c {}

cat $OUTPUT_DIR_V1/*mutated > $OUTPUT_DIR_V1/${ID}_neoepitope_indel_mutation_v1
cat $OUTPUT_DIR_V1/*ref > $OUTPUT_DIR_V1/${ID}_epitope_indel_wildtype_v1

### New version
# running MHCI
awk -v fa=$targeted -v OUTPUT_DIR=$OUTPUT_DIR -v netMHCI=$netMHCpan40 -v pepType="mutated" '{if (substr($1,1,3)=="HLA") print netMHCI " -a " $1 " -BA -l 8,9,10,11 -f " fa " > " OUTPUT_DIR "/netMHCI3.2_" $1 "_" pepType}' ${format_hla_phlat} | xargs -P 8 -n 1 -I {} sh -c {}

awk -v fa=$refPeptide -v OUTPUT_DIR=$OUTPUT_DIR -v netMHCI=$netMHCpan40 -v pepType="ref" '{if (substr($1,1,3)=="HLA") print netMHCI " -a " $1 " -BA -l 8,9,10,11 -f " fa " > " OUTPUT_DIR "/netMHCI3.2_" $1 "_" pepType}' ${format_hla_phlat} | xargs -P 8 -n 1 -I {} sh -c {}

# running MHCII -> DRB 
cat ${format_hla_phlat} | grep "DRB" |  awk -v fa=$targeted -v OUTPUT_DIR=$OUTPUT_DIR -v netMHCII=$netMHCIIpan32 -v pepType="mutated" '{print netMHCII " -a " $1 " -f " fa " > " OUTPUT_DIR "/netMHCII4.0_"$1 "_" pepType}' |  xargs -P 8 -n 1 -I {} sh -c {}
cat ${format_hla_phlat} | grep "DRB" |  awk -v fa=$refPeptide -v OUTPUT_DIR=$OUTPUT_DIR -v netMHCII=$netMHCIIpan32 -v pepType="ref" '{print netMHCII " -a " $1 " -f " fa " > " OUTPUT_DIR "/netMHCII4.0_"$1 "_" pepType}' |  xargs -P 8 -n 1 -I {} sh -c {}

cat $OUTPUT_DIR/*mutated > $OUTPUT_DIR/${ID}_neoepitope_indel_mutation
cat $OUTPUT_DIR/*ref > $OUTPUT_DIR/${ID}_epitope_indel_wildtype

#rm $targeted 
#rm $refPeptide
rm $OUTPUT_DIR/*mutated
rm $OUTPUT_DIR/*ref
rm $OUTPUT_DIR_V1/*mutated
rm $OUTPUT_DIR_V1/*ref

