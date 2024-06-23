#! /bin/bash

#module load perl/5.20.2
#####
# require convert2multiple_transcripts.sh
#####

#if [ -f ${OUT_DIR}/$predictedProtein ];then
#	echo "$predictedProtein already exists."
#	exit 0
#fi

# This script is to generate mutated protein based on VCF files
# input is the VCF file
# output is the mutated protein of fasta format


inputVCF1=VCF_frameshift
inputVCF2=VCF_nonframeshift

echo ">>>> input VCF <<<<<<"
echo `ls $INDEL`


grep -v nonframeshift $INDEL > ${OUT_DIR}/$inputVCF1

echo ">>>> run ANNOVAR <<<<<<"
perl $ANNOVAR/convert2annovar.pl -format vcf4old ${OUT_DIR}/$inputVCF1 > ${OUT_DIR}/$inputVCF1.annovar
echo "perl $ANNOVAR/convert2annovar.pl -format vcf4old ${OUT_DIR}/$inputVCF1 > ${OUT_DIR}/$inputVCF1.annovar"

perl $ANNOVAR/annotate_variation.pl -buildver hg19 ${OUT_DIR}/$inputVCF1.annovar $ANNOVAR/humandb/

sh $PIPELINE_DIR/convert2multiple_transcripts.sh ${OUT_DIR}/$inputVCF1.annovar.exonic_variant_function > ${OUT_DIR}/$inputVCF1.annovar.exonic_variant_function_multiple

perl $ANNOVAR/coding_change.pl -includesnp  ${OUT_DIR}/$inputVCF1.annovar.exonic_variant_function_multiple $ANNOVAR/humandb/hg19_refGene.txt  $ANNOVAR/humandb/hg19_refGeneMrna.fa  | awk ' {if (substr($1,1,1)==">") {if (NF<5) print $1 "," $2 "," $3; else {mutated=$1 "," $2 "," $3 "," $4 "," $9 ":" $6 ":" $11; gsub(")","",mutated);print mutated}} else print}' > ${OUT_DIR}/${predictedProtein}_frameshift


if grep -q nonframeshift $INDEL
then
    grep nonframeshift $INDEL > ${OUT_DIR}/$inputVCF2


    perl $ANNOVAR/convert2annovar.pl -format vcf4old ${OUT_DIR}/$inputVCF2 > ${OUT_DIR}/$inputVCF2.annovar
    echo "perl $ANNOVAR/convert2annovar.pl -format vcf4old ${OUT_DIR}/$inputVCF2 > ${OUT_DIR}/$inputVCF2.annovar"

    perl $ANNOVAR/annotate_variation.pl -buildver hg19 ${OUT_DIR}/$inputVCF2.annovar $ANNOVAR/humandb/

    sh $PIPELINE_DIR/convert2multiple_transcripts.sh ${OUT_DIR}/$inputVCF2.annovar.exonic_variant_function > ${OUT_DIR}/$inputVCF2.annovar.exonic_variant_function_multiple

    perl $ANNOVAR/coding_change.pl -includesnp  ${OUT_DIR}/$inputVCF2.annovar.exonic_variant_function_multiple $ANNOVAR/humandb/hg19_refGene.txt  $ANNOVAR/humandb/hg19_refGeneMrna.fa  | awk ' {if (substr($1,1,1)==">") {if (NF<5) print $1 "," $2 "," $3; else {mutated=$1 "," $2 "," $3 "," $4 "," $9 ":" $6 ":" $11; gsub(")","",mutated);print mutated}} else print}' > ${OUT_DIR}/${predictedProtein}_nonframeshift
    
fi