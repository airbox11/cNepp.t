#! /bin/bash

#module load perl/5.20.2
#####
# require convert2multiple_transcripts.sh
#####


# This script is to generate mutated protein based on VCF files
# input is the VCF file
# output is the mutated protein of fasta format


inputVCF=VCF_nonsynonymous

echo ">>>> input VCF <<<<<<"
echo `ls $SNV`

grep nonsynonymous $SNV > ${OUT_DIR}/$inputVCF

# create ANNOVAR input files from the VCF files
echo ">>>> run ANNOVAR <<<<<<"
perl $ANNOVAR/convert2annovar.pl -format vcf4old ${OUT_DIR}/$inputVCF > ${OUT_DIR}/$inputVCF.annovar
echo "perl $ANNOVAR/convert2annovar.pl -format vcf4old ${OUT_DIR}/$inputVCF > ${OUT_DIR}/$inputVCF.annovar"

# type of variation is annotated
perl $ANNOVAR/annotate_variation.pl -buildver hg19 ${OUT_DIR}/$inputVCF.annovar $ANNOVAR/humandb/

# split single entries (per gene) into multiple entries (per transcript) if multiple transcripts were annotated for a gene
sh $PIPELINE_DIR/convert2multiple_transcripts.sh ${OUT_DIR}/$inputVCF.annovar.exonic_variant_function > ${OUT_DIR}/$inputVCF.annovar.exonic_variant_function_multiple

#wildtype and mutated protein sequences are calculated from the gene annotation
perl $ANNOVAR/coding_change.pl -includesnp  ${OUT_DIR}/$inputVCF.annovar.exonic_variant_function_multiple $ANNOVAR/humandb/hg19_refGene.txt  $ANNOVAR/humandb/hg19_refGeneMrna.fa  | awk ' {if (substr($1,1,1)==">") {if (NF<5) print $1 "," $2 "," $3; else {mutated=$1 "," $2 "," $3 "," $4 "," $9 ":" $6 ":" $11; gsub(")","",mutated);print mutated}} else print}' > ${OUT_DIR}/$predictedProtein
