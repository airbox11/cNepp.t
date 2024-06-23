#! /bin/bash

#####
# require convert2multiple_transcripts.sh
#####


# This script is to generate mutated protein based on VCF files
# input is the VCF file
# output is the mutated protein of fasta format

inputVCF=${VCF}
#inputVCF=${VCF}_nonsynonymous

#grep nonsynonymous $VCF > $inputVCF

perl $ANNOVAR/convert2annovar.pl -format vcf4old $inputVCF > $inputVCF.annovar
echo "perl $ANNOVAR/convert2annovar.pl -format vcf4old $inputVCF > $inputVCF.annovar"

perl $ANNOVAR/annotate_variation.pl -buildver hg19 $inputVCF.annovar $ANNOVAR/humandb/

# choose only nonsynonymous
sh $PIPELINE_DIR/convert2multiple_transcripts.sh $inputVCF.annovar.exonic_variant_function | grep nonsynonymous > $inputVCF.annovar.exonic_variant_function_multiple

## first column of ANNOVAR file has to be the defaul format. Replace with position coordinates is impossible
perl $ANNOVAR/coding_change.pl -includesnp  $inputVCF.annovar.exonic_variant_function_multiple $ANNOVAR/humandb/hg19_refGene.txt  $ANNOVAR/humandb/hg19_refGeneMrna.fa  | awk ' {if (substr($1,1,1)==">") {if (NF<5) print $1 "," $2 "," $3; else {mutated=$1 "," $2 "," $3 "," $4 "," $9 ":" $6 ":" $11; gsub(")","",mutated);print mutated}} else print}' > $predictedProtein


# generate protein NM accession and gene name pairs
awk 'BEGIN{FS="\t";OFS="\t"}{split($3,a,":");print a[2], a[1], $4"_"$5}' $inputVCF.annovar.exonic_variant_function_multiple | sort | uniq | awk '{if (NR==1){print "NM_accessions" "\t" "GENE", "\t" "chrPos"};print $0}'> ${OUTPUT_DIR}/map_NM_geneName


sleep 30s 
