#! /bin/bash
###################
# require summary_MHCI.py and summary_MHCII.py
###################

#OUTPUT_DIR=$1
#PID=$2
#VCF=$3
#extendLen=14
#PIPELINE_DIR=/home/huangz/immunoinfo/immunoinfo_201703

module load python/2.7.9
## reformat netMHC prediction result files for downstream analysis
for file in `ls $OUTPUT_DIR | grep "mut$\|ref$"`
do
        awk '{out=""; for (i=1;i<=NF;i++) out=out"\t"$i; sub("\t","",out);sub("<=","",out);sub("\t\t","\t",out); print out}' $OUTPUT_DIR/$file > $OUTPUT_DIR/${file}_reformat
done


# integrate multiple files (reformat files) for MHC-I type and MHC-II type
# Identify HLA ID
hlaAllele=`ls $OUTPUT_DIR/*mut_reformat | awk '{split($1,a,"_mut_"); print a[1]}' | sort | uniq`
echo '>>>>> hlaAllele <<<<<'
echo $hlaAllele
# Integration
echo '>>>>> run summary_MHCI.py/summary_MHCII.py <<<<<'
echo $hlaAllele | awk -v vcf=VCF_nonsynonymous -v mutPep=$OUTPUT_DIR/${PID}_mut.fa -v refPep=$OUTPUT_DIR/${PID}_ref.fa -v sumI=$PIPELINE_DIR/summary_MHCI.py -v sumII=$PIPELINE_DIR/summary_MHCII.py '{for (i=1;i<=NF;i++) {if (index($i,"HLA-")!=0) print "python", sumI, $i"_mut_reformat", $i"_ref_reformat", vcf, mutPep, refPep, $i"_summary"; else  print "python", sumII, $i"_mut_reformat", $i"_ref_reformat", vcf, mutPep, refPep, $i"_summary" } }' > $OUTPUT_DIR/intermediate.sh
cd $OUTPUT_DIR
while IFS='' read -r line || [[ -n "$line" ]]; do $line ; done < $OUTPUT_DIR/intermediate.sh


#echo $hlaAllele | awk -v vcf=VCF_nonsynonymous -v mutPep=$OUTPUT_DIR/${PID}_mut.fa -v refPep=$OUTPUT_DIR/${PID}_ref.fa -v sumI=$PIPELINE_DIR/summary_MHCI.py -v sumII=$PIPELINE_DIR/summary_MHCII.py '{for (i=1;i<=NF;i++) {if (index($i,"HLA-")!=0) print "python", sumI, $i"_mut_reformat", $i"_ref_reformat", vcf, mutPep, refPep, $i"_summary"; else  print "python", sumII, $i"_mut_reformat", $i"_ref_reformat", vcf, mutPep, refPep, $i"_summary" } }'


# merge epitope prediction results for MHC-I and MHC-II
head -n 1 $VCF > $OUTPUT_DIR/headerVCF
paste $PIPELINE_DIR/header_MHCI $OUTPUT_DIR/headerVCF > $OUTPUT_DIR/headerMHCI
paste $PIPELINE_DIR/header_MHCII $OUTPUT_DIR/headerVCF > $OUTPUT_DIR/headerMHCII
cat $OUTPUT_DIR/headerMHCI $OUTPUT_DIR/netMHCI_*summary > $OUTPUT_DIR/results_${PID}_MHCI_epitopes.tab
cat $OUTPUT_DIR/headerMHCII $OUTPUT_DIR/netMHCII_*summary > $OUTPUT_DIR/results_${PID}_MHCII_epitopes.tab

python $PIPELINE_DIR/filter_for_ident_mut_ref_pep.py $OUTPUT_DIR/results_${PID}_MHCI_epitopes.tab $OUTPUT_DIR/results_${PID}_MHCI_epitopes_filtered_ident.csv
python $PIPELINE_DIR/filter_for_ident_mut_ref_pep.py $OUTPUT_DIR/results_${PID}_MHCII_epitopes.tab $OUTPUT_DIR/results_${PID}_MHCII_epitopes_filtered_ident.csv

# Extract neoeptitopes inhaboring mutated amino acids in core epitopes
# Fix a bug. In the core of "IMAECNA-V", it may not include mutated aa.
#awk -v extendLen=${extendLen} 'BEGIN{FS="\t"} {b=$5; sub(/-/,"",b); if ( ($2+length(b))>=(extendLen+1)  && ($2+$6)<=(extendLen+1) ) print; if (NR==1) print}' $OUTPUT_DIR/results_${PID}_MHCI_epitopes.tab > $OUTPUT_DIR/results_${PID}_MHCI_epitopes_filtered.tab_old


## after header changed
#awk -v extendLen=${extendLen} 'BEGIN{FS="\t"} {b=$4; sub(/-/,"",b); 
#        if (length($35)==(2*extendLen+1)) {if (($1+length(b)-1)>=(extendLen+1)  && ($1+$5)<=(extendLen+1) ) 
#                                                print $0 "\t" extendLen+2-$1} 
#        if (length($35)!=(2*extendLen+1)) {refAA=substr($32,3,1);mutAA=substr($32,length($32),1); 
#                                        if(refAA==substr($36,extendLen+1,1) && mutAA==substr($35,extendLen+1,1)){
#                                                if ( ($1+length(b)-1)>=(extendLen+1)  && ($1+$5)<=(extendLen+1) ) 
#                                                        print $0 "\t" extendLen+2-$1}
#                                        else {resetLen=length($35)-extendLen;
#                                                if ( ($1+length(b)-1)>=resetLen  && ($1+$5)<=resetLen ) print $0 "\t" resetLen+2-$1}}
#        if (NR==1) print $0 "\t" "mutPosInEpitope" }' \
#        $OUTPUT_DIR/results_${PID}_MHCI_epitopes.tab > $OUTPUT_DIR/results_${PID}_MHCI_epitopes_filtered.tab_new
#
#
#awk -v extendLen=${extendLen} 'BEGIN{FS="\t"} {
#        if (length($29)==(2*extendLen+1)) {
#                if ( ($1+1+$5+1+length($6)) >= (extendLen+1) && ($1+1+$5+1)<=(extendLen+1)) print;}
#        else {  refAA=substr($27,3,1);
#                mutAA=substr($27,length($27),1);
#                if(refAA==substr($29,extendLen+1,1) && mutAA==substr($30,extendLen+1,1)){
#                        if ( ($1+1+$5+1+length($6)) >= (extendLen+1) && ($1+1+$5+1)<=(extendLen+1)) print;}
#                else {resetLen=length($29)-extendLen;
#                        if ( ($1+1+$5+1+length($6)) >= resetLen && ($1+1+$5+1)<=resetLen) print;}
#        }               
#        if (NR==1) print}' $OUTPUT_DIR/results_${PID}_MHCII_epitopes.tab >  $OUTPUT_DIR/results_${PID}_MHCII_epitopes_filtered.tab_new



# remove intermediate files
rm $OUTPUT_DIR/*_reformat
rm ${VCF}_nonsynonymous*annovar*

# remove empty files
for file in `ls $OUTPUT_DIR | egrep -v summary`
do
	if [ ! -s "$OUTPUT_DIR/$file" ];then
		rm $OUTPUT_DIR/$file
	fi	
done
