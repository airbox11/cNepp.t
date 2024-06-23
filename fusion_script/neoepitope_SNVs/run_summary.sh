#! /bin/bash
###################
# require summary_MHCI.py and summary_MHCII.py
###################

#OUTPUT_DIR=$1
#PID=$2
#VCF=$3
#extendLen=14
#PIPELINE_DIR=/desktop-home/huangz/immunoinfo/immunoinfo_201703

out1=$OUTPUT_DIR/results_${PID}_MHCI_epitopes_filtered.tab
out2=$OUTPUT_DIR/results_${PID}_MHCII_epitopes_filtered.tab

## reformat netMHC prediction result files for downstream analysis
for file in `ls $OUTPUT_DIR | grep "mut$\|ref$"`
do
        awk '{out=""; for (i=1;i<=NF;i++) out=out"\t"$i; sub("\t","",out);sub("<=","",out);sub("\t\t","\t",out); print out}' $OUTPUT_DIR/$file > $OUTPUT_DIR/${file}_reformat
done


# integrate multiple files (reformat files) for MHC-I type and MHC-II type
# Identify HLA ID
hlaAllele=`ls $OUTPUT_DIR/*mut_reformat | awk '{split($1,a,"_mut_"); print a[1]}' | sort | uniq`
# Integration
###
###  To Check: VCF file merge has to be update if this vcf has not only nonsynonymous mutations
###	 
###

echo $hlaAllele | awk -v vcf=${VCF} -v mutPep=$OUTPUT_DIR/${PID}_mut.fa -v refPep=$OUTPUT_DIR/${PID}_ref.fa -v sumI=$PIPELINE_DIR/summary_MHCI.py -v sumII=$PIPELINE_DIR/summary_MHCII.py '{for (i=1;i<=NF;i++) {if (index($i,"HLA")!=0) print "python", sumI, $i"_mut_reformat", $i"_ref_reformat", vcf, mutPep, refPep, $i"_summary"; else  print "python", sumII, $i"_mut_reformat", $i"_ref_reformat", vcf, mutPep, refPep, $i"_summary" } }' | xargs -P 1 -n 1 -I {} sh -c {}

## print only command line without execution
echo $hlaAllele | awk -v vcf=${VCF} -v mutPep=$OUTPUT_DIR/${PID}_mut.fa -v refPep=$OUTPUT_DIR/${PID}_ref.fa -v sumI=$PIPELINE_DIR/summary_MHCI.py -v sumII=$PIPELINE_DIR/summary_MHCII.py '{for (i=1;i<=NF;i++) {if (index($i,"HLA")!=0) print "python", sumI, $i"_mut_reformat", $i"_ref_reformat", vcf, mutPep, refPep, $i"_summary"; else  print "python", sumII, $i"_mut_reformat", $i"_ref_reformat", vcf, mutPep, refPep, $i"_summary" } }' 

# merge epitope prediction results for MHC-I and MHC-II
#head -n 1 $VCF > $OUTPUT_DIR/headerVCF
grep "#CHROM" $VCF > $OUTPUT_DIR/headerVCF
paste $PIPELINE_DIR/header_MHCI $OUTPUT_DIR/headerVCF > $OUTPUT_DIR/headerMHCI
paste $PIPELINE_DIR/header_MHCII $OUTPUT_DIR/headerVCF > $OUTPUT_DIR/headerMHCII
cat $OUTPUT_DIR/headerMHCI $OUTPUT_DIR/netMHCI_*summary > $OUTPUT_DIR/results_${PID}_MHCI_epitopes.tab
cat $OUTPUT_DIR/headerMHCII $OUTPUT_DIR/netMHCII_*summary > $OUTPUT_DIR/results_${PID}_MHCII_epitopes.tab

# Extract neoeptitopes inhaboring mutated amino acids in core epitopes
# Fix a bug. In the core of "IMAECNA-V", it may not include mutated aa.
#awk -v extendLen=${extendLen} 'BEGIN{FS="\t"} {b=$5; sub(/-/,"",b); if ( ($2+length(b))>=(extendLen+1)  && ($2+$6)<=(extendLen+1) ) print; if (NR==1) print}' $OUTPUT_DIR/results_${PID}_MHCI_epitopes.tab > $OUTPUT_DIR/results_${PID}_MHCI_epitopes_filtered.tab_old

## calculate the mutated position in epitope. 
## filter the core-epitopes which don't include mutated amino acide
## Use core-peptide for CLASS-II
awk -v extendLen=${extendLen} 'BEGIN{FS="\t";OFS="\t"} {b=$3; sub(/-/,"",b);
        mutPos="null"
        if (NR>1){
                for (i=1; i<=length($3); i++) {
                        if (substr($3,i,1)!=substr($18,i,1)){
                                mutPos=i}}} 
        if (mutPos!="null"){
                print $0 "\t" mutPos
        }
        if (NR==1) print $0 "\t" "mutPosInEpitope" }' \
        $OUTPUT_DIR/results_${PID}_MHCI_epitopes.tab > $out1


awk -v extendLen=${extendLen} 'BEGIN{FS="\t";OFS="\t"} {
        mutPos="null"
        if (NR>1){
                for (i=1; i<=length($6); i++) {
                        if (substr($6,i,1)!=substr($18,i,1)){
                                mutPos=i}}} 
        if (mutPos!="null"){
                print $0 "\t" mutPos
        }
        if (NR==1) print $0 "\t" "mutPosInEpitope" }' \
        $OUTPUT_DIR/results_${PID}_MHCII_epitopes.tab >  $out2


# remove intermediate files
rm $OUTPUT_DIR/netMHCI*
mv ${VCF}*annovar* $OUTPUT_DIR/.
rm $OUTPUT_DIR/header*

# remove empty files
for file in `ls $OUTPUT_DIR | egrep -v summary`
do
	if [ ! -s "$OUTPUT_DIR/$file" ];then
		rm $OUTPUT_DIR/$file
	fi	
done
