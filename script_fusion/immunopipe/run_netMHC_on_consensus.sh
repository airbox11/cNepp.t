#! /bin/bash
set -e
# Run netMHC for HLA-A, B, C and DQB1, DRB1 
                    
INPUT1=`cat ${sourceDir}/netMHCI_input.txt`
INPUT2=`cat ${sourceDir}/netMHCII_input.txt`


##read HLAI alleles into array and run netMHC on each one

function filter_seq () {
    script=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/update_netMHCpan/netMHCpan_4.1/filter_short_seq_15_8.py
    python $script $peptideFasta
}

function filter_seq_rm () {
    rm -rf ${peptideFasta}_*
}

filter_seq


while IFS=',' read -ra ADDR; do
    for i in "${ADDR[@]}"; do
        $netMHCpan -BA -f ${peptideFasta}_8  -l 8,9,10,11 -a $i -s > $OUTPUT_DIR/netMHCI_${i}_${pepType}
        # $netMHCpan -BA -f ${peptideFasta}_8  -l 8,9,10,11,12,13 -a $i -s > $OUTPUT_DIR/netMHCI_${i}_${pepType} #27mer
    done
done <<< "$INPUT1"


##read HLAII alleles into array and run netMHC on each one
while IFS=',' read -ra ADDR; do
    for i in "${ADDR[@]}"; do
        $netMHCIIpan -BA -f ${peptideFasta}_15 -a $i -s > $OUTPUT_DIR/netMHCII_${i}_${pepType}
    done
done <<< "$INPUT2"

# filter_seq_rm
