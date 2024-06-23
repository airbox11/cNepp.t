#! /bin/bash

# Run netMHC for HLA-A, B, C and DQB1, DRB1 
                    
INPUT1=`cat ${OUTPUT_DIR}/../netMHCI_input.txt`
INPUT2=`cat ${OUTPUT_DIR}/../netMHCII_input.txt`

##read HLAI alleles into array and run netMHC on each one
while IFS=',' read -ra ADDR; do
    for i in "${ADDR[@]}"; do
        singularity exec -B "${SCRATCHDIR}/${LSB_JOBID}:/tmp" -B "${peptideFasta}:/input/peptideFasta.fa:ro" -B "${OUTPUT_DIR}:/output" \
        "$netMHCpan_SIF" sh -c "/netMHCpan-4.1/netMHCpan -f /input/peptideFasta.fa \
        -l 8,9,10,11 -a $i -s > /output/netMHCI_${i}_${pepType}"
    done
done <<< "$INPUT1"

##read HLAII alleles into array and run netMHC on each one
while IFS=',' read -ra ADDR; do
    for i in "${ADDR[@]}"; do
        singularity exec -B "${SCRATCHDIR}/${LSB_JOBID}:/tmp" -B "${peptideFasta}:/input/peptideFasta.fa:ro" -B "${OUTPUT_DIR}:/output" \
        "$netMHCIIpan_SIF" sh -c "/netMHCIIpan-4.0/netMHCIIpan -f /input/peptideFasta.fa \
        -a $i -s > /output/netMHCII_${i}_${pepType}"
    done
done <<< "$INPUT2"
