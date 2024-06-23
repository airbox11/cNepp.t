#!/bin/bash

mkdir -p ${RESULTS_DIR}/Optitype
cd $RESULTS_DIR/Optitype

FASTQ1=$(ls "${RESULTS_DIR}/"*total_extracted_1.fastq)
FASTQ2=$(ls "${RESULTS_DIR}/"*total_extracted_2.fastq)



# FASTQ1=/omics/groups/OE0422/internal/yanhong/git/optitype/test/exome/NA11995_SRR766010_1_fished.fastq
# FASTQ2=/omics/groups/OE0422/internal/yanhong/git/optitype/test/exome/NA11995_SRR766010_2_fished.fastq

echo "
[mapping]
razers3=/usr/local/bin/razers3 
threads=$LSB_DJOB_NUMPROC

[ilp]
solver=cbc 
threads=$LSB_DJOB_NUMPROC

[behavior]
deletebam=true 
unpaired_weight=0 
use_discordant=false
" > config.ini

singularity run \
	-B "$FASTQ1:/read1.fq:ro" -B "$FASTQ2:/read2.fq:ro" -B "$(pwd)/config.ini:/config.ini:ro" -B "$(pwd):/output" \
	"${OPTITYPE_SIF}" -i /read1.fq /read2.fq -d -c /config.ini -o /output

rm -f config.ini

