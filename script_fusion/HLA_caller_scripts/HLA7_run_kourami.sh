#!/bin/bash

if [ -f "${RESULTS_DIR}/Kourami/"*.result ]; then
	echo "results file already exists."
	exit 0
fi

mkdir "${RESULTS_DIR}/Kourami"
cd "${RESULTS_DIR}/Kourami"

# Align HLA-originating reads to Kourami's HLA panel
FASTQ1=$(ls "${RESULTS_DIR}/"*hla_extracted_1.fastq)
FASTQ2=$(ls "${RESULTS_DIR}/"*hla_extracted_2.fastq)
singularity exec \
	-B "$FASTQ1:/fastq1.gz:ro" -B "$FASTQ2:/fastq2.gz:ro" -B "$(pwd):/output" \
	"${KOURAMI_SIF}" bash -c "bwa mem -t $LSB_DJOB_NUMPROC /kourami/db/All_FINAL_with_Decoy.fa.gz /fastq1.gz /fastq2.gz | samtools view -Sb - > /output/kourami_hla_panel.bam"

# Run Kourami
singularity exec \
	-B "$(pwd)/kourami_hla_panel.bam:/input.bam:ro" -B "$(pwd):/output" \
	"${KOURAMI_SIF}" java -jar /kourami/build/Kourami.jar -d /kourami/db -o "/output/${NAME}" /input.bam

# Cleanup
rm -f *.bam *.fa

