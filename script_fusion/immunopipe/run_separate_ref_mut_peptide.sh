#! /bin/bash

module load python/2.7.9

echo "python ${PIPELINE_DIR}/separate_ref_mut_peptide.py $fastaFile $outputMut $outputRef $extendLen"
python ${PIPELINE_DIR}/separate_ref_mut_peptide.py $fastaFile $outputMut $outputRef $extendLen
