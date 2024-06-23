#! /bin/bash

module load python/2.7.9

echo "python ${PIPELINE_DIR}/separate_mut_ref_indel.py ${fastaFile}_frameshift $outputMut $outputRef $extendLen"
python ${PIPELINE_DIR}/separate_mut_ref_indel.py ${fastaFile}_frameshift $outputMut $outputRef $extendLen

sleep 1m

if [ -e ${fastaFile}_nonframeshift ]
then 
    echo "python ${PIPELINE_DIR}/separate_ref_mut_peptide_indel_nfs.py ${fastaFile}_nonframeshift ${outputMut}_temp ${outputRef}_temp $extendLen"
    python ${PIPELINE_DIR}/separate_ref_mut_peptide_indel_nfs.py ${fastaFile}_nonframeshift ${outputMut}_temp ${outputRef}_temp $extendLen
    cat ${outputMut}_temp >> $outputMut
    cat ${outputRef}_temp >> $outputRef
fi
