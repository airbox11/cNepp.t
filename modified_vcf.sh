#! /bin/env bash

file=$1
header='#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	location	SEQUENCE_CONTEXT	INFO_control(VAF=variant_allele_fraction;TSR=total_variant_supporting_reads_incl_lowqual)	ANNOTATION_control	DBSNP	1K_GENOMES	ANNOVAR_FUNCTION	GENE	EXONIC_CLASSIFICATION	ANNOVAR_TRANSCRIPTS	SEGDUP	CYTOBAND	REPEAT_MASKER	DAC_BLACKLIST	DUKE_EXCLUDED	HISEQDEPTH	SELFCHAIN	MAPABILITY	SIMPLE_TANDEMREPEATS	CONFIDENCE	RECLASSIFICATION	PENALTIES	seqBiasPresent_1	seqingBiasPresent_1	seqBiasPresent_2	seqingBiasPresent_2	Enhancers	CpGislands	TFBScons	ENCODE_DNASE	miRNAs_snoRNAs	miRBase18	COSMIC	miRNAtargets	CgiMountains	phastConsElem20bp	ENCODE_TFBS'
awk '!/^#/' $file | awk -v header="$header" 'BEGIN { print header} {print $0}' > ${file}_add_columnNames


py_script=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/modify_add_geneName.py
/omics/groups/OE0422/internal/yanhong/git/miniconda3/bin/python3.8 $py_script ${file}_add_columnNames

# cat modify_py.tsv
# mv $file ${file}_backup
