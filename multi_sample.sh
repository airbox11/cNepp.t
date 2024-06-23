#!/usr/bin/bash

path1=/icgc/dkfzlsdf/project/D120/immuno_patients_nct/sequencing/exon_sequencing/view-by-pid/IPNCT_EH63
subID1=$1_blood
subID2=EH63_$1

# find  $path1/snv*/*/${subID1}/results*/*somatic*functional_snvs*vcf | xargs  realpath | xargs -I {} echo '{}' | xargs -I {} ln -sf {}  ~/mhc/${subID2}/2_SNVs_based_neoepitope_prediction/
# find  $path1/indel*/*/${subID1}/results*/indel*somatic*functional_indel*vcf | xargs  realpath | xargs -I {} echo '{}' | xargs -I {} ln -sf {}  ~/mhc/${subID2}/4_indel_based_prediction

path2=/icgc/dkfzlsdf/project/D120/immuno_patients_nct/sequencing/rna_sequencing/view-by-pid/IPNCT_EH63
find  $path2/$1/paired/merged-alignment/.merging_2/*merged.mdup.bam | xargs  realpath | xargs -I {} echo '{}' | xargs -I {} ln -sf {}  ~/mhc/${subID2}/3_add_expression
find  $path2/$1/paired/merged-alignment/.merging_2/featureCounts/*tsv | xargs  realpath | xargs -I {} echo '{}' | xargs -I {} ln -sf {}  ~/mhc/${subID2}/3_add_expression

find  $path2/$1/paired/run*/sequence/*fastq.gz | xargs  realpath | xargs -I {} echo '{}' | xargs -I {} ln -sf {}  ~/mhc/${subID2}/3_add_expression


