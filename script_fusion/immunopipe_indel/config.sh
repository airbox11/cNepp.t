#! /bin/bash

# the length of amino acid before/after mutation position. Total length of 29 aa was chosen here.
extendLen=14 #14

# identify your control folder
controlType=(control buffy_coat blood)

#PIPELINE_DIR=/icgc/dkfzlsdf/analysis/G200/pfitzerl/HLA-typing/immunopipe

#data structure
#results_per_pid is used for the snv vcf
#view_by_pid for the fastq - has to be changed for the extracted reads/bams
#RESULTS_PER_PID=/icgc/dkfzlsdf/analysis/D120/kosalogl/Immuno-Patients-NCT/sequencing/exon_sequencing/results_per_pid/
#VIEW_PER_PID=/icgc/dkfzlsdf/analysis/D120/kosalogl/Immuno-Patients-NCT/sequencing/exon_sequencing/view_by_pid/



#soft
ANNOVAR=/icgc/ngs_share/annovar/annovar_Feb2016
netMHCIIpan_SIF=/icgc/dkfzlsdf/analysis/G200/immuno/tools/netMHCIIpan-4.0/netMHCIIpan.sif
netMHCpan_SIF=/icgc/dkfzlsdf/analysis/G200/immuno/tools/netMHCpan-4.1/netMHCpan.sif
phlatRelease=/icgc/dkfzlsdf/analysis/G200/immuno/tools/phlat-release









