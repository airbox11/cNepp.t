#! /bin/bash

# the length of amino acid before/after mutation position. Total length of 29 aa was chosen here.
extendLen=14 #14

# identify your control folder
controlType=(control buffy_coat blood)

#PIPELINE_DIR=/icgc/dkfzlsdf/analysis/G200/pfitzerl/HLA-typing/immunopipe

#soft
ANNOVAR=/icgc/ngs_share/annovar/annovar_Feb2016
netMHCIIpan_SIF=/icgc/dkfzlsdf/analysis/G200/immuno/tools/netMHCIIpan-4.0/netMHCIIpan.sif
netMHCpan_SIF=/icgc/dkfzlsdf/analysis/G200/immuno/tools/netMHCpan-4.1/netMHCpan.sif
phlatRelease=/icgc/dkfzlsdf/analysis/G200/immuno/tools/phlat-release
OPTITYPE_SIF=/icgc/dkfzlsdf/analysis/G200/immuno/tools/optitype/optitype-1.3.1.sif
KOURAMI_SIF=/icgc/dkfzlsdf/analysis/G200/immuno/tools/kourami-0.9.6/kourami.sif

