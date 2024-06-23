#! /bin/bash

# the length of amino acid before/after mutation position. Total length of 29 aa was chosen here.
extendLen=14 # default setting is 14 before and after mutated position

# Your pipeline path including all dependant scripts
# SNVs based neoepitope prediction directory path
PIPELINE_DIR=/icgc/dkfzlsdf/analysis/D120/scripts/neoepitope_SNVs


# Required tools
ANNOVAR=/icgc/ngs_share/annovar/annovar_Feb2016
netMHCIIpan=/lsdf/dkfz/analysis/G200/immuno/tools/netMHCIIpan-3.2/netMHCIIpan
netMHCpan=/lsdf/dkfz/analysis/G200/immuno/tools/netMHCpan-4.0/netMHCpan
