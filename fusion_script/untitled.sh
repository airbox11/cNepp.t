workDir=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline_collection/mhc4.1/H021-RFEX31/9_Fusion_gene_based_neoepitope_identification/2_mainRun
pid=26
log=$workDir/${pid}.log

bsub -N -r -R "rusage[mem=200G]" -J deFuse26 -W 250:00  -n 8 -o deFuse26.log "fusion_pipeline3 H021-RFEX31 deFuse$pid s2 zh"
