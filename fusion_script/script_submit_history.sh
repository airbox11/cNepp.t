




fusion_pipeline S2914Nr1
fusion_pipeline S2975Nr4
fusion_pipeline S3005Nr2


fusion_pipeline H021-RFEX31-metastasis07



fusion_pipeline S2914Nr3
fusion_pipeline S2975Nr2

fusion_pipeline test2 s2


fusion_pipeline3 trcpoc004HLA s1
fusion_pipeline3 trcpoc004HLA s2 &

fusion_pipeline3 H021-RFEX31 s1


356663


fusion_pipeline S2914Nr1 s3 & #
fusion_pipeline S2914Nr3 s3 & # s3 done
fusion_pipeline S2975Nr2 s3 &
fusion_pipeline S2975Nr4 s3 &
fusion_pipeline S3005Nr2 s3 & # s3 done
fusion_pipeline trcpoc004HLA &

fusion_pipeline3 H021-RFEX31 s2 & 


## s3217
fusion_pipeline3 S3217Nr1 s2 & 
fusion_pipeline3 S3217Nr3 s2 & 

fusion_pipeline3 trcpoc004HLA s3 & 
fusion_pipeline3 H021-RFEX31 s2& 
fusion_pipeline3 H021-RFEX31 deFuse21 s2 > deFuse21.log 2>&1 &

fusion_pipeline3 H021-RFEX31 deFuse31 s2 > deFuse31.log 2>&1 &
fusion_pipeline3 H021-RFEX31 deFuse31 s2 > deFuse32.log 2>&1 &

fusion_pipeline3 H021-RFEX31 deFuse32 s2 > deFuse32.log 2>&1 &

fusion_pipeline3 H021-RFEX31 deFuse32 s2 > deFuse32.log 2>&1 &



bsub -N -r -R "rusage[mem=10G]" -J deFuse32 -W 250:00  -n 8 -o deFuse32.log "fusion_pipeline3 H021-RFEX31 deFuse32 s2" 


bsub -N -r -R "rusage[mem=10G]" -J deFuse24 -W 250:00  -n 8 -o deFuse24.log "fusion_pipeline3 H021-RFEX31 deFuse24 s2" 
bsub -N -r -R "rusage[mem=200G]" -J deFuse25 -W 250:00  -n 8 -o deFuse25.log "fusion_pipeline3 H021-RFEX31 deFuse25 s2"
bsub -N -r -R "rusage[mem=200G]" -J deFuse33 -W 250:00  -n 8 -o deFuse33.log "fusion_pipeline3 H021-RFEX31 deFuse33 s2"



bsub -N -r -R "rusage[mem=200G]" -J deFuse26 -W 250:00  -n 8 -o deFuse26.log "fusion_pipeline3 H021-RFEX31 deFuse26 s2 zh"

bsub -r -R "rusage[mem=200G]" -J deFuse26 -W 250:00 "fusion_pipeline3 H021-RFEX31 deFuse28 s2 zh"

fusion_pipeline3 H021-RFEX31 deFuse210 s2 zh 

fusion_pipeline3 H021-RFEX31 deFuse34 s2 yanhong



