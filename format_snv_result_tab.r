
runID <- 'HL50-3'
f0 <- paste('/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline_collection/mhc4.1/', runID, '/2_SNVs_based_neoepitope_prediction', sep = '')
fs <- list.files(path = f0, pattern = '.*tab$', full.names = TRUE)


f1 <- fs[i]
t1 <- read.table(f1, sep = '\t', stringsAsFactors = FALSE, header = FALSE,comment.char = "#")
cols <- c('Mut_pos','MHC','Peptide','Core','Of','Gp','Gl','Ip','Il','Icore','Identity','Score_EL','%Rank_EL','Score_BA','%Rank_BA','Aff(nM)','BindLevel','Ref_pos','MHC','Peptide','Core','Of','Gp','Gl','Ip','Il','Icore','Identity','Score_EL','%Rank_EL','Score_BA','%Rank_BA','Aff(nM)','BindLevel','NM_accessions','nChange','aaChange','typeChange','mutPeptide','refPeptide','#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','bam','SEQUENCE_CONTEXT','INFO_control(VAF=variant_allele_fraction;TSR=total_variant_supporting_reads_incl_lowqual)','ANNOTATION_control','DBSNP','1K_GENOMES','ANNOVAR_FUNCTION','GENE','EXONIC_CLASSIFICATION','ANNOVAR_TRANSCRIPTS','SEGDUP','CYTOBAND','REPEAT_MASKER','DAC_BLACKLIST','DUKE_EXCLUDED','HISEQDEPTH','SELFCHAIN','MAPABILITY','SIMPLE_TANDEMREPEATS','CONFIDENCE','RECLASSIFICATION','PENALTIES','seqBiasPresent_1','seqingBiasPresent_1','seqBiasPresent_2','seqingBiasPresent_2','Enhancers','CpGislands','TFBScons','ENCODE_DNASE','miRNAs_snoRNAs','miRBase18','COSMIC','miRNAtargets','CgiMountains','phastConsElem20bp','ENCODE_TFBS')  
t1[cols[(ncol(t1)+1):length(cols)]] <- NA
colnames(t1) <- cols
f2 <- paste(f1,'.new', sep = '')
write.table(t1, file = f2, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
