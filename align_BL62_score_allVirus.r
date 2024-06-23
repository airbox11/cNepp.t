rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  print('parameters needed for converting to xlsx')
} else{
  workDir <- args[1]
  inputID <- args[2]
  file.bl.score <- args[3]
}

## test ==== ====
test <- FALSE
test <- TRUE
if (test){
  workDir <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/MO79'
  inputID <- 'MHCI_epitopes_RNAseq_netMHCpan4_1.tab'
  file.bl.score <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/MO79/8_chose_neoepitode/blast_score/MHCI_epitopes_RNAseq_netMHCpan4_1.tab_renameCol_loh_query.fa_blscore'
}

main_run <- function(){
  
  ## get max score table ====
  df.maxScore.allLength <- read.table(file.bl.score, header = FALSE, stringsAsFactors = FALSE, sep = '\t')
  colnames(df.maxScore.allLength) <- c('Mutant_peptide','alignBL62_ref_seq','alignBL62_score')
  ## merge to MHC table ====
  file.mhc <- paste0(workDir,'/8_chose_neoepitode/',inputID,'_renameCol_loh_blastpScore' )
  df.mhc <- read.table(file.mhc, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
  col.index <- which(colnames(df.mhc)=='Mutant_peptide')
  
  df.final <- merge(df.mhc, df.maxScore.allLength, by = 'Mutant_peptide', all.x=TRUE)
  df.final2 <- df.final[,c(2:(col.index),(ncol(df.final)-1):ncol(df.final),1,(col.index+1):(ncol(df.final)-2))]
  
  ## write to file ====
  file.output <- paste0(file.mhc, '_alignBL62score')
  write.table(df.final2, file=file.output, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
}

main_run()