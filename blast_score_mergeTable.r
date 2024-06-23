rm(list = ls())
## args ==== ====
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  print('parameters needed for converting to xlsx')
} else{
  workDir <- args[1]
  mhc_input <- args[2]
}

## testing ==== ==== 
test <- TRUE
test <- FALSE
if (test){
  # workDir <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/K26K-P7V4JK-tumor11'
  workDir <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/MO79'
  mhc_input <- 'MHCI_epitopes_RNAseq_netMHCpan4_1.tab_renameCol_loh'
}

main_run <- function () {
  ## table of max score ==== ====
  dir1 <- file.path(workDir,'8_chose_neoepitode/blast_score')
  setwd(dir1)
  file.csv <- paste0(mhc_input,'_output.csv')
  
  df.csv <- read.table(file.csv,header = TRUE, stringsAsFactors = FALSE, sep = '\t')
  df.csv.score <- aggregate(bit_score ~ query_acc_ver, data = df.csv, max)
  colnames(df.csv.score)[2] <- 'blastp_score'
  
  ## table of MHCI
  mhc_input.2 <- paste0('../',mhc_input)
  df.mhc <- read.table(mhc_input.2, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
  col.index <- which(colnames(df.mhc) == 'Mutant_peptide')
  
  ## merge
  df.mhc.score <- merge(df.mhc, df.csv.score, by.x = 'Mutant_peptide', by.y = 'query_acc_ver',all.x = TRUE)
  df.mhc.score2 <- df.mhc.score[,c(2:col.index, ncol(df.mhc.score),1,(col.index+1):(ncol(df.mhc.score)-1))]
  
  ## write table
  write.table(df.mhc.score2, file = paste0('../',mhc_input,'_blastpScore'), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
}
main_run()