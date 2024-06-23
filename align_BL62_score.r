rm(list=ls())

## args ==== ====
library(Biostrings)
library(rlang)
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  print('parameters needed for converting to xlsx')
} else{
  workDir <- args[1]
  mhc_input <- args[2]
}

## test ==== ====
test <- FALSE
test <- TRUE
if (test){
  workDir <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/MO79'
  mhc_input <- 'MHCI_epitopes_RNAseq_netMHCpan4_1.tab_renameCol_loh'
}

aln_matrix   <- "BLOSUM62"
data(list=c(aln_matrix),envir=.GlobalEnv)

## functions ==== ====

parallel_run <- function(qur.mat, ref.mat){
  
  if (nrow(qur.mat) == 1){
    df.maxScore <- get_max_score(qur.mat, ref.mat)
    return(df.maxScore)
  }
  
  ## prpare time, log file
  t1 <- proc.time()
  # parallel.log <- '/home/lyuya/tmpodcf/parallel.log'
  # parallel.log <- ''
  # file.remove(parallel.log, showWarnings = FALSE)
  
  ## set cluster
  num.cores <- nrow(qur.mat)
  if (num.cores > 12) num.cores <- 12
  cl <- parallel::makeCluster(num.cores)
  doParallel::registerDoParallel(cl)
  
  ## export variable
  parallel::clusterExport(cl, varlist = ls(.GlobalEnv), envir = .GlobalEnv)
  
  ## parallel run
  df.maxScore <- parallel::parApply(cl, qur.mat, 1,
                                    function(x)
                                      get_max_score_single(x, ref.mat))
  
  parallel::stopCluster(cl)
  print(paste0("Alignment Time : ", (proc.time() - t1)[3]))
  
  ## convert to data.frame
  df.maxScore.lite <- data.frame(t(df.maxScore))
  colnames(df.maxScore.lite) <-
    c('qur.pep.seq', 'ref.pep.seq', 'score')
  rownames(df.maxScore.lite) <- c(1:nrow(df.maxScore.lite))
  
  return(df.maxScore.lite)
}


main_run <- function(){
  input.dir  <- file.path(workDir, '/8_chose_neoepitode/blast_score')
  output.dir <- file.path(workDir, '/8_chose_neoepitode')
  dir.create(output.dir, showWarnings = FALSE)
  setwd(output.dir)

  ## input reference ====
  # rFile <- '/omics/groups/OE0422/internal/yanhong/git/hex/unipro/hpv_reviewed_sequences_all_subtypes.csv_fasta'
  # rFile <- '/omics/groups/OE0422/internal/yanhong/git/hex/virus_bac_database/sequences_human_virus_freseq_all.fasta'
  # rFile <- '/omics/groups/OE0422/internal/yanhong/git/hex/virus_bac_database/sequences_human_virus_freseq_all.fasta_3'
  # rFile <- '/omics/groups/OE0422/internal/yanhong/git/hex/virus_bac_database/hav_1.fasta'
  file.matrix <- '/omics/groups/OE0422/internal/yanhong/git/hex/virus_bac_database/splice_py/sequences_human_virus_freseq_all.fasta_seq_matrix'


  ## input query ====
  qFile <- file.path(input.dir, mhc_input)
  qFile <- paste0(qFile, '_query.fa')
  # qFile <- '/omics/groups/OE0422/internal/yanhong/git/20220406_paula_hex_test/mart1.fa'
  
  ## get max score table ====
  # align to every length, without window sliding
  df.maxScore.allLength <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df.maxScore.allLength) <- c('qur.pep.seq', 'ref.pep.seq', 'score')
  query.length <-sort(unique(nchar(readLines(qFile)[c(FALSE,TRUE)])))
  
  for (i in query.length){
  # for (i in 9:9){
    ## prepare qur matrix
    qur.mat <- get_matrix_query(qFile, window.width = i)
    
    ## prepare ref matrix
    file.matrix1 <- paste0(file.matrix, '_', i)
    # ref.mat <- read.table(file.matrix1, sep = '\t', stringsAsFactors = FALSE, header = FALSE)
    # saveRDS(ref.mat, file = paste0(file.matrix1, '.rdata'))
    
    # ref.mat <- readRDS(paste0(file.matrix1, '.rdata'))[1:100,]
    ref.mat <- readRDS(paste0(file.matrix1, '.rdata'))
    
    
    ## get max score and seq
    t1 <- proc.time()
    # df.maxScore <- get_max_score(qur.mat, ref.mat)
    df.maxScore.lite <- parallel_run(qur.mat, ref.mat)
    print(paste0("Alignment Time : ", (proc.time() - t1)[3]))
    

    df.maxScore.allLength <- rbind(df.maxScore.allLength, df.maxScore.lite)
  }
  
  colnames(df.maxScore.allLength) <- c('Mutant_peptide','alignBL62_ref_seq','alignBL62_score')
  ## merge to MHC table ====
  file.mhc <- paste0(output.dir,'/',mhc_input,'_blastpScore' )
  df.mhc <- read.table(file.mhc, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
  col.index <- which(colnames(df.mhc)=='Mutant_peptide')
  
  df.final <- merge(df.mhc, df.maxScore.allLength, by = 'Mutant_peptide', all.x=TRUE)
  df.final2 <- df.final[,c(2:(col.index),(ncol(df.final)-1):ncol(df.final),1,(col.index+1):(ncol(df.final)-2))]

  ## write to file ====
  file.output <- paste0(output.dir,'/',mhc_input, '_alignBL62score')
  write.table(df.final2, file=file.output, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
}

##

get_max_score_single <- function(qur.pep, ref.mat){
  # df1 <- data.frame(matrix(ncol = 3, nrow = 0))
  # colnames(df1) <- c('qur.pep.seq', 'ref.pep.seq', 'score')
  
  vec1 <- c()
  seq.score <- align_to_ref_epitopes(qur.pep, ref.mat)
  qur.pep.seq <- paste(qur.pep, collapse='')
  # df1[nrow(df1) + 1,] = list(qur.pep.seq, seq.score[1], seq.score[2])
  
  vec1 <- c(qur.pep.seq, seq.score[1], seq.score[2])
  return(vec1)
}

get_max_score <- function(qur.mat, ref.mat){
  df1 <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df1) <- c('qur.pep.seq', 'ref.pep.seq', 'score')
  
  for (i in 1:nrow(qur.mat)){
    qur.pep <- qur.mat[i,]
    seq.score <- align_to_ref_epitopes(qur.pep, ref.mat)
    
    qur.pep.seq <- paste(qur.pep, collapse='')
    df1[nrow(df1) + 1,] = list(qur.pep.seq, seq.score[1], seq.score[2])
  }
  return(df1)
}



##
get_matrix_query <- function (file.input, window.width){
  aaSet <- readAAStringSet(file.input, format = "fasta")
  aaSet <- aaSet[nchar(aaSet)==window.width]
  
  df1 <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(df1) <- c('pep.name', 'pep.seq')
  for (i in 1:length(aaSet)){
    sequence <-  as.character(aaSet[[i]])
    query.name <- names(aaSet)[i]

    for (j in 1: (nchar(sequence)-window.width+1)){
      seq.sub <- substr(sequence, j, (j+window.width-1))
      pep.name <- paste(query.name, seq.sub, sep='_')
      df1[nrow(df1) + 1,] = list(pep.name, seq.sub)
    }
  }
  vec <- stats::setNames(df1$pep.seq, df1$pep.name)
  aaSet.sub <- Biostrings::AAStringSet(vec)
  mat <- Biostrings::as.matrix(aaSet.sub, use.names = T)
  return(mat)
}




align_to_ref_epitopes <- function(qurSeq, ref, mag = 4) {
  aln_matrix_solid <- get(aln_matrix, envir = .GlobalEnv)
  xPosWt <- get_pos_weights(qurSeq, mag)
  
  vec1 <- apply(ref, 1, function(x) get_score(x, qurSeq, xPosWt, aln_matrix_solid))
  max.score <- max(vec1)
  max.index <- which(vec1==max.score)[[1]]
  refSeq <- paste0(ref[max.index,], collapse = '')

  return(c(refSeq, max.score))
}

get_score <- function (refSeq, qurSeq, xPosWt, aln_matrix_solid) {
  xAlnScore <- 0
  xAlnMat <- diag(aln_matrix_solid[qurSeq, refSeq])
  xAlnAmp <- xPosWt * xAlnMat
  xAlnScore <- sum(xAlnAmp)
  return (xAlnScore)
}

##

get_pos_weights <- function(x, mag = 4) {
  posScore <- c()
  xLen <- length(x)
  xMid <- xLen / 2
  if (xLen %% 2) {
    xMidCeil <- ceiling(xMid)
    xMidScore <- xMidCeil * mag
    posScore <- seq(1, xMidScore, mag)
    posScore <- c(posScore, rev(posScore[1:length(posScore) - 1]))
    xTopFloor <- floor(xLen / 3)
    posScore[1:xTopFloor] <- 1:xTopFloor
    xTail <- xLen - xTopFloor + 1
    posScore[xTail:xLen] <- xTopFloor:1
    
    } else{
    xMidScore <- xMid * mag
    posScore <- seq(1, xMidScore, mag)
    posScore <- c(posScore, rev(posScore))
    xTopFloor <- floor(xLen / 3)
    posScore[1:xTopFloor] <- 1:xTopFloor
    xTail <- xLen - xTopFloor + 1
    posScore[xTail:xLen] <- xTopFloor:1
  }
  return(posScore)
}

## exute main ==== ====
main_run()