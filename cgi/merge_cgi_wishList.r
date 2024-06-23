library(biomaRt)
library(stringr)
options(error=traceback)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  print('parameters needed for converting to xlsx')
} else{
  workDir <- args[1]
}
vcfOnly <- Sys.getenv('vcfOnly')
org.snv <- Sys.getenv('org_snv')
org.ind <- Sys.getenv('org_ind')

print(workDir)
print(vcfOnly)
print(org.snv)
print(org.ind)

## for test ====

  # workDir <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/H021-PAQVNC_meta'
  # vcfOnly <- 'origin'
  # org.snv <- ''
  # org.ind <- ''

## test done

setwd(workDir)
## generic func ====
get_fileName <- function(dir1, pattern){
  l1 <- list.files(dir1)
  file1 <- l1[str_detect(l1, pattern)]
  return(file1)
}

promise_hg38 <- function(df3, snv_indel) {
  if (snv_indel == 'snv') {
    df.promise.hg38 <- read.table(file = org.snv, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  } else {
    df.promise.hg38 <- read.table(file = org.ind, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  }
  
  df.promise.hg38$chrPos <- paste0(df.promise.hg38$CHROM,'_', df.promise.hg38$POS)
  cols1 <- colnames(df.promise.hg38)[str_detect(colnames(df.promise.hg38), pattern = 'ANN|GEN')]
  df.promise.1 <- df.promise.hg38[,colnames(df.promise.hg38) %in% c(cols1,'chrPos')]

  colnames(df3)[which(names(df3)=='chr')] <- 'CHROM'
  colnames(df3)[which(names(df3)=='pos')] <- 'POS'
  df3$chrPos <- paste0(df3$CHROM,'_',df3$POS)
  df.3.promise <- merge(df3,df.promise.1, by='chrPos',all.x = TRUE)
  return(df.3.promise)
}

## prepare cgi  ====
perpare.cgi <- function() {
  file.cgi <- './1_hla_type/cgi/alterations.tsv'
  df.cgi <- read.table(file.cgi, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
  # if(!file.exists(file.cgi)){
  #   df.cgi <- readRDS(file = '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/promise/batch4_new_snv__20240322/cig.fake.rd')
  # }else{
  #   df.cgi <- read.table(file.cgi, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
  # }

  df.cgi$ch.pos <- paste0(df.cgi$chr,'_', df.cgi$pos)
  df.cgi.lite <- df.cgi[,c('CGI.Oncogenic.Summary', 'CGI.Oncogenic.Prediction', 'CGI.Consequence','CGI.Transcript','ch.pos')]
  return(df.cgi.lite)

  ##
  # mart <- useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl")
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
  transcript_ids <- unique(df.cgi$CGI.Transcript)
  res <- getBM(attributes = c('ensembl_transcript_id',
                              'ensembl_gene_id',
                              'external_transcript_name',
                              'external_gene_name'),
               filters = 'ensembl_transcript_id',
               values = transcript_ids,
               useCache = FALSE,
               mart = mart)

  df.cgi$ch.pos <- paste0(df.cgi$chr,'_', df.cgi$pos)
  df.cgi.lite <- df.cgi[,c('CGI.Oncogenic.Summary', 'CGI.Oncogenic.Prediction', 'CGI.Consequence','CGI.Transcript','ch.pos')]
  df.cgi.lite <- df.cgi.lite[!duplicated(df.cgi.lite),]
  df.cgi.lite.res <- merge(df.cgi.lite, res[,c('ensembl_transcript_id','ensembl_gene_id')], by.x = 'CGI.Transcript', by.y = 'ensembl_transcript_id', all.x = TRUE)
  df.cgi.lite.res <- df.cgi.lite.res[!duplicated(df.cgi.lite.res$ensembl_gene_id),]
  df.cgi.lite.res <- df.cgi.lite.res[, -which(colnames(df.cgi.lite.res) %in% c("ensembl_gene_id"))]
  return(df.cgi.lite.res)
}


## snv func ====
merge.snv.run <- function() {

  dir1 <-'8_chose_neoepitode'

  ## snv MHCI
  file.input <- get_fileName(dir1, 'MHCI_.*renameCol_loh$')
  print(file.input)
  if (length(file.input) == 0) {
    file.input <- get_fileName(dir1, 'MHCI_.*renameCol$')
  }
  if (length(file.input) != 1) {
    stop(paste0('=== === No file found: ',workDir,'/',dir1, '/', 'MHCI_.*renameCol$'))
  }
  df1 <- read.table(paste0('./8_chose_neoepitode/', file.input), sep = '\t', stringsAsFactors = FALSE, header = TRUE)
  df1$geneID <- str_match(df1$geneID, pattern = '(ENSG\\d+)\\.?\\d?')[,2]

  df1$ch.pos <- paste0(df1$CHROM, '_', df1$POS)
  df2 <- merge(df1, df.cgi.lite.res, by.x = 'ch.pos', by.y = 'ch.pos', all.x = TRUE)
  df2 <- df2[, -which(colnames(df2) %in% c("ch.pos"))]
  # df2 <- merge(df1, df.cgi.lite.res, by.x = 'geneID', by.y = 'ensembl_gene_id', all.x = TRUE)

  df3 <- merge(df2, df.wish[,!colnames(df.wish) %in% c('external_gene_name')], by.x = 'geneID', by.y = 'ensembl_gene_id', all.x = TRUE)
  df3$wishList[is.na(df3$wishList)] <- FALSE
  cols <- c('geneID','Gene','MHC_allele','Mutant_peptide','aaChange','Epitope_length','Mut_pos_epitope','Score_EL_mut','Rank_EL_Mut','Score_BA_mut','Rank_BA_mut','Aff.nM_mut','BindLevel_mut','Wildtype_peptide','Score_EL_wt','Rank_EL_wt','Score_BA_wt','Rank_BA_wt','Aff.nM_wt','BindLevel_wt','Epi_pos_in_longpep','Mutant_long_peptide','WildType_long_peptide','GenBank_entry_mRNA','TPM','FPKM','Mean','Median','sumReads','numOfBaseExp','freRef','freAlt','expAlt','freAlt.FPKM','dna_freAlt','dna_cov','CHROM','POS','CGI.Oncogenic.Summary','CGI.Oncogenic.Prediction','CGI.Consequence','CGI.Transcript','wishList', 'Hydrophobic_GRAVY', 'ANNOVAR_TRANSCRIPTS')

  for (i in cols[!cols%in%names(df3)]){
    df3[[i]] <- NA
  }
  df3 <- df3[,cols]
  if (vcfOnly == 'promise') {
    df3 <- promise_hg38(df3, 'snv')
  }


  write.table(df3, file = paste0('./8_chose_neoepitode/', file.input,'_wish'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

  ## snv MHCII
  file.input <- get_fileName(dir1, 'MHCII_.*renameCol$')
  if (length(file.input) != 1) {
    print(paste0('=== === No file found: ',workDir,'/',dir1, '/', 'MHCII_.*renameCol$'))
    return()
  }
  file.input <- paste0('./8_chose_neoepitode/', file.input)
  df1 <- read.table(file.input, sep = '\t', stringsAsFactors = FALSE, header = TRUE)
  df1$geneID <- str_match(df1$geneID, pattern = '(ENSG\\d+)\\.?\\d?')[,2]

  df1$ch.pos <- paste0(df1$CHROM, '_', df1$POS)
  df2 <- merge(df1, df.cgi.lite.res, by.x = 'ch.pos', by.y = 'ch.pos', all.x = TRUE)
  df2 <- df2[, -which(colnames(df2) %in% c("ch.pos"))]

  df3 <- merge(df2, df.wish[,!colnames(df.wish) %in% c('external_gene_name')], by.x = 'geneID', by.y = 'ensembl_gene_id', all.x = TRUE)
  df3$wishList[is.na(df3$wishList)] <- FALSE

  cols <- c('geneID','Gene','MHC_allele','Mutant_peptide','aaChange','Epitope_length','Mut_pos_epitope','X9mer_core_mut','Aff.nM_mut','Rank_mut','BindLevel_mut','Wildtype_peptide','X9mer_core_wt','Aff.nM_wt','Rank_wt','BindLevel_wt','WildType_lng_peptide','Mutant_long_peptide','Epi_pos_in_longpep','GenBank_entry_mRNA','TPM','FPKM','Mean','Median','sumReads','numOfBaseExp','freRef','freAlt','expAlt','freAlt.FPKM','dna_freAlt','dna_cov', 'CHROM','POS','CGI.Oncogeic.Summary','CGI.Oncogenic.Prediction','CGI.Consequence','CGI.Transcript','wishList', 'Hydrophobic_GRAVY', 'ANNOVAR_TRANSCRIPTS')

  for (i in cols[!cols%in%names(df3)]){
    df3[[i]] <- NA
  }
  df3 <- df3[,cols]
  if (vcfOnly == 'promise') {
    df3 <- promise_hg38(df3, 'snv')
  }
  write.table(df3, file = paste0(file.input,'_wish'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

}



## indel func ====

merge.indel <- function(dir1, pattern){
  file.input <- get_fileName(dir1, pattern)
  if (length(file.input) != 1) {
    print(paste0('=== === No file found: ',workDir,'/',dir1, '/', pattern))
    return()
  }
  
  df1 <- read.table(paste0(dir1, file.input), sep = '\t', stringsAsFactors = FALSE, header = TRUE)
  if (nrow(df1) == 0) {
    print(paste0('Warning: no valid rows in file:', file.input))
    return()
  }

  df1$ch.pos <- paste0(df1$chr, '_', df1$pos)
  df2 <- merge(df1, df.cgi.lite.res, by.x = 'ch.pos', by.y = 'ch.pos', all.x = TRUE)
  df2 <- df2[, -which(colnames(df2) %in% c("ch.pos"))]
  # df2 <- merge(df1, df.cgi.lite.res, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id', all.x = TRUE)

  df3 <- merge(df2, df.wish[,!colnames(df.wish) %in% c('external_gene_name')], by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id', all.x = TRUE)
  df3$wishList[is.na(df3$wishList)] <- FALSE

  if (vcfOnly == 'promise') {
    df3 <- promise_hg38(df3,'indel')
  }
  write.table(df3, file = paste0(dir1, file.input,'_wish'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
  
}

merge.indel.run <- function() {
  dir1 <-'./4_indel_based_prediction/'
  merge.indel(dir1, 'indel_mutant_MHCI.tsv$')
  merge.indel(dir1, 'indel_mutant_MHCII.tsv$')
  merge.indel(dir1, 'indel_wildType_MHCI.tsv$')
  merge.indel(dir1, 'indel_wildType_MHCII.tsv$')
}

## main run ====
df.cgi.lite.res <- perpare.cgi()
df.wish <- readRDS(file ='/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/cgi/wishList_ensg.rdata')

print('===[ cgi.wish merge: start')
print('=== === code path: /omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/cgi/merge_cgi_wishList.r')
merge.snv.run()
merge.indel.run()
print('===] cgi.wish merge: done')