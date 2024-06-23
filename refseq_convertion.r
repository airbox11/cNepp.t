library(biomaRt)
library(Biostrings)


args <- commandArgs(trailingOnly=TRUE)
workDir <- args[1]
  dir.indel <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/MR80/4_indel_based_prediction/result'
  if (!dir.exists(dir.indel)) {
    dir.create(dir.indel)
  }
  workDir <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/MR80/2_SNVs_based_neoepitope_prediction'
setwd(workDir)

old_mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                    host="grch37.ensembl.org",
                    path="/biomart/martservice",
                    dataset="hsapiens_gene_ensembl"
)



str_sub <- function(str1, start, end, str2){
  left <- substr(str1, 1, start-1)
  right <- substr(str1, end+1, nchar(str1))
  result <- paste0(left, str2, right)
  return(result)
}


## functions


get.pep <- function(seq,j) {
  if (j>10) {
    pep.left <- substring(seq, j-10,j)
  } else {
    pep.left <- substring(seq, 1, j)
  }
  
  if (j+10 <= nchar(seq)) {
    pep.right <- substring(seq, j+1, j+10)
  }else{
    pep.right <- substring(seq, j+1, nchar(seq)) 
  }
  
  pep <- paste0(pep.left, pep.right)
}

## 


file.input <- './snv.tsv'
df1 <- read.table(file.input, sep = '\t', stringsAsFactors = FALSE, header = TRUE)
refseq <- df1$refseq

attributes <- c('refseq_mrna','chromosome_name','coding', 'transcript_start', 'transcript_end')
attributes <- c('refseq_mrna','coding')

gbm <- getBM(
      attributes=attributes,
      filters = 'refseq_mrna', 
      values = refseq,
      mart = old_mart,
      useCache = FALSE
      )
gbm <- gbm[, attributes]

df1$ref.seq <- ''
df1$alt.seq <- ''

pattern1 <- '\\*.*'
replacement1 <- ''

fa.ref <- "./netMHCpan4_1/ref.fa"
fa.alt <- "./netMHCpan4_1/mut.fa"
if (file.exists(fa.ref)) {
  file.remove(fa.ref)
  file.remove(fa.alt)
}

for (i in 1:nrow(df1)) {
  # i <- 1
  refseq <- df1$refseq[i]
  ref <- gbm[gbm$refseq_mrna == refseq,]$coding
  ref.char.true <- substring(ref, df1$start[i], df1$end[i])
  ref.char.df1 <- df1$ref[i]
  
  if (ref.char.df1 != ref.char.true) {
    if (ref.char.df1 != 0){
      stop('something wrong!')
    }
  }
  alt <- str_sub(ref, df1$start[i], df1$end[i], df1$alt[i])
  
  ref.seq <- gsub(pattern1, replacement1, as.character(translate(DNAString(ref), no.init.codon=TRUE)))
  alt.seq <- gsub(pattern1, replacement1, as.character(translate(DNAString(alt), no.init.codon=TRUE)))
  df1$ref.seq[i] <- ref.seq
  df1$alt.seq[i] <- alt.seq
  
  for (j in 1:nchar(ref.seq)) {
    char.ref <- unlist(str_split(ref.seq, ''))[j]
    char.alt <- unlist(str_split(alt.seq, ''))[j]
    
    if (is.na(char.ref) || is.na(char.alt)) {
      break
    }
    
    tryCatch(
      expr = {
        if (char.ref != char.alt) {
          print(df1[i,1:6])
          pep.ref <- get.pep(ref.seq, j)
          pep.alt <- get.pep(alt.seq, j)
          
          # if (df1$ref[i] == 0) {
          #   file.indel.tab <- paste0(dir.indel,'/',df1$gene[i],'.tab')
          #   if (file.exists(file.indel.tab)){
          #     file.remove(file.indel.tab)
          #   }
          #   cat('geneName	aaChange	transID	targetPeptide	refPeptide	mutationType	predictedMut	proteinRef	seqMut	seqRef\n', file = file.indel.tab, append = TRUE)
          #   cat(paste(c(df1$gene[i], '', '', pep.alt, pep.ref, '', '', '', '', '','\n'), collapse = '\t'), file=file.indel.tab,append=TRUE)
          # break
          # }
          
          
          # cat(paste0('> ', df1$gene[i], ', ', df1$refseq[i], ', ', df1$start[i],'-', df1$end[i], ', Wildtype', '\n'), file=fa.ref,append=TRUE)
          cat(paste0('>', df1$gene[i],'l',1,',',df1$refseq[i], ',','WILDTYPE', '\n'), file=fa.ref,append=TRUE)
          cat(paste0(pep.ref, '\n'),file=fa.ref, append=TRUE)
    
          cat(paste0('>', df1$gene[i],'l',1,',',df1$refseq[i], ',', df1$start[i],'-', df1$end[i], ',Mutated',',',char.ref,':',j,':',char.alt,'\n'), file=fa.alt,append=TRUE)
          cat(paste0(pep.alt, '\n'),file=fa.alt, append=TRUE)
          
          break
        }
      },
      error = function(e){ 
        print('debug')
      }
    )
  }
}

# write.table(df1, file = 'snv_sub.tsv', sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
