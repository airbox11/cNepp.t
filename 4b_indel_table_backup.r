library(stringr)
library(biomaRt)

args1 <- commandArgs(trailingOnly = TRUE)
workDir <- args1[1]

### for test:
# workDir <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/promise/output_datasets/S014-2CDKKU/'
# workDir <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/EB72'
### test end
setwd(paste(workDir,'4_indel_based_prediction/result',sep = '/'))

# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
old_mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                    host="grch37.ensembl.org",
                    path="/biomart/martservice",
                    dataset="hsapiens_gene_ensembl"
)

f1 <- list.files('.')
genes <- str_match(f1[str_detect(f1, pattern = '.*.tab$')], pattern = '(.*)\\.tab')[,2]
file_long_pep <- '../indel_long_peptides.tsv'
if (file.exists(file_long_pep)) file.remove(file_long_pep)
write(paste('gene', 'splice_state', 'ref/mut', 'sequence', sep = '\t'),append = TRUE, file = file_long_pep)
## prepare long peps

prepare_long_peps <- function(gene, type){
  pep.long.ref <- readLines(paste(gene,'_refPeptide.fa', sep = ''))[2]
  pep.long.mu <- readLines(paste(gene,'_targeted.fa', sep = ''))[2]
  file_long_pep <<- file_long_pep
  
  m <- 1
  for (m in 1:nchar(pep.long.ref)){
    c.ref <- substring(pep.long.ref,m,m )
    c.mu <- substring(pep.long.mu,m,m )
    if (c.ref != c.mu){
      break
    }
  }
  if (type == 'mutant'){
    write(paste(gene, 'unspliced', 'ref', pep.long.ref, sep = '\t'),append = TRUE, file = file_long_pep)
    write(paste(gene, 'unspliced', 'mut', pep.long.mu , sep = '\t'),append = TRUE, file = file_long_pep)
  }
  
  length.ref <- nchar(pep.long.ref)
  length.mu <- nchar(pep.long.mu)
  if (m <= 10){
    if (length.ref - m >= 10){
      pep.long.ref <- substring(pep.long.ref,1,m+10)
    }
  }else{
    if (length.ref - m <= 10){
      pep.long.ref <- substring(pep.long.ref,m-10,length.ref)
    }else{
      pep.long.ref <- substring(pep.long.ref,m-10,m+10)
    }
    pep.long.mu <- substring(pep.long.mu,m-10,length.mu)
  }
  
  if (type == 'mutant'){
    write(paste(gene, 'spliced', 'ref', pep.long.ref, sep = '\t'),append = TRUE, file = file_long_pep)
    write(paste(gene, 'spliced', 'mut', pep.long.mu , sep = '\t'),append = TRUE, file = file_long_pep)
  }
  
  return(c(pep.long.ref,pep.long.mu))
}

## functions:
list2df <- function(l1){
  # tb <- 0
  for (j in 1:length(l1)){
    if (j == 1){
      tb <- l1[[j]]
    }else{
      (l1[[1]])
      tb <- rbind(tb, l1[[j]])
    }    
  }
  if (length(l1) < 2){
    tb <- rbind(tb, tb)
  }
  tb <- data.frame(tb, stringsAsFactors=FALSE)
  
  return(tb)
}

remove_empty_column <- function(l1,c1){
  for (i in 1:length(l1)){
    l1[[i]] <- l1[[i]][c1]
  }
  return(l1)
}

add_colnames <- function (tb, cols) {
  if (is.numeric(tb)){
    return(tb)
  }
  if (!is.data.frame(tb)){
    stop('debug')
  }
  rownames(tb) <- c(1:nrow(tb))
  colnames(tb) <- cols
  return(tb)
}

write_to_file <- function(tb, type, mhc){
  if (!is.numeric(tb)){
    tb.wt1 <- merge(tb, indelInfo, by.x = 'gene', by.y = 'gene')
    tb.wt1 <- unique(tb.wt1)
    file1 <- paste('../indel_',type,'_',mhc,'.tsv', sep = '')
    write.table(tb.wt1, file = file1, sep = '\t', row.names = FALSE, quote = FALSE)
  }
}

get_tb <- function (type) {
  tb0.1 <- NULL
  tb0.2 <- NULL
  for (i in 1:length(genes)){
    gene <- genes[i]
    peps <- prepare_long_peps(gene, type)
    pep.long.ref <- peps[1]
    pep.long.mu <- peps[2]
    
    if (type == 'mutant'){
      wt <- paste(gene, '_neoepitope_indel_mutation', sep = '')
    }else{
      wt <- paste(gene, '_epitope_indel_wildtype', sep = '')
    }
    if (!file.exists(wt)) next
    
    wt1 <- readLines(wt)
    l1 <-str_split(wt1[str_detect(wt1, pattern = '^(?=.*HLA).*<=')], pattern = '\\s+')
    l2 <-str_split((wt1[str_detect(wt1, pattern = '^(?!.*HLA).*<=')]), pattern = '\\s+')
    
    ## functions  
    check_indel <- function(pep) {
      d.mu <- str_detect(pep.long.mu, pattern = pep)
      d.ref <- str_detect(pep.long.ref, pattern = pep)
      if(type == 'mutant' & d.mu & !(d.ref)){
        return(TRUE)
      }else if(type == 'wildType' & d.ref){
        return(TRUE)
      }else {
        return(FALSE)
      }
    }
    
    mu_ref_classify <- function(tb1){
        peps <- unique(tb1[,2])
        peps2 <- peps[sapply(peps, check_indel,USE.NAMES = FALSE)]
        tb2 <- tb1[tb1[,2] %in% peps2,]
        if (nrow(tb2) < 1){
          return(NULL)
        }
        return(tb2)
    }
    
    ## execution
    if (length(l1) > 0)  {
      l1 <- remove_empty_column(l1, c(3:17,19))
      tb1 <- list2df(l1)
      tb1 <- mu_ref_classify(tb1)
      if (is.null(tb0.1)){
        tb0.1 <- tb1
      }else{
        tb0.1 <- rbind(tb0.1, tb1)
      }
    }
    
    if (length(l2) > 0)  {
      c1 <- c(3,4,8,5,6,7,12,13,14,11,16)
      # l2 <- remove_empty_column(l2, c(3:13))
      l2 <- remove_empty_column(l2, c1)
      tb2 <- list2df(l2)
      tb2 <- mu_ref_classify(tb2)
      
      if (is.null(tb0.2)){
        tb0.2 <- tb2
      }else{
        tb0.2 <- rbind(tb0.2, tb2)
      }
    }

    
  }
  
  cols1 <- c('MHC','Peptide','Core','Of','Gp','Gl','Ip','Il','Icore','gene','Score_EL','%Rank_EL','Score_BA','%Rank_BA','Aff(nM)','BindLevel')
  cols2 <- c('Allele','Peptide','gene','Pos','Core','Core_Rel','1-log50k(aff)','Affinity(nM)','%Rank','Exp_Bind','BindingLevel')
  if (!is.null(tb0.1)){
    tb0.1 <- add_colnames(tb0.1, cols1)
    write_to_file(tb0.1, type, 'MHCI')
  }
  if (!is.null(tb0.2)){
    tb0.2 <- add_colnames(tb0.2, cols2)  
    write_to_file(tb0.2, type, 'MHCII')
  }
  print(dim(tb0.1))
  print(dim(tb0.2))
  print('================  ================')
}


f2 <- list.files('../', pattern = '.*\\.vcf')[1]
f2 <- paste('../', f2, sep = '')
indelInfo <- read.table(f2, stringsAsFactors = FALSE, sep = '\t',
                        comment.char = "",
                        header = TRUE
                        )[,c('POS', 'REF', 'ALT', 'ANNOVAR_FUNCTION', 'GENE')]
indelInfo$GENE <- str_match(indelInfo$GENE, pattern = '([^(]+)\\(?.*')[,2]
colnames(indelInfo) <- c('indel_position','reference','mutation','genomic_location','gene')
tb.geneID <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                   filters = "external_gene_name", 
                   values = as.character(indelInfo$gene),
                   useCache = FALSE,
                   # values = 'AOAH',
                   # verbose = TRUE,
                   # uniqueRows = FALSE,
                   # mart = ensembl)
                   mart = old_mart)
indelInfo <- merge(indelInfo, tb.geneID, by.x = 'gene', by.y = 'external_gene_name', all.x = TRUE)[,c(1,6,2:5)]



## output:


get_tb('wildType')
get_tb('mutant')

