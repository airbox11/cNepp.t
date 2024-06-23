library(stringr)
library(biomaRt)
library(limma)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly=TRUE)
workDir <- args[1]

input_sf <- args[2]
output_sf <- args[3]

input_s2 <- args[4]
output_s3 <- args[5]



  ## test:
  # workDir <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/p28202/3_add_expression'
  # input_sf <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/p28202/2_SNVs_based_neoepitope_prediction/E-2023-28208.quant.sf'
  # output_sf <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/p28202/3_add_expression/fpkm_tpm.featureCounts.tsv'
  # input_s2 <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/p28202/2_SNVs_based_neoepitope_prediction/netMHCpan4_1/results_MHCI_epitopes.tab_splitGenes'
  # output_s3 <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/p28202/3_add_expression/MHCI_epitopes_RNAseq_netMHCpan4_1.tab'

  # workDir <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/42215/3_add_expression'
  # input_sf <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/42215/2_SNVs_based_neoepitope_prediction/E-2023-42215_RNA_71792_Dx.quant.sf'
  # output_sf <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/42215/3_add_expression/fpkm_tpm.featureCounts.tsv'
  # input_s2 <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/42215/2_SNVs_based_neoepitope_prediction/netMHCpan4_1/results_MHCI_epitopes.tab_splitGenes'
  # output_s3 <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/42215/3_add_expression/MHCI_epitopes_RNAseq_netMHCpan4_1.tab'

# workDir <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/p005144/3_add_expression'
# input_sf <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/p005144/2_SNVs_based_neoepitope_prediction/E-2024-05144_RNA_74920_Dx.quant.sf'
# output_sf <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/p005144/3_add_expression/fpkm_tpm.featureCounts.tsv'
# input_s2 <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/p005144/2_SNVs_based_neoepitope_prediction/netMHCpan4_1/results_MHCI_epitopes.tab_splitGenes'
# output_s3 <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/p005144/3_add_expression/MHCI_epitopes_RNAseq_netMHCpan4_1.tab'


setwd(workDir)

## get TPM for genes ====
df1 <- read.table(file=input_sf, sep='\t', stringsAsFactors = FALSE, header = TRUE)
df1 <- df1[,c(1,4)]
colnames(df1) <- c('transID', 'TPM')
df1$transID <- str_match(df1$transID, pattern = '(.*)\\.')[,2]
df1 <- df1[!df1$TPM == 0,]

mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")

if (!file.exists('df.bm.rdata')){
  df.bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "ensembl_transcript_id"),
        filters   = "ensembl_transcript_id",
        values    = df1$transID,
        # values    = df1$transID[1:100], 
        mart      = mart,
        useCache  = FALSE)
  saveRDS(df.bm,file='df.bm.rdata')
}else{
  df.bm <- readRDS(file = 'df.bm.rdata')
}

df2 <- merge(df1, df.bm, by.x='transID', by.y='ensembl_transcript_id', all.x=TRUE)
df2 <- df2[!is.na(df2$ensembl_gene_id),]
  df2.1 <- df2[,c('ensembl_gene_id','hgnc_symbol','TPM')]
  names(df2.1)[names(df2.1) == "hgnc_symbol"] <- "name"
  write.table(df2.1, file=output_sf, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
df4 <- aggregate(TPM~ensembl_gene_id,df2,sum)
  names(df4)[names(df4) == "ensembl_gene_id"] = "geneID"
  df4$FPKM <- NA

  
## get geneIDs ====
df.s2 <- read.table(file=input_s2, sep='\t', stringsAsFactors = FALSE, header = TRUE)

# v1 <- unique(df.s2$GENE)
get_geneID <- function(v1){
  df.bm2 <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", 'chromosome_name'),
                 filters   = "hgnc_symbol",
                 values    = unique(v1), 
                 mart      = mart,
                 useCache  = FALSE)
  chrs <- c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y')
  df.bm3 <- df.bm2[df.bm2$chromosome_name%in%chrs,]
  df.bm3 <- df.bm3[!duplicated(df.bm3),][,c(1,2)]
  colnames(df.bm3) <- c('geneID','GENE')
  return(df.bm3)
}
if ('geneID' %in% colnames(df.s2)){
  df.s2.1 <- df.s2
}else{
  df.bm3 <- get_geneID(df.s2$GENE)
  df.s2.1 <- merge(df.s2, df.bm3, by.x='GENE', by.y='GENE', all.x=TRUE)
}

## merge for TPM ====
df6 <- merge(df.s2.1, df4, by.x='geneID', by.y='geneID', all.x=TRUE)
df6$chrPos <- paste0(df6$CHROM, '_', df6$POS)
  names(df6)[names(df6) == "CHROM"] = "X.CHROM"
  cols1 <- c("sumReads","numOfBaseExp","dna_freAlt","dna_cov","freRef","freAlt","expAlt")
  for (i in cols1){
    df6[[i]] <- NA
  }
  write.table(df6, file=output_s3, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
