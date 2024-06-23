#!/software/r/3.6.2/lib64/R/bin/R
library(biomaRt)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)
workDir <- args[1]

if (length(args) == 2){
  vcfOnly <- args[2]
}
## test ====
  # workDir <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/OE0370_CRC_4251_pdo01_germline01/2_SNVs_based_neoepitope_prediction'

  # vcfOnly <- 'origin'
  # vcfOnly <- 'promise'

## test end

setwd(workDir)
file1 <- '1.vcf'
file.output <- '2.vcf'
df1 <- read.table(file = file1, sep = '\t', stringsAsFactors = FALSE, header = TRUE, fill = TRUE)

if (vcfOnly == 'promise') {
  old_mart <- useEnsembl(biomart = "ensembl", 
                        dataset = "hsapiens_gene_ensembl")
} else if (vcfOnly=='pathology' | vcfOnly == 'origin') {
  old_mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                      host="grch37.ensembl.org",
                      path="/biomart/martservice",
                      dataset="hsapiens_gene_ensembl"
  )
}
attributes <- c("ensembl_gene_id", "hgnc_symbol")
filters <- c("chromosome_name","start","end")

for (i in 1:nrow(df1)){
  chr <- df1[i,]$CHROM
  chr <- str_match(chr, pattern = '(chr)?(.*)')[,3]
  pos <- as.character(df1[i,]$POS)
  values <- list(chromosome=chr,start=pos,end=pos)

  tryCatch(
    expr = {
      all.genes <- getBM(attributes=attributes, filters=filters, values=values, mart=old_mart, useCache = FALSE)
      if (nrow(all.genes)==1) {
        df1[i,]$geneID <- all.genes$ensembl_gene_id
        df1[i,]$GENE <- all.genes$hgnc_symbol
      } else {
        all.genes <- all.genes[order(all.genes$hgnc_symbol, decreasing = TRUE),]
        df1[i,]$geneID <- all.genes[1,]$ensembl_gene_id
        df1[i,]$GENE <- all.genes[1,]$hgnc_symbol
      }
    },
    error = function(e){ 
      print('/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/prepare_vcf.r')
      print('debug: getBM could not find relative gene.')
    }
  )
}
df1.1 <- df1[!is.na(df1$geneID),]
write.table(df1.1, file = file.output, sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
