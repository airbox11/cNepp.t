library(stringr)

args <- commandArgs(trailingOnly=TRUE)
workDir <- args[1]
  # workDir <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/29056/4_indel_based_prediction'
setwd(workDir)


file1 <- list.files(workDir, pattern = 'indel_tmp.vcf')
file.output <- 'indel_somatic.vcf'

df1 <- read.table(file = file1, sep = '\t', stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
colnames(df1)[1:11] <- c('X.CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','CONTROL','TUMOR')
df2 <- df1[,c(colnames(df1)[1:11], "GENE")]

write.table(df2, file = file.output, sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)


