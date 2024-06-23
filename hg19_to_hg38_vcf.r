args <- commandArgs(trailingOnly=TRUE)

file.vcf <- args[1]
file.addRNA1 <- args[2]
file.output <- args[3]

## test
  # workDir <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/promise/batch_process_20230302/result/S014-2NHL7W_T1T2_tumor11/3_add_expression'
  # setwd(workDir)
  # file.vcf <- './3.vcf'
  # file.addRNA1 <- './addRNA1'
  # file.output <- 'addRNA1.5'
## test end.

df1 <- read.table(file = file.vcf, sep = '\t', stringsAsFactors = FALSE, header = TRUE)
df2 <- read.table(file = file.addRNA1, sep = '\t', stringsAsFactors = FALSE, header = TRUE)
df2 <- df2[match(df1$POS, df2$loc),]

df2$chr <- df1$chr_38
df2$loc <- df1$pos_38

write.table(df2, file = file.output, col.names = TRUE, sep = '\t', quote = FALSE, row.names = FALSE)
