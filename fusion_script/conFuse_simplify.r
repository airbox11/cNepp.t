

library('stringr')
workDir <- '/icgc/dkfzlsdf/analysis/D120/yanhong/test_fusion/S2914Nr3/4_confuse'


file1 <- paste(workDir, 'confFuse.tsv', sep = '/')
fileOutput <- paste(workDir, 'confFuse_simplify.tsv', sep = '/')

f1 <- read.table(file1, header = TRUE, stringsAsFactors = FALSE)
f2 <- f1[,c(2,21, 22, 23, 24, 25, 26, 29, 30, 31, 32, 35, 36, 42, 45, 66, 67,71)]
f3 <- f2[f2$confidence_score >= 7,]
f3$geneName <- paste(f3$gene_name1, f3$gene_name2, sep = '_')



for (i in 1:nrow(f3)){
  dna_sequence <- f3$splitr_sequence[i]
  output.fasta <- paste(workDir,'/', i, '_', f3$geneName[i], '.fasta', sep = '')
  command_get_peptide <- paste('curl -s -d "dna_sequence=',
                               dna_sequence,
                               '&output_format=fasta" https://web.expasy.org/cgi-bin/translate/dna2aa.cgi >',
                               output.fasta,
                               sep = ''
  )
  system(command_get_peptide, intern = FALSE) 
}

system()

write.table(f3, file = fileOutput, col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

## get peptide:
for (i in 1:nrow(f3)){
  
  
  dna_sequence <- dna_sequence <- f3$splitr_sequence[i]
  command_get_peptide2 <- paste('Rscript /icgc/dkfzlsdf/analysis/D120/scripts/fusion/peptidePrediction/run_predict_fusion_peptide.R "', 
                                dna_sequence, '"',
                                sep = '')
  system(command_get_peptide2, intern = TRUE)
}

get_peptide <- function(dna_sequence) {
  command_get_peptide2 <- paste('Rscript /icgc/dkfzlsdf/analysis/D120/scripts/fusion/peptidePrediction/run_predict_fusion_peptide.R "', 
                                dna_sequence, '"',
                                sep = '')
  peptide <- system(command_get_peptide2, intern = TRUE)
  return(peptide)
  
}

peptides <- sapply(f3$splitr_sequence, get_peptide, USE.NAMES = FALSE)
?sapply()
f3$splitr_sequence[1]


###

length1 <- nchar(f3$splitr_sequence) -1/3
length(f3$splitr_sequence)

str_match(f3$splitr_sequence, pattern = '([TCGA]+)(\\|)([TCGA]+)')

