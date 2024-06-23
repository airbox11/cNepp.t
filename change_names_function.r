library(stringr)


df.name <- read.table('/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/reName.tsv', header = FALSE, row.names = NULL, sep = '\t', stringsAsFactors = FALSE, quote="")
colnames(df.name) <- c('name_old', 'name_new')
change.names <- function(df1) {
  
  for (name_old in df.name$name_old) {
    if (name_old %in% colnames(df1)){
      new_name <- df.name$name_new[which(df.name$name_old == name_old)]
      colnames(df1)[which(names(df1) == name_old)] <- new_name
      
      detect1 <- str_detect(new_name, pattern = 'Ensembl_Gene_ID')
      if (detect1) {
        df1[,new_name] <- str_match(df1[,new_name], pattern = 'ENSG\\d+')
      }
    }
  }
  if ('Official_Gene_Symbol_1' %in% colnames(df1)) {
    colname.1 <- 'Official_Gene_Symbol_1'
  }
  else{
    colname.1 <- 'Official_Gene_Symbol'
  }
  df1 <- df1[order(df1[,colname.1]),]
  return(df1)
}