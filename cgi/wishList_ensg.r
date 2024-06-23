library(biomaRt)

file.wishList <- '/omics/groups/OE0422/internal/yanhong/documentation/NCT_gene_list_20210511.txt'
df.wish <- read.table(file.wishList, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
gene.wish <- unique(df.wish$GENE)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
df.wish <- getBM(attributes = c('ensembl_gene_id',
                            'external_gene_name'),
             filters = 'external_gene_name',
             values = gene.wish,
             useCache = FALSE,
             mart = mart)
df.wish$wishList <- TRUE

saveRDS(df.wish, file = '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/cgi/wishList_ensg.rdata')
