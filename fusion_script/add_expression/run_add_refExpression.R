### This script is to add reference expression from TCGA
### Furthermore, to calculate the DNA mutant fraction

library(stringr)

args<-commandArgs(TRUE)

if (length(args) != 3){
	args <- c("--help")
}

## Help section
if("--help" %in% args) {
cat("
	Required three arguments:
	arg1: Path to TCGA expression file
	arg2: Path to neoepitope prediction file
	arg3: Path to your output file

	Example:
	Rscript run_add_refExpression.R /to/your/TCGA-BRCA/TCGA-BRCA_expression_addName.tab /to/your/results_HB60_MHCI_epitopes_filtered.tab /to/my/results_HB60_MHCI_epitopes_filtered_TCGA-BRCA.tab
		\n\n")
	q(save="no")
}

rpkm=args[1]
neoepitope=args[2]
output=args[3]

mat <- read.table(rpkm,header=FALSE)
colnames(mat)[1:4]<-c("GeneID","GENE","Mean","Median")

mymat<-mat[,1:4]

mhcPrediction <- read.delim(neoepitope,header=TRUE, sep="\t")

# mutant allele frequence
dp4=str_split_fixed(str_split_fixed(mhcPrediction$INFO, "DP4=", 2)[,2],";",2)[,1]
mutFre=(as.integer(str_split_fixed(dp4,",",4)[,3])+as.integer(str_split_fixed(dp4,",",4)[,4])) / (as.integer(str_split_fixed(dp4,",",4)[,1])+as.integer(str_split_fixed(dp4,",",4)[,2]) + as.integer(str_split_fixed(dp4,",",4)[,3])+as.integer(str_split_fixed(dp4,",",4)[,4]))

mhcPrediction <- cbind(mhcPrediction,mutFre)

mhcPrediction=as.matrix(mhcPrediction)

for (i in 1:nrow(mhcPrediction)){
	geneName=as.character(mhcPrediction[i, "GENE"])
	if (grepl(",",geneName)){
		name=strsplit(geneName,",")[[1]][1]
		mhcPrediction[i, "GENE"]<-name
	}
}


write.table(merge(mhcPrediction,mymat,by="GENE",all.x=TRUE), file=output,quote=FALSE,sep="\t",row.names=FALSE)

