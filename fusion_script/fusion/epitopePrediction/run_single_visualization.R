library("stringr", quietly = TRUE)

## To visualize the predicted epitopes from netMHCpan
source("/icgc/dkfzlsdf/analysis/D120/scripts/netchop/run_netchop.R")

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least two arguments must be supplied (input file and output file) 
        >>> Usage: Rscripts your_script fasta_peptide input_fusion_netMHC_results outputFile \n", call.=FALSE)
} else if (length(args)==3) {
  pepFasta <- args[1]
#  fTab <- args[1]
  netResult <- args[2]
  outputFile <- args[3]
}

### select WB and SB in netMHC results
tmpfile1 <- "__tmp_binding"
command1 <- paste("grep =", netResult, ">", tmpfile1, sep = " ")
system(command1)

# "" refers to any length whitespace as being the delimiter
myBinding<-read.table(tmpfile1,header = F,sep = "",fill=T)
# oder by epitope position in peptide
myBinding <- myBinding[order(myBinding[,1]),]

# read fasta files
library("Biostrings", quietly = TRUE)
myPep <- readAAStringSet(pepFasta)
pepDescription <- names(myPep)
pepSeq <- paste(myPep)

## call netchop, default score 0.5
netchopRes <- funNetchop(pepDescription,pepSeq)
cleavageSite2 <- paste(netchopRes[2,],collapse = "")
netchopRow2 <- cbind(cleavageSite2,"CS 0.5")
# set up threshold 0.8
netchopRes[2, !netchopRes[1, ] >= 0.8] <- "."
cleavageSite1 <- paste(netchopRes[2,],collapse = "")
netchopRow1 <- cbind(cleavageSite1,"CS 0.8")

# length of peptide sequence
pepSeqLen<-nchar(pepSeq)

#
outContent <- cbind(pepSeq,pepDescription)

# neoepitope alignment
for (j in 1:nrow(myBinding)){
	neoepitope <- as.character(myBinding[j,3])
	location<-str_locate(pattern=neoepitope,pepSeq)
	if (! is.na(location[1,1])){
		startPos <- location[1,1]
		endPos <- location[1,2]
		leftSide <- paste(rep("-",startPos-1),collapse = "")
		rightSide <- paste(rep("-",pepSeqLen-endPos),collapse = "")
		neoRecord <- paste(leftSide,neoepitope,rightSide,sep = "")
		affinity <- round(ifelse(grepl("HLA",myBinding[j,2]),myBinding[j,13],myBinding[j,9]), digits=0)
		neoRecordLine <- cbind(neoRecord,paste(as.character(myBinding[j,2]),as.character(affinity),sep=" ,"))
		outContent <- rbind(outContent, neoRecordLine)
	}
}

emptyR <- cbind(" "," ")
outContent <- rbind(netchopRow1, netchopRow2, outContent, emptyR)

write.table(outContent,file=outputFile,sep = "\t",row.names = F, col.names = F, quote = F)

command2=paste("rm ", tmpfile1)
system(command2)

