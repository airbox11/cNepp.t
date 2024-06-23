library("stringr")

## To visualize the predicted epitopes from netMHCpan, if netchop runs well, add this funtion
# source("/icgc/dkfzlsdf/analysis/D120/scripts/netchop/run_netchop.R")


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least two arguments must be supplied (input file and output file) 
        >>> Usage: Rscripts your_script file.tab file_binding outputFile \n", call.=FALSE)
} else if (length(args)==3) {
  fTab <- args[1]
  fBinding <- args[2]
  outputFile <- args[3]
}

# "" refers to any length whitespace as being the delimiter
#myFile<-read.table(fTab,header = T,sep = "\t")
myFile<-read.delim(fTab,header = T,sep = "\t")
myBinding<-read.table(fBinding,header = F,sep = "",fill=T)

myFilermDup <- myFile[!duplicated(myFile[,c("geneName","targetPeptide","refPeptide")]),]
# oder by epitope position in peptide
myBinding <- myBinding[order(myBinding[,1]),]

# empty 
sumOut=c()
for (i in 1:nrow(myFilermDup)){
  tPep <- as.character(myFilermDup[i,c("targetPeptide")])
  rPep <- as.character(myFilermDup[i,c("refPeptide")])
  pepDescription <- as.character(paste(myFilermDup[i,c("geneName")],myFilermDup[i,c("aaChange")],myFilermDup[i,c("mutationType")],sep=","))
  outContent <- rbind(cbind(tPep,pepDescription), cbind(rPep,"Ref"))

# deactivate netchop part
  ## call netchop, default score 0.5
#  netchopRes <- funNetchop(pepDescription,tPep)
#  cleavageSite2 <- paste(netchopRes[2,],collapse = "")
#  netchopRow2 <- cbind(cleavageSite2,"CS 0.5")
  # set up threshold 0.8
#  netchopRes[2, !netchopRes[1, ] >= 0.8] <- "."
#  cleavageSite1 <- paste(netchopRes[2,],collapse = "")
#  netchopRow1 <- cbind(cleavageSite1,"CS 0.8")
  
  targetPepLen <- nchar(tPep)
  for (j in 1:nrow(myBinding)){
    neoepitope <- as.character(myBinding[j,3])
	# find the location of neoepitope in targeted neopeptides
    location<-str_locate(pattern=neoepitope,tPep)
	# "*" is replaced with "X" in netMHCpan, therefore some neoepitope (with "X" end) not able to match to neopeptide
    if (! is.na(location[1,1])){
      startPos <- location[1,1]
      endPos <- location[1,2]
      leftSide <- paste(rep("-",startPos-1),collapse = "")
      rightSide <- paste(rep("-",targetPepLen-endPos),collapse = "")
      neoRecord <- paste(leftSide,neoepitope,rightSide,sep = "")
	  # extract binding affinity from netMHCpan results, consider both class I and class II
	  affinity <- round(ifelse(grepl("HLA",myBinding[j,2]),myBinding[j,13],myBinding[j,9]), digits=0)
      neoRecordLine <- cbind(neoRecord,paste(as.character(myBinding[j,2]),as.character(affinity),sep=" ,"))
      outContent <- rbind(outContent, neoRecordLine)
    }
  }
  emptyR <- cbind(" "," ")
  # deactivate netchop 
  # outContent <- rbind(netchopRow1, netchopRow2, outContent, emptyR)
  outContent <- rbind(outContent, emptyR)
  sumOut <- rbind(outContent, sumOut)
  # convert string in first column to list
  #strContent <- data.frame(do.call(rbind,strsplit(outContent[,1],"")),stringsAsFactors=FALSE)
  #strContent <- cbind(strContent,outContent[,2])
 }

write.table(sumOut,file=outputFile,sep = "\t",row.names = F, col.names = F, quote = F)


