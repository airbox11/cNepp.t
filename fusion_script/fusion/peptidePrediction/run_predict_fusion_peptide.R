
## This R script is to translate the fusion transcripts from confFuse prediction.

lenAA=14

args<-commandArgs(TRUE)

if (length(args) != 1){
    args <- c("--help")
}

## Help section
if("--help" %in% args) {
cat("
    Required one argument:
    arg1: Fustion transcript, breakpoint | must be provided, quote must be used

    Example:	
    Rscript run_predict_fusion_peptide.R \"AAGATGAAGACAATGTTCTCAAAGCATTTACAGTACCTAAAAACAGGTCCCTGGCGAGCCCCTTGCAG|CTGAGGATTTGTGACTGGACCATGAATCAAAACGGCAGGCATTTATACCCCAGTA\"
        \n\n")
    q(save="no")
}

fusionSeq=args[1]

suppressMessages(library("biomaRt"))
suppressMessages(library("gdata"))
suppressMessages(library("Biostrings")) # for using function DNAString # source("https://www.bioconductor.org/biocLite.R"),  biocLite("")
#install.packages("seqinr")
suppressMessages(library("seqinr"))
suppressMessages(library("stringr"))


# reverse complement, for complicated cases
# reverseComplement(DNAString("TAAAGAAGAGGCGGTAATTCTGCGAGCAGGGCTGGCCGCGCTCCTCAGTGTGGTACAGGGCCAT"))

## Break point position
breakPos <- regexpr("\\|",fusionSeq)[1]
selSeq1 <- substr(fusionSeq,breakPos-lenAA*3,breakPos-1)
selSeq2 <- substr(fusionSeq,breakPos+1,breakPos+1+lenAA*3)
predictedAA1<-c2s(translate(s2c(selSeq1)))
#print(selSeq1)
if (nchar(selSeq1) %%3 != 0){
  print("Error in sequence extraction")
}
predictedAA2<-c2s(translate(s2c(selSeq2)))
#predictedAA <- c2s(translate(s2c(paste(selSeq1,selSeq2,sep=""))))
#print(predictedAA)

## This script is to translate the fusion transcripts from confFuse prediction.
#print(predictedAA1)
#print(predictedAA2)

#print(predictedAA)
print(paste(selSeq1,selSeq2,predictedAA1,predictedAA2))




