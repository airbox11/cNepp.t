
# this script is to 
# choose candidates for vaccination

#setwd("/home/kosalogl/Elispots_Summary/selecte_epitope/")

#fileName="results_H021-RK6XZV_MHCI_epitopes_filtered_addRNA.tab_2018"

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3){
    args <- c("--help")
}
## Help section
if("--help" %in% args) {
cat("
    Required three arguments:
    arg1: Path to neoepitope result file
    arg2: Path to your new file
    arg3: Defined ID, e.g. ID-mhcI or ID-mhcII
    
 >>> Usage: Rscripts parse_neoepitope_rna.R inputFile outputFile definedID 
		\n\n")
    q(save="no")
}
# Rscript parse_neoepitope_rna.R /icgc/dkfzlsdf/analysis/D120/immunoinfor/H021/results_per_pid/H021-99G9EH/immuno_201705_tumor03/results_H021-99G9EH_MHCII_epitopes_filtered_add_RNA01.tab test_output_H021-99G9EH_MHCII H021-99G9EH_MHCII

# default output file
fileName <- args[1]
outputName <- args[2]
sampleID <- args[3]

f1 <-  read.delim(fileName,header=TRUE, sep="\t", stringsAsFactors=FALSE)
sampleID <- as.character(sampleID)

if("FPKM" %in% colnames(f1)){
   myDat <- unique(cbind(f1$GENE, f1$TPM,f1$FPKM,f1$ID,nchar(substr(f1$mutPeptide,3,27)),substr(f1$mutPeptide,3,27),substr(f1$refPeptide,3,27),f1$aaChange,f1$expAlt,sampleID,f1$mutPeptide,f1$refPeptide, f1$X.CHROM, f1$POS, f1$sumReads, f1$freAlt, f1$dna_freAlt, f1$dna_cov))
   colnames(myDat) <- c("GENE","TPM","FPKM","snpID","LengthNeoepitope","Mutated(25aa)","WildType(25aa)","aaChange","expAlt","sampleID","mutPeptide","refPeptide", "CHR", "POS", "rna_cov", "rna_freAlt", "dna_freAlt", "dna_cov")
}

if("RPKM" %in% colnames(f1)){
  myDat <- unique(cbind(f1$GENE, f1$TPM,f1$RPKM,f1$ID,nchar(substr(f1$mutPeptide,3,27)),substr(f1$mutPeptide,3,27),substr(f1$refPeptide,3,27),f1$aaChange,f1$expAlt,sampleID,f1$mutPeptide,f1$refPeptide, f1$X.CHROM, f1$POS, f1$sumReads, f1$freAlt, f1$dna_freAlt, f1$dna_cov))
  colnames(myDat) <- c("GENE","TPM","RPKM","snpID","LengthNeoepitope","Mutated(25aa)","WildType(25aa)","aaChange","expAlt","sampleID","mutPeptide","refPeptide", "CHR", "POS", "rna_cov", "rna_freAlt", "dna_freAlt", "dna_cov")
}

write.table(myDat,file=paste0(outputName,"_select",sep=""),quote=F,col.names = T,row.names = F,sep = "\t")



