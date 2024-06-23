
# Input class I and class II neoantigen prediction results
# FPKM or Median column required
# Output results
library(stringr)

# make sure the netchop runs successfully with individual environment
source("/icgc/dkfzlsdf/analysis/D120/scripts/netchop/run_netchop.R")
#source("/home/huangz/scripts/immunoinfo/neoepitope_indels/run_netchop.R")

args<-commandArgs(TRUE)
if (length(args) != 3){
    args <- c("--help")
}


## Help section
if("--help" %in% args) {
cat("
    Required three arguments:
    arg1: Path to MHC-I prediction result file
    arg2: Path to MHC-II prediction result file
    arg3: Path to your output file
	
Example:
Rscript run_visualization_neoepitope_prediction.R classI_file classII_file output_file
        \n\n")
    q(save="no")
}

classI_file <- args[1]
classII_file <- args[2]
outFile <- args[3]

#classI_file  <- "results_AISHA_MHCI_epitopes_filtered_addRNA.tab"
#classII_file <- "results_AISHA_MHCII_epitopes_filtered_addRNA.tab"
#classI_file<- "results_H021-B34BBM_MHCI_epitopes_filtered_addRNA.tab"
#classII_file <- "results_H021-B34BBM_MHCII_epitopes_filtered_addRNA.tab"

classI <- read.delim(classI_file,header = T,sep="\t")
classII <- read.delim(classII_file, header= T,sep="\t")


# rename HLA allele column
colnames(classII)[colnames(classII) == "Allele"] <- "HLA"
colnames(classII)[colnames(classII) == "Mut_seq"] <- "Mut_pos"
colnames(classI)[colnames(classI) == "VARIANTREADS_RNA"] <- "expAlt"
colnames(classII)[colnames(classII) == "VARIANTREADS_RNA"] <- "expAlt"
colnames(classI)[colnames(classI) == "RPKM"] <- "FPKM"
colnames(classII)[colnames(classII) == "RPKM"] <- "FPKM"



 #select columns
# ANNOTATION_RNA instead of expAlt in INFORM
if ("FPKM" %in% colnames(classI)){
  c1 <- classI[,c("mutPeptide","refPeptide","GENE","Peptide","aaChange","Mut_pos","HLA","Aff.nM.","FPKM","expAlt")]
	
}else{
  # for sample without RNA-seq data, taken expression from TCGA
  if ("Median" %in% colnames(classI)){
	c1 <- classI[,c("mutPeptide","refPeptide","GENE","Peptide","aaChange","Mut_pos","HLA","Aff.nM.","Mean","Median")]
	}
}

if ("FPKM" %in% colnames(classII)){
  c2 <- classII[,c("mutPeptide","refPeptide","GENE","Peptide","aaChange","Mut_pos","HLA","Affinity.nM.","FPKM","expAlt")]
}else{ 
  if ("Median" %in% colnames(classI)){
	c2 <- classII[,c("mutPeptide","refPeptide","GENE","Peptide","aaChange","Mut_pos","HLA","Affinity.nM.","Mean","Median")]
	}
}

# unify expression name
colnames(c1)[colnames(c1) == "FPKM"] <- "expression"
colnames(c2)[colnames(c2) == "FPKM"] <- "expression"
colnames(c1)[colnames(c1) == "Mean"] <- "expression"
colnames(c2)[colnames(c2) == "Mean"] <- "expression"
# unify affinity binding in mutant peptide
colnames(c1)[colnames(c1) == "Aff.nM."] <- "Aff_mut"
colnames(c2)[colnames(c2) == "Affinity.nM."] <- "Aff_mut"


# order by expression
mDat <- rbind(c1,c2)
mDat <- mDat[order(-mDat$expression),]

# remove duplicate
mDat <- unique(mDat)
# remove duplicates due to multiple/different aaChange
mDat <- mDat[!duplicated(mDat[,c('mutPeptide','Peptide','HLA')]),]

#TODO: make a dictionary to visulize results.
# define each output unit
neopeptide <- unique(mDat$mutPeptide)

outNeo <- c()
# choose all neoepitope 
for (n in neopeptide){
  sel <- mDat[n == mDat$mutPeptide,]
  # order by neoepitope position
  sel <- sel[order(sel$Mut_pos),]
 
  # create reference peptide row
  outPutRow1 <- cbind(as.character(sel[1, c("refPeptide")]),"Ref")
  #colnames(outPutRow1) <- c("Peptide","Comments")
  
  # create mutated peptide row
  row2Comments <- paste(as.character(sel[1, c("GENE")]),as.character(sel[1, c("aaChange")]), as.character(round(sel[1, c("expression")],digits=1)), as.character(sel[1, c("expAlt")]), sep=",")
  outPutRow2 <- cbind(n,row2Comments)
  #colnames(outPutRow2) <- c("Peptide","Comments")

  # call netChop
  pepDescription<-paste(outPutRow2[1,1],outPutRow2[1,2],sep = ",")
  netchopRes <- funNetchop(pepDescription,n)
  cleavageSite2 <- paste(netchopRes[2,],collapse = "")
  netchopRow2 <- cbind(cleavageSite2,"0.5")
  # set up threshold 0.8
  netchopRes[2, !netchopRes[1, ] >= 0.8] <- "."
  cleavageSite1 <- paste(netchopRes[2,],collapse = "")
  netchopRow1 <- cbind(cleavageSite1,"0.8")
	
  colnames(outPutRow1)=colnames(outPutRow2)=colnames(netchopRow1)=colnames(netchopRow2)=c("Peptide","Comments")
 
  # combination
  outPut <- rbind(netchopRow1,netchopRow2,outPutRow1,outPutRow2)

  for (i in 1:nrow(sel)){
	myPos<-str_locate(n, as.character(sel[i,c("Peptide")]))
	if (! is.na(myPos[1,1])){
		leftSide <- paste(rep("-",myPos[1,1]-1),collapse = "")
		rightSide <- paste(rep("-",nchar(n)-myPos[1,2]),collapse = "")
    	recordNeo <- paste(leftSide,sel[i,c("Peptide")],rightSide,sep = "")
	    recordNeo <- cbind(recordNeo, paste(as.character(sel[i,c("HLA")]),as.character(round(sel[i,c("Aff_mut")]),digits=0),sep=","))
    	colnames(recordNeo) <- c("Peptide","Comments")
	    outPut <- rbind(outPut, recordNeo)
	}
  }
  emptyR <- cbind(" "," ")
  colnames(emptyR) <- c("Peptide","Comments")
  outNeo <- rbind(outNeo,emptyR)
  outNeo <- rbind(outNeo, outPut)
}

write.table(outNeo,file=outFile,sep = "\t",row.names = F, col.names = F, quote = F)

  

