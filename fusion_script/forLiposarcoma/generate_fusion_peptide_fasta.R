## This script takes fusion.tsv as input
# To generate a new file my_modified_fusion.tsv, which add rwo number as ID
# To generate a new folder "fusionPepDir" which contains individual file including fusion peptide, output file name format {sampleID}_{tumorID}_{BreakPointPosition}_{fusionID}
# To calculate recurrent fusion, result file "recurrent_fusion.pdf"


setwd("/icgc/dkfzlsdf/analysis/D120/project/fusion_liposarcoma/scripts/")

inputF<- "/icgc/dkfzlsdf/analysis/D120/project/fusion_liposarcoma/fusions.tsv"

myDat <- read.delim(inputF)

## add numbers as ID
myDat$ID <- seq.int(nrow(myDat))

library(stringr)

sampleID <- paste(str_split_fixed(myDat$filename,"/",14)[,10], str_split_fixed(myDat$filename,"/",14)[,11],myDat$ID, sep="_")

myDat$sampleID <- sampleID

## the breakPos-1 is the fusion-partner1 end position
myDat$breakPos <- str_locate(myDat$peptide_sequence,"[|]")[,1]
myDat$pepLen <- str_length(myDat$peptide_sequence)

# replace | and * with empty ""
myDat$selPep <- gsub("[*]","",gsub("[|]","",myDat$peptide_sequence))

selDat <- myDat[ ! is.na(myDat$breakPos),]

selDat$sampleID <- paste(selDat$sampleID,selDat$breakPos,sep="_")

write.table(selDat,file="my_modified_fusion.tsv", col.names=T, row.names = F, quote = F,sep="\t")

setwd("/icgc/dkfzlsdf/analysis/D120/project/fusion_liposarcoma/fusionPepDir/")

## output each fusion peptide in each file
for (i in 1:nrow(selDat)){
  myLine <- paste(">",selDat[i,c("sampleID")],sep="")
  pLine <- rbind(myLine, selDat[i,c("selPep")])
  write.table(pLine,file=paste(selDat[i,c("sampleID")],".fa",sep=""), col.names=F, row.names = F, quote = F)
}
###
# 
pdf("recurrent_fusions.pdf")
par(mar = c(20, 2, 2, 2))
myFusion <- myDat[!duplicated(myDat[,c(1:3)]),c(2:3)]
mySortFusion <- sort(table(paste(myFusion$gene1,myFusion$gene2)),decreasing = T)
plot(mySortFusion[mySortFusion>=3],las=2, cex.axis =0.7, xlab = "")
title(main="Recurrent fusions")
dev.off()


