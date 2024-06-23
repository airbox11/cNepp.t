
### extract fusion peptides across break points
setwd("/icgc/dkfzlsdf/analysis/D120/project/fusion_liposarcoma/scripts/")

neoDir = "/icgc/dkfzlsdf/analysis/D120/project/fusion_liposarcoma/outputDir/"
setwd(neoDir)
files <- list.files(path=neoDir, pattern = "*visualization")

myOut <- c()
for (file in files){
  myDat <- read.delim(file)
  mySub <- str_split_fixed(myDat[2,2],"_",4)
  lineID<-mySub[3]
  breakPos <- mySub[4]
  print(file)
  if (nrow(myDat)>3){
    eLine <- cbind(lineID,paste(mySub[1],mySub[2],sep = "_"))
    pepCol=""
    for (i in 4:nrow(myDat)){
      epitope <- as.character(myDat[i,1])
      if (substr(epitope,breakPos,breakPos) != "-" & nchar(epitope)>1){
        pep1 <- sub(" ","",paste(myDat[i,2], gsub("-","",epitope),sep = ","))
        pepCol <- paste(pepCol,pep1,sep=";")
      }
    }
    eLine <- cbind(eLine, sub(";","",pepCol))
  }
  myOut <- rbind(myOut,as.data.frame(eLine))
  # n=n-1
  # if(n==0){
  #   break
  #   }
}

colnames(myOut) <- c("ID","sampleID","epitopePredictionList")
write.table(myOut,file="fusion_epitope_report.tab",col.names = T, row.names = F, quote = F, sep = "\t")

f1 <- "/icgc/dkfzlsdf/analysis/D120/project/fusion_liposarcoma/my_modified_fusion.tsv"
myF1 <- read.delim(f1,header =T)

mdat <- merge(myF1, myOut, by="ID")
write.table(mdat,file="final_report.tab",col.names = T, row.names = F, quote = F, sep = "\t")

