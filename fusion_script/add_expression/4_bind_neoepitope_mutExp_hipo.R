 

args <- commandArgs(TRUE)
file1=args[1] # mutant allele expression file
file2=args[2] # neoepitope prediction results, MHC I
file3=args[3] # TPM, FPKM, ensembleID, GeneName
outFile=args[4]	# final output


mutExp=read.delim(file1,header=TRUE)
neoepitope=read.delim(file2,header=TRUE)
exp=read.delim(file3,header=TRUE)

# extract allele expression information
#mat1=mutExp[,c(1,50:65)]
n<-ncol(mutExp)
mat1=mutExp[,c(1,(n-15):n)]

# construct ID for each mutant position
mat2=cbind(paste0(neoepitope$X.CHROM,"_",neoepitope$POS),neoepitope)
colnames(mat2)[1]<-"chrPos"

# merge mutExp and neoepitope files
mDat=merge(mat2,mat1,by="chrPos",all.x=T)


# rename column
if ( ncol(exp) == 23) {
	mydat<-cbind(as.character(exp$name),round(exp$TPM,4),round(exp$FPKM,4))
	colnames(mydat)<-c("GENE","TPM","FPKM")	

}else{colnames(exp)[2]<-"GENE"
	mydat<-exp}

mydat <- as.data.frame(mydat)
##### noted there are same gene symbol, ENSG00000182378.8 from chromosome X and Y, to be explained
# GENE    TPM   FPKM
# 54851 PLCXD1 2.6394 1.7162
# 57243 PLCXD1      0      0

# add gene expression value, for default VCF in DKFZ pipeline
if("GENE" %in% colnames(mDat)){
	mDat2=merge(mydat,mDat,by="GENE",all.y=T)
}else {print("Your neoepitope_results file does not include column name GENE") }



# output file
write.table(mDat2,file=outFile,quote=FALSE,col.names=TRUE,sep="\t",row.names=F)


