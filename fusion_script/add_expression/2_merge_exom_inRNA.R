
# 2017-04-03, Zhiqin Huang
# this scrip is to 
# merge the file which includes the base expression in RNA-seq bam, with mpileup and converting to base txt, namely inRNA
# and the file from SNV calling, namely inExom
# command-line parameters
args <- commandArgs(TRUE)
inRNA=read.delim(args[1],header=TRUE)
#inExom=read.delim(args[2],header=TRUE)
vcf=args[2]

output=args[3]

# read two times the vcf file, first for the columns names, second for the data
tmp.vcf<-readLines(vcf)
tmp.vcf.data<-read.delim(vcf,sep="\t")

# filter for the columns names
tmp.vcf<-tmp.vcf[-(grep("#CHROM",tmp.vcf)+1):-(length(tmp.vcf))]
vcf.names<-unlist(strsplit(tmp.vcf[length(tmp.vcf)],"\t"))
names(tmp.vcf.data)<-vcf.names

inExom<-tmp.vcf.data



#cp=do.call(paste,c(inRNA[c("chr","loc")],sep="_"))

chrPos=paste(inRNA[,1],inRNA[,2],sep="_")

numOfBaseExp=rowSums(inRNA[,4:7]>0)

sumReads=apply(inRNA[,4:7],1,sum)

baseRatio=round(inRNA[,4:7]/sumReads,4)
colnames(baseRatio) <- c("rA","rT","rC","rG")

inRNA=cbind(chrPos,inRNA,baseRatio,sumReads,numOfBaseExp)

colnames(inRNA)[1]<-"chrPos"

inExom[,"chrPos"]=c(paste(inExom[,1],inExom[,2],sep="_"))

#colnames(inExom)[3]="chrPos"

m=merge(inExom,inRNA,by="chrPos")

name=basename(args[1])

#filename=unlist(strsplit(name, "[_]"))[]

write.table(m,file=output,quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")

