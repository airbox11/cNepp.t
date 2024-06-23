# To predict peptides
# This script doesn't work for stoploss frameshift indels.

suppressMessages(library("biomaRt"))
suppressMessages(library("gdata"))
suppressMessages(library("Biostrings")) # for using function DNAString # source("https://www.bioconductor.org/biocLite.R"),  biocLite("")
#install.packages("seqinr")
suppressMessages(library("seqinr"))
suppressMessages(library("stringr"))
#getBM(attributes = c("cdna","hgnc_symbol"),filters = c("refseq_mrna"), values = "NM_001276351", mart=ensembl)
#getBM(attributes = c("peptide","hgnc_symbol"),filters = c("refseq_mrna"), values = "NM_001276351", mart=ensembl)


#grch37 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",dataset="hsapiens_gene_ensembl")
#grch38 <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least three arguments must be supplied (input file).n", call.=FALSE)
} else if (length(args)==3) {
  # default output file
  annotationColumn <- args[1]
  aaLen <- as.integer(args[2])
  outFile <- args[3]
}

# BBS4:NM_001252678:exon7:c.A38T:p.D13V,BBS4:NM_033028:exon8:c.A554T:p.D185V,
# Amino acid length (peptide) = 29 by default
# NM_001252678, 38, A, T, 13, D, V
# accID="NM_001252678"
# codingPos=38
# refDna="A"
# mutDna="T"
# aaPos=13
# refAa="D"
# mutAa="V"
# aaLen=29
# getSeq(accID,codingPos,refDna,mutDna,aaPos,refAa,mutAa,aaLen)

# choose the corresponding annotation version of genome reference. 
grch37 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",dataset="hsapiens_gene_ensembl")
#grch38 <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
MART=grch37

# return a matrix with cdna and peptide sequences. Input is annotation column in VCF and length of peptide you want to have.
## callSeq("HTT:NM_002111:exon48:c.A6509T:p.E2170V,SYNE1:NM_033071:exon138:c.G24833A:p.G8278E,",29)
## callSeq("ZMYM4:ENST00000373297.2:exon26:c.G3865A:p.V1289I,ZMYM4:ENST00000314607.6:exon28:c.G4132A:p.V1378I,",29)
## callSeq("MAMDC2:ENST00000377182.4:exon6:c.658_666del:p.220_222del,",29)
## callSeq("TP53:ENST00000359597.4:exon4:c.454dupC:p.P152fs,",29)
## callSeq("TP53:ENST00000455263.2:exon5:c.518_519del:p.V173fs,",29)
## callSeq("PRR12:ENST00000418929.2:exon6:c.4399_4400insCCTCAG:p.P1467delinsPSA,",29)
## callSeq("NKX2-3:ENST00000344586.7:exon2:c.892_894del:p.298_298del,",29)
## outlier, "NF1:NM_001042492:exon37:c.5108_5110CC,NF1:NM_000267:exon36:c.5045_5047CC," chr17	29653110	.	GCA	CC	7933.16	PASS	AF=0.660221;AO=711;DP=1098;FAO=717;FDP=1086;FDVR=10;FR=.;FRO=369;FSAF=325;FSAR=392;FSRF=142;FSRR=227;FWDB=-0.00885586;FXX=0.0118288;HRUN=1;HS_ONLY=0;LEN=3;MLLD=479.017;OALT=CC;OID=.;OMAPALT=CC;OPOS=29653110;OREF=GCA;PB=0.5;PBP=1;QD=29.2197;RBI=0.00893356;REFB=-0.00559692;REVB=0.00117572;RO=359;SAF=322;SAR=389;SRF=136;SRR=223;SSEN=0;SSEP=0;SSSB=0.0471207;STB=0.523567;STBP=0.032;TYPE=complex;VARB=0.00272484	GT:GQ:DP:FDP:RO:FRO:AO:FAO:AF:SAR:SAF:SRF:SRR:FSAR:FSAF:FSRF:FSRR	0/1:2667:1098:1086:359:369:711:717:0.660221:389:322:136:223:392:325:142:227

callSeq <- function(annotationColumn,aaLen){
  transInfo = strsplit(as.character(annotationColumn),",")[[1]]
  combPep=NULL
  for (item in transInfo){
    transUnit = strsplit(item,":")[[1]]
    geneName <- transUnit[1]
    if(startsWith(transUnit[2],"ENST")){
      accID=substr(transUnit[2],1,regexpr(".",transUnit[2],fixed=T)[1]-1)
    }else{
      accID=transUnit[2]
    }  
    cInfo=transUnit[4] # 4th column
    subInfo=strsplit(as.character(cInfo),"c.")[[1]][2]
    codingPos=str_split(gsub("[A-Z]|[a-z]", "", subInfo),"_")[[1]] # if "_" is there
    # duplication
    if (grepl("dup",subInfo)){
      dupInfo <- str_split(subInfo,"dup")[[1]]
      if (dupInfo[2]==""){
        print("Error: This is duplication for insertion. Pleas check them manually")
      }else{
        insEle <- dupInfo[2]
        type <- "duplication"
      }
    }else if (grepl("ins",subInfo)){
      insInfo <- str_split(subInfo,"ins")[[1]]
      if (insInfo[2]==""){
        print("Error: This is empty insertion. Pleas check them manually")
      }else{
        insEle <- insInfo[2]
        type <- "insertion"
      }
    }else if (grepl("del",subInfo)){
      insEle <- ""
      type <- "deletion"
    }else {
		print(paste("Error: ", item, " This is outlier. Please check them manually"))
		break
	}    
    ## mutate aa information 
    pInfo=strsplit(as.character(transUnit[5]),"p.")[[1]][2]
    aaPos <- as.numeric(str_split(gsub("[A-Z]|[a-z]", "", pInfo),"_")[[1]] [1])
    if (grepl("fs",pInfo)){
      frameshift=TRUE
    }else{
      frameshift=FALSE
    }
    ## getting mutated peptides
	eachPep<-getSeq(accID,codingPos,insEle,type,aaPos,frameshift,aaLen)
	# if peptides could be predicted
	if ( length(eachPep) > 1 ){
		## add information
		eachPep <- cbind(geneName, aaChange=transUnit[5], transID=transUnit[2], eachPep)
		## combinning mutated peptides
		combPep<-rbind(combPep,eachPep)
	}else{
		print(paste(accID, " cannot be found in database"))
	}
  }
  return(unique(combPep))
}

# get the sequence from biomaRt via web server.
getSeq <- function(accID,codingPos,insEle,type,aaPos,frameshift,aaLen){
	insEle <- as.character(insEle)
	#codingPos=as.integer(codingPos)
	#aaPos=as.integer(aaPos)
	#aaLen=as.integer(aaLen)
	#Refseq annotation, "HTT:NM_002111:exon48:c.A6509T:p.E2170V,SYNE1:NM_033071:exon138:c.G24833A:p.G8278E,"
	if(!startsWith(accID,"ENST")){
		try((bmCoding <- getBM(attributes = c("coding","hgnc_symbol","ccds"),filters = c("refseq_mrna"), values = accID, mart=MART)),TRUE)
		try((bmPeptide <- getBM(attributes = c("peptide","hgnc_symbol","ccds"),filters = c("refseq_mrna"), values = accID, mart=MART)),TRUE)
	}else{
		#Ensembl annotation, "MAMDC2:ENST00000377182.4:exon6:c.658_666del:p.220_222del,"
		try((bmCoding <- getBM(attributes = c("coding","hgnc_symbol"),filters = c("ensembl_transcript_id"), values = accID, mart=MART)),TRUE)
		try((bmPeptide <- getBM(attributes = c("peptide","hgnc_symbol"),filters = c("ensembl_transcript_id"), values = accID, mart=MART)),TRUE)
	}  
	sideLen=as.integer(aaLen/2)
	# for coding DNA sequence
	if(exists("bmCoding") && length(bmCoding)>0 && nrow(bmCoding)){     
		#if (refDna==substr(bmCoding$coding,codingPos,codingPos)){
		#  seqRef=bmCoding$coding
		#}else {stop("Error in reference genome");}
		seqRef=bmCoding$coding
		proteinRef=bmPeptide$peptide
		codingPos1 <- as.integer(codingPos[1])    
		if (type=="deletion"){
		  if (is.na(codingPos[2])){
			seqMut <- paste(substr(seqRef,1,codingPos1-1),substr(seqRef,codingPos1+1,nchar(seqRef)),sep="")
		  }else{
			codingPos2 <- as.integer(codingPos[2])
			seqMut <- paste(substr(seqRef,1,codingPos1-1),substr(seqRef,codingPos2+1,nchar(seqRef)),sep="")
		  }
		}
		if (type=="insertion"){
		  if (is.na(codingPos[2])){
			seqMut <- paste(substr(seqRef,1,codingPos1), insEle ,substr(seqRef,codingPos1+1,nchar(seqRef)),sep="")
		  }else{
			codingPos2 <- as.integer(codingPos[2])
			seqMut <- paste(substr(seqRef,1,codingPos1), insEle ,substr(seqRef,codingPos2,nchar(seqRef)),sep="")
		  }
		}
		if (type=="duplication"){
		  if (is.na(codingPos[2])){
			 seqMut <- paste(substr(seqRef,1,codingPos1), insEle ,substr(seqRef,codingPos1+1,nchar(seqRef)),sep="")
			type="insertion" # reclassification
		  }else{
			#codingPos2 <- as.integer(codingPos[2])
			#seqMut <- paste(substr(seqRef,1,codingPos1), insEle ,substr(seqRef,codingPos2,nchar(seqRef)),sep="")
			print("Error: this script doen't not work for segment duplication for insertion, please manually check; return protein reference")
			seqMut=proteinRef
		  }
		}
		# predict peptide based on coding sequence
		# predictedMut<-c2s(translate(DNAString(seqMut),if.fuzzy.codon="error"))
		predictedMut<-c2s(translate(s2c(seqMut)))
		#print(predictedMut)
		# output=cbind(accID,codingPos,refDna,mutDna,aaPos,refAa,mutAa,aaLen,codingRefSeq,codingMutSeq,pepRefSeq,pepMutSeq,predictedRef,predictedMut)
	
		if (frameshift==TRUE){
			if (gregexpr(pattern="\\*",predictedMut)[[1]][1] != -1){
				stopPos <- gregexpr(pattern="\\*",predictedMut)[[1]][1] # exclude asterisk symbol
			}else{
			stopPos <- nchar(predictedMut)  # if no stop sign is found. Search 3' UTR
			# ENST00000508793.1, TP53, p.V157fs, no stop code in cdna found.
			# TP53:ENST00000508793.1:exon5:c.470_473del:p.V157fs
			# callSeq("TP53:ENST00000508793.1:exon5:c.470_473del:p.V157fs,",29)
			type <- paste(type," , stopcode in 3'UTR",sep="")
			# downstream 500bp
			get3utr <- biomaRt::getSequence(id= accID, type="ensembl_transcript_id", seqType="coding_transcript_flank", downstream = 500,mart =  MART)
			seqMut <- paste(seqMut, get3utr$coding_transcript_flank, seq="")
			predictedMut <- c2s(translate(s2c(seqMut)))
			stopPos <- gregexpr(pattern="\\*",predictedMut)[[1]][1] # find the stop code in 3utr
			}
			mutationType <- paste("frameshift", type)
		}else{
			stopPos <- min(aaPos + sideLen,nchar(predictedMut)) # in case right side of the peptide is not long enough
			mutationType <- paste("nonFrameshift", type)
		}

		startPos <- aaPos - sideLen
		#print(seqMut)
		targetPeptide <- substr(predictedMut,startPos,stopPos)
		refPeptide <- substr(proteinRef,startPos,stopPos)
		#output=cbind(targetPeptide, refPeptide, predictedMut, seqMut, bmCoding$coding)
		seqRef<-bmCoding$coding
		output=cbind(targetPeptide, refPeptide, mutationType, predictedMut, proteinRef,seqMut,seqRef)
		#print(output)
		return(output)
	}else{
		return(NULL)
	}
}

results <- callSeq(annotationColumn,aaLen)

write.table(results,file=outFile,quote=F,sep="\t",row.names = F)

