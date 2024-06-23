#! /usr/bin/python

from Bio import SeqIO
import sys

# input fasta file generated from ANNOVAR after renaming title of each protein
#fastaFile="/icgc/dkfzlsdf/analysis/D120/huangz/test_annotation/"
fastaFile=sys.argv[1]
outputMut=sys.argv[2]
outputRef=sys.argv[3]
# peptide length away from mutation point
extendLen=sys.argv[4]
extendLen=int(extendLen)


# output peptide fasta files used for HLA binding affinity prediction
pepWildFasta=open(outputRef,'w+r')
pepMutFasta=open(outputMut,'w+r')
# read the input fasta protein file
handle=open(fastaFile,"rU")
records = list(SeqIO.parse(handle, "fasta"))
handle.close()

for i in range(0,len(records),2):
    try:
        (aaWild,mutPos,aaMut)=records[i+1].id.split(":")
        mutPos=int(mutPos)
    except ValueError:
        print(records[i+1].id + " is not somatic mutation and ignored")
        continue
    wildLen=len(records[i].seq)
    mutLen=len(records[i+1].seq)
    # determin the wild type peptide spanning the mutation position
    if mutPos>extendLen:
        if mutPos+extendLen<wildLen:
            wildPep=records[i].seq[mutPos-extendLen-1:mutPos+extendLen]
        else:
            wildPep=records[i].seq[mutPos-extendLen-1:wildLen]
    else:
        if mutPos+extendLen<wildLen:
            wildPep=records[i].seq[0:mutPos+extendLen]
        else:
            wildPep=records[i].seq[0:wildLen]
    # determin the mutated peptide spanning the mutation position
    if mutPos>extendLen:
        mutPep=records[i+1].seq[mutPos-extendLen-1:mutLen]
    else:
        mutPep=records[i+1].seq[0:mutLen]
    ### check that mutPep does not already exists in new fasta so it is not a duplicate  
    #put cursor for reading back at the start of the file
    pepMutFasta.seek(0)
    if str(mutPep) in pepMutFasta.read():
        next
    else:
        pepWildSeq=str(">" + str(i/2+1) + records[i].id.replace("line","l") + "\n" + wildPep + "\n")
        pepMutSeq=str(">" + str(i/2+1) + records[i+1].id.replace("line","l") + "\n" + mutPep + "\n")
        #print(pepWildSeq + pepMutSeq)
        pepWildFasta.write(pepWildSeq)
        pepMutFasta.write(pepMutSeq)

pepWildFasta.close()
pepMutFasta.close()