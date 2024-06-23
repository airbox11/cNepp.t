#! /usr/bin/python
import os
import sys
from Bio import SeqIO


mutMHCfile=sys.argv[1]
refMHCfile=sys.argv[2]
nonsynonymousFile=sys.argv[3]
mutPeptide=sys.argv[4]
refPeptide=sys.argv[5]
rawResultsFile=sys.argv[6]

#
output=open(rawResultsFile,'w')

# make a dictionary for vcf file
# remove colomn name line because the number of row is part of uniq ID for mutation
# The number of line is used for matching ANNOVAR annotation results
# remove all line start with #
os.system("bash -c 'grep -o ^[^#]* %s > %s_rmHeader'" % (nonsynonymousFile, nonsynonymousFile))
nonsynonymousFile=nonsynonymousFile + "_rmHeader"


# make a dictionary for vcf file
# recall corresponding information in vcf file
vcfDict={}
with open(nonsynonymousFile) as fp2:
	for k, line in enumerate(fp2):
		key=str(k+1)
		vcfDict[key]=line.strip()

os.system("bash -c 'rm %s'" % nonsynonymousFile)

# make a dictionary for mutated 
# recall mutated peptide information in fasta files
handleMut=open(mutPeptide,"rU")
recordsMut = list(SeqIO.parse(handleMut, "fasta"))
handleMut.close()

mutDict={}
for i in range(0,len(recordsMut)):
	arr=recordsMut[i].id.split(",")
	aaChange="p."+arr[4].replace(":","")
	key=arr[0]
	mutDict[key]=arr[1] + "\t" + arr[2] + "\t" + aaChange + "\t" + arr[3] + "\t" + recordsMut[i].seq

# dictionary for reference peptide
# recall reference peptide information in fasta files
handleRef=open(refPeptide,"rU")
recordsRef = list(SeqIO.parse(handleRef, "fasta"))
handleRef.close()

refDict={}
for i in range(0,len(recordsRef)):
        arr=recordsRef[i].id.split(",")
        #lineID=arr[0].split("l")[1]
        #key="line" + str(lineID)
	key=arr[0]
        refDict[key]=recordsRef[i].seq

# recall reference peptide
# dictionary for reference MHC neoepitopes
mhcRefDict={}
for line in open(refMHCfile,'r'):
        line=line.strip()
        arr=line.split("\t")
	# excluding non-relevant lines
        if line[:1]=="#" or len(arr)<2 or "Peptide" in line or "binders:" in line or "Distance" in line:
                continue
        try:
                key=arr[3].split(",")[0] + "_" + arr[0]
                mhcRefDict[key]=line
        except:
                print "Error in line: " + line
                continue

# make a  dictionary for mutated MHC neoepitopes
mhcMutDict={}
for line in open(mutMHCfile,'r'):
        line=line.strip()
        try:
                bindingLevel=line.split("\t")[11]
        except:
                continue
        if bindingLevel=="WB" or bindingLevel=="SB":
                arr=line.split("\t")  # extract line positon for identifying line in file in later steps, max number of input predicted peptides per fasta file is 99999
                key=arr[3].split(",")[0] + "_" + arr[0]
                mhcMutDict[key]=line
                key1=arr[3].split(",")[0]
                lineID=key1.split("l")[1]
                if len(mhcRefDict[key].split("\t"))!=12:
                        data=mhcMutDict[key] + "\t" + mhcRefDict[key] + "\tnonB" + "\t" + mutDict[key1] + "\t" + refDict[key1] + "\t" + vcfDict[lineID] +"\n"
                else:
                        data=mhcMutDict[key] + "\t" + mhcRefDict[key] + "\t" + mutDict[key1] + "\t" + refDict[key1] + "\t" + vcfDict[lineID] +"\n"
                print key1
                output.write(str(data))

output.close()

