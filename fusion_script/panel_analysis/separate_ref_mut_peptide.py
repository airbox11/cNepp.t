#! /usr/bin/python

from Bio import SeqIO
import sys

# input fasta file generated from ANNOVAR after renaming title of each protein
#fastaFile="/icgc/dkfzlsdf/analysis/D120/huangz/test_annotation/"
fastaFile=sys.argv[1]
outputMut=sys.argv[2]
outputRef=sys.argv[3]
# peptide length distance from mutation point position
extendLen=sys.argv[4]
extendLen=int(extendLen)


# output peptide fasta files used for HLA binding affinity prediction
pepWildFasta=open(outputRef,'w')
pepMutFasta=open(outputMut,'w')
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
		if mutPos+extendLen<mutLen:
			mutPep=records[i+1].seq[mutPos-extendLen-1:mutPos+extendLen]
		else:
			mutPep=records[i+1].seq[mutPos-extendLen-1:mutLen]
	else:
		if mutPos+extendLen<wildLen:
			mutPep=records[i+1].seq[0:mutPos+extendLen]
		else:
			mutPep=records[i+1].seq[0:mutLen]
	pepWildSeq=str(">" + str(i/2+1) + records[i].id.replace("line","l") + "\n" + wildPep + "\n")
	pepMutSeq=str(">" + str(i/2+1) + records[i+1].id.replace("line","l") + "\n" + mutPep + "\n")
	

	# put the previous record into buffer
	# filter duplicate mutated peptides, it cannot be sure to remove all duplicates
	if i==0:
		buff_wildPep=wildPep
		buff_mutPep=mutPep
		buff_pepWildSeq=pepWildSeq
		buff_pepMutSeq=pepMutSeq
	else:
		try:
			buff_mutPep
		except NameError:
			buff_wildPep=wildPep
	                buff_mutPep=mutPep
        	        buff_pepWildSeq=pepWildSeq
                	buff_pepMutSeq=pepMutSeq
			continue
		if str(mutPep)!=str(buff_mutPep):
			pepWildFasta.write(buff_pepWildSeq)
			pepMutFasta.write(buff_pepMutSeq)
			buff_wildPep=wildPep
			buff_mutPep=mutPep
	                buff_pepWildSeq=pepWildSeq
	                buff_pepMutSeq=pepMutSeq
		else:
			print("They are the same:" + mutPep)
	# output the last record
	if i==len(records)-2:
		if mutPep==buff_mutPep:
			pepWildFasta.write(buff_pepWildSeq)
			pepMutFasta.write(buff_pepMutSeq)
		else:
			pepWildFasta.write(pepWildSeq)
			pepMutFasta.write(pepMutSeq)


pepWildFasta.close()
pepMutFasta.close()

"""
record_dict=SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/ibios/tbi_cluster/13.1/x86_64/python/python/lib/python2.7/site-packages/Bio/SeqIO/__init__.py", line 673, in to_dict
    raise ValueError("Duplicate key '%s'" % key)
ValueError: Duplicate key 'line1'
"""
