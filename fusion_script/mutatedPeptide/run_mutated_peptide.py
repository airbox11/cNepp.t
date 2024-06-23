#! /user/bin/python
from Bio import SeqIO
import sys
import argparse
## To generate mutated peptide based on ANNOVAR annotation result files
## {prefix}.annovar.exonic_variant_function_multiple

# input fasta file generated from ANNOVAR after renaming title of each protein
#annoFile=sys.argv[1] # snvs_EH63_somatic_functional_snvs_conf_7_to_10.vcf_nonsynonymous.annovar.exonic_variant_function_multiple
#mutFile=sys.argv[2] # EH63_mut.fa
#refFile=sys.argv[3] # EH63_ref.fa
#outFile=sys.argv[4] # EH63_predicted_mutational_peptide.fa

parser = argparse.ArgumentParser()
parser.add_argument('-a', required=True, metavar="annoFile", help="ANNOVAR annotation file")
parser.add_argument('-m', required=True, metavar="fasta", help="mutation fasta file")
parser.add_argument('-r', required=True, metavar="fastq", help="reference fasta file")
parser.add_argument('-o', required=True, metavar="output", help="path to output file")

args = parser.parse_args()

annoFile=args.a 
mutFile=args.m
refFile=args.r
outFile=args.o


# open output file
output=open(outFile,'w')

myD={}
# a key of NM accession is not enough
# it should include mutated position
for line in open(annoFile,'r'):
	arr1=line.split("\t")
	arr2=arr1[2].split(":")
	# include the line poistion, indicating different mutated positions
	key=arr2[1]+arr1[0].split("line")[1]
	value=arr2[0]+":"+arr2[1]+","+arr1[3]+":"+arr1[4]
	myD[key]=value


refHandle=open(refFile,"rU")
refRecords = list(SeqIO.parse(refHandle, "fasta"))
refHandle.close()



mutHandle=open(mutFile,"rU")
mutRecords = list(SeqIO.parse(mutHandle, "fasta"))
mutHandle.close()

for i in range(0,len(mutRecords)):
	arr=mutRecords[i].id.split(",")
	# the same key generated in annoFile
	key=arr[1]+arr[0].split("l")[1]
	substitution=arr[4].replace(":","")
	value=myD[key].split(",")	
	outputLine= str("chr" + value[1] +"\t"+ value[0] +"\t" + substitution +"\t"+ mutRecords[i].seq +"\t"+ refRecords[i].seq + "\n")
	output.write(outputLine)

output.close()



