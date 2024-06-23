#!/bin/env python
import blosum as bl
import os.path
import sys
import pdb


## 
def iterate_character(str1):
	for c in str1:
		yield c

def get_couples (a,b):
	return (f'{a}{b}')


def csv_reader(file_input):
	for ref in open(file_input, "r"):
		yield ref.strip('\n')

def get_max_score (qur, file_input, weight):
	score_max = -1000
	ref_max = ''
	set1 = set()
	count = 0
	for ref in csv_reader(file_input):
		count += 1
		couples = list(map(get_couples, iterate_character(qur), iterate_character(ref) ))
		scores = [ mat[i] for i in couples]
		score = sum(list(map(lambda x,y: x*y, scores, weight)))
		if score > score_max:
			score_max = score
			set1 = set()
		if score == score_max:
			ref_score = f'{ref}\t{score}'
			set1.add(ref_score)
		# if count > 10000:
		# 	break

	return (set1)


##
def genereate_weight (len1, odd, mag = 4):
	a = [*range(1,len1//2*mag,mag)]
	b = [*range(1,len1//3+1)]
	a[:len(b)] = b
	if odd:
		c = [*a[:-1], *a[::-1]]
	else:
		c = [*a,*a[::-1]]
	return(c)

def get_weight (len1):
	if len1%2 == 0:
		weight = genereate_weight(len1, False)
	elif len1%2 == 1:
		len2 = len1 + 1
		weight = genereate_weight(len2, True)

	return(weight)

## 
def select_refs_input (db_type, len1):
	file_input = f'/omics/groups/OE0422/internal/yanhong/git/hex/virus_bac_database/splice_py_all_virus_protein/sequences_human_virus_freseq_all.fasta_seq_matrix_{len1}_unique'
	file_input = f'/omics/groups/OE0422/internal/yanhong/git/hex/virus_bac_database/splice_py_all_virus_protein/all_virus_protein_sequences.fasta_seq_matrix_{len1}_unique_lite_1per'
	file_input = f'/omics/groups/OE0422/internal/yanhong/git/hex/virus_bac_database/splice_py_all_virus_protein/all_virus_protein_sequences.fasta_seq_matrix_{len1}_unique_lite'

	database_virus = f'/omics/groups/OE0422/internal/yanhong/git/hex/virus_bac_database/splice_py_all_virus_protein/all_virus_protein_sequences.fasta_seq_matrix_{len1}_unique_lite_1per'
	database_virus = f'/omics/groups/OE0422/internal/yanhong/git/hex/virus_bac_database/splice_py_all_virus_protein/all_virus_protein_sequences.fasta_seq_matrix_{len1}_unique_lite'
	database_bacteria = f'/omics/groups/OE0422/internal/yanhong/git/hex/virus_bac_database/splice_py_uniProt/all_protien_bacteria.fasta_seq_shortSeq{len1}'
	database_fungi = f'/omics/groups/OE0422/internal/yanhong/git/hex/virus_bac_database/splice_py_uniProt/all_protien_fungi.fasta_seq_shortSeq{len1}'

	if db_type == 'virus':
		return (database_virus)
	elif db_type == 'bacteria':
		return (database_bacteria)
	elif db_type == 'fungi':
		return (database_fungi)

## public variable
qur = sys.argv[1]
file_output = sys.argv[2]
db_type=sys.argv[3]

len1 = len(qur)
mat = bl.BLOSUM(62)

## main:
file_input = select_refs_input(db_type, len1)

weight = get_weight(len1)
ref_score = get_max_score(qur, file_input, weight)

with open(file_output, 'a') as w1:
	for i in ref_score:
		result = f'{qur}\t{i}\t{db_type}\n'
		w1.write(result)

print(f'==== ==== {qur} is done for {db_type} ]]]] ]]]]')

