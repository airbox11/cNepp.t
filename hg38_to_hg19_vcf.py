import sys
from liftover import ChainFile
import pdb

file_input = sys.argv[1]
file_output_19 = sys.argv[2]
file_output_38 = sys.argv[3]


chain = '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/hg38ToHg19.over.chain.gz'
converter = ChainFile(chain, one_based=True)
# chr = '6'
# pos = 160585140

def hg38_to_hg19 (chr, pos): 
	chrpos_19 = converter[chr][pos]
	if len(chrpos_19) > 0:
		chr_new = chrpos_19[0][0].split('chr')[1]
		pos_new = chrpos_19[0][1]
	else:
		chr_new, pos_new = '0','0'
	return(chr_new, pos_new)

with open(file_input, 'r') as r1, open(file_output_19, 'w') as w1, open(file_output_38, 'w') as w2:
	for index, line in enumerate(r1):
		if index == 0:
			line = line.strip()
			w1.write(f'{line}\tchr_38\tpos_38\n')
			w2.write(f'{line}\n')
			continue
		cols = line.strip().split('\t')
		cols_extra_string = '\t'.join(cols[2:])
		chr_38, pos_38 = cols[0], int(cols[1])
		chr_19, pos_19 = hg38_to_hg19(chr_38, pos_38)
		if chr_19 != '0':
			w1.write(f'{chr_19}\t{pos_19}\t{cols_extra_string}\t{chr_38}\t{pos_38}\n')
			w2.write(f'{chr_38}\t{pos_38}\t{cols_extra_string}\n')
