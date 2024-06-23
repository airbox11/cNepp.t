
import sys
import pdb

workDir = sys.argv[1]

file_input        = f'{workDir}/2_SNVs_based_neoepitope_prediction/2.vcf'
file_output_snv   = f'{workDir}/2_SNVs_based_neoepitope_prediction/snv_somatic.vcf'
file_output_indel = f'{workDir}/4_indel_based_prediction/indel_tmp.vcf'

with open(file_input, 'r') as r1, open(file_output_snv, 'w') as w1, open(file_output_indel, 'w') as w2:
	count = 0
	while True:
		count += 1
		line = r1.readline()
		if not line:
			break
		if count == 1:
			w1.writelines(line)
			w2.writelines(line)
			continue

		ref = line.split('\t')[3]
		alt = line.split('\t')[4]
		if len(ref) == 1 and len(alt) == 1:
			w1.writelines(line)
		else:
			w2.writelines(line)



