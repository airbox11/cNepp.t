import sys
import re

file1 = sys.argv[1]
file_output = f'{file1}_addGeneName.vcf'
print(file1)

pattern = '\|([^\|]+)\|ENSG'

# pattern = 'ENSG'
with open(file1, 'r') as f1, open(file_output, 'w') as w1:
	c = 0
	for line in f1:
		# line = '|name|ENSG|name|ENSG|nameS|ENSG'
		line = line.strip()
		# print(line)
		p1 = re.search(pattern, line)
		if p1:
			geneName = p1.group(1)
			# print(geneName)
		else:
			if c== 0:
				geneName = ''
			else: 
				geneName = 'None'
				print('no match.........')

		w1.write(f'{line}\t\t\t\t\t\t\t\t\t{geneName}\n')

		c += 1
		# if c > 3:
		# 	quit()
