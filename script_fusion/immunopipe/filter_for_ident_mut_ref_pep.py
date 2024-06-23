#! /usr/bin/python

import pandas as pd
import sys

rawResultsFile = sys.argv[1]
outFile = sys.argv[2]

#open netMHC output
output_unfiltered=open(rawResultsFile,'r').readlines()

#filter out neoepitopes that are identical to the reference peptide generated for that position 
#(peptide does not include mutated region)
p2 = [ line.split("\t") for line in output_unfiltered ]
df = pd.DataFrame(p2)
df.columns = df.iloc[0]
df = df.iloc[1:]
d={'Peptide': ['Peptide1', 'Peptide2']}
df = df.rename(columns=lambda c: d[c].pop(0) if c in d.keys() else c)
df1 = df.drop_duplicates(subset=['Peptide1'])
df2 = df1.query("Peptide1 != Peptide2")

out_csv = df2.to_csv(sep='\t')
with open(outFile, "w") as text_file:
    text_file.write(out_csv)