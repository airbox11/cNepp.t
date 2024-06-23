#!/bin/bash

AMINO_ACIDS=100 # how many amino to select before and after fusion junction

cat $1 |

# extract fusions with predicted peptides
awk -F '\t' '
	NR == 1 { print $0"\tshort_peptide" } # print header line

	                    # select fusions according to these criteria:
	$6 ~ /^12:/ &&      # first breakpoint on chr12
	$7 ~ /^12:/ &&      # second breakpoint on chr12
	$18 ~ /high/ &&     # high confidence prediction
	$24 != "." &&       # has predicted peptide sequence
	$24 !~ /\*.*\|/ &&  # predicted peptide sequence does not have a stop codon before fusion junction
	$24 !~ /\|\*/ {     # predicted peptide sequence does not have a stop codon right after the fusion junction

		# remove special characters from peptide: [, ], -, and *
		short_peptide = $24
		gsub(/\[|\]|-|\*/, "", short_peptide)
		short_peptide = toupper(short_peptide)

		# select 20 amino acids before and after the junction
		match(short_peptide, /.{1,'$AMINO_ACIDS'}\|.{1,'$AMINO_ACIDS'}/)
		short_peptide = substr(short_peptide, RSTART, RLENGTH)

		# skip peptides with unclear amino acids (marked by "?")
		if (short_peptide !~ /\?/)
			print $0"\t"short_peptide

	}
' |

# count fusion recurrence
awk -F '\t' '
NR == 1 { print $0"\trecurrence" } # print header line
NR > 1{
	# extract PID from path to BAM file
	pid = $1
	sub(".*/view-by-pid/", "", pid)
	sub("/.*", "", pid)

	# count fusions with identical breakpoints in different patients
	if (!duplicate[pid"\t"$6"\t"$7]++) {
		counts[$6"\t"$7]++
		lines[++line_count] = $0
		breakpoints[line_count] = $6"\t"$7
	}
}
END {
	for (line in lines)
		print lines[line]"\t"counts[breakpoints[line]]
}' |
sort -k27,27n -k 6,7

