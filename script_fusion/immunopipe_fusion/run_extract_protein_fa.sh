#!/bin/bash

# tell bash to be verbose and to abort on error
set -x -e -u -o pipefail

# extract fusion peptides from Arriba's output file and convert them to FastA
awk -F '\t' -v FUSION_CONFIDENCE="${FUSION_CONFIDENCE}" '

	# get column names
	NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i }

	NR > 1 &&                                      # skip header
	length($col["peptide_sequence"]) > 1 &&        # skip empty peptides
	match(FUSION_CONFIDENCE, $col["confidence"]) { # skip low-confidence fusions

		peptide = $col["peptide_sequence"]

		# ignore sequences before stop codon (*) in the 5p part; they are probably not translated
		if (index(peptide, "*") < index(peptide, "|"))
			sub(/^[^|]*\*/, "", peptide)

		# ignore lowercase sequences in the 5p part; they are probably UTRs
		if (peptide ~ /^[^|]*[A-Z]/) # are there uppercase sequences (=coding regions) in the 5p part?
			sub(/^[a-z?]*/, "", peptide) # remove lowercase sequences (=UTRs)

		# remove special characters (indicating gaps, insertions, etc.)
		peptide = toupper(peptide)
		gsub(/[^A-Z?]/, "", peptide)

		# Arriba puts question marks when it is unsure about the amino acid => this confuses netMHCpan
		# => remove them by splitting the sequence into the chunks between the question marks
		split(peptide, peptide2, "?")
		for (i in peptide2) {
			print ">fusion"(NR-1)
			print peptide2[i]
		}
	}

' "$FUSION" > "$OUT"

