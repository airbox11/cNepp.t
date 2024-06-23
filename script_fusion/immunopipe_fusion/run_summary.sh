#!/bin/bash



# tell bash to be verbose and to abort on error
set -e -u -o pipefail



# a grep version that only exits on error when there really is an error
# and not just when the search string is not found
function grep0 { grep "$@" || [ $? -lt 2 ]; }



cd "$OUTPUT_DIR"



# reformat netMHCpan output to tab-separated file
for MHC_CLASS in MHCI MHCII; do
	cat "net${MHC_CLASS}_"*"_fusion" | # concatenate results from all HLA alleles
	sed -e 's/<=//' | # remove arrows from table, since they mess up column numbering
	grep0 -v -P '^-+$' | # skip table borders (consisting of lots of dashes)
	grep0 -v 'Number of' | # skip summary lines
	grep0 -v -P '^#' | # skip comment lines
	grep0 -v -P 'Wrong' | # skip lines with 'Wrong input' 
	sed -e 's/^  *//' -e 's/  */\t/g' | # use tab as column separator instead of variable number of blanks
	awk '/^$/{content=0} $3=="Peptide"{content=1} content==1{print}' | # content is between the header "Peptide" and an empty line
	awk -F '\t' -v OFS='\t' '
		NR==1 { cols=NF } # get number of columns from header
		{ for (i=1; i<=cols; i++) if ($i=="") $i="."; print } # fill empty columns with dot
	' |
	awk '!duplicate[$0]++' > "${MHC_CLASS}_reformat" # skip duplicate lines, such as multiple headers
done



# We want to remove peptides that are wild-type.
# Arriba marks all non-references amino acids as lowercase characters, including
# frameshifts, intron retention, non-template bases, SNPs, SNVs and indels.
# So to identify wild-type peptides, all we need to do is check if a peptide
# appears in Arriba's output all uppercase.
# Since we are interested in neoantigens arising from fusions, we should ignore
# non-reference amino acids arising from SNPs, SNVs and indels, however.
# So we treat the latter as normal amino acids by converting them to uppercase.

# extract column containing fusion peptides
awk -F '\t' '
	NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i } # get column names
	NR >  1 { print $col["peptide_sequence"] }
' "$FUSION" |

# convert SNPs and SNVs to uppercase
# assume a single lowercase letter flanked by uppercase letters to be a SNP/SNV
sed -e 's/[A-Z][a-z][A-Z]/\U&/g' -e 's/[A-Z][a-z][A-Z]\|^[a-z][A-Z]\|[A-Z][a-z]$/\U&/g' |

# extract wild-type peptides (i.e., uppercase sequences)
grep0 -o -P '[A-Z]+' > wild-type_peptides

# remove wild-type and non-binding peptides from netMHC output files
for MHC_CLASS in MHCI MHCII; do
	awk -F '\t' '

		# load all wild-type peptides into memory
		FILENAME == "wild-type_peptides" { wild_type[i++] = $0 }

		# extract column names
		FILENAME ~ /reformat/ && FNR == 1 { for (i = 1; i <= NF; i++) col[$i] = i }

		# skip all peptides from netMHCpan that have a wild-type sequence
		FILENAME ~ /reformat/ {
			for (i in wild_type)
				if (match(wild_type[i], $col["Peptide"]))
					next # peptide has a wild-type sequence => skip it
			# if ("BindingLevel" in col && $col["BindingLevel"] != "." || # netMHCIIpan
			#     "BindLevel"    in col && $col["BindLevel"]    != ".")   # netMHCpan
				print
		}
	' wild-type_peptides "${MHC_CLASS}_reformat" > "${MHC_CLASS}_filtered"
done



# join full fusion annotation from Arriba's output file with netMHC output
for MHC_CLASS in MHCI MHCII; do
	cat "$FUSION" |
	awk -F '\t' '
		FNR == 1 { for (i = 1; i <= NF; i++) col[$i] = i } # extract column names
		FILENAME == "/dev/stdin" { fusion["fusion"(FNR-1)] = $0 } # load fusions into memory
		FILENAME != "/dev/stdin" && FNR == 1 { print $0 "\t" fusion["fusion0"] } # print header
		FILENAME != "/dev/stdin" && FNR >  1 { print $0 "\t" fusion[$col["Identity"]] } # print content
	' /dev/stdin "${MHC_CLASS}_filtered" > "results_${PID}_${MHC_CLASS}_epitopes_filtered_ident.tsv"
done



# cleanup
# rm -f  wild-type_peptides *_fusion *_reformat *_filtered