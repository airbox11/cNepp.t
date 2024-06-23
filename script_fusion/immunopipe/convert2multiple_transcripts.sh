#! /bin/bash

# if there exist multiple transcripts for a gene, this script splits the single entry (per gene)
# into multiple entries (per transcript)

variant_function=$1

awk 'BEGIN{FS="\t";OFS="\t"} {split($3,a,","); for(i=1;i<length(a);i++) {$3=a[i];print $0}}' $variant_function


