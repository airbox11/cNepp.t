#! /bin/bash

variant_function=$1

awk 'BEGIN{FS="\t";OFS="\t"} {split($3,a,","); for(i=1;i<length(a);i++) {$3=a[i];print $0}}' $variant_function


