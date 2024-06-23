#! /bin/bash

## To format the PHLAT prediction results for downstream analysis

usage() { echo "Usage: $0 -p phlatFile -o outputFile
    -p: Path to PHLAT hla result file
    -o: Path to output file " 1>&2; exit 1; }

while getopts ":p:o:" opt; do
    case "${opt}" in
        p)
            phlat=${OPTARG}
            ;;
        o)
            OUTPUT_file=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${phlat}" ] || [ -z "${OUTPUT_file}" ] ; then
    usage
fi

cat $phlat | awk '{if (NR>1 && NR<=4) 
                {if (index($2,"*")!=0) 
                    {gsub("*","",$2); 
                    split($2,a,":"); 
                    if(a[2] !~ /[^0-9]/) print "HLA-" a[1] ":" a[2];else print $2}; 
                if (index($3,"*")!=0) 
                {gsub("*","",$3); 
                split($3,b,":"); 
                if(b[2] !~ /[^0-9]/) print "HLA-" b[1] ":" b[2];else print $3 }}; 
            if (NR>4) 
                {if (index($2,"*")!=0) 
                    {split($2,a,":"); 
                    split(a[1],aa,"*"); 
                    if(a[2] !~ /[^0-9]/) print aa[1] "_" aa[2] a[2]} ;  
                if (index($3,"*")!=0) 
                    {split($3,b,":"); 
                    split(b[1],bb,"*"); 
                    if(b[2] !~ /[^0-9]/) print bb[1] "_" bb[2] b[2]}}}' | sort > ${OUTPUT_file}
