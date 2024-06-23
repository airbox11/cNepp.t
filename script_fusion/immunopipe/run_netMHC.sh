#! /bin/bash

# Run netMHC for HLA and DRB 

cat $phlatResultFile | awk '{if (NR>1 && NR<=4) 
				{if (index($2,"*")!=0) 
					{gsub("*","",$2); 
					split($2,a,":"); 
					if(a[2] !~ /[^0-9]/) print "HLA-" a[1] ":" a[2]}; 
				if (index($3,"*")!=0) 
				{gsub("*","",$3); 
				split($3,b,":"); 
				if(b[2] !~ /[^0-9]/) print "HLA-" b[1] ":" b[2]; }}; 
			if (NR>4) 
				{if (index($2,"*")!=0) 
					{split($2,a,":"); 
					split(a[1],aa,"*"); 
					if(a[2] !~ /[^0-9]/) print aa[1] "_" aa[2] a[2]} ;  
				if (index($3,"*")!=0) 
					{split($3,b,":"); 
					split(b[1],bb,"*"); 
					if(b[2] !~ /[^0-9]/) print bb[1] "_" bb[2] b[2]}}}' | sort | uniq | grep 'HLA\|DRB'| awk -v fa=$peptideFasta -v OUTPUT_DIR=$OUTPUT_DIR -v netMHCI=$netMHCpan_SIF -v netMHCII=$netMHCIIpan_SIF -v pepType=$pepType '{if (substr($1,1,3)=="HLA") print "singularity exec -B \"${SCRATCHDIR}/${LSB_JOBID}:/tmp\" -B \"" fa ":/input/peptideFasta.fa:ro\" -B \"" OUTPUT_DIR ":/output\" " netMHCI " sh -c \"/netMHCpan-4.1/netMHCpan -f /input/peptideFasta.fa -a " $1 " -l 8,9,10,11 > /output/netMHCI_" $i "_" pepType"\""; else print "singularity exec -B \"${SCRATCHDIR}/${LSB_JOBID}:/tmp\" -B \"" fa ":/input/peptideFasta.fa:ro\" -B \"" OUTPUT_DIR ":/output\" " netMHCII " sh -c \"/netMHCIIpan-4.0/netMHCIIpan -f /input/peptideFasta.fa -a " $1 " > /output/netMHCII_" $i "_" pepType "\"" }' | xargs -P 8 -n 1 -d '\n' -I {} sh -c {}

#					if(b[2] !~ /[^0-9]/) print bb[1] "_" bb[2] b[2]}}}' | sort | uniq | grep 'HLA\|DRB'| awk -v fa=$peptideFasta -v OUTPUT_DIR=$OUTPUT_DIR -v netMHCI=$netMHCpan -v netMHCII=$netMHCIIpan -v pepType=$pepType '{if (substr($1,1,3)=="HLA") print netMHCI " -tdir " OUTPUT_DIR " -a " $1 " -l 8,9,10,11 -f " fa " > " OUTPUT_DIR "/netMHCI_" $1 "_" pepType; else print netMHCII " -tdir " OUTPUT_DIR " -a " $1 " -f " fa " > " OUTPUT_DIR "/netMHCII_"$1 "_" pepType}' | xargs -P 8 -n 1 -I {} sh -c {}

### 17 Feb,2017. It seems the servers currently cannot generate proper temporary folder. I deleted the option " -tdir " 



