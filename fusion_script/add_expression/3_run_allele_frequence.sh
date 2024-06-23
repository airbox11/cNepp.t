#! /bin/bash

## Be note: X	49098545	0	C	T, Alternative T in DNA, in RNA bp A was identified. This script reports as "na" because A cannot be found.
## This should be improved in future work. A and T, G and C should be biologically equall in RNA-seq level.

file=$1

        awk 'BEGIN {OFS="\t";FS="\t"}
                        {if (NR==1) print $0,"dna_freAlt","dna_cov","freRef","freAlt","expAlt";
                         if (NR>1) { freAlt="na";freRef="na";altExp="no";
                                if ($5==$(NF-10)) {
                                        if ($5=="A") { if ($(NF-9)!=0) { refExp="true"; freRef=$(NF-5) } }
                                        if ($5=="T") { if ($(NF-8)!=0) { refExp="true"; freRef=$(NF-4) } }
                                        if ($5=="C") { if ($(NF-7)!=0) { refExp="true"; freRef=$(NF-3) } }
                                        if ($5=="G") { if ($(NF-6)!=0) { refExp="true"; freRef=$(NF-2) } }
                                        if ($6=="A") { if ($(NF-9)!=0) { altExp="true"; freAlt=$(NF-5) } else { nonAltExp="true" } }
                                        if ($6=="T") { if ($(NF-8)!=0) { altExp="true"; freAlt=$(NF-4) } else { nonAltExp="true" } }
                                        if ($6=="C") { if ($(NF-7)!=0) { altExp="true"; freAlt=$(NF-3) } else { nonAltExp="true" } }
                                        if ($6=="G") { if ($(NF-6)!=0) { altExp="true"; freAlt=$(NF-2) } else { nonAltExp="true" } }
                                        split($0,a,"DP4=");
                                        split(a[2],b,";");
                                        split(b[1],c,",");
                                        if (length(c[1])!=0) dna_cov=c[1]+c[2]+c[3]+c[4]; 
					else {dna_cov="NA";dna_freAlt="NA"}
					if (length(a[1])==0) {dna_cov="NA";dna_freAlt="NA"}
                                        if (dna_cov==0) {print  "Note: Wrong line " $0;}
                                        else { if (dna_cov!="NA") dna_freAlt=(c[3]+c[4])/(dna_cov);}
                                        if (altExp=="true") expAlt="Y";else expAlt="N";
					if (freAlt==1) freRef=0;
					if (freRef==1) freAlt=0;
                                        print $0,dna_freAlt,dna_cov,freRef,freAlt,expAlt;
                                        }
                                else {print "wrong:" $0}
                        }
                }' $file



<<comment
file=$1

        awk 'BEGIN {OFS="\t";FS="\t"}
                        {if (NR==1) print $0,"dna_freAlt","dna_cov","freRef","freAlt","expAlt";
			 if (NR>1) { 
                                if ($5=$50) {
                                        if ($5=="A") { if ($51!=0) { refExp="true"; freRef=$55 } }
                                        if ($5=="T") { if ($52!=0) { refExp="true"; freRef=$56 } }
                                        if ($5=="C") { if ($53!=0) { refExp="true"; freRef=$57 } }
                                        if ($5=="G") { if ($54!=0) { refExp="true"; freRef=$58 } }
                                        if ($6=="A") { if ($51!=0) { altExp="true"; freAlt=$55 } else { nonAltExp="true" } }
                                        if ($6=="T") { if ($52!=0) { altExp="true"; freAlt=$56 } else { nonAltExp="true" } }
                                        if ($6=="C") { if ($53!=0) { altExp="true"; freAlt=$57 } else { nonAltExp="true" } }
                                        if ($6=="G") { if ($54!=0) { altExp="true"; freAlt=$58 } else { nonAltExp="true" } }
                                        split($9,a,"DP4=");
                                        split(a[2],b,";");
                                        split(b[1],c,",");
                                        dna_cov=c[1]+c[2]+c[3]+c[4]
                                        if (dna_cov==0) { print "Note: Wrong line " $0;}
                                        dna_freAlt=(c[3]+c[4])/(dna_cov);
					if (altExp=="true") expAlt="Y";else expAlt="N";
					print $0,dna_freAlt,dna_cov,freRef,freAlt,expAlt;
                                        }
                                else {print "wrong"}
                                freAlt="na";freRef="na";altExp="no"
                        }
                }' $file

comment
#print $2,$3,$4,$5,dna_freAlt,dna_cov,$34,$7,$8,freRef,freAlt,$49,$50,$51,$52,$57,$58;
#if (NR=1) print $0,"dna_freAlt","dna_cov","freRef","freAlt";
