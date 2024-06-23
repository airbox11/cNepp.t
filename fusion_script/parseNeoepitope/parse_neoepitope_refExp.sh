# this script to select vaccinated candidates
# for patients without RNA seq samples

usage() { echo "Usage: $0 -n neoepitopeFile -t mhcType -a affinity -e expression -b wtBind
    -n: path to neoepitope result file
    -t: 1 or 2 indicating mhc-I or mhc-II
    -a: affinity threshold, e.g. 2000.
    -e: minimum gene expression value
    -b: wild type binding level.  nonB indicates wildtype epitope are non binding level. Must be nonB or empty
>>> Note: Absolute path needed" 1>&2; exit 1; }

while getopts ":n:t:a:e:b:" o; do
    case "${o}" in
        n)
            neoFile=${OPTARG}
            ;;
        t)
            mhcType=${OPTARG}
            ;;
        a)
            affnM=${OPTARG}
            ;;
        e)
            expV=${OPTARG}
            ;;
        b)
            wtBind=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ ! -f "${neoFile}" ] || [ -z "${mhcType}" ] || [ -z "${affnM}" ] || [  -z "${expV}" ] ; then
    usage
fi

### For samples without RNA-seq data
# Class I
# results_GM60_MHCI_epitopes_filtered_addExp.tab
if [ "$mhcType" = "1" ];
then
	awk 'BEGIN{OFS="\t";FS="\t"} {if (NR==1) print $1,$2,$3,$4,$14,$16,$29,$31,$34,$36,$37,$(NF-4),$(NF-3),$(NF-2),$(NF-1),$NF,"lengthVaccinatedPeptide"}' $neoFile
	if [ "$wtBind" = "nonB" ]; 
	then		
		awk -v affnM1=$affnM -v expV1=$expV 'BEGIN{OFS="\t"} {if (NR>1 && $14<int(affnM1) && $NF>expV1) print $1,$2,$3,$4,$14,$16,$29,$31,$34,substr($36,3,25),substr($37,3,25),$(NF-4),$(NF-3),$(NF-2),$(NF-1),$NF, length(substr($36,3,25))}' $neoFile | grep "nonB" | sort -nrk16,16
	else
		awk -v affnM1=$affnM -v expV1=$expV 'BEGIN{OFS="\t"} {if (NR>1 && $14<int(affnM1) && $NF>expV1) print $1,$2,$3,$4,$14,$16,$29,$31,$34,substr($36,3,25),substr($37,3,25),$(NF-4),$(NF-3),$(NF-2),$(NF-1),$NF, length(substr($36,3,25))}' $neoFile | sort -nrk16,16
	fi
fi
# Class II


if [ "$mhcType" = "2" ];
then
	awk 'BEGIN{FS="\t";OFS="\t"} {if (NR==1) print $1,$2,$3,$4,$10,$13,$22,$25,$28,$30, $31, $(NF-4),$(NF-3),$(NF-2),$(NF-1),$NF, "lengthVaccinatedPeptide"}' $neoFile
	if [ "$wtBind" = "nonB" ];
	then
		awk -v affnM1=$affnM -v expV1=$expV 'BEGIN{FS="\t";OFS="\t"} {if (NR>1 && $10<int(affnM1) && $NF>expV1) print $1,$2,$3,$4,$10,$13,$22,$25,$28,substr($30,3,25), substr($31,3,25),$(NF-4), $(NF-3),$(NF-2),$(NF-1),$NF, length(substr($31,3,25))}' $neoFile | grep "nonB" | sort -nrk15,15
	else
		awk -v affnM1=$affnM -v expV1=$expV 'BEGIN{FS="\t";OFS="\t"} {if (NR>1 && $10<int(affnM1) && $NF>expV1) print $1,$2,$3,$4,$10,$13,$22,$25,$28,substr($30,3,25), substr($31,3,25),$(NF-4), $(NF-3),$(NF-2),$(NF-1),$NF, length(substr($31,3,25))}' $neoFile | sort -nrk15,15
	fi
fi
