#! /bin/bash

usage() { echo "Usage: $0 -i PID -o defuse_result_folder
    -i: Sample ID, same as in previous defuse run
    -o: Path to defuse result folder
>>> Note: Absolute path needed" 1>&2; exit 1; }



while getopts ":i:o:" opt; do
    case "${opt}" in
        i)
            sampleID=${OPTARG}
            ;;
        o)
            defuseFolder=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${sampleID}" ] || [ -z "${defuseFolder}" ] ; then
    usage
fi


module load perl/5.24.1
module load R/3.4.0


# defuseFolder=/icgc/dkfzlsdf/analysis/D120/yanhong/test_fusion/test2
# sampleID=libraryNam

perl /icgc/dkfzlsdf/analysis/D120/tools/deFuse/defuse-0.6.1/scripts/annotate_fusions.pl -c /icgc/dkfzlsdf/analysis/D120/tools/deFuse/defuse-0.6.1/scripts/config.txt -o $defuseFolder -n $sampleID > $defuseFolder/annotations

perl /icgc/dkfzlsdf/analysis/D120/tools/deFuse/defuse-0.6.1/scripts/coallate_fusions.pl -c /icgc/dkfzlsdf/analysis/D120/tools/deFuse/defuse-0.6.1/scripts/config.txt -o $defuseFolder -l $defuseFolder/annotations > $defuseFolder/results.tsv

Rscript /icgc/dkfzlsdf/analysis/D120/tools/deFuse/defuse-0.6.1/scripts/run_adaboost.R /lsdf/dkfz/analysis/D120/tools/deFuse/defuse-0.6.1/data/controls.txt $defuseFolder/results.tsv $defuseFolder/results.classify.tsv

perl /icgc/dkfzlsdf/analysis/D120/tools/deFuse/defuse-0.6.1/scripts/filter.pl probability '> 0.5' < $defuseFolder/results.classify.tsv >  $defuseFolder/results.filtered.tsv


#######
######
## because module load cannot work properly in cluster, the following step are manually done

# perl /icgc/dkfzlsdf/analysis/D120/tools/deFuse/defuse-0.6.1/scripts/annotate_fusions.pl -c /icgc/dkfzlsdf/analysis/D120/tools/deFuse/defuse-0.6.1/scripts/config.txt -o /icgc/dkfzlsdf/analysis/D120/huangz/ACC_Yang_Nature_medicine/results_per_pid/ACC_M9_SRR8752373_mRNA/deFuse -n ACC_M9_SRR8752373_mRNA > /icgc/dkfzlsdf/analysis/D120/huangz/ACC_Yang_Nature_medicine/results_per_pid/ACC_M9_SRR8752373_mRNA/deFuse/annotations.tmp

# mv /icgc/dkfzlsdf/analysis/D120/huangz/ACC_Yang_Nature_medicine/results_per_pid/ACC_M9_SRR8752373_mRNA/deFuse/annotations.tmp /icgc/dkfzlsdf/analysis/D120/huangz/ACC_Yang_Nature_medicine/results_per_pid/ACC_M9_SRR8752373_mRNA/deFuse/annotations

# perl /icgc/dkfzlsdf/analysis/D120/tools/deFuse/defuse-0.6.1/scripts/coallate_fusions.pl -c /icgc/dkfzlsdf/analysis/D120/tools/deFuse/defuse-0.6.1/scripts/config.txt -o /icgc/dkfzlsdf/analysis/D120/huangz/ACC_Yang_Nature_medicine/results_per_pid/ACC_M9_SRR8752373_mRNA/deFuse -l /icgc/dkfzlsdf/analysis/D120/huangz/ACC_Yang_Nature_medicine/results_per_pid/ACC_M9_SRR8752373_mRNA/deFuse/annotations > /icgc/dkfzlsdf/analysis/D120/huangz/ACC_Yang_Nature_medicine/results_per_pid/ACC_M9_SRR8752373_mRNA/deFuse/results.tsv


# Rscript /icgc/dkfzlsdf/analysis/D120/tools/deFuse/defuse-0.6.1/scripts/run_adaboost.R /lsdf/dkfz/analysis/D120/tools/deFuse/defuse-0.6.1/data/controls.txt /icgc/dkfzlsdf/analysis/D120/huangz/ACC_Yang_Nature_medicine/results_per_pid/ACC_M9_SRR8752373_mRNA/deFuse/results.tsv /icgc/dkfzlsdf/analysis/D120/huangz/ACC_Yang_Nature_medicine/results_per_pid/ACC_M9_SRR8752373_mRNA/deFuse/results.classify.tsv


# perl /icgc/dkfzlsdf/analysis/D120/tools/deFuse/defuse-0.6.1/scripts/filter.pl probability '> 0.5' < /icgc/dkfzlsdf/analysis/D120/huangz/ACC_Yang_Nature_medicine/results_per_pid/ACC_M9_SRR8752373_mRNA/deFuse/results.classify.tsv >  /icgc/dkfzlsdf/analysis/D120/huangz/ACC_Yang_Nature_medicine/results_per_pid/ACC_M9_SRR8752373_mRNA/deFuse/results.filtered.tsv
