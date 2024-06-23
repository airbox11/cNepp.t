workDir=$1



# cp ${workDir}/1.1_hla_typing/phlat/*_HLA.sum ${workDir}/1.1_hla_typing/

input=${workDir}/1.1_hla_typing/$phlat_dir/*_HLA.sum
output=${workDir}/1.2_format_hla_result/format_hla_phlat
script=/icgc/dkfzlsdf/analysis/D120/scripts/HLA_typing/run_format_phlat.sh
sh $script -p $input -o $output

